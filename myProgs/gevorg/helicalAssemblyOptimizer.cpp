#include <stdio.h>
#include <iostream>
#include <fstream>
#include <time.h>

#include "Atom.h"
#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "MslTools.h"
#include "release.h"
#include "Transforms.h"
#include "OptimalRMSDCalculator.h"
#include "Frame.h"
#include "OptionParser.h"
#include "glab/ProximitySearch.h"
#include <map>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_multimin.h>
#include <math.h>

using namespace std;
using namespace MSL;

void assert(bool cond, string msg) {
  if (!cond) {
    cout << msg << endl;
    exit(-1);
  }
}

class options {
  public:
    int numUnits, numCycles, numIters;
    string type, pdbFile, motifsFile, filteredMotifsFile, outPdbFile;
};

string option(string opt, string mes, int w, int p1, int p2) {
  // first print the name of the option
  string text(p1, ' ');
  text += opt;
  if (p2 > text.size()) text += string(p2 - text.size(), ' ');

  // next print the description text
  int i = 0, k, L = text.size(), n;
  while (i < mes.size()) {
    k = mes.find_first_of(" ", i);
    if (k == string::npos) k = mes.size();
    n = k - i;
    if ((L + n >= w) && (L > 0)) { text += "\n" + string(p2, ' '); L = p2; }
    text += mes.substr(i, n) + " ";
    L += n + 1;
    i = k+1;
  }
  return text;
}

void usage() {
  int w = 80, p1 = 3, p2 = p1+12;
  cout << endl << option("", "Builds cool assemblies using motifs. Options:", w, 0, 0) << endl;
  cout << option("--pdb", "input PDB file.", w, p1, p2) << endl;
  cout << option("--motifs", "file with a list of motif splits to consider for scoring interfaces.", w, p1, p2) << endl;
  cout << option("--cycles", "number of independent optimization cycles to run.", w, p1, p2) << endl;
  cout << option("--iters", "number of iterations per cycle.", w, p1, p2) << endl;
  cout << option("--opdb", "output PDB file. Name will be modified to include cycle number.", w, p1, p2) << endl;
  cout << option("--type", "optimization type. Currently 'helix' or 'circle' are possible.", w, p1, p2) << endl;
  cout << option("--units", "optional: number of units per turn: interpreted as the starting value for helices and fixed for circles.", w, p1, p2) << endl;
  cout << option("--outMotifs", "optional: file name where surviving motif splits will be written. If provided, will quit after writing this file.", w, p1, p2) << endl;
}

void parseCommandLine(int argc, char* argv[], options& opt) {
  vector<string> required, optional;
  required.push_back("pdb");
  required.push_back("motifs");
  required.push_back("iters");
  required.push_back("cycles");
  required.push_back("opdb");
  required.push_back("type");
  optional.push_back("units");
  optional.push_back("outMotifs");

  OptionParser op;
  op.readArgv(argc, argv);
  op.setRequired(required);
  op.setAllowed(optional);

  if (op.countOptions() == 0) { usage(); exit(0); }

  opt.pdbFile = op.getString("pdb");
  assert(!op.fail(), "PDB file not specified!");
  opt.motifsFile = op.getString("motifs");
  assert(!op.fail(), "Motifs file not specified!");
  opt.outPdbFile = op.getString("opdb");
  assert(!op.fail(), "Output PDB file name not specified!");
  opt.numIters = op.getInt("iters");
  assert(!op.fail(), "Number of iterations not specified!");
  opt.numCycles = op.getInt("cycles");
  assert(!op.fail(), "Number of cycles not specified!");
  opt.type = op.getString("type");
  assert(!op.fail(), "Optimization type not specified!");
  assert((opt.type == "helix") || (opt.type == "circle"), "Unknown optimization type specified!");
 
  opt.numUnits = op.getInt("units");
  if(op.fail()) { opt.numUnits = 6; }
  opt.filteredMotifsFile = op.getString("outMotifs");
  if(op.fail()) opt.filteredMotifsFile = "";
}


bool fileExists(const string filename) {
  struct stat buffer ;
  if (stat(filename.c_str(), &buffer) == 0) return true;
  return false;
}

string parentPath(string path, string& fileName) {
  int pos = path.rfind("/");
  if (pos == string::npos) {
    fileName = path;
    return "./";
  } else if (pos == 0) {
    fileName = path.substr(1);
    return "/";
  } else {
    fileName = path.substr(pos+1);
    return path.substr(0, pos);
  }
}

double min(double a, double b) {
  if (a < b) return a;
  return b;
}

double max(double a, double b) {
	if (a > b) return a;
	return b;
}

void openFileCPP(ofstream& ofs, const string filename, std::ios_base::openmode mode = ios_base::out) {
	ofs.open(filename.c_str(), mode);
	if (ofs.fail()) {
		printf("error: failed to open file '%s'\n", filename.c_str());
		exit(-1);
	}
}

template <class T>
string stringifyVector(vector<T>& vec, string del = " ") {
	stringstream ss;
	for (int i = 0; i < vec.size(); i++) {
		ss << vec[i];
		if (i != vec.size()-1) ss << del;
	}
	return ss.str();
}

double calculateRMSDCutoff(vector<int> L) {
  double rmsdMax = 0.8;
  int L0 = 10;
  double a = exp(-1./L0);
  int N = 0, n;
  double c = 0;

  // disjoint segments are counted as independent, so their correlation
  // with respect to each other is zero
  for (int i = 0; i < L.size(); i++) {
    N += L[i];
    n = L[i];
    c = c + (a/(1-a))*(n-1) - pow((a/(1-a)), 2)*(1 - pow(a, n-1));
  }
  double df = N*(1 - (2.0/(N*(N-1)))*c);

  return rmsdMax/sqrt(N/df);
}

void splitBackboneByResidue(AtomPointerVector& atoms, vector<vector<int> >& residues) {
  int k;
  // split atoms by residue
  residues.clear();
  for (int i = 0; i < atoms.size(); i++) {
    if ( (i == 0) || (atoms[i]->getResidueNumber() != atoms[i-1]->getResidueNumber()) || 
         (atoms[i]->getChainId().compare(atoms[i-1]->getChainId()) != 0) ||
         (atoms[i]->getResidueIcode().compare(atoms[i-1]->getResidueIcode()) != 0) ) {
      residues.push_back(vector<int>(4, -1));
    }
    if (atoms[i]->getName().compare("N") == 0) k = 0;
    else if (atoms[i]->getName().compare("CA") == 0) k = 1;
    else if (atoms[i]->getName().compare("C") == 0) k = 2;
    else if (atoms[i]->getName().compare("O") == 0) k = 3;
    else continue;
    residues.back()[k] = i;
  }
}


bool copyBackbone(AtomContainer& from, AtomPointerVector& to, vector<vector<int> >* residues = NULL) {
	string an;
	for (int i = 0; i < from.atomSize(); i++) {
		an = from[i].getName();
		if (!an.compare("N") || !an.compare("CA") || !an.compare("C") || !an.compare("O")) {
			to.push_back(new Atom(from[i]));
		}
		if (!an.compare("N")) to.back()->setRadius(1.6);
		if (!an.compare("CA")) to.back()->setRadius(2.365);
		if (!an.compare("C")) to.back()->setRadius(2.1);
		if (!an.compare("O")) to.back()->setRadius(1.6);
	}

	if (residues != NULL) {
    splitBackboneByResidue(to, *residues);
    // make sure all necessary backbone atoms were found
    for (int i = 0; i < residues->size(); i++) {
      for (int j = 0; j < (*residues)[i].size(); j++) {
        if ((*residues)[i][j] == -1) {
          cout << "error: not all heavy backbone atoms found by copyBackbone in the following:\n" << from << endl;
          cout << "specifically, backbone atom " << j+1 << " in residue " << i+1 << " ended up not found\n";
          return false;
        }
      }
    }
	}
	return true;
}

void copyResidueAtoms(AtomPointerVector& unit, int beg, int len, vector<vector<int> >& residues, vector<Atom*>& sub) {
  int k = 0;
  for (int i = beg; i < beg + len; i++) {
    vector<int>& res = residues[i];
    for (int j = 0; j < res.size(); j++) {
      sub[k] = unit[res[j]]; k++;
    }
  }
}

void concatenateResidueAtoms(AtomPointerVector& unit, int beg, int len, vector<vector<int> >& residues, vector<Atom*>& sub) {
  for (int i = beg; i < beg + len; i++) {
    vector<int>& res = residues[i];
    for (int j = 0; j < res.size(); j++) {
      sub.push_back(unit[res[j]]);
    }
  }
}

void concatenateResidueAtoms(AtomPointerVector& unit, vector<int>& indices, vector<vector<int> >& residues, vector<Atom*>& sub) {
  for (int i = 0; i < indices.size(); i++) {
    vector<int>& res = residues[indices[i]];
    for (int j = 0; j < res.size(); j++) {
      sub.push_back(unit[res[j]]);
    }
  }
}

class MotifSplit {
	private:
		vector<vector<int> > parts;   // indices in the motif that correspond to the first and second "subs"
		vector<vector<int> > algn;    // possible alignment locations for the two subs in the target
		vector<vector<double> > rmsd; // RMSDs corresponding to each plausible alignment of each sub
		vector<vector<int> > match;   // matching alignments (those that pass the joint subA and subB cutoff)
		vector<double> rCut;          // RMSD cutoffs for counting a match for each sub

	public:
		MotifSplit() {}
		MotifSplit(vector<int>& _partA, vector<int>& _partB);
		void addAlgnPos(int pidx, int targetIdx, double _rmsd);
		void setRMSDCut(int pidx, double cut) { rCut[pidx] = cut; }
		int numAlgnPos(int pidx) { return algn[pidx].size(); }
		int getAlgnPos(int pidx, int algnIdx) { return algn[pidx][algnIdx]; }
    double getAlgnRMSD(int pidx, int algnIdx) { return rmsd[pidx][algnIdx]; }
		vector<int>& getAlgnPositions(int pidx) { return algn[pidx]; }
		vector<int>& getSub(int pidx) { return parts[pidx]; }
		int length(int pidx) { return parts[pidx].size(); }
		double getRMSDCut(int pidx) { return rCut[pidx]; }
};

MotifSplit::MotifSplit(vector<int>& _partA, vector<int>& _partB) {
	parts.push_back(_partA);
	parts.push_back(_partB);
	rCut.resize(2, 0);
	algn.resize(2);
	rmsd.resize(2);
	match.resize(2);
}

void MotifSplit::addAlgnPos(int pidx, int targetIdx, double _rmsd) {
	algn[pidx].push_back(targetIdx);
	rmsd[pidx].push_back(_rmsd);
}

class Motif {
	private:
		AtomPointerVector atoms;
		vector<vector<int> > residues;
		vector<MotifSplit> splits;
		double rCut;
		string name, pdbfile;

	public:
		Motif() { rCut = -1; }
    ~Motif();
		void addSplit(vector<int>& partA, vector<int>& partB);
    void removeSplit(int _si);
		bool assignRMSD();
		bool readAtoms(string pdbFile);
		void setName(string _name) { name = _name; }
		void setPDB(string _pdbfile) { pdbfile = _pdbfile; }
		int numSplits() { return splits.size(); }
		int length() { return atoms.size(); }
		int numResidues() { return residues.size(); }
		string getName() { return name; }
		string getPDB() { return pdbfile; }
		double getRMSDCut() { return rCut; }
		MotifSplit& getSplit(int _si) { return splits[_si]; }
		AtomPointerVector getSub(int _si, int partID);
		AtomPointerVector& getAtoms() { return atoms; }
};

Motif::~Motif() {
  for (int i = 0; i < atoms.size(); i++) delete(atoms[i]);
}

bool Motif::readAtoms(string pdbFile) {
	// clear anything previously
	for (int i = 0; i < atoms.size(); i++) delete(atoms[i]);
	atoms.clear();

  // read atoms, copy backbone, and split by residue
  bool tmpFile = false;
  AtomContainer t;
  if (!fileExists(pdbFile)) {
    // try seeing if the parent directory has been archived
    string pdbFileName;
    string tarDir = parentPath(pdbFile, pdbFileName) + ".tar.gz";
    assert(fileExists(tarDir), "could not find '" + pdbFile + "' or its archived parent directory.");
    pdbFile = "/tmp/mot" + MslTools::intToString((int) getpid()) + ".pdb";
    string cmd = "tar -xzOf " + tarDir + " *" + pdbFileName + " > " + pdbFile;
    system(cmd.c_str());
    tmpFile = true;
  }
	assert(t.readPdb(pdbFile), "could not read unit PDB file '" + pdbFile + "'");
  if (tmpFile) remove(pdbFile.c_str());
	bool suc = copyBackbone(t, atoms, &residues);
  if (suc) {
    // map backbone atoms in each residue
    vector<map<string, Atom*> > resBB(residues.size());
    for (int i = 0; i < residues.size(); i++) {
      for (int j = 0; j < residues[i].size(); j++) {
        Atom* a = atoms[residues[i][j]];
        resBB[i][a->getName()] = a;
      }
    }
    return true;
  } else {
    return false;
  }
}

void Motif::addSplit(vector<int>& partA, vector<int>& partB) {
	// make sure residue indices for the split are not out of range
	for (int k = 0; k < 2; k++) {
		vector<int>& part = (k == 0 ? partA : partB);
		for (int j = 0; j < part.size(); j++) {
			if ((part[j] < 0) || (part[j] >= residues.size())) {
				cout << "error: split indices out of range for motif " << name << ", which was found to have " << residues.size() << " residues" << endl;
				exit(-1);
			}
		}
	}
	splits.push_back(MotifSplit(partA, partB));

  // helper map for backbone atoms in each residue
  vector<map<string, Atom*> > resBB(residues.size());
  for (int i = 0; i < residues.size(); i++) {
    for (int j = 0; j < residues[i].size(); j++) {
      Atom* a = atoms[residues[i][j]];
      resBB[i][a->getName()] = a;
    }
  }
  // calculate the RMSD cutoff for the segments of this split
  MotifSplit& sp = splits.back();
  vector<int> segsT;       // lengths of segments in both subs
  for (int k = 0; k < 2; k++) {
    vector<int> segs;    // lengths of segments in this subs
    vector<int> breaks;
    vector<int>& sub = sp.getSub(k);
    for (int j = 0; j < sub.size(); j++) {
      if ((j == sub.size() - 1) || (resBB[sub[j]]["C"]->distance(*(resBB[sub[j+1]]["N"])) > 2.0)) {
        breaks.push_back(j);
      }
    }
    segs.push_back(breaks[0] + 1);
    int len = segs.back();
    for (int j = 1; j < breaks.size(); j++) {
      segs.push_back(breaks[j] - len + 1);
      len += segs.back();
    }
    sp.setRMSDCut(k, calculateRMSDCutoff(segs));
    segsT.insert(segsT.end(), segs.begin(), segs.end());
  }
  // now do the combined RMSD cutoff
  if (splits.size() == 1) rCut = calculateRMSDCutoff(segsT);
  else assert(fabs(getRMSDCut() - calculateRMSDCutoff(segsT)) < 0.0001, "for segment '" + name + "', different splits give different RMSD cutoffs (probably have different segment lengths)!");
}

void Motif::removeSplit(int _si) {
  assert((_si >= 0) && (_si < splits.size()), "Motif::removeSplit split index " + MslTools::intToString(_si) + " out of range.");
  splits.erase(splits.begin() + _si);
}

bool Motif::assignRMSD() {
	// map backbone atoms in each residue
	vector<map<string, Atom*> > resBB(residues.size());
	for (int i = 0; i < residues.size(); i++) {
		for (int j = 0; j < residues[i].size(); j++) {
			Atom* a = atoms[residues[i][j]];
			resBB[i][a->getName()] = a;
		}
	}
	// we assume that the two subs are already disjoint, so
	// we just need to count the number of segments in each
	for (int i = 0; i < splits.size(); i++) {
		vector<int> segsT;       // lengths of segments in both subs
		for (int k = 0; k < 2; k++) {
			vector<int> segs;    // lengths of segments in this subs
			vector<int> breaks;
			vector<int>& sub = splits[i].getSub(k);
			for (int j = 0; j < sub.size(); j++) {
				if ((j == sub.size() - 1) || (resBB[sub[j]]["C"]->distance(*(resBB[sub[j+1]]["N"])) > 2.0)) {
					breaks.push_back(j);
				}
			}
			segs.push_back(breaks[0] + 1);
			int len = segs.back();
			for (int j = 1; j < breaks.size(); j++) {
				segs.push_back(breaks[j] - len + 1);
				len += segs.back();
			}
			splits[i].setRMSDCut(k, calculateRMSDCutoff(segs));
			segsT.insert(segsT.end(), segs.begin(), segs.end());
		}
		// now do the combined RMSD cutoff
		if (i == 0) { rCut = calculateRMSDCutoff(segsT); }
		if (fabs(getRMSDCut() - calculateRMSDCutoff(segsT)) > 0.0001) {
			cout << "for segment '" << name << "', different splits appear to produce different RMSD cutoffs (probably have different segment lengths)!" << endl;
			return false;
		}
	}
	return true;
}

AtomPointerVector Motif::getSub(int _si, int pidx) {
	vector<int>& subIndices = splits[_si].getSub(pidx);
	AtomPointerVector sub;
	concatenateResidueAtoms(atoms, subIndices, residues, sub);
	return sub;
}

class ProblemDef {
	public:
		vector<Motif*>* motifs;
		AtomPointerVector* unit;
		vector<vector<int> > residues;
		bool scoreClashes;
    int numUnitsToGenerate;
    double circRadius;
    map<string, double> fixedParams;
};

class solution {
  public:
    solution() { R = w = P = a = b = g = score = 999.99; }
    double R, w, P, a, b, g, score;
    vector<AtomPointerVector> units;
};

void transformUnit(vector<AtomPointerVector>& res, double xrot, double yrot, double zrot, double xtran, double ytran, double ztran) {
	Transforms tr;
	CartesianPoint trans(xtran, ytran, ztran);
	
	for (int i = 0; i < res.size(); i++) {
		for (int j = 0; j < res[i].size(); j++) {
			tr.Xrotate(*(res[i][j]), xrot);
			tr.Yrotate(*(res[i][j]), yrot);
			tr.Zrotate(*(res[i][j]), zrot);
			tr.translate(*(res[i][j]), trans);
		}
	}
}

void setChainId(AtomPointerVector& atoms, string cid) {
	for (int i = 0 ; i < atoms.size(); i++) atoms[i]->setChainId(cid);
}

int numberOfAtomPairsWithin(AtomPointerVector& unitA, AtomPointerVector& unitB, vector<vector<int> >& residues, double dcut) {
  int c = 0;
  for (int riA = 0; riA < residues.size(); riA++) {
    for (int aiA = 0; aiA < residues[riA].size(); aiA++) {
      int a = residues[riA][aiA];
      for (int riB = 0; riB < residues.size(); riB++) {
        for (int aiB = 0; aiB < residues[riB].size(); aiB++) {
          int b = residues[riB][aiB];
          if (unitA[a]->distance(*(unitB[b])) < dcut) c++;
        }
      }
    }
  }
  return c;
}

// evaluate the designability of a pair of units
double evaluatePose(vector<Motif*>& motifs, AtomPointerVector& unit1, AtomPointerVector& unit2, vector<vector<int> >& residues, bool scoreClashes = false) {
  // NOTE: could consider a quick check with ProximitySearch to see if these units are even worth looking into
  OptimalRMSDCalculator orc;
  double score = 0;
  bool suc = true;
  // look for all combinations of plausible sub alignments and see
  // which motifs fully pass under their respective RMSD cutoffs
  for (int i = 0; i < motifs.size(); i++) {
    Motif& m = *(motifs[i]);
    for (int j = 0; j < m.numSplits(); j++) {
      MotifSplit& sp = m.getSplit(j);
      vector<int>& algnA = sp.getAlgnPositions(0);
      vector<int>& algnB = sp.getAlgnPositions(1);
      vector<Atom*> motA = m.getSub(j, 0);
      vector<Atom*> motB = m.getSub(j, 1);
      vector<Atom*> motAB(motA); motAB.insert(motAB.end(), motB.begin(), motB.end());
      int lenA = sp.length(0);
      int lenB = sp.length(1);
      // each motif can fit by aligning sub A onto unit 1 and sub B onto unit 2, or sub A onto unit 2 and sub B onto unit 1
      for (int order = 0; order < 2; order++) {
        AtomPointerVector& unitA = (order ? unit1 : unit2);
        AtomPointerVector& unitB = (order ? unit2 : unit1);
        // iterate over options of the first sub
        for (int a = 0; a < algnA.size(); a++) {
          vector<Atom*> regA(motA.size(), NULL);
          copyResidueAtoms(unitA, algnA[a], lenA, residues, regA);
          // iterate over options of the second sub
          for (int b = 0; b < algnB.size(); b++) {
            vector<Atom*> regB(motB.size(), NULL);
            copyResidueAtoms(unitB, algnB[b], lenB, residues, regB);
            vector<Atom*> regAB(regA.size() + regB.size(), NULL);
            int ii = 0;
            for (int jj = 0; jj < regA.size(); ii++, jj++) regAB[ii] = regA[jj];
            for (int jj = 0; jj < regB.size(); ii++, jj++) regAB[ii] = regB[jj];
            // consider this combination of alignments
            double rmsd = orc.bestRMSD(regAB, motAB, &suc);
            assert(suc, "RMSD calculation failed!");
//            if (rmsd < m.getRMSDCut()) {
              // calculate the contribution of this match
//              score -= m.numResidues() / (rmsd + 0.1) - m.numResidues() / (m.getRMSDCut() + 0.1);
//              score -= m.numResidues() * exp(-max(rmsd - m.getRMSDCut(), 0));
//              score -= m.numResidues() * exp(-pow(max(rmsd - m.getRMSDCut(), 0), 2));
              score -= m.numResidues() / (1 + exp(10*(rmsd - m.getRMSDCut())));
//if (rmsd < 3) printf("RMSD = %f, cut = %f, name = %s, score now is %e\n", rmsd, m.getRMSDCut(), m.getName().c_str(), score);
//            }
          }
        }
      }
    }
  }
  
  if (scoreClashes) {
    double f = 0.8;
    ProximitySearch ps(unit1, 6.0);
    for (int i = 0; i < unit2.size(); i++) {
      Atom& a2 = *(unit2[i]);
      vector<int> closeOnes;
      ps.pointsWithin(a2.getCoor(), 0, 6.0, closeOnes);
      for (int j = 0; j < closeOnes.size(); j++) {
        Atom& a1 = *(unit1[closeOnes[j]]);
        if (a1.distance(a2) < f*(a1.getRadius() + a2.getRadius())) {
          score += 50*(f - a1.distance(a2)/(a1.getRadius() + a2.getRadius()));
        }
      }
    }
  }
  return score;
}

void cloneAtoms(AtomPointerVector& from, AtomPointerVector& to) {
  to.deletePointers();
  for (int i = 0; i < from.size(); i++) to.push_back(new Atom(*(from[i])));
}

// rad -- radius of a sphere circumscribing the unit
// If N is non-negative, will build the first N units
// If N is negative, will only preserve the first unit and all whose circumscribed spheres are within potential contacting distance.
void buildHelicalAssembly(AtomPointerVector& U, vector<AtomPointerVector>& units, double R, double w, double P, double alpha, double beta, double gamma, int N, double rad = -1) {
  double dZ = (P*w/360); // rise per unit
  w = w*M_PI/180;
  double dcut = (2*rad + 6.5);
  double dcut2 = dcut*dcut;
  Transforms tr;
  Frame O, Fo, Fi;
  // We shall define a frame associated with each unit in the helix. The X-axis in this frame
  // will be the vector from the helix axis to the unit center of mass, in the XY laboratory plane.
  // The Y-axis in this frame will be the tangent to the helical curve in the current point; note,
  // this is by definition orthogonal to the X-axis. The Z-axis will then be the cross product of X and Y.
  // X = [R*cos(wr*t) R*sin(wr*t), 0] (normalized)
  // Y = [-R*wr*sin(wr*t), R*wr*cos(wr*t), dZ]
  // Z = cross(X, Y)
  double b1 = 1/sqrt(1 + pow(dZ/(R*w), 2));
  double b2 = 1/sqrt(1 + pow((R*w)/dZ, 2));
  CartesianPoint Xo(1, 0, 0);
  CartesianPoint Yo(0, b1, b2);
  CartesianPoint Zo = Xo.cross(Yo);
  CoordAxes axo(Xo, Yo, Zo);
  Fo.computeFrameFromAxes(axo);
  CoordAxes xyz(CartesianPoint(1, 0, 0), CartesianPoint(0, 1, 0), CartesianPoint(0, 0, 1));
  O.computeFrameFromAxes(xyz);
  
  // rotate the original unit by a, b, g
  AtomPointerVector unit;
  cloneAtoms(U, unit);
  tr.translate(unit, -(unit.getGeometricCenter()));
  tr.Xrotate(unit, alpha);
  tr.Yrotate(unit, beta);
  tr.Zrotate(unit, gamma);
  
  // now generate clones around the helical curve
  for (int t = 0; t < ((N < 0) ?  min(100, ceil(fabs(dcut/dZ))) : N); t++) {
    double xo = R*cos(w*t); double yo = R*sin(w*t); double zo = dZ*t;
    if ((N < 0) && (t > 0)) {
      if (2*pow(R, 2)*(1 - cos(w*t)) + pow(dZ*t, 2) > dcut2) continue;
    }
    AtomPointerVector newUnit;
    cloneAtoms(unit, newUnit);

    // rotate and translate into helical frame
    CartesianPoint Xi(cos(w*t), sin(w*t), 0);
    CartesianPoint Yi(-b1*sin(w*t), b1*cos(w*t), b2);
    CartesianPoint Zi = Xi.cross(Yi);
    CoordAxes axi(Xi, Yi, Zi);
    Fi.computeFrameFromAxes(axi);
    Frame::transformAtoms(newUnit, Fo, O);
    Frame::transformAtoms(newUnit, O, Fi);
    tr.translate(newUnit, CartesianPoint(xo, yo, zo));

    // destructor of AtomPointerVector does not destroy the atoms, so when push_back makes a copy of the object,
    // atoms should persist even though they are not deep-copied
    units.push_back(newUnit);
  }

  // free space
  for (int i = 0; i < unit.size(); i++) delete(unit[i]);
}

void combineUnits(vector<AtomPointerVector>& units, AtomContainer& con) {
  char segid[10];
  for (int i = 0; i < units.size(); i++) {
    AtomPointerVector& unit = units[i];
    for (int j = 0; j < unit.size(); j++) {
      sprintf(segid, "%s%03d", unit[j]->getChainId().c_str(), i);
      Atom na(*(unit[j]));
      na.setChainId(segid);
      na.setSegID(segid);
      con.addAtom(na);
    }
  }
}

void freeUnits(vector<AtomPointerVector>& units) {
  for (int i = 0; i < units.size(); i++) {
    for (int j = 0; j < units[i].size(); j++) { delete(units[i][j]); }
  }
}

double evaluateHelix(vector<Motif*>& motifs, AtomPointerVector& U, vector<vector<int> >& residues, double R, double w, double P, double alpha, double beta, double gamma, int N, double rad, bool scoreClashes = false) {
  double score = 0;
  vector<AtomPointerVector> units;
  buildHelicalAssembly(U, units, R, w, P, alpha, beta, gamma, N, rad);
  int c = 0;
  for (int i = 1; i < units.size(); i++) {
    score += evaluatePose(motifs, units[0], units[i], residues, scoreClashes);
    c++;
  }
  if (c > 0) score /= c;

  // free memory
  freeUnits(units);

  return score;
}

double evaluateOptHelix(const gsl_vector* v, void* params) {
	Transforms tr;
	ProblemDef* p = (ProblemDef*) params;
	double R = gsl_vector_get(v, 0);
	double w = gsl_vector_get(v, 1);
	double P = gsl_vector_get(v, 2);
  double alpha = gsl_vector_get(v, 3);
  double beta = gsl_vector_get(v, 4);
  double gamma = gsl_vector_get(v, 5);

  return evaluateHelix(*(p->motifs), *(p->unit), p->residues, R, w, P, alpha, beta, gamma, p->numUnitsToGenerate, p->circRadius, p->scoreClashes);
}

double evaluateOptCircle(const gsl_vector* v, void* params) {
  Transforms tr;
  ProblemDef* p = (ProblemDef*) params;
  double R = gsl_vector_get(v, 0);
  double w = p->fixedParams["angle"];
  double P = 0;
  double alpha = gsl_vector_get(v, 1);
  double beta = gsl_vector_get(v, 2);
  double gamma = gsl_vector_get(v, 3);

  return evaluateHelix(*(p->motifs), *(p->unit), p->residues, R, w, P, alpha, beta, gamma, p->numUnitsToGenerate, p->circRadius, p->scoreClashes);
}

bool findProductiveMotifSplitAlignments(Motif& m, int splitIndex, AtomPointerVector& U, vector<vector<int> >& resU) {
  OptimalRMSDCalculator orc;
  bool suc;
  MotifSplit& sp = m.getSplit(splitIndex);
  for (int k = 0; k < 2; k++) {
    AtomPointerVector sub = m.getSub(splitIndex, k);
    for (int l = 0; l < (int) resU.size() - (int) sp.length(k) + 1; l++) {
      AtomPointerVector reg;
      concatenateResidueAtoms(U, l, sp.length(k), resU, reg);
      double rmsd = orc.bestRMSD(sub, reg, &suc);
      assert(suc, "RMSD calculation failed!");
      if (rmsd < sp.getRMSDCut(k)) {
        sp.addAlgnPos(k, l, rmsd);
      }
    }
    if (sp.numAlgnPos(k) == 0) return false;
  }
  return true;
}

// radius of sphere containing the atoms
double computeRadius(AtomPointerVector& A) {
  CartesianPoint cent = A.getGeometricCenter();
  double R = 0;
  for (int i = 0; i < A.size(); i++) {
    R = max(R, cent.distance((A[i])->getCoor()));
  }
  return R;
}

// radius of circle circumscribing the XY projection of atoms, after applying the given rotations around X, Y, and Z, respectively
double computeRadiusXY(AtomPointerVector& A, double alpha = 0, double beta = 0, double gamma = 0) {
  // apply the transformations
  Transforms tr;
  AtomPointerVector unit;
  cloneAtoms(A, unit);
  tr.translate(unit, -(unit.getGeometricCenter()));
  tr.Xrotate(unit, alpha);
  tr.Yrotate(unit, beta);
  tr.Zrotate(unit, gamma);

  // compute tha radius in XY plane
  CartesianPoint cent = A.getGeometricCenter();
  double R = 0;
  for (int i = 0; i < A.size(); i++) {
    R = max(R, pow((A[i])->getX() - cent.getX(), 2) + pow((A[i])->getY() - cent.getY(), 2));
  }
  return sqrt(R);
}

void parseMotifs(vector<Motif*>& motifs, string motifsFile, AtomPointerVector& U, vector<vector<int> >& resU, string filteredMotifsFile) {
  map<string, Motif*> M;
  Motif* m;
  vector<string> lines;
  assert(MslTools::readTextFile(lines, motifsFile), "could not read '" + motifsFile + "'");
  ofstream ofs;
  if (!filteredMotifsFile.empty()) openFileCPP(ofs, filteredMotifsFile.c_str());
  int nsplits = 0;
  for (int i = 0; i < lines.size(); i++) {
    if (i % 100 == 0) printf("reading %d/%d...\n", i+1, (int) lines.size());
    vector<string> entries = MslTools::tokenize(lines[i], " ");
    assert(entries.size() == 4, "could not parse line '" + lines[0] + "' from motif list file");
    bool prev = (M.find(entries[0]) != M.end());
    if (prev) {
      m = M[entries[0]];
    } else {
      m = new Motif();
      M[entries[0]] = m;
      if (!(m->readAtoms(entries[1]))) {
        cout << "Skipping motif " << entries[0] << " ..." << endl;
        M.erase(entries[0]);
        delete(m);
        continue;
      }
      m->setName(entries[0]);
      m->setPDB(entries[1]);
    }
    vector<string> partAs = MslTools::tokenize(entries[2], ",");
    vector<string> partBs = MslTools::tokenize(entries[3], ",");
    vector<int> partA(partAs.size()), partB(partBs.size());
    for (int j = 0; j < partAs.size(); j++) partA[j] = MslTools::toInt(partAs[j]);
    for (int j = 0; j < partBs.size(); j++) partB[j] = MslTools::toInt(partBs[j]);
    m->addSplit(partA, partB);
    if (!findProductiveMotifSplitAlignments(*m, m->numSplits() - 1, U, resU)) {
      m->removeSplit(m->numSplits() - 1);
    } else {
      nsplits++;
      if (!filteredMotifsFile.empty()) {
        MotifSplit& sp = m->getSplit(m->numSplits() - 1);
        ofs << m->getName() << " " << m->getPDB() << " ";
        ofs << stringifyVector(sp.getSub(0), ",") << " " << stringifyVector(sp.getSub(1), ",") << endl;
      }
    }
    // if the point is to filter only, don't need to keep track of motifs once productive splits are recorded
    if (!filteredMotifsFile.empty()) { M.erase(m->getName()); delete(m); }
  }
  if (!filteredMotifsFile.empty()) { ofs.close(); printf("Filtering done, %d splits survived. Exiting...\n", nsplits); exit(0); }
  lines.clear();

  // remove any motifs with no remaining splits
  for (map<string, Motif*>::iterator it = M.begin(); it != M.end(); ++it) {
    if (it->second->numSplits() > 0) motifs.push_back(it->second);
  }
  printf("%d motifs and %d motif splits survived\n", (int) motifs.size(), nsplits);
}

vector<solution> optimizeHelix(vector<Motif*>& motifs, AtomPointerVector& U, vector<vector<int> >& resU, int n, int Ncyc, int Niters, string opdbf = "") {
  ofstream ofs;
  double rad = computeRadius(U); // get radius of circumscribed sphere

  srand(time(NULL));
  vector<solution> sols(Ncyc);
  for (int c = 0; c < Ncyc; c++) {
    printf("--> Starting iteration %d\n", c+1);
    // initial parameters (base everything on the initial rough number of units per helical turn)
    double Ri = rad*n/M_PI/2;
    double wi = 360/n;
    double Pi = 2*rad + 6;
    double ai = rand() % 360;
    double bi = rand() % 360;
    double gi = rand() % 360;
    // problem parameters
    ProblemDef prob;
    prob.motifs = &motifs;
    prob.unit = &U;
    prob.residues = resU;
    prob.circRadius = rad;
    // cost function
    gsl_multimin_function myfunc;
    myfunc.n = 6;
    myfunc.f = &evaluateOptHelix;
    myfunc.params = (void*) &prob;
    // starting point and initial step size
    gsl_vector *x = gsl_vector_alloc (myfunc.n);
    gsl_vector_set(x, 0, Ri);
    gsl_vector_set(x, 1, wi);
    gsl_vector_set(x, 2, Pi);
    gsl_vector_set(x, 3, ai);
    gsl_vector_set(x, 4, bi);
    gsl_vector_set(x, 5, gi);
    gsl_vector *ss = gsl_vector_alloc (myfunc.n);
    gsl_vector_set(ss, 0, 3);
    gsl_vector_set(ss, 1, 5);
    gsl_vector_set(ss, 2, 3);
    gsl_vector_set(ss, 3, 5);
    gsl_vector_set(ss, 4, 5);
    gsl_vector_set(ss, 5, 5);
    gsl_vector *bestx = gsl_vector_alloc (myfunc.n); gsl_vector_set_all(bestx, 0);
    double bestScore = 0;
    // the loop
    for (int k = 0; k < 1; k++) {
      prob.scoreClashes = true;
      prob.numUnitsToGenerate = -1; // generate all that are within possible contact distance of the original
      // optimizer
      const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
      gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc (T, myfunc.n);
      gsl_multimin_fminimizer_set (minimizer, &myfunc, x, ss);
      size_t iter = 0; int status; double sz;
      do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(minimizer);
        if (status) break;
        sz = gsl_multimin_fminimizer_size (minimizer);
        status = gsl_multimin_test_size (sz, 1e-2);
        if (status == GSL_SUCCESS) {
          printf ("converged to minimum at\n");
        }
        printf ("%5d %10.3e %10.3e %10.3e %10.3e f() = %7.3f size = %.3f\n", (int) iter, gsl_vector_get(minimizer->x, 0), gsl_vector_get(minimizer->x, 1), gsl_vector_get(minimizer->x, 2), gsl_vector_get(minimizer->x, 3), minimizer->fval, sz);
      } while (status == GSL_CONTINUE && iter < Niters);

      // on next iteration start from the best point from this iteration
      gsl_vector_memcpy(x, minimizer->x);

      // update best score (with clashes) if necessary
      prob.scoreClashes = true; prob.numUnitsToGenerate = -1;
      double tempBest = evaluateOptHelix(x, &prob);
      if ((k == 0) || (tempBest < bestScore)) {
        bestScore = tempBest;
        gsl_vector_memcpy(bestx, x);
      }

      // clean old minimizer so can start over from previous point
      gsl_multimin_fminimizer_free(minimizer);
    }
    
    // output best structure
    printf("best score = %e\n", bestScore);
    sols[c].R = gsl_vector_get(bestx, 0);
    sols[c].w = gsl_vector_get(bestx, 1);
    sols[c].P = gsl_vector_get(bestx, 2);
    sols[c].a = gsl_vector_get(bestx, 3);
    sols[c].b = gsl_vector_get(bestx, 4);
    sols[c].g = gsl_vector_get(bestx, 5);
    sols[c].score = bestScore;
    vector<AtomPointerVector> units;
    buildHelicalAssembly(U, units, sols[c].R, sols[c].w, sols[c].P, sols[c].a, sols[c].b, sols[c].g, int(ceil(min(100, 360/sols[c].w) + 1)));
    if (!opdbf.empty()) {
      AtomContainer assembly;
      combineUnits(units, assembly);
      string opdbfc = MslTools::pathRoot(opdbf) + MslTools::intToString(c) + ".pdb";
      assembly.writePdb(opdbfc);
      openFileCPP(ofs, opdbfc, ios::app);
      ofs << "REMARK score = " << sols[c].score << ", R = " << sols[c].R << ", w = " << sols[c].w << ", P = " << sols[c].P << ", alpha = " << sols[c].a << ", beta = " << sols[c].b << ", gamma = " << sols[c].g << endl;
      ofs.close();
      freeUnits(units);
    } else {
      sols[c].units = units;
    }

    // clean up
    gsl_vector_free(x);
    gsl_vector_free(bestx);
    gsl_vector_free(ss);
  }
  return sols;
}

vector<solution> optimizeCircle(vector<Motif*>& motifs, AtomPointerVector& U, vector<vector<int> >& resU, int n, int Ncyc, int Niters, string opdbf = "") {
  ofstream ofs;

  srand(time(NULL));
  vector<solution> sols(Ncyc);
  for (int c = 0; c < Ncyc; c++) {
    printf("--> Starting iteration %d\n", c+1);
    // initial parameters based on rough size of unit
    double ai = rand() % 360;
    double bi = rand() % 360;
    double gi = rand() % 360;
    double rad = computeRadiusXY(U, ai, bi, gi); // get radius of circumscribing circle in XY plane
    double Ri = 2*rad/sqrt(2*(1 - cos(2*M_PI/n))); // you get this by writing the theorem of cosines for the anglular turn between adjacent units (plus some room between adjacent units)
    // problem parameters
    ProblemDef prob;
    prob.motifs = &motifs;
    prob.unit = &U;
    prob.residues = resU;
    prob.circRadius = rad;
    prob.fixedParams["angle"] = 360/n;
    // cost function
    gsl_multimin_function myfunc;
    myfunc.n = 4;
    myfunc.f = &evaluateOptCircle;
    myfunc.params = (void*) &prob;
    // starting point and initial step size
    gsl_vector *x = gsl_vector_alloc (myfunc.n);
    gsl_vector_set(x, 0, Ri);
    gsl_vector_set(x, 1, ai);
    gsl_vector_set(x, 2, bi);
    gsl_vector_set(x, 3, gi);
    gsl_vector *ss = gsl_vector_alloc (myfunc.n);
    gsl_vector_set(ss, 0, 3);
    gsl_vector_set(ss, 1, 10);
    gsl_vector_set(ss, 2, 10);
    gsl_vector_set(ss, 3, 10);
    gsl_vector *bestx = gsl_vector_alloc (myfunc.n); gsl_vector_set_all(bestx, 0);
    double bestScore = 0;
    // the loop
    for (int k = 0; k < 1; k++) {
      prob.scoreClashes = true;
      prob.numUnitsToGenerate = 2; // since working on a circle, will only need to build two adjacent units
      // optimizer
      const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
      gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc (T, myfunc.n);
      gsl_multimin_fminimizer_set (minimizer, &myfunc, x, ss);
      size_t iter = 0; int status; double sz;
      do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(minimizer);
        if (status) break;
        sz = gsl_multimin_fminimizer_size (minimizer);
        status = gsl_multimin_test_size (sz, 1e-2);
        if (status == GSL_SUCCESS) {
          printf ("converged to minimum at\n");
        }
        printf ("%5d %10.3e %10.3e %10.3e %10.3e f() = %7.3f size = %.3f\n", (int) iter, gsl_vector_get(minimizer->x, 0), gsl_vector_get(minimizer->x, 1), gsl_vector_get(minimizer->x, 2), gsl_vector_get(minimizer->x, 3), minimizer->fval, sz);
      } while (status == GSL_CONTINUE && iter < Niters);

      // on next iteration start from the best point from this iteration
      gsl_vector_memcpy(x, minimizer->x);

      // update best score if necessary
      double tempBest = evaluateOptCircle(x, &prob);
      if ((k == 0) || (tempBest < bestScore)) {
        bestScore = tempBest;
        gsl_vector_memcpy(bestx, x);
      }

      // clean old minimizer so can start over from previous point
      gsl_multimin_fminimizer_free(minimizer);
    }
    
    // output best structure
    printf("best score = %e\n", bestScore);
    sols[c].R = gsl_vector_get(bestx, 0);
    sols[c].a = gsl_vector_get(bestx, 1);
    sols[c].b = gsl_vector_get(bestx, 2);
    sols[c].g = gsl_vector_get(bestx, 3);
    sols[c].score = bestScore;
    vector<AtomPointerVector> units;
    buildHelicalAssembly(U, units, sols[c].R, prob.fixedParams["angle"], 0, sols[c].a, sols[c].b, sols[c].g, n);
    if (!opdbf.empty()) {
      AtomContainer assembly;
      combineUnits(units, assembly);
      string opdbfc = MslTools::pathRoot(opdbf) + MslTools::intToString(c) + ".pdb";
      assembly.writePdb(opdbfc);
      openFileCPP(ofs, opdbfc, ios::app);
      ofs << "REMARK score = " << sols[c].score << ", R = " << sols[c].R << ", alpha = " << sols[c].a << ", beta = " << sols[c].b << ", gamma = " << sols[c].g << endl;
      ofs.close();
      freeUnits(units);
    } else {
      sols[c].units = units;
    }

    // clean up
    gsl_vector_free(x);
    gsl_vector_free(bestx);
    gsl_vector_free(ss);
  }
  return sols;
}

int main(int argc, char *argv[]) {
  Transforms tr;

  // parse parameters
	options opt;
  parseCommandLine(argc, argv, opt);

	// read unit, copy backbone, split by residue, and translate to origin
	AtomContainer t;
	AtomPointerVector U;
	vector<vector<int> > resU;
	assert(t.readPdb(opt.pdbFile), "could not read unit 1");
	assert(copyBackbone(t, U, &resU), "could not copy backbone");
  tr.translate(U, -U.getGeometricCenter());
  t.removeAllAtoms();

  // aling unit principal axes to the laboratory frame
  Frame O, F;
  CoordAxes xyz(CartesianPoint(1, 0, 0), CartesianPoint(0, 1, 0), CartesianPoint(0, 0, 1));
  O.computeFrameFromAxes(xyz);
  F.computeFrameFromPCA(U);
  Frame::transformAtoms(U, F, O);
  tr.translate(U, -U.getGeometricCenter());

// // test building and combining
// vector<AtomPointerVector> unitsTest;
// buildHelicalAssembly(U, unitsTest, 50.8971, 30.9748, 57.1866, 29.9424, 6.43647, 29.4324, 12+1);
// AtomContainer assemblyTest;
// combineUnits(unitsTest, assemblyTest);
// assemblyTest.writePdb("ass.pdb");
// exit(0);
  
  // parse the list of motif splits
  vector<Motif*> motifs;
  parseMotifs(motifs, opt.motifsFile, U, resU, opt.filteredMotifsFile);

//  printf("example score = %e\n", evaluateHelix(motifs, U, resU, 50, 30, 50, 15, 15, 15, 12+1, computeRadius(U), true));

  // --- optimize
  vector<solution> sols;
  if (opt.type == "helix") {
    vector<solution> sols = optimizeHelix(motifs, U, resU, opt.numUnits, opt.numCycles, opt.numIters, opt.outPdbFile);
  } else if (opt.type == "circle") {
    vector<solution> sols = optimizeCircle(motifs, U, resU, opt.numUnits, opt.numCycles, opt.numIters, opt.outPdbFile);
  }

  // DONE: 1. write out best parameters as a remark in the output PDB
  // DONE: 2. fix problem--assembly not generated correctly
  // 3. output motifs that matched below cutoff
  // 3.1 perhaps make the iteration for a given motif a separate function
  // 3.2 then use the same function in a function for outputting matching motifs
  // DONE: 4. modularize better, organize classes, make optimization functions for different cases
  
 
}  




