#include <stdio.h>
#include <iostream>
#include <fstream>

#include "Atom.h"
#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "MslTools.h"
#include "release.h"
#include "Transforms.h"
#include "OptimalRMSDCalculator.h"
#include "Frame.h"
#include "glab/ProximitySearch.h"
#include <map>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_multimin.h>
#include <math.h>

using namespace std;
using namespace MSL;

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

void openFileCPP(ofstream& ofs, const char* filename, std::ios_base::openmode mode = ios_base::out) {
	ofs.open(filename, mode);
	if (ofs.fail()) {
		printf("error: failed to open file '%s'\n", filename);
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

void assert(bool cond, string msg) {
	if (!cond) {
		cout << msg << endl;
		exit(-1);
	}
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

void concatenateResidueAtoms(AtomPointerVector& unit, int beg, int len, vector<vector<int> >& residues, AtomPointerVector& sub) {
  for (int i = beg; i < beg + len; i++) {
    vector<int>& res = residues[i];
    for (int j = 0; j < res.size(); j++) {
      sub.push_back(unit[res[j]]);
    }
  }
}

void concatenateResidueAtoms(AtomPointerVector& unit, vector<int>& indices, vector<vector<int> >& residues, AtomPointerVector& sub) {
  for (int i = 0; i < indices.size(); i++) {
    vector<int>& res = residues[indices[i]];
    for (int j = 0; j < res.size(); j++) {
      sub.push_back(unit[res[j]]);
    }
  }
}

void concatenateResidueAtoms(vector<AtomPointerVector>& residues, AtomPointerVector& sub) {
  for (int i = 0; i < residues.size(); i++) {
    for (int j = 0; j < residues[i].size(); j++) {
      sub.push_back(residues[i][j]);
    }
  }
}

void copyResidues(AtomPointerVector& unit, int beg, int len, vector<vector<int> >& residues, vector<AtomPointerVector>& into, vector<int>& intoIndices) {
  for (int i = 0; i < intoIndices.size(); i++) {
    into[intoIndices[i]].clear();
  }
  for (int i = 0; i < len; i++) {
    vector<int>& res = residues[i + beg];
    for (int j = 0; j < res.size(); j++) {
      into[intoIndices[i]].push_back(unit[res[j]]);
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
		int numRes;
		double rCut;
		string name, pdbfile;

	public:
		Motif() { rCut = -1; }
    ~Motif();
		void addSplit(vector<int>& partA, vector<int>& partB);
		bool assignRMSD();
		bool readAtoms(string pdbFile);
		void setName(string _name) { name = _name; }
		void setPDB(string _pdbfile) { pdbfile = _pdbfile; }
		int numSplits() { return splits.size(); }
		int length() { return atoms.size(); }
		int numResidues() { return numRes; }
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
	return copyBackbone(t, atoms, &residues);
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
	numRes = residues.size();
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
		bool scoreScrew;
    double circRadius;
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
      AtomPointerVector& motif = m.getAtoms();
      vector<int>& algnA = sp.getAlgnPositions(0);
      vector<int>& algnB = sp.getAlgnPositions(1);
      vector<int>& subIndsA = sp.getSub(0);
      vector<int>& subIndsB = sp.getSub(1);
      int lenA = sp.length(0);
      int lenB = sp.length(1);
      // each motif can fit by aligning sub A onto unit 1 and sub B onto unit 2, or sub A onto unit 2 and sub B onto unit 1
      vector<AtomPointerVector> residuesAB(lenA + lenB);
      for (int order = 0; order < 2; order++) {
        AtomPointerVector& unitA = (order ? unit1 : unit2);
        AtomPointerVector& unitB = (order ? unit2 : unit1);
        // iterate over options of the first sub
        for (int a = 0; a < algnA.size(); a++) {
          AtomPointerVector regAB;
          copyResidues(unitA, algnA[a], lenA, residues, residuesAB, subIndsA);
          // iterate over options of the second sub
          for (int b = 0; b < algnB.size(); b++) {
            copyResidues(unitB, algnB[b], lenB, residues, residuesAB, subIndsB);
            AtomPointerVector regAB;
            concatenateResidueAtoms(residuesAB, regAB);
            // consider this combination of alignments
            double rmsd = orc.bestRMSD(regAB, motif, &suc);
            assert(suc, "RMSD calculation failed!");
//            if (rmsd < m.getRMSDCut()) {
              // calculate the contribution of this match
//              score -= m.numResidues() / (rmsd + 0.1) - m.numResidues() / (m.getRMSDCut() + 0.1);
//              score -= m.numResidues() * exp(-max(rmsd - m.getRMSDCut(), 0));
//              score -= m.numResidues() * exp(-pow(max(rmsd - m.getRMSDCut(), 0), 2));
              score -= m.numResidues() / (1 + exp(10*(rmsd - m.getRMSDCut())));
//            }
          }
        }
      }
    }
  }
  
  if (scoreClashes) {
    ProximitySearch ps(unit1, 6.0);
    for (int i = 0; i < unit2.size(); i++) {
      Atom& a2 = *(unit2[i]);
      vector<int> closeOnes;
      ps.pointsWithin(a2.getCoor(), 0, 6.0, closeOnes);
      for (int j = 0; j < closeOnes.size(); j++) {
        Atom& a1 = *(unit1[closeOnes[j]]);
        if (a1.distance(a2) < 0.5*(a1.getRadius() + a2.getRadius())) {
          score += 50;
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
// If rad is not specified, will build the first N units, otherwise will only preserve the first unit and all whose circumscribed spheres are within potential contacting distance.
void buildHelicalAssembly(AtomPointerVector& U, vector<AtomPointerVector>& units, double R, double w, double P, double alpha, double beta, double gamma, int N, double rad = -1) {
  double dZ = (P*w/360); // rise per unit
  double wr = w*M_PI/180;
  double dcut = (2*rad + 6.5);
  double dcut2 = dcut*dcut;
  Transforms tr;
  Frame Fo, Fi;
  // We shall define a frame associated with each unit in the helix. The X-axis in this frame
  // will be the vector from the helix axis to the unit center of mass, in the XY laboratory plane.
  // The Y-axis in this frame will be the tangent to the helical curve in the current point; note,
  // this is by definition orthogonal to the X-axis. The Z-axis will then be the cross product of X and Y.
  // X = [R*cos(wr*t) R*sin(wr*t), 0] (normalized)
  // Y = [-R*wr*sin(wr*t), R*wr*cos(wr*t), dZ]
  // Z = cross(X, Y)
  double b1 = 1/sqrt(1 + pow(dZ/R*wr, 2));
  double b2 = 1/sqrt(1 + pow(R*wr/dZ, 2));
  CartesianPoint Xo(1, 0, 0);
  CartesianPoint Yo(0, b1, b2);
  CartesianPoint Zo = Xo.cross(Yo);
  CoordAxes axo(Xo, Yo, Zo);
  Fo.computeFrameFromAxes(axo);
  
  // rotate the original unit by a, b, g
  AtomPointerVector unit;
  cloneAtoms(U, unit);
  tr.translate(unit, -(unit.getGeometricCenter()));
  tr.Xrotate(unit, alpha);
  tr.Yrotate(unit, beta);
  tr.Zrotate(unit, gamma);

  // now generate clones around the helical curve
  for (int t = 0; t < ((rad < 0) ? N : min(100, ceil(abs(dcut/dZ)))); t++) {
    double xo = R*cos(wr*t); double yo = R*sin(wr*t); double zo = dZ*t;
    if ((rad >= 0) && (t > 0)) {
      if (2*pow(R, 2)*(1 - cos(wr*t)) + pow(dZ*t, 2) > dcut2) continue;
    }
    AtomPointerVector newUnit;
    cloneAtoms(unit, newUnit);

    // rotate and translate into helical frame
    CartesianPoint Xi(cos(wr*t), sin(wr*t), 0);
    CartesianPoint Yi(-b1*sin(wr*t), b1*cos(wr*t), b2);
    CartesianPoint Zi = Xi.cross(Yi);
    CoordAxes axi(Xi, Yi, Zi);
    Fi.computeFrameFromAxes(axi);
    Frame::transformAtoms(newUnit, Fo, Fi);
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
      sprintf(segid, "%s%03d", unit[i]->getChainId().c_str(), i);
      Atom na(*(unit[j]));
      na.setChainId(segid);
      na.setSegID(segid);
      con.addAtom(na);
    }
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

printf("to evaluate [%f, %f, %f, %f, %f, %f] built %d units, SCORE = %f\n", R, w, P, alpha, beta, gamma, (int) units.size(), score);
// if (score < 1) {
//   printf("SCORE is %f\n", score);
//   AtomContainer assembly;
//   combineUnits(units, assembly);
//   assembly.writePdb("ass.pdb");
//   exit(-1);
// }

  // free memory
  for (int i = 0; i < units.size(); i++) {
    for (int j = 0; j < units[i].size(); j++) { delete(units[i][j]); }
  }

  return score;
}

double evaluateOpt(const gsl_vector* v, void* params) {
	Transforms tr;
	ProblemDef* p = (ProblemDef*) params;
	double R = gsl_vector_get(v, 0);
	double w = gsl_vector_get(v, 1);
	double P = gsl_vector_get(v, 2);
  double alpha = gsl_vector_get(v, 3);
  double beta = gsl_vector_get(v, 4);
  double gamma = gsl_vector_get(v, 5);

	if (p->scoreScrew) {
    return evaluateHelix(*(p->motifs), *(p->unit), p->residues, R, w, P, alpha, beta, gamma, -1, p->circRadius, p->scoreClashes);
	} else {
    return evaluateHelix(*(p->motifs), *(p->unit), p->residues, R, w, P, alpha, beta, gamma, 2, p->circRadius, p->scoreClashes);
	}

}

double computeRadius(AtomPointerVector& A) {
  CartesianPoint cent = A.getGeometricCenter();
  double R = 0;
  for (int i = 0; i < A.size(); i++) {
    R = max(R, cent.distance((A[i])->getCoor()));
  }
  return R;
}

int main(int argc, char *argv[]) {
	OptimalRMSDCalculator orc;
	Transforms tr;
	ofstream ofs;
	RandomNumberGenerator rng;
	rng.setTimeBasedSeed();
	bool suc = true;

  // parse parameters
	if (argc < 5) {
		printf("Usage: [unit.pdb] [motifList.list] [N iterations] [pdbOut] {filtered_motif_list.out}\n");
    printf("list example: /u/gevorg/work/proteinUniverse/assemblies/optimizerInput1.list\n\n");
		exit(-1);
	}
  string ipdbf(argv[1]);
  string motifsFile(argv[2]);
  int Niters = MslTools::toInt(argv[3]); assert(Niters > 0, "number of iterations has to be positive!");
  string opdbf(argv[4]);
  string filteredMotifsFile;
  if (argc > 5) filteredMotifsFile = argv[5];

	// read unit, copy backbone, split by residue, and translate to origin
	AtomContainer t;
	AtomPointerVector U;
	vector<vector<int> > resU;
	assert(t.readPdb(ipdbf), "could not read unit 1");
	assert(copyBackbone(t, U, &resU), "could not copy backbone");
  tr.translate(U, -U.getGeometricCenter());
  setChainId(U, "A");
  t.removeAllAtoms();

  // aling unit principal axes to the laboratory frame
  Frame O, F;
  CoordAxes xyz(CartesianPoint(1, 0, 0), CartesianPoint(0, 1, 0), CartesianPoint(0, 0, 1));
  O.computeFrameFromAxes(xyz);
  F.computeFrameFromPCA(U);
  Frame::transformAtoms(U, F, O);
  tr.translate(U, -U.getGeometricCenter());
  double rad = computeRadius(U); // get radius of circumscribed sphere

// // test building and combining
// vector<AtomPointerVector> unitsTest;
// buildHelicalAssembly(U, unitsTest, 50, 30, 50, 15, 15, 15, 12+1);
// AtomContainer assemblyTest;
// combineUnits(unitsTest, assemblyTest);
// assemblyTest.writePdb("ass.pdb");
// exit(0);
  
	// parse the list of motif splits
	map<string, Motif> M;
	vector<string> names;
	vector<string> lines;
	MslTools::readTextFile(lines, motifsFile);
	for (int i = 0; i < lines.size(); i++) {
		if (i % 100 == 0) printf("reading %d/%d...\n", i+1, (int) lines.size());
		vector<string> entries = MslTools::tokenize(lines[i], " ");
		assert(entries.size() == 4, "could not parse line '" + lines[0] + "' from motif list file");
		bool prev = (M.find(entries[0]) != M.end());
		Motif& m = M[entries[0]];
		if (!prev) {
			if (!(m.readAtoms(entries[1]))) {
				cout << "Skipping motif " << entries[0] << " ..." << endl;
				M.erase(entries[0]);
				continue;
			}
			m.setName(entries[0]);
			m.setPDB(entries[1]);
			names.push_back(m.getName());
		}
		vector<string> partAs = MslTools::tokenize(entries[2], ",");
		vector<string> partBs = MslTools::tokenize(entries[3], ",");
		vector<int> partA(partAs.size()), partB(partBs.size());
		for (int j = 0; j < partAs.size(); j++) partA[j] = MslTools::toInt(partAs[j]);
		for (int j = 0; j < partBs.size(); j++) partB[j] = MslTools::toInt(partBs[j]);
		m.addSplit(partA, partB);
	}
	lines.clear();

	// set RMSD cutoffs for each motif and for each sub of each motif.
	// note that for subs it matters if they are contiguous or have
	// multiple disjoint units, so look for breaks
	for (int i = 0; i < names.size(); i++) {
		if (!M[names[i]].assignRMSD()) {
			M.erase(names[i]);
		}
	}
	names.clear();
	for (map<string, Motif>::iterator it = M.begin(); it != M.end(); ++it) {
		names.push_back(it->first);
	}
	printf("%d motifs survived RMSD assignment\n", (int) names.size());

	// try out all motif subs in all positions of unit, to determine all the
	// plausible alignment positions
	for (int i = 0; i < names.size(); i++) {
		Motif& m = M[names[i]];
		for (int j = 0; j < m.numSplits(); j++) {
			MotifSplit& sp = m.getSplit(j);
			bool matches = true;
			for (int k = 0; k < 2; k++) {
				AtomPointerVector sub = m.getSub(j, k);
				for (int l = 0; l < (int) resU.size() - (int) sp.length(k) + 1; l++) {
					AtomPointerVector reg;
          concatenateResidueAtoms(U, l, sp.length(k), resU, reg);
					double rmsd = orc.bestRMSD(sub, reg, &suc);
					assert(suc, "RMSD calculation failed!");
					if (rmsd < sp.getRMSDCut(k)) {
						sp.addAlgnPos(k, l, rmsd);
					}
				}
				if (sp.numAlgnPos(k) == 0) { matches = false; break; }
			}
			if (!matches) {
				M.erase(names[i]);
				break;
			}
		}
	}
	names.clear();
	printf("%d motif splits survived after filtering for individual alignments...\n", (int) M.size());
	vector<Motif*> motifs;
	if (!filteredMotifsFile.empty()) openFileCPP(ofs, filteredMotifsFile.c_str());
	for (map<string, Motif>::iterator it = M.begin(); it != M.end(); ++it) {
		motifs.push_back(&(it->second));
		if (!filteredMotifsFile.empty()) {
			Motif& m = it->second;
			for (int j = 0; j < m.numSplits(); j++) {
				MotifSplit& sp = m.getSplit(j);
				ofs << m.getName() << " " << m.getPDB() << " ";
				ofs << stringifyVector(sp.getSub(0), ",") << " " << stringifyVector(sp.getSub(1), ",") << endl;
			}
		}
	}
	if (!filteredMotifsFile.empty()) ofs.close();

  // one evaluation
  printf("initial score = %f\n", evaluateHelix(motifs, U, resU, 50, 30, 50, 15, 15, 15, 12+1, rad, true));

	// --- optimization loop
	// problem parameters
	ProblemDef prob;
  prob.motifs = &motifs;
  prob.unit = &U;
  prob.residues = resU;
  prob.circRadius = rad;
	// cost function
	gsl_multimin_function myfunc;
	myfunc.n = 6;
	myfunc.f = &evaluateOpt;
	myfunc.params = (void*) &prob;
	// starting point and initial step size
	gsl_vector *x = gsl_vector_alloc (myfunc.n);
  gsl_vector_set(x, 0, 30);
  gsl_vector_set(x, 1, 20);
  gsl_vector_set(x, 2, 50);
  gsl_vector_set(x, 3, 0);
  gsl_vector_set(x, 4, 0);
  gsl_vector_set(x, 5, 0);
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
	for (int k = 0; k < 3*4; k++) {
// 		prob.scoreClashes = k % 3;
// 		prob.scoreScrew = (k % 3) == 2;
    prob.scoreClashes = true;
    prob.scoreScrew = true;
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
		prob.scoreClashes = true; prob.scoreScrew = true;
		double tempBest = evaluateOpt(x, &prob);
		if ((k == 0) || (tempBest < bestScore)) {
			bestScore = tempBest;
			gsl_vector_memcpy(bestx, x);
		}

		// clean old minimizer so can start over from previous point
		gsl_multimin_fminimizer_free(minimizer);
	}
	
	// output best solution
	prob.scoreClashes = true;
	prob.scoreScrew = true;
	evaluateOpt(bestx, &prob);
	printf("best score = %f\n", bestScore);

  // output best structure
  double R = gsl_vector_get(bestx, 0);
  double w = gsl_vector_get(bestx, 1);
  double P = gsl_vector_get(bestx, 2);
  double a = gsl_vector_get(bestx, 3);
  double b = gsl_vector_get(bestx, 4);
  double g = gsl_vector_get(bestx, 5);
  vector<AtomPointerVector> units;
  buildHelicalAssembly(U, units, R, w, P, a, b, g, int(min(100, 360/w) + 1));
  AtomContainer assembly;
  combineUnits(units, assembly);
  assembly.writePdb(opdbf);
  
	// clean up
	gsl_vector_free(x);
	gsl_vector_free(bestx);
	gsl_vector_free(ss);
 
}  




