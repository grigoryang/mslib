/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/
#include "CCD.h"
#include "CartesianPoint.h"
#include "CartesianGeometry.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "PDBWriter.h"
#include "AtomVector.h"

CCD::CCD(){
}

CCD::CCD(string _BBQTableForBackboneAtoms){
	bbqT.openReader(_BBQTableForBackboneAtoms);
}

CCD::~CCD(){

}

string CCD::localSample(AtomVector &_av,int numFragments, int maxDegree){

	// Take one atom off each end for BBQ purposes
	Atom *forBBQ_N = new Atom(_av(_av.size()-1));
	Atom *forBBQ_C = new Atom(_av(0));

	_av.erase(_av.begin());
	_av.erase(_av.end()-1);

	AtomVector forBBQ;
	forBBQ.push_back(forBBQ_N);
	forBBQ.push_back(forBBQ_C);



	// Store output PDB string
	stringstream ss;

	RandomNumberGenerator rng;
	rng.setRNGTimeBasedSeed();
	PDBWriter pout;


	_av.saveCoor("pre");
	Atom *fixedAtom = NULL;
	bool reversed = false;
	for (uint d = 0; d < numFragments; d++){
		reversed = false;

		// Figure out which way to order the loop
		if (rng.getRandomDouble() > 0.5){
			std::reverse(_av.begin(), _av.end());
			reversed = true;
		}

		// Fixed atom is last Atom.
		fixedAtom = new Atom(_av(_av.size()-1));

		// Now move the _av
		Transforms t;
		for (uint i = 1; i < _av.size()-2;i++){
		
			CartesianPoint axis = _av(i+1).getCoor()-_av(i).getCoor();
			axis = axis.getUnit();
			axis += _av(i).getCoor();
		        double angle = rng.getRandomIntLimit(maxDegree);
			for (uint j=i+2; j < _av.size();j++){
				t.rotate(_av(j),angle, axis, _av(i).getCoor());
			}

		
		}

		closeFragment(_av,*fixedAtom);

		if (reversed){
			std::reverse(_av.begin(), _av.end());
		}

		Chain c;
		c.addAtoms(forBBQ);
		c.addAtoms(_av);
		
		bbqT.fillInMissingBBAtoms(c);
		
		
		ss << "MODEL "<<endl;
		stringstream tmp;
		pout.open(tmp);
		pout.write(c.getAtoms());
		pout.close();

		ss << tmp.str()<< "ENDMDL\n";


		_av.applySavedCoor("pre");
		delete(fixedAtom);
		fixedAtom = NULL;
	}

	

	return ss.str();
}

void CCD::closeFragment(AtomVector &_av, Atom &_fixedEnd){
	

	int numIterations = 0;
	bool converged = false;
	double convergedDist = MslTools::doubleMax;
	while (true) {
		converged = false;

		for (uint i = _av.size()-3; i > 1;i--){
			//cout << "Pivot index: "<<i<<" "<<_av(i).getResidueNumber()<<endl;

			// Get angle of rotation
			double angle = getMinimumAngle(_av, i,_fixedEnd);
			if (angle == MslTools::doubleMax) continue;
			//cout << "ANGLE: "<<angle*180/M_PI<<endl;

			// Rotate the fragment from pivotIndex+1 to end
			rotateFragment(_av, i, angle);
		
			double dist  = _fixedEnd.distance(_av(_av.size()-1));
			//cout << "\n***Dist: "<<dist<<endl<<endl;;
			if (dist < 0.02){
				convergedDist = dist;
				converged = true;
				break;
			}
		}
		if (converged) break;

		for (uint i = 0; i < _av.size()-2;i++){
			//cout << "Pivot index: "<<i<<" "<<_av(i).getResidueNumber()<<endl;

			// Get angle of rotation
			double angle = getMinimumAngle(_av, i,_fixedEnd);
			if (angle == MslTools::doubleMax) continue;
			//cout << "ANGLE: "<<angle*180/M_PI<<endl;

			// Rotate the fragment from pivotIndex+1 to end
			rotateFragment(_av, i, angle);
		
			double dist  = _fixedEnd.distance(_av(_av.size()-1));
			//cout << "\n***Dist: "<<dist<<endl<<endl;;
			if (dist < 0.02){
				convergedDist = dist;
				converged = true;
				break;
			}
		}

		numIterations++;
		if (converged) break;
		if (numIterations > 1000000) break;
	}

	if (converged) {
		fprintf(stdout, "Fragment closure in %6d\n", numIterations);
	}
}




double CCD::getMinimumAngle(AtomVector &_av, int _indexOfPivot, Atom &_fixedEnd){
	
	// Get Pid and Pih
	CartesianPoint pid  = (_fixedEnd.getCoor()            - _av(_indexOfPivot).getCoor());
	CartesianPoint pih  = (_av(_av.size()-1).getCoor()    - _av(_indexOfPivot).getCoor());
	
	CartesianPoint pivotBondVector =  (_av(_indexOfPivot+1).getCoor() - _av(_indexOfPivot).getCoor());
	pivotBondVector = pivotBondVector.getUnit();

	// Compute Ks
	double k1 = (pid * pivotBondVector);
	k1       *= (pih * pivotBondVector);
	double k2 = (pid * pih);
	CartesianPoint cross = pivotBondVector.cross(pih);
	double k3 = (pid * cross);

//	cout << "Ks: "<<k1<<" "<<k2<<" "<<k3<<endl;
//	cout << "Pid: "<<pid.toString()<<endl;
//	cout << "Pih: "<<pih.toString()<<endl;
//	cout << "z: "<<pivotBondVector.toString()<<endl;

	// Compute 1st derivative at 0.
	double psi = atan(k3 / (k2 - k1));

	// Compute second derivative
	double secondDerviative = (k1-k2)*cos(psi) - k3* sin(psi);
	//cout << "Second derviative: "<<secondDerviative<<endl;
	if (secondDerviative >= 0) {
		cout << "Failed second derivate test"<<endl;
		return MslTools::doubleMax;
	}


	double distToMaxOfPsi   = k1 * (1 - cos(psi))      + k2 * cos(psi)      + k3 * sin(psi);
	double distToMaxOfPsiPi = k1 * (1 - cos(psi+M_PI)) + k2 * cos(psi+M_PI) + k3 * sin(psi+M_PI);
	

	if (distToMaxOfPsi > distToMaxOfPsiPi){
		return psi;
	}

	return psi+M_PI;
	
}


void CCD::rotateFragment(AtomVector &_av, int _indexOfPivot, double _angleOfRotation){


	// Axis of rotation
	CartesianPoint axisOfRotation = (_av(_indexOfPivot+1).getCoor() - _av(_indexOfPivot).getCoor());
	axisOfRotation = axisOfRotation.getUnit();
	axisOfRotation += _av(_indexOfPivot).getCoor();

	// Transform each atom downstream.
	Transforms t;
	for (uint i = _indexOfPivot+2; i < _av.size();i++){
		t.rotate(_av(i),_angleOfRotation, axisOfRotation,_av(_indexOfPivot).getCoor());
	}
	
}