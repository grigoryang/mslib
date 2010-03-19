/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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
#include <iostream>
#include <cstdlib>

#include "System.h"
#include "MslTools.h"
#include "Transforms.h"


using namespace std;
using namespace MSL;

/*******************************************************************
 *  This program illustrates how to use the AtomPointerVector
 *  object, which is the object used for communication by a number
 *  of objects that manipulate atoms.
 *  
 *  The AtomPointerVector derives from inheritance from the stl
 *  vector<Atom*>, maintaining its basic functions and adding a few
 *  more.
 *******************************************************************/

int main(int argc, char *argv[]) {

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to use the AtomPointerVector (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;

	// Let's create a AtomPointerVector from scratch.  The atom contructor can take
	// an atomId (chain,resnum,res name, atom name) and the coordinates
	AtomPointerVector atoms;
	atoms.push_back(new Atom("A,1,ALA,N",    2.143,  1.328,  0.000));
	atoms.push_back(new Atom("A,1,ALA,CA",   1.539,  0.000,  0.000));
	atoms.push_back(new Atom("A,1,ALA,CB",   2.095, -0.791,  1.207));
	atoms.push_back(new Atom("A,1,ALA,C",    0.000,  0.000,  0.000));
	atoms.push_back(new Atom("A,1,ALA,O",   -0.661,  1.037,  0.000));
	atoms.push_back(new Atom("A,2,ILE,N",   -0.612, -1.210, -0.000));
	atoms.push_back(new Atom("A,2,ILE,CA",  -2.052, -1.462, -0.000));
	atoms.push_back(new Atom("A,2,ILE,CB",  -2.790, -0.764, -1.194));
	atoms.push_back(new Atom("A,2,ILE,CG2", -2.141, -1.025, -2.572));
	atoms.push_back(new Atom("A,2,ILE,CG1", -4.319, -1.018, -1.213));
	atoms.push_back(new Atom("A,2,ILE,CD1", -5.054, -0.338, -2.380));
	atoms.push_back(new Atom("A,2,ILE,C",   -2.221, -2.971, -0.000));
	atoms.push_back(new Atom("A,2,ILE,O",   -1.239, -3.712, -0.000));
	atoms.push_back(new Atom("A,3,ALA,N",   -3.474, -3.466, -0.000));
	atoms.push_back(new Atom("A,3,ALA,CA",  -3.791, -4.877, -0.000));
	atoms.push_back(new Atom("A,3,ALA,CB",  -3.085, -5.538,  1.207));
	atoms.push_back(new Atom("A,3,ALA,C",   -5.297, -5.192, -0.000));
	atoms.push_back(new Atom("A,3,ALA,OT1", -6.104, -4.223,  0.000));
	atoms.push_back(new Atom("A,3,ALA,OT2", -5.649, -6.401, -0.000));
	
	// one of the added functions to a regular vector<Atom*> is a << operator,
	// which print all atoms
	cout << "Print the AtomPointerVector using the << operator" << endl;
	cout << atoms << endl;
	cout << endl;
	cout << "===============================" << endl;

	cout << "Print all atoms by looping over all elements and using the << operator on the Atom" << endl;
	for (unsigned int i=0; i<atoms.size(); i++) {
		cout << *(atoms[i]) << endl; // the [] returns the pointer
		//cout << atoms(i) << endl; // alternative way, the () operator returns and Atom object (by reference)
	}
	cout << endl;
	cout << "===============================" << endl;

	/*************************************************************************
	 *  The AtomPointerVector is mainly used to "communicate" between objects
	 *  that deal with atoms.
	 *
	 *  For example, an AtomPointerVector can be fed to create a System
	 *************************************************************************/
	System sys(atoms);
	cout << "Print the System created from the AtomPointerVector" << endl;
	cout << sys << endl;

	/*************************************************************************
	 *  Note, the System makes internal copies of the atoms, their memory address is
	 *  different from that of the original pointers in the AtomPointerVector
	 *************************************************************************/
	// get atoms[0] atomId and get the corrispoing atom from the System
	string atomId = atoms[0]->getAtomId(); 
	Atom * pAtom = &(sys.getAtom(atomId)); 
	// print both atoms and their address, show they are different
	cout << *pAtom << " (" << pAtom << ")" << endl;
	cout << *(atoms[0]) << " (" << atoms[0] << ")" << endl;

	cout << endl;
	cout << "===============================" << endl;

	/*************************************************************************
	 *  Here is an example on how objects communicate using the AtomPointerVector.
	 *
	 *  The atom pointers from the System can get obtained using the .getAtomPointers
	 *  function.
	 * 
	 *  The Transforms objects takes an AtomPointerVector and can apply a translation.
	 *
	 *  Here we translate all atoms in the System by (3.0, 0.0, 0.0).
	 *************************************************************************/
	Transforms tr;

	tr.translate(sys.getAtomPointers(), CartesianPoint(3.0, 0.0, 0.0));
	/*************************************************************************
	 *  Note: we could have done the same more explicitely this way
	 *  	AtomPointerVector sysAtms = sys.getAtomPointers();
	 *      tr.translate(sysAtms, CartesianPoint(3.0, 0.0, 0.0);
	 *************************************************************************/

	cout << "Print the System atoms after the trasnlation by (3, 0, 0)" << endl;
	cout << sys.getAtomPointers() << endl;

	// since we manually created new Atom pointers, before closing we need to delete the allocated memory
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete (*k);
		*k = NULL;
	}
	atoms.clear();

	return 0;
}
