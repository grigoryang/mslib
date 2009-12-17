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

#ifndef FRAME_H
#define FRAME_H


// MSL Includes
#include "Line.h"
#include "AtomVector.h"
#include "PrincipleComponentAnalysis.h"
#include "CoordAxes.h"

// STL Includes
#include <iostream>
#include <fstream>
#include <map>

// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#endif

// Namespaces
using namespace std;

class Frame {
public:
    // Constructors
    Frame();
    Frame(const Frame &_frame);
    ~Frame();

    // Operators
    void operator=(const Frame & _frame);
    friend ostream & operator<<(ostream &_os, Frame * _frame);
    friend ofstream & operator<<(ofstream &_of, Frame * _frame);
    Line & operator[](string _n);

    // Member Functions
    map<string, Line> getLines() const;
    PrincipleComponentAnalysis getPCA() const;

    // Create a Frame somehow..
    void computeFrameFromPCA(AtomVector &_atoms);
    void computeFrameFrom3Atoms(Atom &_at1, Atom &_at2, Atom &_at3);
    void computeFrameFromAxes(CoordAxes &_axes);
    bool computeFrameFromFunctionalGroup(Residue &_res); // return false if no frame computed
    bool computeFrameFrom3AtomNames(Residue &_res, string & atom1, string & atom2, string & atom3); // return false if no frame computed
    void computeFrameFrom2Lines(Line &_Z, Line &_X);

    // Transformation Matrix between frames
    Matrix getBasisTransformMatrix(Frame &_frame);
    static Matrix getBasisTransformMatrix(Frame &_fromFrame, Frame &_toFrame);

    // Transform atoms using local to global frame transformation matrix
    void transformToFromGlobalBasis(AtomVector &_atoms, bool bToGlobal);
    void transformToGlobalBasis(AtomVector &_atoms);
    void transformFromGlobalBasis(AtomVector &_atoms);
    static void transformAtoms(AtomVector &_atoms, Frame &_fromFrame, Frame &_toFrame);


    void setName(string _name);
    string getName() const;

    string toString();


    double distanceToFrame(Frame &_frame);
    Matrix anglesBetweenFrame(Frame &_frame);

    CartesianPoint getCenter() const;

    // Geometric functions
    void translate(const CartesianPoint & vec);

private:

    // Copy used in assignment " = " operator
    void copy(const Frame & _frame);


    string name;
    map<string, Line> lines;
    PrincipleComponentAnalysis pca;
    CartesianPoint center;



    // BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__
public:

    void save_checkpoint(string filename) const {
        std::ofstream fout(filename.c_str());
        boost::archive::text_oarchive oa(fout);
        oa << (*this);
    }

    void load_checkpoint(string filename) {
        std::ifstream fin(filename.c_str(), std::ios::binary);
        boost::archive::text_iarchive ia(fin);
        ia >> (*this);
    }
private:

    friend class boost::serialization::access;

    template<class Archive> void serialize(Archive & ar, const unsigned int version) {
        ar & name;
        ar & lines;
        ar & pca;
        ar & center;
    }
#endif

};
#endif

// INLINES

inline ostream & operator<<(ostream &_os, Frame * _frame) {
    _os << _frame->toString();
    return _os;
}

inline ofstream & operator<<(ofstream &_of, Frame * _frame) {
    _of << _frame->toString();
    return _of;
}

inline CartesianPoint Frame::getCenter() const {
    return center;
}

inline string Frame::getName() const {
    return name;
}
