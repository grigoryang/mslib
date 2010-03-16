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


#ifndef CHARMMSYSTEMBUILDER_H
#define CHARMMSYSTEMBUILDER_H

#include <iostream>
#include <vector>

#include "System.h"
#include "CharmmTopologyResidue.h"
#include "CharmmTopologyReader.h"
#include "CharmmParameterReader.h"
#include "CharmmEEF1ParameterReader.h"
#include "PolymerSequence.h"
#include "CharmmVdwInteraction.h"
#include "CharmmBondInteraction.h"
#include "CharmmElectrostaticInteraction.h"
#include "CharmmUreyBradleyInteraction.h"
#include "CharmmAngleInteraction.h"
#include "CharmmDihedralInteraction.h"
#include "CharmmImproperInteraction.h"
#include "CharmmEEF1Interaction.h"
#include "CharmmEEF1RefInteraction.h"


namespace MSL { 
class CharmmSystemBuilder {
	public:
		CharmmSystemBuilder();
		CharmmSystemBuilder(std::string _topologyFile, std::string _parameterFile, std::string _solvationFile="");
		CharmmSystemBuilder(const CharmmSystemBuilder & _sysBuild);
		~CharmmSystemBuilder();

		void operator=(const CharmmSystemBuilder & _sysBuild);

		bool readTopology(std::string _topologyFile);
		bool readParameters(std::string _parameterFile);
		bool readSolvation(std::string _solvationFile);
		void setSolvent(std::string _solvent);

		bool buildSystem(System & _system, const PolymerSequence & _sequence);
		bool updateNonBonded(System & _system, double _ctonnb=0.0, double _ctofnb=0.0, double _cutnb=0.0);
	//	bool updateSolvation(System & _system, std::string _solvent, double _ctonnb, double _ctofnb, double _cutnb);
		
		bool getBuildNonBondedInteractions();
		void setBuildNonBondedInteractions(bool _flag);

		CharmmTopologyReader * getCharmmTopologyReader();
		CharmmParameterReader * getCharmmParameterReader();
		CharmmEEF1ParameterReader * getCharmmEEF1ParameterReader();

		void setVdwRescalingFactor(double _factor);
		double setVdwRescalingFactor() const;

		// rescaling of the 1-4 electostatic interactions.  should be 1 for charmm 22
		// and 0.6 for charmm 19
		void setElec14factor(double _e14);
		double getElec14factor() const;

		// the dielectric constant
		void setDielectricConstant(double _diel);
		double getDielectricConstant() const;

		// use a distance dependent dielectric
		void setUseRdielectric(bool _flag);
		bool getUseRdielectric() const;

	private:
		void setup();
		void copy(const CharmmSystemBuilder & _sysBuild);
		void deletePointers();
		void getAtomPointersFromMulti(std::string _name, std::vector<Atom*> & _out, std::vector<CharmmTopologyResidue*> & _position, std::vector<std::map<std::string, Atom*> > & _atomMap);
		std::vector<Atom*> getAtomPointers(std::string _name, std::vector<std::vector<std::vector<CharmmTopologyResidue*> > >::iterator & _chItr, std::vector<std::vector<CharmmTopologyResidue*> >::iterator & _posItr, std::vector<CharmmTopologyResidue*>::iterator & _idItr);

		std::vector<std::vector<std::vector<CharmmTopologyResidue*> > > polymerDefi;
		std::vector<std::vector<std::vector<std::map<std::string, Atom*> > > > atomMap;

		CharmmTopologyReader * pTopReader;
		CharmmParameterReader * pParReader;
		CharmmEEF1ParameterReader * pEEF1ParReader;

		bool buildNonBondedInteractions;

		double vdwRescalingFactor;
		
		double elec14factor;
		double dielectricConstant;
		bool useRdielectric;
		bool useSolvation;
		std::string solvent;

};

inline bool CharmmSystemBuilder::readTopology(std::string _topologyFile) {
	pTopReader->reset();
	if (!pTopReader->open(_topologyFile)) {
		return false;
	}
	bool out = false;
	if (pTopReader->read()) {
		out = true;
	}
	pTopReader->close();
	return out;
}
inline bool CharmmSystemBuilder::readParameters(std::string _parameterFile) {
	pParReader->reset();
	if (!pParReader->open(_parameterFile)) {
		return false;
	} 
	bool out = false;
	if (pParReader->read()) {
		out = true;
	}
	pParReader->close();
	return out;
}
inline bool CharmmSystemBuilder::readSolvation(std::string _solvationFile) {
	useSolvation = true;
	pEEF1ParReader->reset();
	if (!pEEF1ParReader->open(_solvationFile)) {
		return false;
	} 
	bool out = false;
	if (pEEF1ParReader->read()) {
		out = true;
	}
	pEEF1ParReader->close();
	useRdielectric = true;

	return out;
}

inline bool CharmmSystemBuilder::getBuildNonBondedInteractions()  { return buildNonBondedInteractions; }
inline void CharmmSystemBuilder::setBuildNonBondedInteractions(bool _flag) { buildNonBondedInteractions = _flag; }

inline CharmmTopologyReader  * CharmmSystemBuilder::getCharmmTopologyReader(){return pTopReader; }
inline CharmmParameterReader * CharmmSystemBuilder::getCharmmParameterReader(){return pParReader; }
inline CharmmEEF1ParameterReader * CharmmSystemBuilder::getCharmmEEF1ParameterReader(){return pEEF1ParReader; }
inline void CharmmSystemBuilder::setVdwRescalingFactor(double _factor) {vdwRescalingFactor = _factor;}
inline double CharmmSystemBuilder::setVdwRescalingFactor() const {return vdwRescalingFactor;}
inline void CharmmSystemBuilder::setElec14factor(double _e14) {elec14factor = _e14;}
inline double CharmmSystemBuilder::getElec14factor() const {return elec14factor;}
inline void CharmmSystemBuilder::setDielectricConstant(double _diel) {dielectricConstant = _diel;}
inline double CharmmSystemBuilder::getDielectricConstant() const {return dielectricConstant;}
inline void CharmmSystemBuilder::setUseRdielectric(bool _flag) {useRdielectric = _flag;}
inline bool CharmmSystemBuilder::getUseRdielectric() const {return useRdielectric;}

}

#endif
