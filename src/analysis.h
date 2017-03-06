/*=========================================================================
 * 
 *  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers
 * 
 *  Copyright (c) 2017 U. Lehnert
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * =========================================================================*/

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "bunch.h"
#include "particle.h"
#include "vector.h"
#include "undulator.h"
#include <tuple>
using namespace std;

// a class for data analysis and emittance calculations

class Analysis
{

  public:

	Analysis(Bunch *bunch);
	
	//evaluate average positions and momentum at every time step
	//write to a file named trajectory.txt
	//write sigma of values to stdx.txt
	void avg();
	void FindStd();
	double avgGamma();
	void FindEmittance();
  private:
	Bunch *B;
};

#endif







