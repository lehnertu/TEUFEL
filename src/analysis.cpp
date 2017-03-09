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

#include "bunch.h"
#include "particle.h"
#include "analysis.h"
#include "global.h"
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <tuple>
#include <externalfield.h>
#include "omp.h"
#include <iostream>
using namespace std;

Analysis::Analysis(Bunch *bunch)
{
  //create a copy of bunch
  B=bunch;
  
}

void Analysis::avg()
{
	std::ofstream Out2("avg.txt", std::ofstream::out);
  	std::ofstream Out3("stdx.txt", std::ofstream::out);
	int NOTS = B->getNOTS();
	int NOP = B->getNOP();
	for (int i=0;i<NOTS;i++)
  	{
	  double sumx=0.0;
	  double sumy=0.0;
	  double sumz=0.0;
	  double sumpx=0.0;
	  double sumpy=0.0;
	  double sumpz = 0.0;
	  double sum2x=0.0;
	  double sum2y=0.0;
	  double sum2z=0.0;
	  double sum2px=0.0;
	  double sum2py=0.0;
	  double sum2pz=0.0;
	  double stdx = 0.0;
	  double stdy = 0.0;
	  double stdz = 0.0;
	  double stdpx=0.0;
	  
#pragma omp parallel for reduction(+:sumx,sumy,sumz,sumpx,sumpy,sumpz) shared(i)
	  for (int k=0; k<NOP;k++)
		{
			Vector XP =(B->b[k].TrajPoint(i));
			Vector XP1=(B->b[k].TrajMomentum(i));
			sumx += XP.x;
			sumy += XP.y;
			sumz += XP.z;
			sumpx += XP1.x;
			sumpy += XP1.y;
			sumpz +=XP1.z;
		}
		
	  double NP =(double)NOP;
	  double avgx = sumx/NP;
	  double avgy = sumy/NP;
	  double avgz = sumz/NP;
	  double avgpx = sumpx/NP;
	  double avgpy = sumpy/NP;
	  double avgpz = sumpz/NP;
	  Out2<<avgz<<"\t"<<avgx<<"\t"<<avgy<<"\t"<<avgz<<"\t"<<avgpx<<"\t"<<avgpy<<"\t"<<avgpz<<"\n";
#pragma omp parallel for reduction(+:sum2x,sum2y,sum2z,sum2px,sum2py,sum2pz) 
	  for (int k=0; k<NOP;k++)
		{
			Vector XP =(B->b[k].TrajPoint(i));
			Vector XP1=(B->b[k].TrajMomentum(i));
			sum2x += pow((XP.x-avgx),2);
			sum2y += pow((XP.y-avgy),2);
			sum2z += pow((XP.z-avgz),2);
			sum2px += pow(XP1.x-avgpx,2);
			sum2py += pow(XP1.y-avgpy,2);
			sum2pz += pow(XP1.z-avgpz,2);
		}
	  stdx = sqrt(sum2x/NP);
	  stdy = sqrt(sum2y/NP);
	  stdz = sqrt(sum2z/NP);
	  stdpx = sqrt(sum2px/NP);
	  double stdpy = sqrt(sum2py/NP);
	  double stdpz = sqrt(sum2pz/NP);
	  Out3<<avgz<<"\t"<<stdx<<"\t"<<stdy<<"\t"<<stdz<<"\t"<<stdpx<<"\t"<<stdpy<<"\t"<<stdpz<<"\n";
 	 }

Out2.close();
Out3.close();
}


//dump all the trajectory to a file named "trajectory.txt"
void Analysis::DumpTrajectory()
{

	std::ofstream Out1("trajectory.txt", std::ofstream::out);
	int NOTS = B->getNOTS();
	int NOP = B->getNOP();
	for (int k = 0;k<NOTS;k++)
	{
	  for (int i = 0;i<NOP;i++)
		{
			Vector XP = B->b[i].TrajPoint(k);
			Vector XP1 = B->b[i].TrajMomentum(k);
			double t = B->b[i].TrajTime(k);
			Out1<<t<<"\t"<<XP.x<<"\t"<<XP.y<<"\t"<<XP.z<<"\t"<<XP1.x<<"\t"<<XP1.y<<"\t"<<XP1.z<<"\n";
		}
	}

Out1.close();
}

void Analysis::BunchFactor(double lambdau, double lambdar)
{
	std::ofstream Out1("Bfactor.txt", std::ofstream::out);
	double ku = 2.0*Pi/lambdau;
	double kr = 2.0*Pi/lambdar;
	int NOTS = B->getNOTS();
	int NOP = B->getNOP();
	for (int k = 0;k<NOTS;k++)
	{
	  
	  double sumz = 0;
	  double breal=0.0;
	  double bimag=0.0;
#pragma omp parallel for reduction(+:sumz,breal,bimag) private(ku,kr)
	  for (int i = 0;i<NOP;i++)
		{
			Vector XP = B->b[i].TrajPoint(k);
			sumz += XP.z;
			breal+=cos((ku+kr)*XP.z);
			bimag+=sin((ku+kr)*XP.z);
		}
		double avgz = sumz/(double)NOP;
		breal = breal/(double)NOP;
		bimag = bimag/(double)NOP;
		double b = sqrt(breal*breal + bimag*bimag);
		Out1<<avgz<<"\t"<<b<<"\n";
	}


Out1.close();
}

