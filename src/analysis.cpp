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

using namespace std;

Analysis::Analysis(Bunch *bunch)
{
  //create a copy of bunch
  B=bunch;
  
}

void Analysis::avg()
{
	std::ofstream Out2("trajectory.txt", std::ofstream::out);
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
	  double sum2pz = 0.0;
	  double stdx = 0.0;
	  double stdy = 0.0;
	  double stdz = 0.0;
	  double stdpx=0.0;
	  double stdpy=0.0;
	  double stdpz=0.0;
#pragma omp parallel for reduction(+:sumx,sumy,sumz,sumpx,sumpy,sumpz,sum2x,sum2y,sum2z,sum2px,sum2py,sum2pz)
	  for (int k=0; k<NOP;k++)
		{
			Vector XP =(B->b[k].TrajPoint(i));
			Vector XP1=(B->b[k].TrajMomentum(i));
			sumx += XP.x;
			sumy += XP.y;
			sumz += XP.z;
			sumpx += XP1.x;
			sumpy += XP1.y;
			sumpz += XP1.z;
			sum2x += pow(XP.x,2);
			sum2y += pow(XP.y,2);
			sum2z += pow(XP.z,2);
			sum2px += pow(XP1.x,2);
			sum2py += pow(XP1.y,2);
			sum2pz += pow(XP1.z,2);
		}
		
		double tt=B->b[0].TrajTime(i);
		Out2<<tt<<"\t"<<sumx/(double)NOP<<"\t"<<sumy/(double)NOP<<"\t"<<sumz/(double)NOP<<"\t"<<sumpx/(double)NOP<<"\t"<<sumpy/(double)NOP<<"\t"<<sumpz/(double)NOP<<"\n";
		stdx = sqrt(sum2x/(double)NOP - pow(sumx/(double)NOP,2));
		stdy = sqrt(sum2y/(double)NOP - pow(sumy/(double)NOP,2));
		stdz = sqrt(sum2z/(double)NOP - pow(sumz/(double)NOP,2));
		stdpx = sqrt(sum2px/(double)NOP - pow(sumpx/(double)NOP,2));
		stdpy = sqrt(sum2py/(double)NOP - pow(sumpy/(double)NOP,2));
		stdpz = sqrt(sum2pz/(double)NOP - pow(sumpz/(double)NOP,2));
		Out3<<tt<<"\t"<<stdx<<"\t"<<stdy<<"\t"<<stdz<<"\t"<<stdpx<<"\t"<<stdpy<<"\t"<<stdpz<<"\n";
  }

Out2.close();
Out3.close();
}

