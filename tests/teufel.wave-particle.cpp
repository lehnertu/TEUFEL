/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    Homogeneous magnetic field test case

  Copyright (c) 2017 U. Lehnert

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

/*!
    \brief Plane Wave - Electron Interaction

    @author Ulf Lehnert
    @author Vipul Joshi
    @date 10.2.2017
    @file teufel.magnet.cpp
    
    This test case simulates the interaction between a plane wave and a electron.
    The electron performs oscillations along x axis due to electric field of the wave.
    Except Ex, the other electric field components are zero.
    Due to oscillations of electron, it emits radiation.
    The radiation emitted is observed at a screen placed at a distance 'a';
    a is chosen so that the near field effects are negligible and only the far field radiation is observed.
    The equation for radiation field emitted along x axis and z axis reduce to following form:
    \f[
	E_{x} = 0
    \f]
    \f[
	\E_{z} = \frac{q \dot{\beta}}{4 \pi \epsilon_{0} c a}
    \f]

    The program computes the trajectory of the electron which starts moving under the force field of plane wave.
    The power density and wavelength of the plane wave chosen is 10e-4 Watts/cm2, 10e-4 m (3Thz) respectively.
    The radiation emitted by the motion of electron is captured at a screen placed and compared with the known values.
    Following points are then checked:
    @li the arrival time of the radiation is correct. Must lie in 1% of Expected Arrival Time i.e. Robs/c
    @li the peak value of radiation oberved along z axis compares with analytical formula
    @li the radiation at any arbitrary point along x axis is zero
    @li the overall displacement along z axis is less than a micron   
    @return The number of errors encountered in the above list of checks is reported.
    
 */


#include "particle.h"
#include "undulator.h"
#include <iostream>
#include <stdlib.h>
#include "global.h"
#include <vector>
#include <random>
#include <fstream>
#include "externalfield.h"
#include "vector.h"
#include "bunch.h"
#include "wave.h"
#include "gen_grid.h"
using namespace std;

int main()
{
	// create an empty bunch and add a particle resting at origin to it.
	Bunch *BB= new Bunch();
	ChargedParticle *e1 = new ChargedParticle(-1,1,Vector(0,0,0.00),Vector(0,0,0),0);
	BB->AddParticles(e1);


	//define plane wave and lattice 
	Lattice *FEL=new Lattice();
	PlaneWave *Wave = new PlaneWave(0.0001,0.00001);
	FEL->addElement(Wave);

	// define the tracking parameters i.e. the time steps, total number of time steps
	double dt = 2.5e-15;
	int NOTS = 20000;
	BB->Track_Vay(NOTS, dt, FEL,1);

	//define the obervation point (Robs), time window (to observe radiation) and the number of points in the time window
	double a = 3.0;
	double time_begin = 1.0005988e-8;
	double time_end =1.0058711004e-8;
	int numPoints=10000;
	double dt1 = (time_end-time_begin)/10000.0;
	
	// get the peak electric field value for the plane wave and find maximum value particle acceleration i.e \dot{\beta}
	double Epeak1 = Wave->EPeak();
	double Omega = Wave->Omega();
	double beta_dot = ElementaryCharge*Epeak1/(m_e*SpeedOfLight);
	double scale = ElementaryCharge/(4.0*Pi*EpsNull);
	double ExpectedEpeak = scale*beta_dot/(SpeedOfLight*a);
	double ExpectedArrivalTime = a/3.0e8;
	
	cout<<"Expected Peak Value of Radiation Emitted: "<<ExpectedEpeak<<endl;
	Vector Epeak = Vector(0,0,0);
	int counter = 0;
	double arrivalTime=0;
	for (int i=0;i<10000;i++)
	{
		tuple<Vector,Vector>Fields=BB -> RadiationField(Vector(0,0,a),time_begin+i*dt1);
		if(Epeak.norm() < get<0>(Fields).norm())
		{
			Epeak = get<0>(Fields);
		}
		
		//arrival of the radiation pulse is assumed to be when x Efield reaches 1% of the Expected Epeak 
		if(fabs(Epeak.x)>0.1*ExpectedEpeak && counter==0)
		{
			arrivalTime = time_begin+i*dt;
			counter +=1;
		}

	}

	//find the final displacement in z of the particle and compare it to 0.1 micron --> Ideally it should be zero
	double FinalZ = BB->getTrajPoint(0,NOTS-1).z;

	//find radiation at observation point (3,0,0) i.e. axial field after coloumb part has been removed. Must be zero. Time is a random point in the 'Time Window'
	tuple<Vector,Vector>Fields=BB -> RadiationField(Vector(a,0,0),1.006e-8);
	double E_Axial = get<0>(Fields).norm();
	double CoulombField = scale/pow(3.0,2);
	int errors=0;
	if (fabs(ExpectedEpeak-fabs(Epeak.x)) > ExpectedEpeak*0.01)
	{
	errors++;
	printf("Epeak= %12.9g V/m - \033[1;31m test failed!\033[0m\n", Epeak.x);
    	}
	else
	{
	printf("Epeak= %12.9g V/m - \033[1;32m OK\033[0m\n", Epeak.x);
    	}

	if (fabs(ExpectedArrivalTime-arrivalTime) > ExpectedArrivalTime*0.01)
	{
	errors++;
	printf("Arrival Time= %12.9g s - \033[1;31m test failed!\033[0m\n", arrivalTime);
    	}
	else
	{
	printf("Arrival Time= %12.9g s - \033[1;32m OK\033[0m\n", arrivalTime);
    	}

	if (fabs(FinalZ) > 0.1e-6)
	{
	errors++;
	printf("final z= %12.9g m - \033[1;31m test failed!\033[0m\n", FinalZ);
    	}
	else
	{
	printf("final z= %12.9g m - \033[1;32m OK\033[0m\n", FinalZ);
    	}

	if (fabs(E_Axial-CoulombField) > (ExpectedEpeak)*0.001)
	{
	errors++;
	printf("Axial Field Strength = %12.9g V/m - \033[1;31m test failed!\033[0m\n", E_Axial-CoulombField);
    	}
	else
	{
	printf("Axial Field Strength = %12.9g V/m - \033[1;32m OK\033[0m\n", E_Axial-CoulombField);
    	}

	return errors;
}
