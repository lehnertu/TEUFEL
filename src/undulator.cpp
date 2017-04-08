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

#include "undulator.h"
#include "global.h"

#include <math.h>

Undulator::Undulator(   double B,                         // peak field [T]
                        double lambda,                    // undulator period [m]
                        int    N                          // number of undulator periods
             )
{
  BPeak = B;
  LambdaU = lambda;
  NPeriods = N;
  Krms=LambdaU*SpeedOfLight*BPeak/(2.0*Pi*mecsquared)/sqrt(2.0);
  ky = 2.0*Pi/LambdaU;
  kz = 2.0*Pi/LambdaU;
}

double Undulator::GetBPeak()
{
  return BPeak;
}

double Undulator::GetLambdaU()
{
  return LambdaU;
}

int Undulator::GetNPeriods()
{
  return NPeriods;
}

double Undulator::GetKpeak()
{
  return Krms*sqrt(2.0);
}

double Undulator::GetKrms()
{
  return Krms;
}

ElMagField Undulator::LocalField(double t, Vector X)
{
    Vector E = Vector(0.0 ,0.0 ,0.0);
    Vector B=Vector (0,0,0);
    double z1=-GetLambdaU()/2;
    double z2=GetLambdaU()/2;
    B.x=0.0;
    if(X.z<z1)
    {
	B.y=0;
	B.z=0;
    }
    if(z1<=X.z&&X.z<z2)
    {
	B.y=((X.z-z1)/(z2-z1))*BPeak*sin(kz*X.z)*cosh(kz*X.y);
	B.z=((X.z-z1)/(z2-z1))*BPeak*cos(kz*X.z)*sin(kz*X.y);
    }
    if(z2<=X.z&&X.z<LambdaU*NPeriods-z2)
    {
	B.y=BPeak*sin(kz*X.z)*cosh(kz*X.y);
	B.z=BPeak*cos(kz*X.z)*sinh(kz*X.y);
    }
    if(LambdaU*NPeriods-z2<=X.z&&X.z<=LambdaU*NPeriods-z1)
    {
	B.y=((LambdaU*NPeriods-z1-X.z)/(z2-z1))*BPeak*sin(kz*X.z)*cosh(kz*X.y);
	B.z=((LambdaU*NPeriods-z1-X.z)/(z2-z1))*BPeak*cos(kz*X.z)*sinh(kz*X.y);
    }
    return ElMagField(E,B);
}
