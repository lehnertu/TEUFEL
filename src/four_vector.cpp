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

#include "four_vector.h"
#include "global.h"
#include <math.h>
#include <iostream>

using namespace std;

Four_Vector::Four_Vector()
{
}

Four_Vector::Four_Vector(double t0,double x0, double y0, double z0)
{
  
  t =SpeedOfLight*t0;
  
  x=x0;
  
  y=y0;
  
  z=z0;
}

Four_Vector Four_Vector::operator+ (Four_Vector v2)
{
  Four_Vector temp;
  temp.t = t +v2.t;
  temp.x = x + v2.x;
  temp.y = y + v2.y;
  temp.z = z + v2.z;
  return (temp);
}

Four_Vector& Four_Vector::operator+= (Four_Vector v2)
{
  t+=v2.t;
  x+=v2.x;
  y+=v2.y;
  z+=v2.z;
  return (*this);
}

Four_Vector Four_Vector::operator- (Four_Vector v2)
{
  Four_Vector temp;
  temp.t = t - v2.t;
  temp.x = x - v2.x;
  temp.y = y - v2.y;
  temp.z = z - v2.z;
  return (temp);
}

Four_Vector Four_Vector::operator- ()
{
  Four_Vector temp;
  temp.t = -t;
  temp.x = -x;
  temp.y = -y;
  temp.z = -z;
  return (temp);
}

Four_Vector Four_Vector::operator* (double factor)
{
  Four_Vector temp;
  temp.t = t*factor;
  temp.x = x*factor;
  temp.y = y*factor;
  temp.z = z*factor;
  return (temp);
}

Four_Vector& Four_Vector::operator*= (double factor)
{
  t*=factor;
  x*=factor;
  y*=factor;
  z*=factor;
  return (*this);
}

Four_Vector Four_Vector::operator/ (double factor)
{
  Four_Vector temp;
  temp.t = t/factor;
  temp.x = x/factor;
  temp.y = y/factor;
  temp.z = z/factor;
  return (temp);
}

Four_Vector& Four_Vector::operator/= (double factor)
{
  t/=factor;
  x/=factor;
  y/=factor;
  z/=factor;
  return (*this);
}

double Four_Vector::norm()
{
  return(sqrt(t*t+x*x+y*y+z*z));
}

double Four_Vector::abs2nd()
{
  return(t*t+x*x+y*y+z*z);
}

void Four_Vector::normalize()
{
  double abs = norm();
  t /=abs;
  x /= abs;
  y /= abs;
  z /= abs;
}

double Four_Vector:: dot(Four_Vector a, Four_Vector b)
{
  double temp;
  temp = a.t * b.t;
  temp += -a.x * b.x;
  temp += -a.y * b.y;
  temp += -a.z * b.z;
  return (temp);
}

Four_Vector Four_Vector::Four_Momentum(double E, Vector p)
{
  Four_Vector temp;
  temp.t = E;
  p = p*SpeedOfLight;
  temp.x = p.x;
  temp.y = p.y;
  temp.z = p.z;
  return temp;
}

double Four_Vector::TransferMatrix(int i, int j,Vector S1)
{
	
	double gammabeta=S1.norm();
	double gamma = sqrt(gammabeta*gammabeta+1);
	double betax = S1.x/gamma;
	double betay = S1.y/gamma;
	double betaz = S1.z/gamma;
	double beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
	double L[4][4];
	L[0][0] = gamma;
	L[0][1] = -gamma*betax;
	L[0][2]	= -gamma*betay;
	L[0][3] = -gamma*betaz;
	L[1][0] = -gamma*betax;
	L[1][1] = 1+(gamma-1)*betax*betax/(beta*beta);
	L[1][2] = (gamma-1)*betax*betay/(beta*beta);
	L[1][3] = (gamma-1)*betax*betaz/(beta*beta);
	L[2][0] = -gamma*betay;
	L[2][1] = (gamma-1)*betax*betay/(beta*beta);
	L[2][2] = 1+(gamma-1)*betay*betay/(beta*beta);
	L[2][3] = (gamma-1)*betaz*betay/(beta*beta);
	L[3][0] = -gamma*betaz;
	L[3][1] = (gamma-1)*betax*betaz/(beta*beta);
	L[3][2] = (gamma-1)*betaz*betay/(beta*beta);
	L[3][3] = 1+(gamma-1)*betaz*betaz/(beta*beta);
	return (L[i][j]);
}

double Four_Vector::InverseTransferMatrix(int i, int j ,Vector S1)
{
	//S1 is a vector having gamma beta of the moving frame;
	//so that vector S1=Vector(gammabetax, gammabetay,gammabetaz)
	
	double gammabeta=S1.norm();
	double gamma = sqrt(gammabeta*gammabeta+1);
	double betax = S1.x/gamma;
	double betay = S1.y/gamma;
	double betaz = S1.z/gamma;
	double beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
	double L[4][4];
	L[0][0] = gamma;
	L[0][1] = gamma*betax;
	L[0][2]	= gamma*betay;
	L[0][3] = gamma*betaz;
	L[1][0] = gamma*betax;
	L[1][1] = 1+(gamma-1)*betax*betax/(beta*beta);
	L[1][2] = -(gamma-1)*betax*betay/(beta*beta);
	L[1][3] = -(gamma-1)*betax*betaz/(beta*beta);
	L[2][0] = gamma*betay;
	L[2][1] = -(gamma-1)*betax*betay/(beta*beta);
	L[2][2] = 1+(gamma-1)*betay*betay/(beta*beta);
	L[2][3] = -(gamma-1)*betaz*betay/(beta*beta);
	L[3][0] = gamma*betaz;
	L[3][1] = -(gamma-1)*betax*betaz/(beta*beta);
	L[3][2] = -(gamma-1)*betaz*betay/(beta*beta);
	L[3][3] = 1+(gamma-1)*betaz*betaz/(beta*beta);
	
	return L[i][j];
}

Four_Vector Four_Vector:: LorentzTransform( Vector S1)
{
	//generate the transformation matrix for the given relative motion of frame S1
	// its velocities (expected gamma beta)
	
	double B[4];
	B[0]=SpeedOfLight*t;
	B[1]=x;
	B[2]=y;
	B[3]=z;
	double A[4]={0,0,0,0};
	for (int i=0;i<4;i++)
	{
		for (int k=0;k<4;k++)	
		{
			A[i]+=TransferMatrix(i,k,S1)*B[k];
			cout<<TransferMatrix(i,k,S1)<<"\t";
		}
		cout<<"\n";

	}
	Four_Vector Inverted = Four_Vector(A[0]/SpeedOfLight,A[1],A[2],A[3]);
	return Inverted;
}

Four_Vector Four_Vector:: InverseLorentzTransform( Vector S1)
{
	//generate the transformation matrix for the given relative motion of frame S1
	// its velocities (expected gamma beta)
	double B[4];
	B[0]=SpeedOfLight*t;
	B[1]=x;
	B[2]=y;
	B[3]=z;
	double A[4]={0,0,0,0};
	for (int i=0;i<4;i++)
	{
		for (int k=0;k<4;k++)	
		{
			A[i]+=InverseTransferMatrix(i,k,S1)*B[k];
		}
	

	}
	Four_Vector Inverted = Four_Vector(A[0]/SpeedOfLight,A[1],A[2],A[3]);
	return Inverted;
}
