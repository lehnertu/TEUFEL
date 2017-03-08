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
#include <vector>
#include <tuple>

using namespace std;

Four_Vector::Four_Vector()
{
}


// creates a four vector with coordinates (time, x,y,z)
Four_Vector::Four_Vector(double t0,double x0, double y0, double z0)
{
  
  t =t0;
  
  x=x0;
  
  y=y0;
  
  z=z0;
}

// operators that work with Four vectors


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
  return(sqrt(SpeedOfLight*SpeedOfLight*t*t-(x*x+y*y+z*z)));
}

double Four_Vector::abs2nd()
{
  return(SpeedOfLight*SpeedOfLight*t*t-(x*x+y*y+z*z));
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


// methods for four vector transformations
// methods for EM field transformations





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
	vector<vector<double>> L=TransferMatrix(S1);
	
	for (int i=0;i<4;i++)
	{
		double value = 0;
		for (int k=0;k<4;k++)	
		{
			value+=L[i][k]*B[k];
		}

		A[i] = value;
		

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
	vector<vector<double>> L=InverseTransferMatrix(S1);
	
	for (int i=0;i<4;i++)
	{

		double value = 0;
		for (int k=0;k<4;k++)	
		{
			value+=L[i][k]*B[k];			
		}
		
		A[i] = value;
	}
	Four_Vector Inverted = Four_Vector(A[0]/SpeedOfLight,A[1],A[2],A[3]);
	return Inverted;
}
tuple<Vector,Vector> EMTransform(Vector S1,Vector EField,Vector BField)
{
	vector<vector<double>> L(4,vector<double>(4));
	vector<vector<double>> FPrime(4,vector<double>(4));
	vector<vector<double>> L1(4,vector<double>(4));
	vector<vector<double>> LT(4,vector<double>(4));
	FPrime = EMTensor(EField,BField);
	L = TransferMatrix(S1);
	L1= Multiply(L,FPrime);
	LT = TransposeMatrix(L);
	FPrime = Multiply(L1,LT);
	Vector E = Vector(FPrime[0][1],FPrime[0][2],FPrime[0][3])*SpeedOfLight;
	Vector B = Vector(FPrime[2][3],-FPrime[1][3],FPrime[1][2]);
	return make_tuple(E,B);		

}

tuple<Vector,Vector>EMInverseTransform(Vector S1,Vector EField,Vector BField)
{
	vector<vector<double>> L(4,vector<double>(4));
	vector<vector<double>> F(4,vector<double>(4));
	vector<vector<double>> L1(4,vector<double>(4));
	vector<vector<double>> LT(4,vector<double>(4));
	F = InverseEMTensor(EField,BField);
	L = TransferMatrix(S1);
	L1 = InvertMatrix(L);
	L1= Multiply(L1,F);
	LT = TransposeMatrix(L);
	LT = InvertMatrix(LT);
	F = Multiply(L1,LT);
	for (int i=0;i<4;i++)
	{
		for (int k=0;k<4;k++)	
		{
			cout<<F[i][k]<<"\t";
		}	
		cout<<"\n";
	}
	Vector E = Vector(F[0][1],F[0][2],F[0][3]);
	Vector B = Vector(F[2][3],-F[1][3],F[1][2]);
	return make_tuple(E,B);		

}


//generates the transformation matrix
//for transformation from lab frame to frame S1
// Vector S1 gives the momentum (gammabetax,gammabetay,gammabetaz)
// tranformation operation is defined as matrix operation
// X1 = L X
//where X1 is the four vector in S1 frame
// L is the Transfer matrix
// X is the four vector in S frame (lab)
vector<vector<double>> TransferMatrix(Vector S1)
{
	
	double gammabeta=S1.norm();
	double gamma = sqrt(gammabeta*gammabeta+1);
	double betax = S1.x/gamma;
	double betay = S1.y/gamma;
	double betaz = S1.z/gamma;
	//cout<<betax<<"\t"<<betay<<"\t"<<betaz<<"\n";
	double beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
	vector<vector<double>> L(4,vector<double>(4));
	L[0][0] = gamma;
	L[0][1] = -gamma*betax;
	L[0][2]	= -gamma*betay;
	L[0][3] = -gamma*betaz;
	L[1][0] = -gamma*betax;
	L[1][1] = 1.0+(gamma-1.0)*betax*betax/(beta*beta);
	L[1][2] = (gamma-1.0)*betax*betay/(beta*beta);
	L[1][3] = (gamma-1.0)*betax*betaz/(beta*beta);
	L[2][0] = -gamma*betay;
	L[2][1] = (gamma-1.0)*betax*betay/(beta*beta);
	L[2][2] = 1.0+(gamma-1.0)*betay*betay/(beta*beta);
	L[2][3] = (gamma-1.0)*betaz*betay/(beta*beta);
	L[3][0] = -gamma*betaz;
	L[3][1] = (gamma-1.0)*betax*betaz/(beta*beta);
	L[3][2] = (gamma-1.0)*betaz*betay/(beta*beta);
	L[3][3] = 1.0+(gamma-1.0)*betaz*betaz/(beta*beta);
	
	return L;
}


//generates the inverse transformation matrix
//for transformation of a four vector X1 from frame S1 to S
// Matrix operation is defined as
// X = L_inverse X1
// Vector S1 gives the description of a moving frame
vector<vector<double>> InverseTransferMatrix(Vector S1)
{
	//S1 is a vector having gamma beta of the moving frame;
	//so that vector S1=Vector(gammabetax, gammabetay,gammabetaz)
	
	double gammabeta=S1.norm();
	double gamma = sqrt(gammabeta*gammabeta+1);
	double betax = S1.x/gamma;
	double betay = S1.y/gamma;
	double betaz = S1.z/gamma;
	double beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
	vector<vector<double>> L(4,vector<double>(4));
	L[0][0] = gamma;
	L[0][1] = gamma*betax;
	L[0][2]	= gamma*betay;
	L[0][3] = gamma*betaz;
	L[1][0] = gamma*betax;
	L[1][1] = 1+(gamma-1)*betax*betax/(beta*beta);
	L[1][2] = (gamma-1)*betax*betay/(beta*beta);
	L[1][3] = (gamma-1)*betax*betaz/(beta*beta);
	L[2][0] = gamma*betay;
	L[2][1] = (gamma-1)*betax*betay/(beta*beta);
	L[2][2] = 1+(gamma-1)*betay*betay/(beta*beta);
	L[2][3] = (gamma-1)*betaz*betay/(beta*beta);
	L[3][0] = gamma*betaz;
	L[3][1] = (gamma-1)*betax*betaz/(beta*beta);
	L[3][2] = (gamma-1)*betaz*betay/(beta*beta);
	L[3][3] = 1+(gamma-1)*betaz*betaz/(beta*beta);
		
	return L;
}


// generates 4-4 matrix for EM Field transformation
// from the electric and magnetic field vectors

/* 0	-Ex	-Ey	-Ez
   
   Ex	0	-Bz	By

   Ey	Bz	0	Bx

   Ez	-By	Bx	0

*/
vector<vector<double>>EMTensor(Vector Efield,Vector BField)
{
	vector<vector<double>> L(4,vector<double>(4));
	L[0][0] = 0.0;
	L[0][1] = Efield.x;
	L[0][2]	= Efield.y;
	L[0][3] = Efield.z;
	L[1][0] = -Efield.x;
	L[1][1] = 0.0;
	L[1][2] = BField.z;
	L[1][3] = -BField.y;
	L[2][0] = -Efield.y;
	L[2][1] = -BField.z;
	L[2][2] = 0.0;
	L[2][3] = BField.x;
	L[3][0] = -Efield.z;
	L[3][1] = BField.y;
	L[3][2] = -BField.x;
	L[3][3] = 0.0;
	
	return L;

}
vector<vector<double>>InverseEMTensor(Vector Efield,Vector BField)
{
	vector<vector<double>> L(4,vector<double>(4));
	L[0][0] = 0.0;
	L[0][1] = Efield.x/SpeedOfLight;
	L[0][2]	= Efield.y/SpeedOfLight;
	L[0][3] = Efield.z/SpeedOfLight;
	L[1][0] = -Efield.x/SpeedOfLight;
	L[1][1] = 0.0;
	L[1][2] = BField.z;
	L[1][3] = -BField.y;
	L[2][0] = -Efield.y/SpeedOfLight;
	L[2][1] = -BField.z;
	L[2][2] = 0.0;
	L[2][3] = BField.x;
	L[3][0] = -Efield.z/SpeedOfLight;
	L[3][1] = BField.y;
	L[3][2] = -BField.x;
	L[3][3] = 0.0;
	
	return L;

}


// matrix operations


// 
vector<vector<double>> TransposeMatrix(vector<vector<double>> S1)
{
	vector<vector<double>> L(4,vector<double>(4));
	// do the transpose
	for (int i=0;i<4;i++)
	{
		for (int k=0;k<4;k++)	
		{
			L[k][i] = S1[i][k];
		}	
	}
	
	return L;
}

//invert a 4 by 4 matrix
// S1 is a 4-4 matrix
vector<vector<double>>InvertMatrix(vector<vector<double>> S1)
{
	vector<vector<double>>L(4,vector<double>(4));
	double a11,a12,a13,a14,a22,a23,a21,a24,a31,a32,a33,a34,a41,a42,a43,a44;
	double b11,b12,b13,b14,b22,b23,b21,b24,b31,b32,b33,b34,b41,b42,b43,b44;
	a11 = S1[0][0];
	a12 = S1[0][1];
	a13 = S1[0][2];
	a14 = S1[0][3];
	a21 = S1[1][0];
	a22 = S1[1][1];
	a23 = S1[1][2];
	a24 = S1[1][3];
	a31 = S1[2][0];
	a32 = S1[2][1];
	a33 = S1[2][2];
	a34 = S1[2][3];
	a41 = S1[3][0];
	a42 = S1[3][1];
	a43 = S1[3][2];
	a44 = S1[3][3];

	double det= a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43+a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41+a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42+a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41-a11*a22*a34*a43- a11*a23*a32*a44-a11*a24*a33*a42-a12*a21*a33*a44-a12*a23*a34*a41 -a12*a24*a31*a43-a13*a21*a34*a42-a13*a22*a31*a44-a13*a24*a32*a41-a14*a21*a32*a43-a14*a22*a33*a41-a14*a23*a31*a42;

	if (det != 0.0)
	{
		b11 = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42;
		b12 = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43;
		b13 = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42;
		b14 = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33;
		b21 = a21*a34*a43 + a23*a31*a44 + a24*a33*a41 - a21*a33*a44 - a23*a34*a41 - a24*a31*a43;
		b22 = a11*a33*a44 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 - a14*a33*a41;

		b23 = a11*a24*a43 + a13*a21*a44 + a14*a23*a41 - a11*a23*a44 - a13*a24*a41 - a14*a21*a43;
		b24 = a11*a23*a34 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 - a14*a23*a31;
		b31 = a21*a32*a44 + a22*a34*a41+ a24*a31*a42 - a21*a34*a42 - a22*a31*a44 - a24*a32*a41;
		b32 = a11*a34*a42 + a12*a31*a44 + a14*a32*a41 - a11*a32*a44 - a12*a34*a41 - a14*a31*a42;
		b33 = a11*a22*a44 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 - a14*a22*a41;
		b34 = a11*a24*a32 + a12*a21*a34 + a14*a22*a31 - a11*a22*a34 - a12*a24*a31 - a14*a21*a32;
		b41 = a21*a33*a42 + a22*a31*a43 + a23*a32*a41 - a21*a32*a43 - a22*a33*a41 - a23*a31*a42;
		b42 = a11*a32*a43 + a12*a33*a41+ a13*a31*a42 - a11*a33*a42 - a12*a31*a43 - a13*a32*a41;
		b43 = a11*a23*a42 + a12*a21*a43 + a13*a22*a41- a11*a22*a43 - a12*a23*a41 - a13*a21*a42;
		b44 = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31;

	

	L[0][0] = (1/det)*b11;
	L[0][1] = (1/det)*b12;
	L[0][2]	= (1/det)*b13;
	L[0][3] = (1/det)*b14;
	L[1][0] = (1/det)*b21;
	L[1][1] = (1/det)*b22;
	L[1][2] = (1/det)*b23;
	L[1][3] = (1/det)*b24;
	L[2][0] = (1/det)*b31;
	L[2][1] = (1/det)*b32;
	L[2][2] = (1/det)*b33;
	L[2][3] = (1/det)*b34;
	L[3][0] = (1/det)*b41;
	L[3][1] = (1/det)*b42;
	L[3][2] = (1/det)*b43;
	L[3][3] = (1/det)*b44;
	return L;
	}

	else
	{
		cout<<"Matrix Non Invertible - det A = 0"<<endl;
	}	

	
	
}

//Multiply 4 by 4 m
vector<vector<double>>Multiply(vector<vector<double>> S1,vector<vector<double>> S2)
{

	vector<vector<double>> L(4,vector<double>(4));
	for (int i=0;i<4;i++)
	{
		
		for (int k=0;k<4;k++)	
		{
			double value;
			value = 0;
			for (int in = 0;in<4;in++)
			{
				value=value+S1[i][in]*S2[in][k];
			};
			L[i][k] = value;
		};	
		
	}
	
	return L;
}

