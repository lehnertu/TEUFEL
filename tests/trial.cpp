#include "bunch.h"
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "omp.h"
#include <math.h>
using namespace std;
int main()
{
  std::ofstream Out ("test.txt", std::ofstream::out);
  const char *filename= "BeamProfile.txt";
  Bunch *BB= new Bunch();
  int NOP=1000;
  int NOTS=2000;
  double dt=3.6e-12;
  BB= new Bunch(filename, NOP, -100, 100);
  Undulator *undu =new Undulator(0.534,0.037,54);
  Lattice *field =new Lattice();
  field->addElement(undu);
  BB->Track_Euler(NOTS,dt,field);
  int FT=pow(2,12);
  tuple<Vector,Vector> Field[FT];
  double time_begin=1.03398e-8;
  double time_end=1.03506e-8;
  double dt1;
#pragma omp parallel for shared(time_end,time_begin, dt) 
  for (int i=0;i<FT;i++)
	{	
		dt1 =(time_end-time_begin)/(double)FT;
		Vector Robs=Vector(0.0,0.0,3.0);
		Field[i]=BB->RadiationField(Robs,time_begin+i*dt1);
	}
  for (int i=0;i<FT;i++)
	{
		Vector E=get<0>(Field[i]);
		Vector B=get<1>(Field[i]);
		Vector Poynting=cross(E,B);
		
		Out<<time_begin+i*dt1<<"\t"<<get<0>(Field[i]).x<<"\t"<<Poynting.norm()<<"\n";
	}
  Out.close();
  return 0;
}
