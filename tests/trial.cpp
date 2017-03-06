#include "bunch.h"
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "omp.h"
#include <math.h>
#include "analysis.h"
#include "gen_grid.h"
#include "global.h"

using namespace std;
int main()
{
  std::ofstream Out("test.txt", std::ofstream::out);
  std::ofstream Out1("test1.txt", std::ofstream::out);
  const char *filename= "BeamProfile.txt";
  Bunch *BB;
  int NOP=1000;
  int NOTS=3000;
  double dt=1.4e-12;
  //BB=new Bunch();
  BB= new Bunch(filename, NOP, -1000, 1000);
  Undulator *undu =new Undulator(0.1,0.045,25);
  Lattice *field =new Lattice();
  field->addElement(undu);
  BB->Track_Vay(NOTS,dt,field);
  int FT=3000;
  tuple<Vector,Vector> Field[FT];
  double time_begin=1.03406e-8;
  double time_end=1.03501e-8;
  Vector Robs=Vector(0.0,0.0,3.0);
  double dt1 =(time_end-time_begin)/FT;
  tuple<Vector,Vector> Field1[FT];
  GenGrid Grid;
  Grid.GenPlanarGrid(0.01,0.01,3.0,150);
  double Area=Grid.GetMeshArea();
  Vector Ar=Vector(0.0,0.0,1.0)*Area;
  double t_max;
  tuple<Vector,Vector> Grid1[Grid.GetGridSize()];
  Analysis(BB).avg();
 
#pragma omp parallel for shared(time_begin, dt, Robs) 
  for (int i=0;i<FT;i++)
	{	
		Field[i]=BB->RadiationField(Robs,time_begin+i*dt1);
	}
  for (int i=0;i<FT;i++)
	{
		Vector E=get<0>(Field[i]);
		Vector B=get<1>(Field[i]);
		Vector Poynting1=cross(E,B);
		
		Out<<time_begin+i*dt1<<"\t"<<E.x<<"\t"<<B.y<<"\n";
	}
#pragma omp parallel for 
  for (int i=0;i<Grid.GetGridSize();i++)
	{	
		Grid1[i]=BB->RadiationField(Grid.GetGridPoints(i),1.034545e-8);
	}
  for (int i=0;i<Grid.GetGridSize();i++)
	{
		Vector E=get<0>(Grid1[i]);
		Vector B=get<1>(Grid1[i]);
		Vector Poynting = cross(E,B);
		Vector XP=Grid.GetGridPoints(i);
		double Power=dot(Poynting,Ar);		
		Out1<<XP.x<<"\t"<<XP.y<<"\t"<<XP.z<<"\t"<<Power<<"\n";
	}

  Out.close();  
  Out1.close();
  
  return 0;
}
