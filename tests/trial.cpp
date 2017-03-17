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
  //const char *filename= "TEUFEL-EXAMPLE.txt";
  const char *filename= "BeamProfile.txt";
  Bunch *BB;
  int NOP=100;
  int NOTS=1700;
  double dt=3.2692307692307693e-12;
  //BB=new Bunch();
  BB= new Bunch(filename, NOP, -12500, 12500);
  Undulator *undu =new Undulator(0.3653,0.045,34);
  Lattice *field =new Lattice();
  field->addElement(undu);
  BB->Track_Vay(NOTS,dt,field);
  int FT=pow(2,16);
  tuple<Vector,Vector> Field[FT];
  double time_begin=5.48886e-9;
  double time_end=5.52795e-9;
  Vector Robs=Vector(0.0,0.0,4.0);
  double dt1 =(time_end-time_begin)/FT;
  tuple<Vector,Vector> Field1[FT];
  GenGrid Grid;
  Grid.GenPlanarGrid(0.006,0.006,4,150);
  double Area=Grid.GetMeshArea();
  Vector Ar=Vector(0.0,0.0,1.0)*Area;
  double t_max;
  tuple<Vector,Vector> Grid1[Grid.GetGridSize()];
  Analysis an = Analysis(BB);
  an.avg();
  an.DumpTrajectory();
  double lambdau = undu->GetLambdaU();
  double K = 0.934*lambdau*100*(undu->GetBPeak());
  double gamma =8.511/0.511;
  double lambdar = lambdau*(1+K*K/2.0)/(2.0*gamma*gamma);
  double freq = SpeedOfLight/lambdar;
  cout<<freq<<endl;
  an.BunchFactor(lambdau,lambdar);
  int output = BB->WriteSDDSTrajectory();
  int output2= BB->WriteSDDSTime();
  Bunch *BB2;
  BB2 = new Bunch(BB);
  for (int i = 0;i<NOP;i++)
	{
		
		BB2->b[i].MirrorY(0.014);
		//cout<<BB->b[i].TrajPoint(10).y+BB2->b[i].TrajPoint(10).y<<endl;
	}		
#pragma omp parallel for shared(time_begin, dt, Robs) 
  for (int i=0;i<FT;i++)
	{	
		
		Field[i]=BB->RadiationField(Robs,time_begin+i*dt1);
		Field1[i]=BB2->RadiationField(Robs,time_begin+i*dt1);
	}

  double Eavg = 0;
  for (int i=0;i<FT;i++)
	{
		Vector E=get<0>(Field[i])+get<0>(Field1[i]);
		Vector B=get<1>(Field[i])+get<1>(Field1[i]);
		Vector Poynting1=cross(E,B);		
		Out<<time_begin+i*dt1<<"\t"<<E.x<<"\t"<<Poynting1.z<<"\n";
	}
  Eavg = Eavg/FT;
#pragma omp parallel for 
  for (int i=0;i<Grid.GetGridSize();i++)
	{	
		Grid1[i]=BB->RadiationField(Grid.GetGridPoints(i),1.43464e-8);
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
