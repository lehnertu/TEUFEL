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
  const char *filename= "TEUFEL-EXAMPLE.txt";
  //const char *filename= "BeamProfile.txt";
  Bunch *BB;
  int NOP=999;
  int NOTS=1700;
  double dt=3.2692307692307693e-12;
  //BB=new Bunch();
  BB= new Bunch(filename, NOP, -12500, 12500);
  ChargedParticle *e = new ChargedParticle(-12500,12500,Vector(0,0,0),Vector(0,0,16),0);
  cout<<"Now trying to add particles to a bunch"<<endl;
  //BB->AddParticles(e);
  Undulator *undu =new Undulator(0.3653,0.045,34);
  Lattice *field =new Lattice();
  field->addElement(undu);
  cout<<"beginning tracking of particles"<<endl;
  BB->Track_Vay(NOTS,dt,field);
  int FT=pow(2,16);
  tuple<Vector,Vector> Field[FT];
  double time_begin=1.33e-8;
  double time_end=1.35e-8;
  Vector Robs=Vector(0.0,0.0,4.0);
  double dt1 =(time_end-time_begin)/(double)FT;
  tuple<Vector,Vector> Field1[FT];
  GenGrid Grid;
  Grid.GenPlanarGrid(0.006,0.006,4,150);
  double Area=Grid.GetMeshArea();
  Vector Ar=Vector(0.0,0.0,1.0)*Area;
  double t_max;
  tuple<Vector,Vector> Grid1[Grid.GetGridSize()];
 // Analysis an = Analysis(BB);
 // an.avg();
 // an.DumpTrajectory();
  double lambdau = undu->GetLambdaU();
  double K = 0.934*lambdau*100*(undu->GetBPeak());
  double gamma =8.511/0.511;
  double lambdar = lambdau*(1+K*K/2.0)/(2.0*gamma*gamma);
  double freq = SpeedOfLight/lambdar;
  //cout<<freq<<endl;
  //an.BunchFactor(lambdau,lambdar);
  //ChargedParticle *e2 = new ChargedParticle(e);
  cout<<"Now trying to copy a bunch"<<endl;
  Bunch *BB2= new Bunch(BB);
  BB2->MirrorY(0.012);
  /*for (int i=0;i<BB2->getNOP();i++)
	{

		BB2->b[i]->MirrorY(0.012);
	}*/
 
#pragma omp parallel for shared(time_begin, dt1, Robs) 
  for (int i=0;i<FT;i++)
	{	
		
		//Field[i]=BB->RadiationField(Robs,time_begin+i*dt1);
		Field1[i]=BB2->RadiationField(Robs,time_begin+i*dt1);
	}

  
  for (int i=0;i<FT;i++)
	{
		Vector E=get<0>(Field1[i]);//+get<0>(Field[i]);
		Vector B=get<1>(Field1[i]);//+get<1>(Field[i]);
		Vector Poynting1=cross(E,B);		
		Out<<time_begin+i*dt1<<"\t"<<E.x<<"\t"<<Poynting1.z<<"\n";
	}
  
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

  BB->JoinBunch(BB2);
  Out.close();  
  Out1.close();
  int output = BB->WriteSDDSTrajectory();
  cout<<"Write Trajectory return value: "<<output<<endl;
  int output2= BB->WriteSDDSTime();
  cout<<"Write Time-Trajectory return value: "<<output2<<endl;	
  return 0;
}
