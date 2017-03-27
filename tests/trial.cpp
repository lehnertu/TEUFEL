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
  int NOP=100000;
  int NOTS=1700;
  double dt=3.2692307692307693e-12;
  //BB=new Bunch();
  BB= new Bunch(filename, NOP, -12500, 12500);
  //ChargedParticle *e = new ChargedParticle(-12500,12500,Vector(0,0,0),Vector(0,0,16),0);
  //cout<<"Now trying to add particles to a bunch"<<endl;
  //BB->AddParticles(e);
  Undulator *undu =new Undulator(0.1,0.045,34);
  Lattice *field =new Lattice();
  field->addElement(undu);
  BB->Track_Vay(NOTS,dt,field,0);
  int FT=pow(2,14);
  tuple<Vector,Vector> Field[FT];
  double time_begin=1.331e-8;
  double time_end=1.341e-8;
  Vector Robs=Vector(0.0,0.0,4.0);
  double dt1 =(time_end-time_begin)/(double)FT;
  tuple<Vector,Vector> Field1[FT];
  GenGrid Grid;
  Grid.GenPlanarGrid(0.006,0.006,4,150);
  double Area=Grid.GetMeshArea();
  Vector Ar=Vector(0.0,0.0,1.0)*Area;
  //double t_max;
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
#pragma omp parallel for 
  for (int i=0;i<Grid.GetGridSize();i++)
	{	
		Grid1[i]=BB->RadiationField(Grid.GetGridPoints(i),1.335428478e-8);
	}
  for (int i=0;i<Grid.GetGridSize();i++)
	{
		Vector E=get<0>(Grid1[i]);
		Vector B=get<1>(Grid1[i])/MuNull;
		Vector Poynting = cross(E,B);
		Vector XP=Grid.GetGridPoints(i);
		double Power=dot(Poynting,Ar);		
		Out1<<XP.x<<"\t"<<XP.y<<"\t"<<XP.z<<"\t"<<Poynting.norm()<<"\n";
	}

  //BB->JoinBunch(BB2);
  Out.close();  
  Out1.close();
  int output = BB->WriteSDDSTrajectory();
  cout<<"Write Trajectory return value: "<<output<<endl;
  int output2= BB->WriteSDDSTime();
  cout<<"Write Time-Trajectory return value: "<<output2<<endl;	
  int output3 = BB->WriteSDDSRadiation(Robs, time_begin, time_end, FT);
  cout<<"Write Radiation return value: "<<output2<<endl;
  BB->~Bunch();
  //e->~ChargedParticle();	
  return 0;
}
