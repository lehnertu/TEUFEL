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
  //const char *filename= "BeamProfile.txt";
  Bunch *BB;
  //int NOP=200;
  int NOTS=34000;
  double dt=3.2692307692307693e-13/2.0;
  BB=new Bunch();
  //BB= new Bunch(filename, NOP, -12500, 12500);
  ChargedParticle *e = new ChargedParticle(-12500,12500,Vector(0,0,0),Vector(0,0,16),0);
  ChargedParticle *e2 = new ChargedParticle(-12500,12500,Vector(0,0,-0.000553/2.0),Vector(0,0,16),0);
  //cout<<"Now trying to add particles to a bunch"<<endl;
  BB->AddParticles(e);
  BB -> AddParticles(e2);
  Undulator *undu =new Undulator(0.765,0.045,33);
  Lattice *field =new Lattice();
  field->addElement(undu);
  BB->Track_Vay(NOTS,dt,field,0);
  Bunch *BB2 = new Bunch(BB);
  BB2 -> MirrorY(0.0128);
  //Bunch *BB3 = new Bunch(BB);
  //BB3 -> MirrorY(-0.0128);
  //BB2 -> JoinBunch(BB3);
  //Bunch *BB4  = new Bunch(BB);
  //BB4 -> MirrorY(0.0336);
  //BB2 -> JoinBunch(BB4);
  //BB -> JoinBunch(BB2);  
  int FT=pow(2,14);
  tuple<Vector,Vector> Field[FT];
  double time_begin=1.334e-8;
  double time_end=1.336e-8;
  Vector Robs=Vector(0.0,0.0,4.0);
  double dt1 =(time_end-time_begin)/(double)FT;
  tuple<Vector,Vector> Field1[FT];
  GenGrid Grid;
  Grid.GenPlanarGrid(0.006,0.006,4,150);
  double Area=Grid.GetMeshArea();
  Vector Ar=Vector(0.0,0.0,1.0)*Area;
  tuple<Vector,Vector> Grid1[Grid.GetGridSize()];
  double lambdau = undu->GetLambdaU();
  double K = 0.934*lambdau*100*(undu->GetBPeak());
  double gamma =8.511/0.511;
  double lambdar = lambdau*(1+K*K/2.0)/(2.0*gamma*gamma);
  double freq = SpeedOfLight/lambdar;
  Grid.GenXZGrid(Vector(0,0,3.993335),0.02,0.0,0.001,2048);
  Grid.SDDSRadiationAtGrid(BB,1.335005e-8,time_end, 1, "radiation@grid.sdds");
  Out.close();  
  Out1.close();
  //int output = BB->WriteSDDSTrajectory();
  //cout<<"Write Trajectory return value: "<<output<<endl;
  //int output2= BB->WriteSDDSTime();
  //cout<<"Write Time-Trajectory return value: "<<output2<<endl;	
 // int output3 = BB->WriteSDDSRadiation(Robs, time_begin, time_end, FT);
  //cout<<"Write Radiation return value: "<<output2<<endl;
  BB->~Bunch();
 // BB2->~Bunch();
  //e->~ChargedParticle();	
  return 0;
}
