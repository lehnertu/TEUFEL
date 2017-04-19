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
  const char *filename= "TEUFEL-EXAMPLE.txt";
  //const char *filename= "BeamProfile.txt";
  Bunch *BB;
  int NOP=999;
  int NOTS=1700;
  double dt=3.2692307692307693e-12;
  //BB=new Bunch();
  BB= new Bunch(filename, NOP, -12500, 12500);
  //BB2 -> Translate(Vector(0,0,-0.00053));
  /*BB -> JoinBunch(BB2);
  ChargedParticle *e = new ChargedParticle(-12500,12500,Vector(0,0,0),Vector(0,0,16),0);
  ChargedParticle *e2 = new ChargedParticle(-12500,12500,Vector(0,0,-0.000553),Vector(0,0,16),0);
  ChargedParticle *e3 = new ChargedParticle(-12500,12500,Vector(0,0,-0.001),Vector(0,0,16),0);
  ChargedParticle *e4 = new ChargedParticle(-12500,12500,Vector(0,0,-0.001553),Vector(0,0,16),0);
  //cout<<"Now trying to add particles to a bunch"<<endl;
  BB->AddParticles(e);
  BB->AddParticles(e3);
  BB2 -> AddParticles(e2);
  BB2 -> AddParticles(e4);
  BB ->AddParticles(BB2->PointerToParticle(0));
  BB -> AddParticles(BB2->PointerToParticle(1));
  //BB ->JoinBunch(BB);*/
  Undulator *undu =new Undulator(0.1,0.045,33);
  double lambdau = undu->GetLambdaU();
  double K = 0.934*lambdau*100*(undu->GetBPeak());
  double gamma =8.511/0.511;
  double lambdar = lambdau*(1+K*K/2.0)/(2.0*gamma*gamma);
  double freq = SpeedOfLight/lambdar;
  Lattice *field =new Lattice();
  field->addElement(undu);
  Bunch *BB2 = new Bunch(BB);
  Bunch *BB3 = new Bunch(BB);
  Bunch *BB4 = new Bunch(BB);
  BB2 ->Translate(Vector(0,0,-SpeedOfLight/freq));
  BB3 ->Translate(Vector(0,0,-2.0*SpeedOfLight/freq));
  BB4 ->Translate(Vector(0,0,-3.0*SpeedOfLight/freq));
  BB -> JoinBunch(BB2);
  BB -> JoinBunch(BB3);
  BB -> JoinBunch(BB4);
  /*for (int i=0; i<999;i++)
  {
	
	BB -> AddParticles(BB3->PointerToParticle(i));
  }

  for (int i=0; i<999;i++)
  {
	
	BB -> AddParticles(BB2->PointerToParticle(i));
  }*/
  BB->Track_Vay(NOTS,dt,field,0);

  //int output = BB->WriteSDDSTrajectory();
  //cout<<"Write Trajectory return value: "<<output<<endl;
  int FT=pow(2,14);
  tuple<Vector,Vector> Field[FT];
  double time_begin=1.334e-8/2.0;
  double time_end=1.336e-8/2.0;
  Vector Robs=Vector(0.0,0.0,4.0/2.0);
  double dt1 =(time_end-time_begin)/(double)FT;
  GenGrid Grid;  
  Grid.GenXZGrid(Vector(0,0,4.0/2.0),0.1,0.0,80*SpeedOfLight/freq,512);
  Grid.SDDSRadiationAtGrid(BB,time_begin,time_end,5, "radiation@grid.sdds");
  //cout<<freq<<endl;
  //an.BunchFactor(lambdau,lambdar);
  //ChargedParticle *e2 = new ChargedParticle(e);
  //BB->JoinBunch(BB2);
  
  //int output2= BB->WriteSDDSTime();
  //cout<<"Write Time-Trajectory return value: "<<output2<<endl;	
  //int output3 = BB->WriteSDDSRadiation(Robs, time_begin, time_end, FT);
  //cout<<"Write Radiation return value: "<<output2<<endl;
  BB->~Bunch();
 // BB2->~Bunch();
  //e->~ChargedParticle();	
  return 0;
}
