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
#include <vector.h>
using namespace std;
int main()
{
  const char *filename= "TEUFEL-EXAMPLE.txt";
  //const char *filename= "BeamProfile.txt";
  Bunch *BB;//=new Bunch();
  int NOP=2000;
  int NOTS=2200/4;
  double dt=1.00e-11;
  double QTot  = 1.0*15e-12;
  double QperMacro=QTot/(NOP*ElementaryCharge);
  /*for(int i=0;i<1;i++)
  {
	BB ->AddParticles(new ChargedParticle(-1,1,Vector(0,0,-i*100.0e-15*3.0e8/100.0),Vector(0,0.01,46.96673189823875),0));
  }*/
  //ChargedParticle *e =  new ChargedParticle(-4.37e8,4.37e8,Vector(0,0,-0.0001),Vector(0.1,0.1,48),0);
 //ChargedParticle *e =  new ChargedParticle(-4.37e7,4.37e7,Vector(0,0,-0.0001),Vector(0.002,0.002,8.511/0.511),0);
 //ChargedParticle *e1 =  new ChargedParticle(-4.37e7,4.37e7,Vector(0,0,-0.0001-SpeedOfLight/1.212e+12),Vector(0.002,0.002,8.511/0.511),0);
  //BB -> AddParticles(e);
  //Bunch *B2  = new Bunch();
  //BB -> AddParticles(e1);
 // BB -> AddParticles(e1);
  /*ChargedParticle *e1 =  new ChargedParticle(-1,1,Vector(0,0,-0.000108305),Vector(0,0.0021,46.96673189823875),-0.000108305/SpeedOfLight);
  ChargedParticle *e2 =  new ChargedParticle(-1,1,Vector(0,0,-2*0.000108305),Vector(0,0.0022,46.96673189823875),-0.000108305/SpeedOfLight);
  ChargedParticle *e3 =  new ChargedParticle(-1,1,Vector(0,0,-3*0.000108305),Vector(0,0.002,46.96673189823875),-0.000108305/SpeedOfLight);
  
  BB -> AddParticles(e1);
  BB -> AddParticles(e2);
  BB -> AddParticles(e3);*/
  cout<<QperMacro<<endl;
  BB= new Bunch(filename, NOP, -(int)QperMacro,(int)QperMacro);
  //Bunch *BB2 = new Bunch(BB);
  //BB2 -> JoinBunch(BB);
  //Undulator *undu =new Undulator(0.1,0.3,8);
  Undulator *undu =new Undulator(0.45,0.050,30);
  Lattice *field =new Lattice();
  field->addElement(undu);
  BB->Track_Vay(NOTS,dt,field,1);
  vector<ChargedParticle*>Mirror_e1;
  vector<ChargedParticle*>Mirror_e2;
  //B2 -> AddParticles(e1);
/* for(int i=0;i<40;i++)
  {
	Mirror_e1.push_back(new ChargedParticle(e));
	Mirror_e2.push_back(new ChargedParticle(e));
	Mirror_e1[i] -> MirrorY(i*0.010);
	Mirror_e2[i] ->MirrorY(-i*0.010);
	
  }
  for(int i=0;i<40;i++)
  {
	BB -> AddParticles(Mirror_e1[i]);
	BB -> AddParticles(Mirror_e2[i]);
  }*/
  /*Bunch *Copy_Bunch1 = new Bunch(BB);
  Bunch *Copy_Bunch2 = new Bunch(BB);
  Bunch *Copy_Bunch3 = new Bunch(BB);
  Bunch *Copy_Bunch4 = new Bunch(BB);
  Bunch *Copy_Bunch5 = new Bunch(BB);
  Bunch *Copy_Bunch6 = new Bunch(BB);
  Bunch *Copy_Bunch7 = new Bunch(BB);
  Bunch *Copy_Bunch8 = new Bunch(BB);
  Copy_Bunch1 -> MirrorY(0.05);
  Copy_Bunch2 -> MirrorY(2*0.05);
  Copy_Bunch3 -> MirrorY(3*0.05);
  Copy_Bunch4 -> MirrorY(4*0.05);
  Copy_Bunch5 -> MirrorY(-0.05);
  Copy_Bunch6 -> MirrorY(-2*0.05);
  Copy_Bunch7 -> MirrorY(-3*0.05);
  Copy_Bunch8 -> MirrorY(-4*0.05);
  BB -> JoinBunch(Copy_Bunch1);
  BB -> JoinBunch(Copy_Bunch2);
  BB -> JoinBunch(Copy_Bunch3);
  BB -> JoinBunch(Copy_Bunch4);
  BB -> JoinBunch(Copy_Bunch5);
  BB -> JoinBunch(Copy_Bunch6);
  BB -> JoinBunch(Copy_Bunch7);
  BB -> JoinBunch(Copy_Bunch8);*/
  double lambdau = undu->GetLambdaU();
  double K = 0.934*lambdau*100.0*(undu->GetBPeak());
  double gamma =BB->getTrajMomentum(0,0).norm();
  double lambdar = lambdau*(1+K*K/2.0)/(2.0*gamma*gamma);
  double freq = SpeedOfLight/lambdar;
  cout<<freq<<endl;
  int FT=pow(2,10);
  pair<Vector,Vector> Field[FT];
  double time_begin=9.00567e-9*7.1/2.7;
  double time_end=9.02e-9*7.1/2.7;
  Vector Robs=Vector(0.0,0.0,7.1);
  double dt1 =(time_end-time_begin)/(double)FT;
  //GenGrid Grid;
  GenGrid Grid2;
  Grid2.GenPlanarGrid(0.06/4.0,0.07/4.0,7.1,32); 
  //Grid.GenXZGrid(Vector(0,0,2.7),0.035,0.0,2.0*27.3e-4,4*256);
  //Grid.SDDSRadiationAtGrid(BB,time_begin,time_end,20, "radiation@grid.sdds"); //1
  //Grid2.SDDSRadiationAtGrid(BB,time_begin,time_end,20, "radiation@xygrid.sdds"); //2
  //cout<<SpeedOfLight/freq<<endl;
  //Grid2.SDDSFilteredRadiationAtGrid(freq,"FilteredRadiation.sdds", BB, time_begin, time_end,FT ); //3
  BB -> WriteSDDSTrajectory(); //4
  int output2= BB->WriteSDDSTime(); //5
  cout<<"Write Time-Trajectory return value: "<<output2<<endl;	
  int output3 = BB->WriteSDDSRadiation(Robs, time_begin, time_end, FT); //6
  cout<<"Write Radiation return value: "<<output2<<endl;
  BB->~Bunch();
  return 0;
}
