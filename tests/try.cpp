#include "particle.h"
#include "undulator.h"
#include <iostream>
#include <stdlib.h>
#include "global.h"
#include <vector>
#include <random>
#include <fstream>
#include "externalfield.h"
#include "vector.h"
#include "bunch.h"
using namespace std;

int main()
{
	ofstream Out;
	Out.open("StepVay_Trial.txt");
	const char *filename= "BeamProfile.txt";
	int NOP=10;
  	int NOTS=1700;
	double dt=3.2692307692307693e-12;
	Bunch *BB= new Bunch(filename, NOP, -1, 1);
	Undulator *undu = new Undulator(0.534,0.037,54);
	Lattice *FEL = new Lattice();
	FEL -> addElement(undu);
	BB->Track_Vay(NOTS, dt, FEL);
	ChargedParticle *electron = new ChargedParticle(-1,1,Vector(0,0,0),Vector(0,0,16),0);
	electron->TrackVay(NOTS,dt,Vector(0,0,0),Vector(0,0,16),FEL);
	ChargedParticle *copy_electron = new ChargedParticle(electron);
	copy_electron->MirrorY(0.012);
	for (int i=0; i<1700;i++)
	{
		Vector XP = electron->TrajPoint(i);
		Vector XP2 =electron->TrajMomentum(i);
		Vector XP3 =electron->TrajAccel(i);
		double t =copy_electron->TrajTime(i);
		Out<<t<<"\t"<<XP.x<<"\t"<<XP.y<<"\t"<<XP.z<<"\t"<<XP2.x<<"\t"<<XP2.y<<"\t"<<XP2.z<<"\t"<<XP3.x<<"\t"<<XP3.y<<"\t"<<XP3.z<<"\n";
	}

	Out.close();
	return 1;
}
