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
	int NOP=3;
  	int NOTS=1700;
	double dt=3.2692307692307693e-12;
	Bunch *BB= new Bunch(filename, NOP, -1, 1);
	Undulator *undu = new Undulator(0.534,0.037,54);
	Lattice *FEL = new Lattice();
	FEL -> addElement(undu);
	BB->Track_Vay(NOTS, dt, FEL);
	
	for (int i=0; i<1700;i++)
	{
		Vector XP = BB->b[0].TrajPoint(i);
		Vector XP2 =BB->b[0].TrajMomentum(i);
		Vector XP3 =BB->b[0].TrajAccel(i);
		double t =BB->b[0].TrajTime(i);
		Out<<t<<"\t"<<XP.x<<"\t"<<XP.y<<"\t"<<XP.z<<"\t"<<XP2.x<<"\t"<<XP2.y<<"\t"<<XP2.z<<"\t"<<XP3.x<<"\t"<<XP3.y<<"\t"<<XP3.z<<"\n";
	}

	Out.close();
	return 1;
}
