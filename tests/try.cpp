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
	Bunch *BB= new Bunch();
	ChargedParticle *e1 = new ChargedParticle(-1,1,Vector(0,0,0.05),Vector(0,0,0),0);
	ChargedParticle *e2 = new ChargedParticle(-1,1,Vector(0,0,-0.05),Vector(0,0,0),0);
	BB->AddParticles(e1);
	BB->AddParticles(e2);
	double dt = 2e-8;
	int NOTS = 100000;
	Lattice *FEL=new Lattice();
	Undulator *undu = new Undulator(0,0,0);
	FEL->addElement(undu);
	cout<<BB->getNOP()<<endl;
	BB->Track_Vay(NOTS, dt, FEL,1);
	BB->WriteSDDSTrajectory();
	return 1;
}
