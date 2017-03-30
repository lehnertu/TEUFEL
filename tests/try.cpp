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
	double Energy = 8; //8Mev
	double gamma = 8.511/0.511;
	double gammabeta=sqrt(pow(gamma,2)-1);
	ChargedParticle *e1 = new ChargedParticle(-1,1,Vector(0,0,0.00),Vector(0,0,0),0);
	ChargedParticle *e2 = new ChargedParticle(-1,1,Vector(0,0,0.001),Vector(0,0,0),0);
	BB->AddParticles(e1);
	BB->AddParticles(e2);
	double dt = 2.5e-9;
	int NOTS = 20000;
	Lattice *FEL=new Lattice();
	Undulator *undu = new Undulator(0.0,0.045,10);
	FEL->addElement(undu);
	cout<<BB->getNOP()<<endl;
	BB->Track_Vay(NOTS, dt, FEL,1);
	BB->WriteSDDSTrajectory();
	double time_begin = 2.001e-9;
	double time_end = 2.002e-9;
	BB->WriteSDDSRadiation(Vector(0,0,0.6),time_begin,time_end,100000);
	return 1;
}
