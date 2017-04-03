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
#include "wave.h"
#include "gen_grid.h"
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
	PlaneWave *Wave = new PlaneWave(0.0001,0.00001);
	BB->AddParticles(e1);
	double dt = 2.5e-15;
	int NOTS = 20000;
	Lattice *FEL=new Lattice();
	FEL->addElement(Wave);
	int numPoints=100000;
	tuple<Vector,Vector>Field[numPoints];
	Vector Obs[numPoints];
	BB->Track_Vay(NOTS, dt, FEL,1);	
	BB->WriteSDDSTrajectory();
	double time_begin = 1e-8;
	double time_end =1.008004e-8;
	double dt1 = (time_end-time_begin)/10000.0;
	double Epeak1 = Wave->EPeak();
	double Omega = Wave->Omega();
	double a = 3.0;
	double ExpectedEpeak = pow(ElementaryCharge,2)*Epeak1/(4.0*Pi*EpsNull*a*m_e*pow(SpeedOfLight,2));
	cout<<ExpectedEpeak<<endl;
	Vector Epeak = Vector(0,0,0);
	for (int i=0;i<10000;i++)
	{
		tuple<Vector,Vector>Fields=BB -> RadiationField(Vector(0,0,a),time_begin+i*dt1);
		if(Epeak.norm() < get<0>(Fields).norm())
		{
			Epeak = get<0>(Fields);
		}
	}

	cout<< "Error in computation: "<<ExpectedEpeak - Epeak.x<<endl;
	BB->WriteSDDSRadiation(Vector(0,0,a),time_begin,time_end,10000);
	BB -> WriteSDDSTime();
	Out.close();
	return 1;
}
