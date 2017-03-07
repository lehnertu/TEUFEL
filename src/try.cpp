#include "four_vector.h"
#include "four_vector.cpp"
#include <iostream>
#include <stdlib.h>
#include "global.h"
#include <vector>
#include <random>
using namespace std;

int main()
{
   Four_Vector Tree = Four_Vector(10.0e-6,0.0,0.0,0.0);
   Four_Vector Pole = Four_Vector(10.0e-6,20000.0,110.0,1130.0);
   double betax = 0.5;
   double betay = 0.04;
   double betaz = 0.4;
   double beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
   double gamma = sqrt(1.0/(1.0-beta*beta));
   double gbx=gamma*betax;
   double gby=gamma*betay;
   double gbz=gamma*betaz;
   Vector S1 = Vector(gbx,gby,gbz);
   Four_Vector Y = Pole.LorentzTransform(S1);
   Four_Vector Z = Y.InverseLorentzTransform(S1);
   //cout<<Y.t<<"\t"<<Y.x<<"\t"<<Y.y<<"\t"<<Y.z<<"\n"<<endl;
   //cout<<Z.t<<"\t"<<Z.x<<"\t"<<Z.y<<"\t"<<Z.z<<"\n"<<endl;
   cout<<gamma<<endl;
   Vector E = Vector(1.0,1.7,1.90);
   Vector B = Vector(1.5,1.1,1.4);
   tuple<Vector,Vector> EE = EMTransform(S1,E,B);
   Vector E1  = get<0>(EE);
   Vector B1 = get<1>(EE);
   tuple<Vector,Vector> EE1= EMInverseTransform(S1,E1,B1);
   E1  = get<0>(EE1);
   B1 = get<1>(EE1);
   cout<<E1.x<<"\t"<<E1.y<<"\t"<<E1.z<<"\n"<<endl;
   cout<<B1.x<<"\t"<<B1.y<<"\t"<<B1.z<<"\n"<<endl;
   
   return 1;
}
