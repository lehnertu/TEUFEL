#include "four_vector.h"
#include "four_vector.cpp"
#include <iostream>
#include <stdlib.h>
#include "global.h"
#include <vector>

using namespace std;

int main()
{
   Four_Vector Tree = Four_Vector(10.0e-6,0.0,0.0,0.0);
   Four_Vector Pole = Four_Vector(10.0e-6,20000.0,0.0,0.0);
   double beta = 0.5;
   double gamma = sqrt(1.0/(1.0-0.5*0.5));
   double gb=gamma*beta;
   Vector S1 = Vector(gb,0.0,0.0);
   Four_Vector Y = Tree.LorentzTransform(S1);
   Four_Vector Z = Pole.LorentzTransform(S1);
   //cout<<Y.t<<"\t"<<Y.x<<"\t"<<Y.y<<"\t"<<Y.z<<"\n"<<endl;
   //cout<<Z.t<<"\t"<<Z.x<<"\t"<<Z.y<<"\t"<<Z.z<<"\n"<<endl;
   //Four_Vector A = Four_Vector.Four_Momentum(30.35,Vector(0.0,0.0,30.35));
   Vector E = Vector(1.0,1.0,0.0);
   Vector B = Vector(0.0,0.0,0.0);
   tuple<Vector,Vector> EE = EMTransform(S1,E,B);
   EE = EMInverseTransform(S1,E,B);
   E  = get<0>(EE);
   B = get<1>(EE);
   cout<<E.x<<"\t"<<E.y<<"\t"<<E.z<<"\n"<<endl;
   cout<<B.x<<"\t"<<B.y<<"\t"<<B.z<<"\n"<<endl;
   return 1;
}
