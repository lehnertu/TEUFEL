#include "four_vector.h"
#include "four_vector.cpp"
#include <iostream>
#include <stdlib.h>
#include "global.h"
#include "vector.h"
using namespace std;

int main()
{
   Four_Vector Tree = Four_Vector(10.0e-6,0.0,0.0,0.0);
   Four_Vector Pole = Four_Vector(10.0e-6,20000.0,0.0,0.0);
   double beta = 0.5;
   double gamma = sqrt(1.0/(1.0-0.5*0.5));
   cout <<gamma*beta<<endl;
   double gb=gamma*beta;
   Vector S1 = Vector(gb,0.0,0.0);
   Four_Vector Y = Tree.LorentzTransform(S1);
   Four_Vector Z = Y.InverseLorentzTransform(S1);
   cout<<Y.t<<"\t"<<Y.x<<"\t"<<Y.y<<"\t"<<Y.z<<"\n"<<endl;
   cout<<Z.t<<"\t"<<Z.x<<"\t"<<Z.y<<"\t"<<Z.z<<"\n"<<endl;
   return 1;



}
