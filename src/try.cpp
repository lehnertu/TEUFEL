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

   //test case1 : muon
   // lifetime in its frame: 26ns
   // its beta wrt earth 0.95
   Four_Vector Tree = Four_Vector(26.0e-9,0.0,0.0,0.0);
   Four_Vector Pole = Four_Vector(10.0e-6,20000.0,110.0,1130.0);
   double betax = 0.95;
   double betay = 0.0;
   double betaz = 0.0;
   double beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
   double gamma = sqrt(1.0/(1.0-beta*beta));
   double gbx=gamma*betax;
   double gby=gamma*betay;
   double gbz=gamma*betaz;
   Vector S1 = Vector(gbx,gby,gbz);
   cout<<"************************************************************\n\n";
   cout<<"\t Test Case of Muon's Lifetime and Distace Travelled before decaying\n\n";
   // Use inverse transform to see its lifetime in earth's frame
   Four_Vector Y = Tree.InverseLorentzTransform(S1);
   //Four_Vector Z = Y.InverseLorentzTransform(S1);
   cout<<"Lifetime of muon in its frame: 26ns, Speed = 0.95c\n";
   cout<<"Four Vectors for muon in earth's frame:(lifetime, distance travelled)\n ";
   cout<<Y.t<<"\t"<<Y.x<<"\t"<<Y.y<<"\t"<<Y.z<<"\n"<<endl;
   cout<<"************************************************************\n\n";




   cout<<"************************************************************\n\n";
   cout<<"\t Time Dilation for a Slow Moving Object (Moving Clocks Run Slow)\n\n";
   //cout<<Z.t<<"\t"<<Z.x<<"\t"<<Z.y<<"\t"<<Z.z<<"\n"<<endl;
   Tree = Four_Vector(3.6e3,0.0,0.0,0.0);
   betax = 1.3333333333333334e-06;
   betay = 0.0;
   betaz = 0.0;
   beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
   gamma = sqrt(1.0/(1.0-beta*beta));
   gbx=gamma*betax;
   gby=gamma*betay;
   gbz=gamma*betaz;
   S1 = Vector(gbx,gby,gbz);
   Y = Tree.InverseLorentzTransform(S1);
   cout<<"Event time interval in moving frame: "<<Tree.t<<endl;
   cout<<"Event time interval in earth's frame: "<<Y.t<<endl;
   cout<<"Difference in Time interval in two frames: "<<Y.t-Tree.t<<endl;
   cout<<"\n************************************************************\n\n";




   cout<<"************************************************************\n\n";
   cout<<"\t Length Contraction and Rotation\n\n";
   //cout<<Z.t<<"\t"<<Z.x<<"\t"<<Z.y<<"\t"<<Z.z<<"\n"<<endl;
   double L0 = 100; //length of rod in its own frame
   double theta = Pi/4.0; //angle rod makes with xprime(in its frame)
   Tree = Four_Vector(0.0,L0*cos(theta),L0*sin(theta),0.0);
   betax = 0.5;  //rod moves with veloctiy 0.5c
   betay = 0.0;
   betaz = 0.0;
   beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
   gamma = sqrt(1.0/(1.0-beta*beta));
   gbx=gamma*betax;
   gby=gamma*betay;
   gbz=gamma*betaz;
   S1 = Vector(gbx,gby,gbz);
   Y = Tree.InverseLorentzTransform(S1);
   cout<<"Length of rod in its frame: "<<L0<<endl;
   cout<<"Length of Rod in earth's frame: "<<sqrt(Y.x*Y.x+Y.y+Y.y+Y.z+Y.z)<<endl;
   cout<<"Angle made by Rod in earth's frame with x: "<<atan(Y.y/Y.x)*180.0/Pi<<endl;




   cout<<"\n************************************************************\n\n";
   cout<<"\t Lorentz Velocity Transformation\n\n";
   // four velocity is just a special case of four vector
   // its components are d_dtau(ct,Vector R)
   // thus it becomes FourVector(gamma,gamma*Vector V)
   // 
   // let a particle have speed 0.5c
   // let another particle have speed 0.8c
   // we see first particle from second particle
   // and find the velocities
   // gamma of particle = 1.1547
   Tree = Four_Vector(1.1547,0.5*1.1547*SpeedOfLight,0.0,0.0);
   betax = 0.8;  //rod moves with veloctiy 0.5c
   betay = 0.0;
   betaz = 0.0;
   beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
   gamma = sqrt(1.0/(1.0-beta*beta));
   gbx=gamma*betax;
   gby=gamma*betay;
   gbz=gamma*betaz;
   S1 = Vector(gbx,gby,gbz);
   Y = Tree.LorentzTransform(S1);
   // find the new gammabeta
   double gammabeta_newframe = sqrt(pow(Y.x/SpeedOfLight,2)+pow(Y.y/SpeedOfLight,2)+pow(Y.z/SpeedOfLight,2));
   double gamma2 = sqrt(1+gammabeta_newframe*gammabeta_newframe);
   cout<<"Speed of particle in earth's frame: "<<Tree.x/(gamma2*SpeedOfLight)<<endl;
   cout<<"Speed of particle in a particle's frame moving with (0.8c): "<<Y.x/(1.1547*SpeedOfLight)<<endl;



   cout<<"\n************************************************************\n\n";
   cout<<"\t Lorentz Velocity Transformations for Two Components\n\n";
   // here we see the transformation of veloctiy of a spacecraft
   // from earth spaceship A has velocities (gamma1,0,gamma1*beta1*SpeedOfLight,0.0)
   // from earth spaceship B has velocities (gamma,gamma*beta*SpeeOfLight,0.0,0.0)
   // we compute the velocity as seen by B (of A)
   double gamma1 = sqrt(1.0/(1.0-0.9*0.9));
   Tree = Four_Vector(gamma1,0.0,-0.9*gamma1*SpeedOfLight,0.0);
   betax = 0.9;  //rod moves with veloctiy 0.5c
   betay = 0.0;
   betaz = 0.0;
   beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
   gamma = sqrt(1.0/(1.0-beta*beta));
   gbx=gamma*betax;
   gby=gamma*betay;
   gbz=gamma*betaz;
   S1 = Vector(gbx,gby,gbz);
   Y = Tree.LorentzTransform(S1);
   // the transformed quantity is a four vector(gamma,Vector of gamma*beta*c)
   // so we need to find the new gammabeta
   gammabeta_newframe = sqrt(pow(Y.x/SpeedOfLight,2)+pow(Y.y/SpeedOfLight,2)+pow(Y.z/SpeedOfLight,2));
   gamma2 = sqrt(1+gammabeta_newframe*gammabeta_newframe);
   cout<<"Velocity of ship in earth's frame: "<<Tree.x/(gamma1*SpeedOfLight)<<" "<<Tree.y/(gamma1*SpeedOfLight)<<" "<<Tree.z/(gamma1*SpeedOfLight)<<endl;
   cout<<"Velocity of ship in a ships's frame moving with 0.9c: "<<Y.x/(gamma2*SpeedOfLight)<<" "<<Y.y/(gamma2*SpeedOfLight)<<" "<<Y.z/(gamma2*SpeedOfLight)<<endl;





   cout<<"\n************************************************************\n\n";
   cout<<"\t Lorentz Transformation of Earth's Magnetic Fields\n\n";

   // here we transfom earth's magnetic i.e 30microT
   // assume it to be z direction
   // transform it to a frame moving with speed of airplane
   // speed ~ 900Km/hr
   // beta = 8.3e-7
   Vector E = Vector(0.0,0.0,0.0);
   Vector B = Vector(0.0,0.0,30.0e-6);
   betax = 8.3e-7;  //rod moves with veloctiy 0.5c
   betay = 0.0;
   betaz = 0.0;
   beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
   gamma = sqrt(1.0/(1.0-beta*beta));
   gbx=gamma*betax;
   gby=gamma*betay;
   gbz=gamma*betaz;
   S1 = Vector(gbx,gby,gbz);
   tuple<Vector,Vector> EE = EMTransform(S1,E,B);
   Vector E1  = get<0>(EE);
   Vector B1 = get<1>(EE);
   cout<<"Electric Field observed: "<<E1.x<<" "<<E1.y<<" "<<E1.z<<"\n"<<endl;
   cout<<"Magnetic Field observed: "<<B1.x<<"\t"<<B1.y<<"\t"<<B1.z<<"\n"<<endl;


   cout<<"\n************************************************************\n\n";
   cout<<"\t Lorentz Transformation of Dipole Fields\n\n";

   // here we transfom dipoles magnetic i.e 1T
   // assume it to be z direction
   // transform it to a frame moving with speed of 100keV Hydrogen atom
   // speed ~ 4.38e6
   // beta = 0.0146
   E = Vector(0.0,0.0,0.0);
   B = Vector(0.0,0.0,1);
   betax = 0.0146;  //rod moves with veloctiy 0.5c
   betay = 0.0;
   betaz = 0.0;
   beta  = sqrt(betax*betax+betay*betay+betaz*betaz);
   gamma = sqrt(1.0/(1.0-beta*beta));
   gbx=gamma*betax;
   gby=gamma*betay;
   gbz=gamma*betaz;
   S1 = Vector(gbx,gby,gbz);
   EE = EMTransform(S1,E,B);
   E1  = get<0>(EE);
   B1 = get<1>(EE);
   cout<<"Electric Field observed in particle's frame: "<<E1.x<<" "<<E1.y<<" "<<E1.z<<"\n"<<endl;
   cout<<"Magnetic Field observed in particle's frame: "<<B1.x<<"\t"<<B1.y<<"\t"<<B1.z<<"\n"<<endl;


   cout<<"\n************************************************************\n\n";
   cout<<"\t Inverse Lorentz Transformation of Dipole Fields\n\n";

   // here we inverse transfom dipoles magnetic i.e 1T from above example
   
   tuple<Vector,Vector> EE1= EMInverseTransform(S1,E1,B1);
   E1  = get<0>(EE1);
   B1 = get<1>(EE1);
   cout<<"Electric Field in lab frame: "<<E1.x<<"\t"<<E1.y<<"\t"<<E1.z<<"\n"<<endl;
   cout<<"Magnetic Field in lab frame: "<<B1.x<<"\t"<<B1.y<<"\t"<<B1.z<<"\n"<<endl;

   cout<<"\n************************************************************\n\n";
   cout<<"\t Matrix Multiply [4,4]&[4,1]\n\n";
   vector<vector<double>>L(4,vector<double>(4));
   vector<vector<double>>L1(4,vector<double>(1));
   cout<<"First Matrix is: \n\n";
   for(int i = 0;i<4;i++)
    	{
		for (int k=0;k<4;k++)
		{
			L[i][k] = rand()%10;
			cout<<L[i][k]<<"\t";
		}
		cout<<"\n";

	}

    cout<<"Second Matrix is: \n\n";
    for(int i = 0;i<4;i++)
    	{
		for (int k=0;k<1;k++)
		{
			L1[i][k] = rand()%10;
			cout<<L1[i][k]<<"\t";
		}
		cout<<"\n";

	}
   
   vector<vector<double>>L2(4,vector<double>(1));
   L2=Multiply(L,L1);
   //cout<<L2.size()<<"\t"<<L2[0].size()<<"\n";

   cout<<"Result Matrix is: \n\n";
   Four_Vector D = Four_Vector(L2[0][0],L2[1][0],L2[2][0],L2[3][0]);
   cout<<D.t<<"\t"<<D.x<<"\t"<<D.y<<"\t"<<D.z<<"\n"<<endl;
   //delete L2;
   return 1;
}
