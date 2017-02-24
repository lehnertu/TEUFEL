#include "bunch.h"
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <iostream>

using namespace std;
int main()
{
  std::ofstream Out ("test.txt", std::ofstream::out);
  const char *filename= "Bunch.txt";
  Bunch *BB= new Bunch();
  int NOP=2;
  int NOTS=6000;
  double dt=1e-9;
  //BB= new Bunch(filename, NOP, -1, 1);
  Undulator *undu =new Undulator(0,0.037,54);
  Lattice *field =new Lattice();
  field->addElement(undu);
  BB->Track_Euler(NOTS,dt,field);
  for(int k=0;k<NOP;k++)
  {
  	for (int i=0;i<NOTS;i++)
  	{
		Vector XP=BB->b[k].TrajPoint(i);
		Out<<XP.x<<"\t"<<XP.y<<"\t"<<XP.z<<"\n";
  	}
  }
  cout<<(BB->b[0].TrajMomentum(NOTS-1)).norm()*9.1e-31*3e8*3e8<<endl;
  cout<<(BB->b[0].TrajPoint(NOTS-1)).norm()<<endl;
  Out.close();
  return 0;
}
