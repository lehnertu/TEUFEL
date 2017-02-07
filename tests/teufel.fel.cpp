/************************************************************************/
/*                                                                      */
/*  TEUFEL - THz Emission from Undulators and Free-Electron Lasers      */
/*                                                                      */
/*  written by  U.Lehnert                               12/2016         */
/*                                                                      */
/************************************************************************/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <QtGui>
#include <QtXml/QDomDocument>

#include <boost/math/special_functions/erf.hpp>

#include "global.h"
#include "particle.h"
#include "wave.h"

// <bunch> parameters
int NOP;                        // number of (real) particles
int NOTS;                       // number of time steps
double dt;                      // time step
double x0;                      // start x position of the particles
double SIGMA_Z;                 // width of the gaussian distribution for z
double px0;                     // start x momentum of the particles
double pz0;                     // start z momentum of the particles

// <undulator> parameters
double B_PEAK;                  // peak magnetic field
double LAMBDA_U;                // undulator period length
int NPERIODS;                   // number of undulator periods

// <wave> parameters
double INTENSITY;               // incident plane wave intensity [W/cm²]
double LAMBDA;                  // wavelength

// <waveguide> parameters
int NOM;                        // number of mirror particles per side for waveguide simulation
double WG_GAP;                  // the full gap height of the waveguide

void parse_input()
{
  QFile file("./wgfel.in");
  QDomDocument doc("WGFEL");
  QString s;
  if( !file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      printf("error opening wgfel.in\n");
      exit(-1);
    }
  if( !doc.setContent( &file ) )
    {
      printf("error reading (XML) wgfel.in\n");
      exit(-1);
    };
  file.close();
  // interpret data
  QDomElement root = doc.documentElement();
  if( root.tagName() != "wgfel" )
    {
      printf("wrong XML file type - didn't find root <wgfwl>\n");
      exit(-1);
    };
  // read <bunch> parameters
  QDomElement bunch = root.firstChildElement ("bunch");
  if (bunch.tagName() != "bunch")
    {
      printf("didn't find node <bunch>\n");
      exit(-1);
    };
  s = bunch.attribute("nop","1");
  NOP = s.toInt();
  s = bunch.attribute("timesteps","10");
  NOTS = s.toInt();
  s = bunch.attribute("dt","1e-12");
  dt = s.toDouble();
  s = bunch.attribute("x0","0.0");
  x0 = s.toDouble();
  s = bunch.attribute("sigz","0.001");
  SIGMA_Z = s.toDouble();
  s = bunch.attribute("px0","0.0");
  px0 = s.toDouble();
  s = bunch.attribute("pz0","10.0");
  pz0 = s.toDouble();
  // read <undulator> parameters
  QDomElement undu = root.firstChildElement ("undulator");
  if (undu.tagName() != "undulator")
    {
      printf("didn't find node <undulator>\n");
      exit(-1);
    };
  s = undu.attribute("bpeak","0.0");
  B_PEAK = s.toDouble();
  s = undu.attribute("period","0.1");
  LAMBDA_U = s.toDouble();
  s = undu.attribute("nop","10");
  NPERIODS = s.toInt();
  // read <wave> parameters
  QDomElement wv = root.firstChildElement ("wave");
  if (wv.tagName() != "wave")
    {
      printf("didn't find node <wave>\n");
      exit(-1);
    };
  s = wv.attribute("wattsqcm","1e18");
  INTENSITY = s.toDouble();
  s = wv.attribute("lambda","1e-6");
  LAMBDA = s.toDouble();
  // read <waveguide> parameters
  QDomElement wg = root.firstChildElement ("waveguide");
  if (wg.tagName() != "waveguide")
    {
      printf("didn't find node <waveguide>\n");
      exit(-1);
    };
  s = wg.attribute("nom","0");
  NOM = s.toInt();
  s = wg.attribute("h","0.001");
  WG_GAP = s.toDouble();
}

int main ()
{
  FILE *dump;
  int d;                // integer buffer used for file writing
  Vector X;

  printf("\nWGFEL Version 0.03 U.Lehnert 10/2014\n\n");
  parse_input();

  Lattice *FEL = new Lattice;

  // define the undulator
  Undulator *U1 = new Undulator(B_PEAK, LAMBDA_U, NPERIODS);
  // Undulator *U1 = new Undulator(0.37, 0.11, 40);
  FEL->addElement(U1);
  printf("WGFEL - undulator defined:\n");
  printf("peak field %9.3f T\n",U1->GetBPeak());
  printf("period     %9.6f m\n",U1->GetLambdaU());
  printf("%d periods\n",U1->GetNPeriods());
  printf("K(rms)     %9.3f\n\n",U1->GetKrms());

  // define plane wave
  PlaneWave *w1 = new PlaneWave(LAMBDA, INTENSITY);
  FEL->addElement(w1);
  printf("WGFEL - plane wave defined:\n");
  printf("omega =  %9.6g Hz\n",w1->Omega());
  printf("electric peak field = %9.6g V/m\n",w1->EPeak());
  printf("magnetic peak field = %9.6g T\n",w1->BPeak());
  printf("intensity a0 = %9.6g\n",w1->A0());
  printf("\n");

  ChargedParticle *bunch = new ChargedParticle[NOP];

  for (int i=0; i<NOP; i++)
    {
      // initial position of the particle
      // gaussian distribution along the z axis
      double z0 = SIGMA_Z*boost::math::erf_inv(1.96*((double)i+0.5)/NOP - 0.98);
      Vector X0 = Vector(x0, 0.0, z0);
      // initial momentum of the particle
      Vector P0 = Vector(px0, 0.0, pz0);
      // bunch[i].TrackEuler(NOTS, dt, X0, P0, w1);
      bunch[i].TrackVay(NOTS, dt, X0, P0, FEL);
    }

  // dump the particle trajectories into a binary file
  dump = fopen("particles.bin","wb");
  // first write the number of particles
  d = NOP;
  fwrite(&d, sizeof(int), 1, dump);
  // write the number of time steps
  d = NOTS+1;
  fwrite(&d, sizeof(int), 1, dump);
  // write the traces
  // it is transposed while writing because z (the 1-step index)
  // ist used as the outer loop index
  for (int ip=0; ip<NOP; ip++)
    for (int it=0; it<NOTS+1; it++)
      {
        Vector XP = bunch[ip].TrajPoint(it);
        // fwrite expects a pointer to the value to be written
        fwrite(&(XP.x), sizeof(double), 1, dump);
        fwrite(&(XP.y), sizeof(double), 1, dump);
        fwrite(&(XP.z), sizeof(double), 1, dump);
         // write the momentum in addition
        XP = bunch[ip].TrajMomentum(it);
        fwrite(&(XP.x), sizeof(double), 1, dump);
        fwrite(&(XP.y), sizeof(double), 1, dump);
        fwrite(&(XP.z), sizeof(double), 1, dump);
     };
  fclose(dump);

  return(0);
}
