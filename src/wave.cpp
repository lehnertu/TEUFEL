/************************************************************************/
/*                                                                      */
/*  TEUFEL - THz Emission from Undulators and Free-Electron Lasers      */
/*                                                                      */
/*  written by  U.Lehnert                               12/2016         */
/*                                                                      */
/************************************************************************/

#include "wave.h"
#include "global.h"

#include <math.h>

PlaneWave::PlaneWave( double lambda,
                      double intensity )
  {
    Lambda = lambda;
    omega = 2.0*Pi*SpeedOfLight/Lambda;
    Intensity = 1.0e4 * intensity;      // scale from W/cm² to W/m²
    Epeak = sqrt( 2.0*Intensity / (SpeedOfLight*EpsNull) );
    Bpeak = Epeak / SpeedOfLight;
    a0 = Epeak * SpeedOfLight / mecsquared / omega;
    k = Vector(0.0 ,0.0 ,omega/SpeedOfLight);
  };

double PlaneWave::Omega()
{
  return omega;
}

double PlaneWave::EPeak()
{
  return Epeak;
}

double PlaneWave::BPeak()
{
  return Bpeak;
}

double PlaneWave::A0()
{
  return a0;
}

Vector PlaneWave::ElementLocalEField(double t, Vector X)
  {
    Vector E = Vector(Epeak*cos(omega*t-dot(k,X)), 0.0, 0.0);
    return(E);
  }

Vector PlaneWave::ElementLocalBField(double t, Vector X)
  {
    Vector B = Vector(0.0, Bpeak*cos(omega*t-dot(k,X)), 0.0);
    return(B);
  };
