/************************************************************************/
/*                                                                      */
/*  TEUFEL - THz Emission from Undulators and Free-Electron Lasers      */
/*                                                                      */
/*  written by  U.Lehnert                               12/2016         */
/*                                                                      */
/************************************************************************/

#include "undulator.h"
#include "global.h"

#include <math.h>

Undulator::Undulator(   double B,                         // peak field [T]
                        double lambda,                    // undulator period [m]
                        int    N                          // number of undulator periods
             )
{
  BPeak = B;
  LambdaU = lambda;
  NPeriods = N;
  Krms=LambdaU*SpeedOfLight*BPeak/(2.0*Pi*mecsquared)/sqrt(2.0);
  ky = 2.0*Pi/LambdaU;
  kz = 2.0*Pi/LambdaU;
}

double Undulator::GetBPeak()
{
  return BPeak;
}

double Undulator::GetLambdaU()
{
  return LambdaU;
}

int Undulator::GetNPeriods()
{
  return NPeriods;
}

double Undulator::GetKpeak()
{
  return Krms*sqrt(2.0);
}

double Undulator::GetKrms()
{
  return Krms;
}

Vector Undulator::ElementLocalEField(double t, Vector X)
  {
    return Vector(0.0 ,0.0 ,0.0);
  }

Vector Undulator::ElementLocalBField(double t, Vector X)
  {
    Vector B;
    B.x=0.0;
    if ((X.z>=0.0) && (X.z<=LambdaU*NPeriods)) {
      B.y=BPeak*cosh(ky*X.y)*sin(kz*X.z);
      B.z=BPeak*sinh(ky*X.y)*cos(kz*X.z);
      }
    else {
      B.y=0.0;
      B.z=0.0;
      };
    return(B);
  };
