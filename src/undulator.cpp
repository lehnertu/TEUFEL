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
    Vector B=Vector (0,0,0);
    double z1=-GetLambdaU()/2;
    double z2=GetLambdaU()/2;
    B.x=0.0;
    if(X.z<z1)
    {
        B.y=0;
        B.z=0;
    }
    if(z1<=X.z&&X.z<z2)
    {
        B.y=((X.z-z1)/(z2-z1))*BPeak*sin(kz*X.z)*cosh(kz*X.y);
        B.z=((X.z-z1)/(z2-z1))*BPeak*cos(kz*X.z)*sin(kz*X.y);
    }
    if(z2<=X.z&&X.z<LambdaU*NPeriods-z2)
    {
        B.y=BPeak*sin(kz*X.z)*cosh(kz*X.y);
        B.z=BPeak*cos(kz*X.z)*sinh(kz*X.y);
    }
    if(LambdaU*NPeriods-z2<=X.z&&X.z<=LambdaU*NPeriods-z1)
    {
        B.y=((LambdaU*NPeriods-z1-X.z)/(z2-z1))*BPeak*sin(kz*X.z)*cosh(kz*X.y);
        B.z=((LambdaU*NPeriods-z1-X.z)/(z2-z1))*BPeak*cos(kz*X.z)*sinh(kz*X.y);
    }

    return(B);
  };
