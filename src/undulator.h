/************************************************************************/
/*                                                                      */
/*  TEUFEL - THz Emission from Undulators and Free-Electron Lasers      */
/*                                                                      */
/*  written by  U.Lehnert                               12/2016         */
/*                                                                      */
/************************************************************************/

#ifndef UNDULATOR_H
#define UNDULATOR_H

#include "vector.h"
#include "externalfield.h"

using namespace std;

// a class for a planar undulator

// the beam direction is along the z axis starting at z=0
// the main magnetic field component is in y direction
// so the undulator deflects in x direction
// the field profile is flat in x direction (no focusssing)

class Undulator : public ExternalField
{

  public:

    // constructor
    Undulator(double B,                         // peak field [T]
              double lambda,                    // undulator period [m]
              int    N                          // number of undulator periods
             );

    double  GetBPeak();
    double  GetLambdaU();
    int     GetNPeriods();
    double  GetKpeak();
    double  GetKrms();

  private:

    Vector ElementLocalEField(double t, Vector X);
    Vector ElementLocalBField(double t, Vector X);

    double  BPeak;                              // peak field [T]
    double  LambdaU;                            // undulator period [m]
    int     NPeriods;                           // number of undulator periods
    double  Krms;                               // undulator parameter
    double  ky, kz;                             // undulator periodicity [1/m]
};

#endif