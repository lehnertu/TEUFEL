/************************************************************************/
/*                                                                      */
/*  TEUFEL - THz Emission from Undulators and Free-Electron Lasers      */
/*                                                                      */
/*  written by  U.Lehnert                               12/2016         */
/*                                                                      */
/************************************************************************/

#ifndef CAVITY_H
#define CAVITY_H

#include "vector.h"
#include "externalfield.h"

using namespace std;

// a class for a pillbox accelerator cavity

// the cavity is aligned with the z axis
// reference position is the cavity entrance wall
// the mode number is the number of half-wave periods
// in longitudinal direction present inside the cavity - TM(01n)

// the beam direction is along the z axis starting at z=0

class PillboxCavity : public ExternalField
{

  public:

    // constructor
    PillboxCavity( double freq,                        // cavity frequency [Hz]
                   double grad,                        // peak electric field gradient on axis [V/m]
                   double length,                      // length of the cavity [m]
                   double rcell,                       // radius of the cavity [m]
                   int    nlong );                     // mode number in longitudinal direction

  private:

    double Freq;                                // cavity frequency [Hz]
    double Grad;                                // peak electric field gradient on axis [V/m]
    double Length;                              // length of a (half) cell [m]
    double R;                                   // radius of the cavity [m]
    double R2;                                  // the square of it
    int    Nlong;                               // number of (half) cells

    double kz;                                  // cavity periodicity [1/m]
    double omega;                               // cavity frequency [Hz]

    Vector ElementLocalEField(double t, Vector X);
    Vector ElementLocalBField(double t, Vector X);

};

#endif