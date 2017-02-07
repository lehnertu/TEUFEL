/************************************************************************/
/*                                                                      */
/*  TEUFEL - THz Emission from Undulators and Free-Electron Lasers      */
/*                                                                      */
/*  written by  U.Lehnert                               12/2016         */
/*                                                                      */
/************************************************************************/

#ifndef WAVE_H
#define WAVE_H

#include "vector.h"
#include "externalfield.h"

using namespace std;

// a class for an electromagnetic plane wave in vacuum
// the wave moves along the z axis
// the polarization direction is along the x axis

class PlaneWave : public ExternalField
{

  public:

    // constructor
    PlaneWave( double lambda,           // wavelength [m]
               double intensity );      // power density [W/cm²]

    double Omega();                     // report the frequency [Hz]
    double EPeak();                     // report the maximum electric field strength [V/m]
    double BPeak();                     // report the maximum magnetic field strength [T]
    double A0();                        // report the dimensionless intensity

  private:

    double Lambda;                              // wavelength [m]
    double omega;                               // frequency [Hz]
    double Intensity;                           // power density (averaged pointing vector) [W/m²]
    double Epeak;                               // peak electric field [V/m]
    double Bpeak;                               // peak magnetic field [T]
    double a0;                                  // dimensionless intensity
    Vector k;                                   // wave vector [1/m]

    Vector ElementLocalEField(double t, Vector X);
    Vector ElementLocalBField(double t, Vector X);

};

#endif