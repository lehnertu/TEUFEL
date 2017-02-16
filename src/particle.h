/************************************************************************/
/*                                                                      */
/*  TEUFEL - THz Emission from Undulators and Free-Electron Lasers      */
/*                                                                      */
/*  written by  U.Lehnert                               12/2016         */
/*                                                                      */
/************************************************************************/

#ifndef PARTICLE_H
#define PARTICLE_H

#include "externalfield.h"
#include "vector.h"
#include "undulator.h"

using namespace std;

// a class for charged particle trajectories

class ChargedParticle
{

  public:

    // default constructor
    // charge set to electron charge
    ChargedParticle();

    // copy constructor
    ChargedParticle(const ChargedParticle *Part);

    // destructor
    ~ChargedParticle();

    // return the number of points in the track
    int GetNP();

    // return a position from the trajectory
    double TrajTime(int step);
    Vector TrajPoint(int step);
    Vector TrajMomentum(int step);

    // track the particle through a given field
    // use the most simple Euler algorithm
    void TrackEuler
      (  int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
         Lattice *field);

    // track the particle through a given field
    // using a leap-frog algorithm
    void TrackLF
      (  int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
         Lattice *field);

    // track the particle through a given field
    // using the Vay algorithm
    void TrackVay
      (  int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
         Lattice *field);

    // translate a given particle trajectory
    void Translate(Vector R);

    // mirroring a given particle trajectory
    // on a plane y=MirrorY
    // (this also reverses the charge)
    void MirrorY(double MirrorY);

    // electric field radiated by the particle
    // at a given observation point at time t in lab frame
    // retardation is properly accounted for
    Vector RetardedEField(double time, Vector ObservationPoint);

  private:

    int NP;			// number of trajectory points
    int Charge;			// charge in units of ElementaryCharge
    int Mass;			// mass in unit of the electron rest mass
    double *Time;		// time in lab-frame [s]
    Vector *X;			// position in lab frame [m]
    Vector *P;                  // momentum in lab frame : c p = beta gamma mc²
                                // dimensionless in units of mc²
    Vector *A;                  // acceleration in lab frame a = d/dt(p*c)
				// in unit of 1/s (scaled by mc²)
};

#endif
