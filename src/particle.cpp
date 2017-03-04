/*=========================================================================
 * 
 *  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers
 * 
 *  Copyright (c) 2017 U. Lehnert
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * =========================================================================*/

#include "particle.h"
#include "global.h"
#include <iostream>
#include <math.h>
#include <assert.h>

ChargedParticle::ChargedParticle()
{
  // no trajectory points - by default 
  // no memory allocated, all pointers are zero
  Time = 0;
  X = 0;
  P = 0;
  A = 0;
}

ChargedParticle::ChargedParticle(const ChargedParticle *Part)
{
  Charge = Part->Charge;
  Mass = Part->Mass;
  NP = Part->NP;
  if (NP>0)
    {
      Time = new double[NP];
      X = new Vector[NP];
      P = new Vector[NP];
      A = new Vector[NP];
      for (int i=0; i<NP; i++)
        {
          Time[i] = Part->Time[i];
          X[i] = Part->X[i];
          P[i] = Part->P[i];
          A[i] = Part->A[i];
        };
    } else {
      Time = 0;
      X = 0;
      P = 0;
      A = 0;
    }
}

ChargedParticle::~ChargedParticle()
{
  if (Time!=0) delete[] Time;
  if (X!=0) delete[] X;
  if (P!=0) delete[] P;
  if (A!=0) delete[] A;
}

int ChargedParticle::GetNP()
{
  return NP;
}

void ChargedParticle::setCharge(int charge)
{
  Charge=charge;
}

int ChargedParticle::getCharge()
{
  return Charge;
}

void ChargedParticle::setTrajPoint(int stepnumber,Vector Position)
{
  assert(0<=stepnumber &&stepnumber<=NP);
  X[stepnumber]=Position;
}

void ChargedParticle::setTrajMomentum(int stepnumber,Vector Momentum)
{
  assert(0<=stepnumber &&stepnumber<=NP);
  P[stepnumber]=Momentum;
}

void ChargedParticle:: setTrajAcceleration(int stepnumber,Vector Accel)
{
  assert(stepnumber<=NP);
  A[stepnumber]=Accel;
}

void ChargedParticle:: setTrajTime(int stepnumber, double time)
{
  assert(0<=stepnumber && stepnumber<=NP);
  Time[stepnumber]=time;
}

//set the number of trajectory points the particle can have
//initialise the trajectory arrays
void ChargedParticle::setNP(int NOTS)
{
  NP=NOTS;
  X=new Vector[NP];
  X[0]=getInitialPosition();
  P=new Vector[NP];
  P[0]=getInitialMomentum();
  A=new Vector[NP];
  Time=new double[NP];
  Time[0]=getInitialTime();
}


void ChargedParticle::setMass(int mass)
{
  Mass=mass;
}

int ChargedParticle::getMass()
{
  return Mass;
}

void ChargedParticle::setParticleSize(double L)
{
   Radius=L;
}

double ChargedParticle::getParticleSize()
{
   return Radius;

}

void ChargedParticle::setParticleID(int ParticleID)
{
  ID=ParticleID;
}

int ChargedParticle::getParticleID()
{
  return ID;
}
void ChargedParticle::setInitialPosition(Vector InitialPosition)
{
  X0=InitialPosition;
}

void ChargedParticle::setInitialMomentum(Vector InitialMomentum)
{
  P0=InitialMomentum;
}

Vector ChargedParticle::getInitialPosition()
{
  return X0;
}

void ChargedParticle::setInitialTime(double t)
{
  T0=t;
}

double ChargedParticle::getInitialTime()
{
  return T0;
}

Vector ChargedParticle::getInitialMomentum()
{
  return P0;
}

double ChargedParticle::TrajTime(int step)
{
  double t = 0.0;
  if (step>=0 && step<NP)
    t = Time[step];
  return t;
}

Vector ChargedParticle::TrajPoint(int step)
{
  Vector p = Vector(0.0, 0.0, 0.0);
  if (step>=0 && step<NP)
    p = X[step];
  return p;
}

Vector ChargedParticle::TrajMomentum(int step)
{
  Vector p = Vector(0.0, 0.0, 0.0);
  if (step>=0 && step<NP)
    p = P[step];
  return p;
}

void ChargedParticle::TrackEuler(
         int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
         Lattice *field )
{
  NP = Nstep+1;
  if (Time!=0) delete[] Time;
  Time = new double[NP];
  if (X!=0) delete[] X;
  X = new Vector[NP];
  if (P!=0) delete[] P;
  P = new Vector[NP];
  if (A!=0) delete[] A;
  A = new Vector[NP];
  Time[0] = 0.0;
  X[0] = X0;
  P[0] = P0;
  double qm = Charge*InvRestMass/Mass;     // charge over mass
  for (int i=0; i<NP-1; i++)
  {
    // compute derivatives of the particle motion
    // used to solve the equation of motion inside an undulator
    Vector p = P[i];
    double betagamma2 = p.abs2nd();
    double gamma = sqrt(betagamma2 + 1.0);
    Vector beta = p / gamma;
    Vector efield = field->EField(Time[i], X[i]);
    Vector bfield = field->BField(Time[i], X[i]);
    Vector force = cross(beta, bfield) + efield/SpeedOfLight;
    Vector dX_dt = beta * SpeedOfLight;
    Vector dP_dt = force * qm;
    // store acceleration
    A[i] = dP_dt;
    // integrator step
    Time[i+1] = Time[i] + tstep;
    X[i+1] = X[i] + dX_dt * tstep;
    P[i+1] = P[i] + dP_dt * tstep;
  };
}

void ChargedParticle::TrackLF(
         int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
         Lattice *field )
{
  NP = Nstep+1;
  if (Time!=0) delete[] Time;
  Time = new double[NP];
  if (X!=0) delete[] X;
  X = new Vector[NP];
  if (P!=0) delete[] P;
  P = new Vector[NP];
  if (A!=0) delete[] A;
  A = new Vector[NP];
  // The algorithm defines velocities (momenta) at integer time steps
  // positions and field are computed at half-step points.
  // We store all quantities at the half-step points.
  double qm = Charge*InvRestMass/Mass;		// charge over mass
  double t_h = 0.0;             		// time at the half-step point
  Vector x_h = X0;              		// position at the half-step point
  Vector p_h = P0;              		// momentum at the half-step point
  double gamma_h = sqrt(p_h.abs2nd() + 1.0);
  Vector beta_h = p_h / gamma_h;
  Vector E_h = field->EField(t_h, x_h);
  Vector B_h = field->BField(t_h, x_h);
  Vector dp_dt = (cross(beta_h, B_h) + E_h/SpeedOfLight) * qm;
  A[0] = dp_dt;
  // The leap-frog algorithm starts one half-step "before" the initial values
  // velocity p^(i+1) to be computed at the integer time step i+1
  // we store this initial velocity in p_i1 which is copied into p_i when we actually do the time step
  Vector p_i1 = p_h - dp_dt * 0.5 * tstep;
  double gamma_i1 = sqrt(p_i1.abs2nd() + 1.0);
  Vector beta_i1 = p_i1 / gamma_i1;
  for (int i=0; i<NP; i++)
  {
    // velocity u^i at the integer time step i
    // has been computed in the last time step
    Vector p_i = p_i1;
    Vector beta_i = beta_i1;
    // compute the velocity change over the integer step
    E_h = field->EField(t_h, x_h);
    B_h = field->BField(t_h, x_h);
    // this is wrong, one should use beta at the half-step point which we don't know
    dp_dt = (cross(beta_i, B_h) + E_h/SpeedOfLight) * qm;
    p_i1 = p_i + dp_dt * tstep;
    gamma_i1 = sqrt(p_i1.abs2nd() + 1.0);
    beta_i1 = p_i1 / gamma_i1;
    // store all half-step quantities
    Time[i] = t_h;
    X[i] = x_h;
    P[i] = (p_i + p_i1)*0.5;
    A[i] = dp_dt;
    // compute the change in position and time
    // for the next step
    x_h += beta_i1*SpeedOfLight*tstep;
    t_h += tstep;
  };
}

void ChargedParticle::TrackVay(
         int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
         Lattice *field )
{
  NP = Nstep+1;
  if (Time!=0) delete[] Time;
  Time = new double[NP];
  if (X!=0) delete[] X;
  X = new Vector[NP];
  if (P!=0) delete[] P;
  P = new Vector[NP];
  if (A!=0) delete[] A;
  A = new Vector[NP];
  // The Vay algorithm defines velocities (momenta) at integer time steps
  // positions and field are computed at half-step points.
  // We store all quantities at the half-step points.
  // The stored velocity is beta_h * gamma_h = u_h / c
  // u_h := u^(i+1/2)
  double qm = Charge*InvRestMass/Mass;		// charge over mass
  double t2 = 0.5 * tstep;              	// half time step
  double qmt2 = qm*t2;
  double t_h = 0.0;             // time at the half-step point
  Vector x_h = X0;              // position at the half-step point
  Vector p_h = P0;              // momentum at the half-step point
  double gamma_h = sqrt(p_h.abs2nd() + 1.0);
  Vector beta_h = p_h / gamma_h;
  Vector E_h = field->EField(t_h, x_h);
  Vector B_h = field->BField(t_h, x_h);
  Vector dp_dt = (cross(beta_h, B_h) + E_h/SpeedOfLight) * qm;
  A[0] = dp_dt;
  // The leap-frog algorithm starts one half-step "before" the initial values
  // velocity p^(i+1) to be computed at the integer time step i+1
  // we store this initial velocity in p_i1 which is copied into p_i when we actually do the time step
  Vector p_i1 = p_h - dp_dt * 0.5 * tstep;
  double gamma_i1 = sqrt(p_i1.abs2nd() + 1.0);
  Vector beta_i1 = p_i1 / gamma_i1;
  for (int i=0; i<NP; i++)
  {
    // velocity u^i at the integer time step i
    // has been computed in the last time step
    Vector p_i = p_i1;
    Vector beta_i = beta_i1;
    // compute the velocity change over the integer step
    E_h = field->EField(t_h, x_h);
    B_h = field->BField(t_h, x_h);
    dp_dt = (cross(beta_i, B_h) + E_h/SpeedOfLight) * qm;
    p_h = p_i + dp_dt * t2;
    Vector p_prime = p_h + E_h/SpeedOfLight * qmt2;
    double gamma_prime = sqrt(p_prime.abs2nd() + 1.0);
    Vector tau = B_h * qmt2;
    double u_star = dot(p_prime,tau);
    double tau2nd = tau.abs2nd();
    double sigma = gamma_prime*gamma_prime - tau2nd;
    gamma_i1 = sqrt(0.5*(sigma+sqrt(sigma*sigma+4.0*(tau2nd+u_star*u_star))));
    Vector t = tau / gamma_i1;
    p_i1 = (p_prime + t*dot(p_prime,t) + cross(p_prime,t))/(1+t.abs2nd());
    beta_i1 = p_i1 / gamma_i1;
    // store all half-step quantities
    Time[i] = t_h;
    X[i] = x_h;
    P[i] = p_h;
    A[i] = dp_dt;
    // compute the change in position and time
    // for the next step
    x_h += beta_i1*SpeedOfLight*tstep;
    t_h += tstep;
  };
}

void ChargedParticle::Translate(Vector R)
{
  for (int i=0; i<NP; i++)
    {
      X[i] = X[i] + R;
    };
}

void ChargedParticle::MirrorY(double MirrorY)
{
  Charge = - Charge;
  for (int i=0; i<NP; i++)
    {
      X[i].y = 2.0*MirrorY - X[i].y;
      P[i].y = -P[i].y;
      A[i].y = -A[i].y;
    };
}

tuple<Vector,Vector> ChargedParticle::RetardedEField(double time, Vector ObservationPoint)
{
  Vector EField = Vector(0.0, 0.0, 0.0);
  Vector BField = Vector(0.0,0.0,0.0);
  double scale=(Charge*ElementaryCharge/(4.0*Pi*8.85e-12));
  int i1 = 0;                                   // index of the first trajectory point
  Vector RVec = ObservationPoint - X[i1];
  double R = RVec.norm();
  double t1 = Time[i1] + R / SpeedOfLight;      // retarded observation time

  int i2 = NP-1;                                // index of the last trajectory point
  RVec = ObservationPoint - X[i2];
  R = RVec.norm();
  double t2 = Time[i2] + R / SpeedOfLight;      // retarded observation time

  // the field is different from zero only if the observation
  // time is within the possible retarded time interval
  if ((time>=t1) && (time<=t2))
  {
    // reduce the interval until the trajectory segment is found
    while (i2-i1>1)
    {
      int i = (i2+i1)/2;
      RVec = ObservationPoint - X[i];
      R = RVec.norm();
      double t = Time[i] + R / SpeedOfLight;    // retarded observation time
      if (t < time)
      {
        i1 = i;
        t1 = t;
      }
      else
      {
        i2 = i;
        t2 = t;
      }
    }
    // interpolate the source point within the interval
    // interpolation could be improved using higher-order terms
    double frac = (time-t1)/(t2-t1);
    Vector SourceX = X[i1]*(1.0-frac) + X[i2]*frac;
    Vector SourceBeta = P[i1]*(1.0-frac) + P[i2]*frac;
    double betagamma2 = SourceBeta.abs2nd();
    double gamma = sqrt(betagamma2 + 1.0);
    SourceBeta/=gamma;
    Vector SourceBetaPrime = A[i1]*(1.0-frac) + A[i2]*frac;
    SourceBetaPrime/=gamma;
    // now compute the field emitted from an interpolated source point
    RVec = ObservationPoint - SourceX;
    R = RVec.norm();
    Vector N = RVec;
    N.normalize();
    double bn3rd = pow(1.0-dot(SourceBeta,N),3.0);
    // velocity term
    EField += (N-SourceBeta)/(R*R*bn3rd)/(gamma*gamma);
    // acceleration term
    EField += cross(N,cross(N-SourceBeta,SourceBetaPrime))/(R*bn3rd)/SpeedOfLight;
    
    EField=EField*scale;
    BField =cross(N/SpeedOfLight,EField);
  } 
  return make_tuple(EField,BField);
}



//find the radiation field
// because of particle identified by ParticleID1 --source
// on Particle identified by ParticleID2 --observer
//Particle1D1 != ParticleID2
//stepnumber is the iteration step number --to identify the last known trajectory point
//time is the lab time or iteration time 
//Observation Point is the position of the particle identified by ParticleID2

tuple<Vector,Vector> ChargedParticle::InteractionField(int ParticleID2,int stepnumber,double time, Vector ObservationPoint)
{

  Vector EField = Vector(0.0, 0.0, 0.0);
  Vector BField = Vector(0.0,0.0,0.0);
  int ParticleID1 = getParticleID();
  double scale=(Charge*ElementaryCharge/(4.0*Pi*8.85e-12)); 
  if(ParticleID1 != ParticleID2)
  {
	  
	  
	  int i1 = 0;                                   // index of the first trajectory point
	  Vector RVec = ObservationPoint - X[i1];
	  double R = RVec.norm();
	  double t1 = Time[i1] + R / SpeedOfLight;      // retarded observation time

	  int i2 = stepnumber;                                // index of the last known trajectory point
	  RVec = ObservationPoint - X[i2];
	  R = RVec.norm();
	  double t2 = Time[i2] + R / SpeedOfLight;      // retarded observation time

	  // the field is different from zero only if the observation
	  // time is within the possible retarded time interval
	  if ((time>=t1) && (time<=t2))
	  {
	    // reduce the interval until the trajectory segment is found
	    while (i2-i1>1)
	    {
	      
	      int i = (i2+i1)/2;
	      RVec = ObservationPoint - X[i];
	      R = RVec.norm();
	      double t = Time[i] + R / SpeedOfLight;    // retarded observation time
	      if (t < time)
	      {
		i1 = i;
		t1 = t;
	      }
	      else
	      {
		i2 = i;
		t2 = t;
	      }
	    }
	    // interpolate the source point within the interval
	    // interpolation could be improved using higher-order terms
	    double frac = (time-t1)/(t2-t1);
	    Vector SourceX = X[i1]*(1.0-frac) + X[i2]*frac;
	    Vector SourceBeta = P[i1]*(1.0-frac) + P[i2]*frac;
	    double betagamma2 = SourceBeta.abs2nd();
	    double gamma = sqrt(betagamma2 + 1.0);
	    SourceBeta/=gamma;
	    Vector SourceBetaPrime = A[i1]*(1.0-frac) + A[i2]*frac;
	    SourceBetaPrime/=gamma;
	    // now compute the field emitted from an interpolated source point
	    RVec = ObservationPoint - SourceX;
	    R = RVec.norm();
	    //check whether the particles are overlapping
	    //if overlapping; set the distance between them
	    //but it can also be set according to need
	    if(R<Radius)
	    {
		//cout<<"Warning: Particles Overlap\n";
		R=2*Radius;
	    }
	    Vector N = RVec;
	    N.normalize();
	    double bn3rd = pow(1.0-dot(SourceBeta,N),3.0);
	    // velocity term
	    EField += (N-SourceBeta)/(R*R*bn3rd)/(gamma*gamma);
	    // acceleration term
	    EField += cross(N,cross(N-SourceBeta,SourceBetaPrime))/(R*bn3rd)/SpeedOfLight;
	    EField=EField*scale;
	    BField =cross(N/SpeedOfLight,EField);
	  }
	  else
	  {
	    Vector RVec=ObservationPoint-X[stepnumber];
	    if (RVec.norm()>Radius)
	    {    
	      EField=RVec/(pow(RVec.norm(),3.0));
	      EField=EField*scale;
	    }
	    else
	    {
	      //cout<<"Warning: Particles Overlap\n";
	      RVec.normalize();
	      //separate the particles by a diameter
	      EField=RVec/pow(2.0*Radius,2.0);
	      EField=EField*scale;
	     }
	    
	  }
  }
  return make_tuple(EField,BField);
}
