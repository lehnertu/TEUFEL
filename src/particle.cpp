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
#include <utility>
ChargedParticle::ChargedParticle()
{
	Charge = -1;
	Mass = 1;
	x0=Vector(0,0,0);
	p0=Vector(0,0,0);
	t0=0;
  	
	
}


ChargedParticle::ChargedParticle(int charge, int mass, Vector X0, Vector P0, double T0)
{
	Charge = charge;
	Mass = mass;
	x0 = X0;
	p0 = P0;
	t0 = T0;
	
	

}
ChargedParticle::ChargedParticle(const ChargedParticle *Part)
{
  	Charge = Part->Charge;
  	Mass = Part->Mass;
  	NOTS = Part->NOTS;
	t0 = Part->t0;
      	x0 = Part->x0;
      	p0= Part->p0;
	
	if(NOTS>0)
	{
   		Time = Part->Time;
       		X = Part->X;
   		P = Part->P;
   		A = Part->A;
   		
	}
	
   
      
}

ChargedParticle::~ChargedParticle()
{
  
}

void ChargedParticle :: init( int TrajLength)
{
	
	NOTS = TrajLength;		// set the number of time steps equal to trajectory length
	X.push_back(x0);
	P.push_back(p0);
	Time.push_back(t0);
	
	
}

void ChargedParticle::TrackEuler(
         int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
         Lattice *field )
{
  NOTS = Nstep+1;
  Time.clear();
  X.clear();
  P.clear();
  A.clear();
  Time.push_back(t0);
  X.push_back(X0);
  P.push_back(P0);
  double qm = Charge/Mass * InvRestMass ;    // charge over mass
  double betagamma2 = P0.abs2nd();
  double gamma = sqrt(betagamma2 + 1.0);
  Vector beta = P0 / gamma;
  pair<Vector,Vector>Fields = field->Field(t0,x0);
  Vector efield = Fields.first;
  Vector bfield = Fields.second;
  Vector force = cross(beta, bfield) + efield / SpeedOfLight;
  Vector dP_dt = force * qm;
    // store acceleration
  A.push_back(dP_dt);
    for (int i = 0; i < NOTS - 1; i++)
    {
	// compute derivatives of the particle motion
	// used to solve the equation of motion inside an undulator
	Vector p = P[i];
	betagamma2 = p.abs2nd();
	gamma = sqrt(betagamma2 + 1.0);
	beta = p / gamma;
	pair<Vector,Vector>Fields = field->Field(Time[i],X[i]);
  	efield = Fields.first;
 	bfield = Fields.second;
  	force = cross(beta, bfield) + efield / SpeedOfLight;
	Vector dX_dt = beta * SpeedOfLight;
	dP_dt = force * qm;
	// store the next trajectory point
	Time.push_back(Time[i] + tstep);
	X.push_back(X[i] + dX_dt * tstep);
	P.push_back(P[i] + dP_dt * tstep);
	A.push_back(dP_dt);
    };
}

void ChargedParticle::TrackLF(
         int Nstep,            // number of timesteps
         double tstep,         // time step size
         Vector X0,            // initial position
         Vector P0,            // initial momentum
         Lattice *field )
{
  NOTS=Nstep+1;
  Time.clear();
  X.clear();
  P.clear();
  A.clear();
  Time.push_back(0.0);
  X.push_back(X0);
  P.push_back(P0);
  // The algorithm defines velocities (momenta) at integer time steps
  // positions and field are computed at half-step points.
  // We store all quantities at the half-step points.
  double qm = Charge*InvRestMass/Mass;		// charge over mass
  double t_h = 0.0;             		// time at the half-step point
  Vector x_h = X0;              		// position at the half-step point
  Vector p_h = P0;              		// momentum at the half-step point
  double gamma_h = sqrt(p_h.abs2nd() + 1.0);
  Vector beta_h = p_h / gamma_h;
  pair<Vector,Vector>Fields = field->Field(t_h,x_h);
  Vector E_h = (Fields.first);
  Vector B_h = (Fields.second);
  Vector dp_dt = (cross(beta_h, B_h) + E_h/SpeedOfLight) * qm;
  A.push_back(dp_dt);
  // The leap-frog algorithm starts one half-step "before" the initial values
  // velocity p^(i+1) to be computed at the integer time step i+1
  // we store this initial velocity in p_i1 which is copied into p_i when we actually do the time step
  Vector p_i1 = p_h - dp_dt * 0.5 * tstep;
  double gamma_i1 = sqrt(p_i1.abs2nd() + 1.0);
  Vector beta_i1 = p_i1 / gamma_i1;
  for (int i=0; i<NOTS; i++)
  {
    // velocity u^i at the integer time step i
    // has been computed in the last time step
    Vector p_i = p_i1;
    Vector beta_i = beta_i1;
    // compute the velocity change over the integer step
    Fields = field->Field(t_h,x_h);
    Vector E_h = (Fields.first);
    Vector B_h = (Fields.second);
    // this is wrong, one should use beta at the half-step point which we don't know
    dp_dt = (cross(beta_i, B_h) + E_h/SpeedOfLight) * qm;
    p_i1 = p_i + dp_dt * tstep;
    gamma_i1 = sqrt(p_i1.abs2nd() + 1.0);
    beta_i1 = p_i1 / gamma_i1;
    // store all half-step quantities
    Time.push_back(t_h);
    X.push_back(x_h);
    P.push_back((p_i + p_i1)*0.5);
    A.push_back(dp_dt);
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
    NOTS = Nstep + 1;

    Time.clear();
    
    X.clear();
    
    P.clear();
    
    A.clear();
    
    // The Vay algorithm defines velocities (momenta) at integer time steps
    // positions and field are computed at half-step points.
    // We store all quantities at the half-step points.
    // The stored velocity is beta_h * gamma_h = u_h / c
    // u_h := u^(i+1/2)
    double qm = Charge/Mass * InvRestMass;      // charge over mass
    double t2 = 0.5 * tstep;                    // half time step
    double qmt2 = qm * t2;
    double t_h = 0.0;    // time at the half-step point
    Vector x_h = X0;     // position at the half-step point
    Vector p_h = P0;     // momentum at the half-step point
    double gamma_h = sqrt(p_h.abs2nd() + 1.0);
    Vector beta_h = p_h / gamma_h;
    pair<Vector,Vector>Fields = field->Field(t_h,x_h);
    Vector E_h = (Fields.first);
    Vector B_h = (Fields.second);
    Vector dp_dt = (cross(beta_h, B_h) + E_h/SpeedOfLight) * qm;
    A.push_back(dp_dt);
    // The leap-frog algorithm starts one half-step "before" the initial values
    // velocity p^(i+1) to be computed at the integer time step i+1
    // we store this initial velocity in p_i1 which is copied into p_i when we actually do the time step
    Vector p_i1 = p_h - dp_dt * 0.5 * tstep;
    double gamma_i1 = sqrt(p_i1.abs2nd() + 1.0);
    Vector beta_i1 = p_i1 / gamma_i1;
    for (int i = 0; i < NOTS; i++)
    {
	// velocity u^i at the integer time step i
	// has been computed in the last time step
	Vector p_i = p_i1;
	Vector beta_i = beta_i1;
	// compute the velocity change over the integer step
	pair<Vector,Vector>Fields = field->Field(t_h,x_h);
        Vector E_h = (Fields.first);
        Vector B_h = (Fields.second);
	dp_dt = (cross(beta_i, B_h) + E_h / SpeedOfLight) * qm;
	p_h = p_i + dp_dt * t2;
	Vector p_prime = p_h + E_h / SpeedOfLight * qmt2;
	double gamma_prime = sqrt(p_prime.abs2nd() + 1.0);
	Vector tau = B_h * qmt2;
	double u_star = dot(p_prime, tau);
	double tau2nd = tau.abs2nd();
	double sigma = gamma_prime * gamma_prime - tau2nd;
	// eq. (11)
	gamma_i1 = sqrt(0.5 * (sigma + sqrt(sigma * sigma + 4.0 * (tau2nd + u_star * u_star))));
	Vector t = tau / gamma_i1;
	// eq. (12)
	p_i1 = (p_prime + t * dot(p_prime, t) + cross(p_prime, t)) / (1 + t.abs2nd());
	beta_i1 = p_i1 / gamma_i1;
	// store all half-step quantities
	Time.push_back(t_h);
	X.push_back(x_h);
	P.push_back(p_h);
	A.push_back(dp_dt);
	// compute the change in position and time
	// for the next step
	x_h += beta_i1 * SpeedOfLight * tstep;
	t_h += tstep;
    };
}
void ChargedParticle::InitVay(double tstep,Lattice *Elfield)
{
	  Tstep = tstep;
	  field = Elfield;
	  qm = Charge*InvRestMass/Mass;		// charge over mass
	  t2 = 0.5 * tstep;              	// half time step
	  qmt2 = qm*t2;
	  t_h = Time[0];            // time at the half-step pointupload to a new 
	  x_h = X[0];              // position at the half-step point
	  p_h = P[0];              // momentum at the half-step point
	  gamma_h = sqrt(p_h.abs2nd() + 1.0);
	  beta_h = p_h / gamma_h;
	  pair<Vector,Vector>Fields = field->Field(t_h,x_h);
          Vector E_h = (Fields.first);
          Vector B_h = (Fields.second);
	  Vector dp_dt = (cross(beta_h, B_h) + E_h/SpeedOfLight) * qm;
	  A.push_back(dp_dt);
	  // The leap-frog algorithm starts one half-step "before" the initial values
	  // velocity p^(i+1) to be computed at the integer time step i+1
	  // we store this initial velocity in p_i1 which is copied into p_i when we actually do the time step
	  p_i1 = p_h - dp_dt * 0.5 * tstep;
	  gamma_i1 = sqrt(p_i1.abs2nd() + 1.0);
	  beta_i1 = p_i1 / gamma_i1;
	  
}
void ChargedParticle::StepVay(
	 Vector EField,	// extra fields that can be superimposed over lattice fields--ex: ambient field
	 Vector BField)	
{
  
    // velocity u^i at the integer time step i
    // has been computed in the last time step
    Vector p_i = p_i1;
    Vector beta_i = beta_i1;
    // compute the velocity change over the integer step
    pair<Vector,Vector>Fields = field->Field(t_h,x_h);
    Vector E_h = (Fields.first)+EField;
    Vector B_h = (Fields.second)+BField;
    Vector dp_dt = (cross(beta_i, B_h) + E_h/SpeedOfLight) * qm;
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
    counter+=1;
    // store all half-step quantities
    Time.push_back(t_h);
    X.push_back(x_h);
    P.push_back(p_h);
    A.push_back(dp_dt);
    // compute the change in position and time
    // for the next step
    x_h += beta_i1*SpeedOfLight*Tstep;
    t_h += Tstep;
    
}

pair<Vector,Vector> ChargedParticle::RetardedEField(double time, Vector ObservationPoint)
{
  Vector EField = Vector(0.0, 0.0, 0.0);
  Vector BField = Vector(0.0,0.0,0.0);
 
  int i1 = 0;                                   // index of the first trajectory point
  Vector RVec = ObservationPoint - X[i1];
  double R = RVec.norm();
  double t1 = Time[i1] + R / SpeedOfLight;      // retarded observation time

  int i2 = counter;                                // index of the last trajectory point
  RVec = ObservationPoint - X[i2];
  R = RVec.norm();
  double t2 = Time[i2] + R / SpeedOfLight;      // retarded observation time
  
  // the field is different from zero only if the observation
  // time is within the possible retarded time interval
  if ((time>=t1) && (time<t2))
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
    
    EField=EField*Charge*Coulomb;
    BField =cross(N/SpeedOfLight,EField);
  } 
  
  /*else
 {
	Vector R = ObservationPoint - X[counter];
	Vector UnitR = R/R.norm();
	if (R.norm()<=1.0e-6)
	{
		EField=Vector(0,0,0);
	}
	else
	{
		EField = UnitR*Charge*Coulomb/(R.abs2nd());
		
	}
	BField = Vector(0,0,0);
  }*/
  return make_pair(EField,BField);
}

void ChargedParticle::Translate(Vector R)
{
  if(NOTS !=0)
  {
	
  	for (int i=0; i<NOTS; i++)
    	{
      		X[i] = X[i] + R;
		
		
    	};
  }

  else
  {
	x0 = x0+R;
	
  }
}


void ChargedParticle::MirrorY(double MirrorY)
{
  Charge = - getCharge();
  Mass = getMass();
  if(NOTS !=0)
  {
  	for (int i=0; i<NOTS; i++)
    	{
      		X[i].y = 2.0*MirrorY - X[i].y;
      		P[i].y = -P[i].y;
      		A[i].y = -A[i].y;
    	}

  }
  
}

double ChargedParticle::TrajTime(int i)
{
	if(NOTS != 0)
	{
		return Time[i];
	}

	else
	{
		return t0;
	}
}

Vector ChargedParticle::TrajPoint(int i)
{
	if(NOTS != 0)
	{
		return X[i];
	}

	else
	{
		return x0;
	}
	
}


Vector ChargedParticle::TrajAccel(int i)
{
	if(NOTS != 0)
	{
		return A[i];
	}

	else
	{
		return Vector(0,0,0);
	}
	
}
Vector ChargedParticle::TrajMomentum(int  i)
{
	if(NOTS != 0)
	{
		return P[i];
	}

	else
	{
		return p0;
	}
	
}
int ChargedParticle::getNOTS()
{
  return NOTS;
}


int ChargedParticle::getCharge()
{
  return Charge;
}

int ChargedParticle::getMass()
{
  return Mass;
}





