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

ChargedParticle::ChargedParticle()
{
    Charge = -1.0;
    Mass = 1.0;
    dt = 0.0;
    // no trajectory points
    NP = 0;
}

ChargedParticle::ChargedParticle(double charge, double mass)
{
    Charge = charge;
    Mass = mass;
    dt = 0.0;
    // no trajectory points
    NP = 0;
}

ChargedParticle::ChargedParticle(const ChargedParticle *origin)
{
    Charge = origin->Charge;
    Mass = origin->Mass;
    NP = origin->NP;
    if (NP > 0)
    {
	Time = origin->Time;
	X = origin->X;
	P = origin->P;
	A = origin->A;
    }
    dt = origin->dt;
    qm = origin->qm;
    qmt2 = origin->qmt2;
    VY_p_i1 = origin->VY_p_i1;
    VY_gamma_i1 = origin->VY_gamma_i1;
}

ChargedParticle::~ChargedParticle()
{
}

double ChargedParticle::getTime()
{
    if (NP>0)
	return Time.back();
    else
	return 0.0;
}

Vector ChargedParticle::getPosition()
{
    if (NP>0)
	return X.back();
    else
	return Vector(0.0, 0.0, 0.0);
}

Vector ChargedParticle::getMomentum()
{
    if (NP>0)
	return P.back();
    else
	return Vector(0.0, 0.0, 0.0);
}

double ChargedParticle::TrajTime(int step)
{
    double t = 0.0;
    if (step >= 0 && step < NP)
	t = Time[step];
    return t;
}

Vector ChargedParticle::TrajPoint(int step)
{
    Vector p = Vector(0.0, 0.0, 0.0);
    if (step >= 0 && step < NP)
	p = X[step];
    return p;
}

Vector ChargedParticle::TrajMomentum(int step)
{
    Vector p = Vector(0.0, 0.0, 0.0);
    if (step >= 0 && step < NP)
	p = P[step];
    return p;
}

Vector ChargedParticle::TrajAccel(int step)
{
    Vector a = Vector(0.0, 0.0, 0.0);
    if (step >= 0 && step < NP)
	a = A[step];
    return a;
}

void ChargedParticle::initTrajectory(double t, Vector x, Vector p, Vector a)
{
    NP = 1;
    Time.clear();
    X.clear();
    P.clear();
    A.clear();
    Time.push_back(t);
    X.push_back(x);
    P.push_back(p);
    A.push_back(a);
}

/*!
 * Implementation details:
 * ----------------------
 * 
 * The trajectory data (position, momentum, acceleration)
 * is stored for the half-integer time steps of the algorithm.
 * 
 * This initialization routine sets the position and velocity of the first
 * point to the given start values and computes the velocity at the end of the
 * time step (half way to the next storage point). The given inital momentum is
 * treated as \f$\mathbf{u}^{i+1/2}\f$ of the algorithm. The computed momentum at the end of
 * the time step is stored in the ChargedParticle::VY_p_i1 variable for use in subsequent time steps.
 * (Note, that instead of \f$\mathbf{u}\f$ we are using \f$\mathbf{p}=\mathbf{u}/c\f$).
 * In addition, we store the relativistic factor belonging to this momentum ChargedParticle::VY_gamma_i1,
 * the charge over mass ratio ChargedParticle::qm
 * and charge over mass divided by the double time step ChargedParticle::qmt2.
 * These values as well as the time step ChargedParticle::dt are constant during the tracking.
 */
void ChargedParticle::InitVay(
	   double tstep,
	   GeneralField* field)
{
    // check prerequisites
    if (NP<1)
	throw("ChargedParticle::InitVay() - must initalize trajectory before tracking");
    // preset the fixed values
    dt = tstep;
    qm = Charge/Mass * InvRestMass;
    qmt2 = qm * 0.5 * dt;
    // compute the fields at the center position of the step (which is stored in the trajectory)
    ElMagField EB_h = field->Field(Time.back(), X.back());
    Vector E_h = EB_h.E();
    Vector B_h = EB_h.B();
    // perform the second half-step to compute the
    // velocity at the end of the time step
    Vector p_prime = P.back() + E_h / SpeedOfLight * qmt2;
    double gamma2_prime = p_prime.abs2nd() + 1.0;
    Vector tau = B_h * qmt2;
    double u_star = dot(p_prime, tau);
    double tau2 = tau.abs2nd();
    double sigma = gamma2_prime - tau2;
    // eq. (11)
    VY_gamma_i1 = sqrt(0.5 * (sigma + sqrt(sigma * sigma + 4.0 * (tau2 + u_star * u_star))));
    Vector t = tau / VY_gamma_i1;
    // eq. (12)
    VY_p_i1 = (p_prime + t * dot(p_prime, t) + cross(p_prime, t)) / (1 + t.abs2nd());
    double gamma_h = sqrt(P.back().abs2nd() + 1.0);
    Vector beta_h = P.back() / gamma_h;
    A.back() = (cross(beta_h, B_h) + E_h / SpeedOfLight) * qm;
}

/*!
 * Implementation details:
 * ----------------------
 * 
 * The trajectory data is stored for the half-integer time steps of the algorithm.
 * The time step takes the velocity at the end point of the last step (at the intermediate
 * point between two storage points) and computes the new helf-step (center) position
 * of the trajectory. This is where the fields are sampled. Then the velocity for the next
 * step is computed. (Note, that instead of \f$\mathbf{u}\f$ we are using \f$\mathbf{p}=\mathbf{u}/c\f$).
 * 
 * This stepper assumes that InitVay() has been called to perform the frist tracking
 * step. This sets a number of variables which carry data from one tracking step to the next.
 * These are : ChargedParticle::dt ChargedParticle::qm ChargedParticle::qmt2 ChargedParticle::VY_p_i1
 * 
 * The time step and particle properties can not change during the tracking. Only the field
 * is allowed to change. Therefore, the field reference is given for every tracking step.
 */
void ChargedParticle::StepVay(GeneralField* field)
{
    // the momentum vector at the beginning of the step is known from the previous step
    Vector p_i = VY_p_i1;
    Vector beta_i = p_i / VY_gamma_i1;
    // use eq. (3) to compute the new center position
    double t_h = Time.back() + dt;
    Vector x_h = X.back() + beta_i * SpeedOfLight * dt;
    // compute the fields at the center position of the step
    ElMagField EB_h = field->Field(t_h, x_h);
    Vector E_h = EB_h.E();
    Vector B_h = EB_h.B();
    // perform the first half-step to compute the
    // velocity at the center (half) point using eq. (13)
    Vector p_h = VY_p_i1 + (E_h/SpeedOfLight + cross(beta_i, B_h)) * qmt2;
    double gamma_h = sqrt(p_h.abs2nd() + 1.0);
    Vector beta_h = p_h / gamma_h;
    // perform the second half-step to compute the
    // velocity at the end of the time step
    Vector p_prime = p_h + E_h / SpeedOfLight * qmt2;
    double gamma2_prime = p_prime.abs2nd() + 1.0;
    Vector tau = B_h * qmt2;
    double u_star = dot(p_prime, tau);
    double tau2 = tau.abs2nd();
    double sigma = gamma2_prime - tau2;
    // eq. (11)
    VY_gamma_i1 = sqrt(0.5 * (sigma + sqrt(sigma * sigma + 4.0 * (tau2 + u_star * u_star))));
    Vector t = tau / VY_gamma_i1;
    // eq. (12)
    VY_p_i1 = (p_prime + t * dot(p_prime, t) + cross(p_prime, t)) / (1 + t.abs2nd());
    // store the trajectory data
    Time.push_back(t_h);
    X.push_back(x_h);
    P.push_back(p_h);
    A.push_back((cross(beta_h, B_h) + E_h / SpeedOfLight) * qm);
    // keep track of the number of stored trajectory points
    NP++;
}

void ChargedParticle::CoordinatesAtTime(double time, Vector *position, Vector *momentum)
{
    // return zero if there is no trajectory
    if (NP<1)
    {
	*position = Vector(1.0, 1.0, 1.0);
	*momentum = Vector(0.0, 0.0, 0.0);
	return;
    };
    int i1 = 0;    // index of the first trajectory point
    double t1 = Time[i1];
    int i2 = NP - 1;    // index of the last trajectory point
    double t2 = Time[i2];
    if (time>=t1 && time<=t2)
    {
	// reduce the interval until the trajectory segment is found
	while (i2 - i1 > 1)
	{
	    int i = (i2 + i1) / 2;
	    double t = Time[i];
	    if (t < time)
	    {
		i1 = i;
		t1 = t;
	    } else {
		i2 = i;
		t2 = t;
	    }
	};
	if (i2==i1)
	{
	    *position = X[i1];
	    *momentum = P[i1];
	} else {
	    // interpolate the coordinates within the interval
	    // interpolation could be improved using higher-order terms
	    double frac = (time - t1) / (t2 - t1);
	    *position = X[i1] * (1.0 - frac) + X[i2] * frac;
	    *momentum = P[i1] * (1.0 - frac) + P[i2] * frac;
	};
    } else {
	*position = Vector(2.0, 2.0, 2.0);
	*momentum = Vector(0.0, 0.0, 0.0);
    }
}

double ChargedParticle::RetardedTime(int index,
				     Vector ObservationPoint)
{
    Vector RVec = ObservationPoint - X[index];
    double R = RVec.norm();
    return(Time[index] + R / SpeedOfLight);
}

ElMagField ChargedParticle::RetardedField(int index,
					  Vector ObservationPoint)
{
    Vector SourceX = X[index];
    Vector SourceBeta = P[index];
    Vector SourceBetaPrime = A[index];
    double betagamma2 = SourceBeta.abs2nd();
    double gamma2 = betagamma2 + 1.0;
    double gamma = sqrt(gamma2);
    SourceBeta /= gamma;
    SourceBetaPrime /= gamma;
    // now compute the distance and direction to the observer
    Vector RVec = ObservationPoint - SourceX;
    double R = RVec.norm();
    Vector N = RVec;
    N.normalize();
    // now compute the radiated field
    double scale = Charge*ElementaryCharge/(4.0*Pi*EpsNull);
    double bn3rd = pow(1.0 - dot(SourceBeta, N), 3.0);
    // velocity term
    Vector EField = (N - SourceBeta) / (R*R*bn3rd*gamma2);
    // acceleration term
    EField += cross(N, cross(N - SourceBeta, SourceBetaPrime)) / (R*bn3rd*SpeedOfLight);
    EField *= scale;
    Vector BField = cross(N,EField) / SpeedOfLight;
    return ElMagField(EField, BField);
}

void ChargedParticle::integrateFieldTrace(
    Vector ObservationPoint,
    double t0,
    double dt,
    int nots,
    std::vector<ElMagField> *ObservationField)
{
    cout << "ChargedParticle::integrateFieldTrace()" << endl;
    double tmax=t0+dt*nots;
    // if there is no trajectory of at least one step we do nothing
    if (NP<2) return;
    // loop index i over the time steps of the trajectory
    for (int i=0; i<NP-2; i++)
    {
	double ts1 = RetardedTime(i,ObservationPoint);
	double ts2 = RetardedTime(i+1,ObservationPoint);
	double dts = ts2-ts1;
	ElMagField field;
	if ((ts2>t0) && (ts1<tmax))
	    // else - time step is completely outside observation range
	{
	    ElMagField f1 = RetardedField(i,ObservationPoint);
	    ElMagField f2 = RetardedField(i+1,ObservationPoint);
	    if (ts1<=t0)
		// if (ts2<=tmax) time step is entering the observation range
		// else time step is covering the entire observation range
		// both cases are handled with the same code
	    {
		// observation bucket in which the step end falls
		int idx2 = floor((ts2-t0)/dt);
		if (idx2>nots) idx2=nots;
		// handle all fully covered buckets
		for (int idx=0; idx<idx2; idx++)
		{
		    double t_center = t0 + (idx+0.5)*dt;
		    field = f1*((ts2-t_center)/dts) + f2*((t_center-ts1)/dts);
		    ObservationField->at(idx) += field;
		    cout << " idx=" << idx;
		};
		// handle the last (partially covered) bucket
		if (idx2<nots)
		{
		    double t_start = t0+idx2*dt;
		    ElMagField f_start = f1*((ts2-t_start)/dts) + f2*((t_start-ts1)/dts);
		    field = (f_start+f2)*0.5*((ts2-t_start)/dt);
		    ObservationField->at(idx2) += field;
		    cout << " idx2=" << idx2;
		};
	    }
	    else
		// if (ts2<=tmax) time step is fully inside observation range
		// else time step is leaving the observation range
		// both cases are handled with the same code
	    {
		// observation bucket in which the step start falls
		int idx1 = floor((ts1-t0)/dt);
		// observation bucket in which the step end falls
		int idx2 = floor((ts2-t0)/dt);
		if (idx2>nots) idx2=nots;
		if (idx1==idx2)
		    // time step if fully within one bucket
		{
		    field = (f1+f2)*0.5*(dts/dt);
		    ObservationField->at(idx1) += field;
		    cout << " idx12=" << idx1;
		}
		else
		{
		    // handle the first (partially covered) bucket
		    double t_end = t0+(idx1+1)*dt;
		    ElMagField f_end = f1*((ts2-t_end)/dts) + f2*((t_end-ts1)/dts);
		    field = (f1+f_end)*0.5*((t_end-ts1)/dt);
		    ObservationField->at(idx1) += field;
		    cout << " idx1=" << idx1;
		    // handle all fully covered buckets
		    for (int idx=idx1+1; idx<idx2; idx++)
		    {
			double t_center = t0 + (idx+0.5)*dt;
			field = f1*((ts2-t_center)/dts) + f2*((t_center-ts1)/dts);
			ObservationField->at(idx) += field;
			cout << " idx=" << idx;
		    };
		    // handle the last (partially covered) bucket
		    if (idx2<nots)
		    {
			double t_start = t0+idx2*dt;
			ElMagField f_start = f1*((ts2-t_start)/dts) + f2*((t_start-ts1)/dts);
			field = (f_start+f2)*0.5*((ts2-t_start)/dt);
			ObservationField->at(idx2) += field;
			cout << " idx2=" << idx2;
		    };
		}
	    };
	};
    };
    // cout << endl;
}

