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
    PreviousTime = 0.0;
    Time = 0.0;
    PreviousX = Vector(0.0,0.0,0.0);
    X = Vector(0.0,0.0,0.0);
    PreviousP = Vector(0.0,0.0,0.0);
    P = Vector(0.0,0.0,0.0);
    PreviousA = Vector(0.0,0.0,0.0);
    A = Vector(0.0,0.0,0.0);
    dt = 0.0;
}

ChargedParticle::ChargedParticle(double charge, double mass)
{
    Charge = charge;
    Mass = mass;
    PreviousTime = 0.0;
    Time = 0.0;
    PreviousX = Vector(0.0,0.0,0.0);
    X = Vector(0.0,0.0,0.0);
    PreviousP = Vector(0.0,0.0,0.0);
    P = Vector(0.0,0.0,0.0);
    PreviousA = Vector(0.0,0.0,0.0);
    A = Vector(0.0,0.0,0.0);
    dt = 0.0;
}

ChargedParticle::ChargedParticle(const ChargedParticle *origin)
{
    Charge = origin->Charge;
    Mass = origin->Mass;
    PreviousTime = origin->PreviousTime;
    Time = origin->Time;
    PreviousX = origin->PreviousX;
    X = origin->X;
    PreviousP = origin->PreviousP;
    P = origin->P;
    PreviousA = origin->PreviousA;
    A = origin->A;
    dt = origin->dt;
    qm = origin->qm;
    qmt2 = origin->qmt2;
    VY_p_i1 = origin->VY_p_i1;
    VY_gamma_i1 = origin->VY_gamma_i1;
}

ChargedParticle::~ChargedParticle()
{
}

int ChargedParticle::getCharge()
{
    return Charge;
}

double ChargedParticle::getTime()
{
    return Time;
}

void ChargedParticle::setTime(double t)
{
    Time = t;
}
    
Vector ChargedParticle::getPosition()
{
    return X;
}

void ChargedParticle::setPosition(Vector x)
{
    X = x;
}

Vector ChargedParticle::getMomentum()
{
    return P;
}

void ChargedParticle::setMomentum(Vector p)
{
    P = p;
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
    // clear the one-before trajectory information
    PreviousTime = Time;
    PreviousX = X;
    PreviousP = P;
    PreviousA = A;
    // preset the fixed values
    dt = tstep;
    qm = Charge/Mass * InvRestMass;
    qmt2 = qm * 0.5 * dt;
    // compute the fields at the center position of the step (which is given)
    ElMagField EB_h = field->Field(Time, X);
    Vector E_h = EB_h.E();
    Vector B_h = EB_h.B();
    // perform the second half-step to compute the
    // velocity at the end of the time step
    Vector p_prime = P + E_h / SpeedOfLight * qmt2;
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
    double gamma_h = sqrt(P.abs2nd() + 1.0);
    Vector beta_h = P / gamma_h;
    A = (cross(beta_h, B_h) + E_h / SpeedOfLight) * qm;
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
    // keep the one-before trajectory information
    PreviousTime = Time;
    PreviousX = X;
    PreviousP = P;
    PreviousA = A;
    // the momentum vector at the beginning of the step is known from the previous step
    Vector p_i = VY_p_i1;
    Vector beta_i = p_i / VY_gamma_i1;
    // use eq. (3) to compute the new center position
    Time = PreviousTime + dt;
    X = X + beta_i * SpeedOfLight * dt;
    // compute the fields at the center position of the step
    ElMagField EB_h = field->Field(Time, X);
    Vector E_h = EB_h.E();
    Vector B_h = EB_h.B();
    // perform the first half-step to compute the
    // velocity at the center (half) point using eq. (13)
    P = VY_p_i1 + (E_h/SpeedOfLight + cross(beta_i, B_h)) * qmt2;
    double gamma_h = sqrt(P.abs2nd() + 1.0);
    Vector beta_h = P / gamma_h;
    // perform the second half-step to compute the
    // velocity at the end of the time step
    Vector p_prime = P + E_h / SpeedOfLight * qmt2;
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
    A = (cross(beta_h, B_h) + E_h / SpeedOfLight) * qm;
}

double ChargedParticle::PreviousRetardedTime(Vector ObservationPoint)
{
    Vector RVec = ObservationPoint - PreviousX;
    double R = RVec.norm();
    return(PreviousTime + R / SpeedOfLight);
}

double ChargedParticle::RetardedTime(Vector ObservationPoint)
{
    Vector RVec = ObservationPoint - X;
    double R = RVec.norm();
    return(Time + R / SpeedOfLight);
}

ElMagField ChargedParticle::PreviousRetardedField(Vector ObservationPoint)
{
    Vector SourceX = PreviousX;
    Vector SourceBeta = PreviousP;
    Vector SourceBetaPrime = PreviousA;
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

ElMagField ChargedParticle::RetardedField(Vector ObservationPoint)
{
    Vector SourceX = X;
    Vector SourceBeta = P;
    Vector SourceBetaPrime = A;
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

ElMagObs ChargedParticle::PreviousObservation(Vector ObservationPoint)
{
    double t = PreviousRetardedTime(ObservationPoint);
    ElMagField field = PreviousRetardedField(ObservationPoint);
    ElMagObs obs = ElMagObs(t, field.E(), field.B());
    return obs;
}

ElMagObs ChargedParticle::Observation(Vector ObservationPoint)
{
    double t = RetardedTime(ObservationPoint);
    ElMagField field = RetardedField(ObservationPoint);
    ElMagObs obs = ElMagObs(t, field.E(), field.B());
    return obs;
}


























