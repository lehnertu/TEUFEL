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

ChargedParticle::ChargedParticle(double *buffer)
{
    Charge = buffer[0];
    Mass = buffer[1];
    dt = 0.0;
    double t = buffer[2];
    Vector x = Vector(buffer[3],buffer[4],buffer[5]);
    Vector p = Vector(buffer[6],buffer[7],buffer[8]);
    Vector a = Vector(buffer[9],buffer[10],buffer[11]);
    initTrajectory(t, x, p, a);
}

ChargedParticle::ChargedParticle(double *buffer, int nTraj)
{
    NP = nTraj;
    double *b = buffer;
    Charge = *b++;
    Mass = *b++;
    dt = *b++;
    qm = *b++;
    qmt2 = *b++;
    double vx = *b++;
    double vy = *b++;
    double vz = *b++;
    VY_p_i1 = Vector(vx,vy,vz);
    VY_gamma_i1 = *b++;
    // now restore the trajectory data
    Time.clear();
    X.clear();
    P.clear();
    A.clear();
    for (int n=0; n<NP; n++)
    {
        double t = *b++;
        Time.push_back(t);
        double xx = *b++;
        double xy = *b++;
        double xz = *b++;
        X.push_back(Vector(xx,xy,xz));
        double px = *b++;
        double py = *b++;
        double pz = *b++;
        P.push_back(Vector(px,py,pz));
        double ax = *b++;
        double ay = *b++;
        double az = *b++;
        A.push_back(Vector(ax,ay,az));
    }
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

void ChargedParticle::serialize(double *buffer)
{
    double *b = buffer;
    *b++ = Charge;
    *b++ = Mass;
    if (NP>0)
    {
        *b++ = Time.back();
        Vector x = X.back();
        Vector p = P.back();
        Vector a = A.back();
        *b++ = x.x;
        *b++ = x.y;
        *b++ = x.z;
        *b++ = p.x;
        *b++ = p.y;
        *b++ = p.z;
        *b++ = a.x;
        *b++ = a.y;
        *b++ = a.z;
    } else {
        *b++ = 0.0;
        *b++ = 0.0;
        *b++ = 0.0;
        *b++ = 0.0;
        *b++ = 0.0;
        *b++ = 0.0;
        *b++ = 0.0;
        *b++ = 0.0;
        *b++ = 0.0;
        *b++ = 0.0;
    }
}

int ChargedParticle::TrajBufSize()
{
    return(9+10*NP);
}

void ChargedParticle::serializeTraj(double *buffer)
{
    double *b = buffer;
    *b++ = Charge;
    *b++ = Mass;
    *b++ = dt;
    *b++ = qm;
    *b++ = qmt2;
    *b++ = VY_p_i1.x;
    *b++ = VY_p_i1.y;
    *b++ = VY_p_i1.z;
    *b++ = VY_gamma_i1;
    if (NP>0)
        for (int n=0; n<NP; n++)
        {
            *b++ = Time[n];
            Vector x = X[n];
            Vector p = P[n];
            Vector a = A[n];
            *b++ = x.x;
            *b++ = x.y;
            *b++ = x.z;
            *b++ = p.x;
            *b++ = p.y;
            *b++ = p.z;
            *b++ = a.x;
            *b++ = a.y;
            *b++ = a.z;
        }
}

double ChargedParticle::TrajTime(int step)
{
    double t = 0.0;
    if (step >= 0 && step < NP)
        t = Time[step];
    else
        throw(IOexception("ChargedParticle::TrajTime - index out of range."));
    return t;
}

Vector ChargedParticle::TrajPoint(int step)
{
    Vector p = Vector(0.0, 0.0, 0.0);
    if (step >= 0 && step < NP)
        p = X[step];
    else
        throw(IOexception("ChargedParticle::TrajPoint - index out of range."));
    return p;
}

Vector ChargedParticle::TrajMomentum(int step)
{
    Vector p = Vector(0.0, 0.0, 0.0);
    if (step >= 0 && step < NP)
        p = P[step];
    else
        throw(IOexception("ChargedParticle::TrajMomentum - index out of range."));
    return p;
}

Vector ChargedParticle::TrajAccel(int step)
{
    Vector a = Vector(0.0, 0.0, 0.0);
    if (step >= 0 && step < NP)
        a = A[step];
    else
        throw(IOexception("ChargedParticle::TrajAccel - index out of range."));
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

void ChargedParticle::Shift(Vector delta, double f_charge)
{
    for (int i=0; i<NP-1; i++)
    {
        Vector position = X[i];
        X[i] = position + delta;
    }
    Charge *= f_charge;
    qm *= f_charge;
    qmt2 *= f_charge;
}

void ChargedParticle::Mirror(Vector origin, Vector normal, double f_charge)
{
    // make sure the normal vector has unit length
    normal.normalize();
    for (int i=0; i<NP-1; i++)
    {
        // normal distance to mirror plane is added 2x to position
        Vector position = X[i];
        Vector dx = normal*dot(origin-position,normal);
        X[i] = position + dx*2.0;
        // normal momentum is reversed
        Vector momentum = P[i];
        Vector dp = normal*dot(momentum,normal);
        P[i] = momentum - dp*2.0;
        // normal acceleration is reversed
        Vector acc = A[i];
        Vector da = normal*dot(acc,normal);
        A[i] = acc - da*2.0;
    }
    Charge *= f_charge;
    qm *= f_charge;
    qmt2 *= f_charge;
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
        throw(RANGEexception("ChargedParticle::InitVay() - must initalize trajectory before tracking"));
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
        *position = Vector(0.0, 0.0, 0.0);
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
        *position = Vector(0.0, 0.0, 0.0);
        *momentum = Vector(0.0, 0.0, 0.0);
    }
}

double ChargedParticle::RetardedTime(int index,
				     Vector ObservationPoint)
{
    double t = TrajTime(index);
    Vector p = TrajPoint(index);
    Vector RVec = ObservationPoint - p;
    double R = RVec.norm();
    return(t + R / SpeedOfLight);
}

ElMagField ChargedParticle::RetardedField(int index,
					  Vector ObservationPoint)
{
    Vector SourceX = TrajPoint(index);
    Vector SourceBeta = TrajMomentum(index);
    Vector SourceBetaPrime = TrajAccel(index);
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

/*!
 * Implementation details:
 * ----------------------
 * 
 * At first we compute relative retardation times at the observation point.
 * This is the time difference between the observation time and
 * the arrival time of a signal from a trajectory point at the
 * observation point. Positive means signal arrives early
 * negative means signal arrives late. If the signals from both ends of the trajectory
 * arrive early, no field can be generated. If both are late the trajectory
 * will be extrapolated to negative times. In the intermediate case a bisection
 * of the trajectory interval is done until the trajectory interval has been reduced
 * to a single time step.
 */
ElMagField ChargedParticle::RetardedField(double obs_time, Vector obs_pos)
{
    // check prerequisites
    if (NP<1) throw(RANGEexception("ChargedParticle::LightConePosition() - must initalize trajectory before calling"));
    // consider the first point of the trajectory.
    int i1 = 0;
    Vector RVec1 = obs_pos - X[i1];
    double R1 = RVec1.norm();
    double t1 = obs_time - (Time[i1] + R1 / SpeedOfLight);
    // consider the last point of the trajectory.
    int i2 = NP-1;
    Vector RVec2 = obs_pos - X[i2];
    double R2 = RVec2.norm();
    double t2 = obs_time - (Time[i2] + R2 / SpeedOfLight);
    // initialize fields as zero
    Vector EField;
    Vector BField;
    // If t2 is positive no signal can be generated
    if (t2<=0.0)
    {
        // these are the quantities we need from the trajectory
        double SourceT;
        Vector SourceX;
        Vector SourceBeta;
        Vector SourceBetaPrime;
        // reletivistic factors to be computed in every case
        double betagamma2, gamma2, gamma;
        if (t1<=0.0)
        {
            // printf("ChargedParticle::RetardedField - extrapolating from zero\n");
            // one has to extrapolate to negative trajectory times
            // if the is only one trajectory point available (i1==i2) we will be here as well
            Vector Rc = (obs_pos - X[0]) / SpeedOfLight;
            SourceBeta = P[0];
            betagamma2 = SourceBeta.abs2nd();
            gamma2 = betagamma2 + 1.0;
            gamma = sqrt(gamma2);
            SourceBeta/=gamma;
            // compute the time t at emission of the signal (field)
            // t is always negative
            double t0 = obs_time;
            double bR = dot(SourceBeta,Rc);
            double bb = dot(SourceBeta,SourceBeta);
            double RR = dot(Rc,Rc);
            SourceT = (t0-bR-sqrt((t0-bR)*(t0-bR)-(1-bb)*(t0*t0-RR)))/(1-bb);
            // compute the trajectory point at emission
            SourceX = X[0] + SourceBeta*SpeedOfLight*SourceT;
            SourceBetaPrime = Vector(0.0,0.0,0.0);
        }
        else
        {
            // printf("ChargedParticle::RetardedField - interpolating within trajectory\n");
            // Bisect the trajectory interval until only a single
            // step is left for interpolation (i2-i2 == 1).
            // There always must be t2<0 and t1>0.
            while (i2-i1 > 1)
            {
                int iMid = (i2+i1)/2;
                Vector RVecMid = obs_pos - X[iMid];
                double RMid = RVecMid.norm();
                double tMid = obs_time - (Time[iMid] + RMid / SpeedOfLight);
                if (tMid > 0.0)
                {
                    i1 = iMid;
                    t1 = tMid;
                }
                else
                {
                    i2 = iMid;
                    t2 = tMid;
                }
            }
            double frac = t1/(t1-t2);
            // printf("i1=%d  i2=%d  frac=%12.9f\n",i1,i2,frac);
            // Simple linear interpolation is not enough.
            // Do a refinement until the interpolaed position changes
            // by less than a 1.0e-5 fraction of the trajectory step.
            //! @todo Linear interpolation fails, as does recursive refinement.<br>
            //! One should do ananalytic solution as it is done for the back-extrapolation.
            //! The present code works at suffucuently large distances from the trajectory.
            double t1 = Time[i1];
            double t2 = Time[i2];
            Vector X1 = X[i1];
            Vector X2 = X[i2];
            SourceT = t1*(1.0-frac) + t2*frac;
            SourceX = X1*(1.0-frac) + X2*frac;
            double R = (obs_pos-SourceX).norm();
            // how much off is the time
            double dt = SourceT + R / SpeedOfLight - obs_time;
            // how does dt change with frac
            // the first refinement step uses the trajectory time step for computing the derivative
            double r1 = (obs_pos-X1).norm();
            double r2 = (obs_pos-X2).norm();
            double dtdfrac = (t2-t1) + (r2-r1)/SpeedOfLight;
            double dfrac = dt/dtdfrac;
            // printf("frac= %12.9f R=%9.6g  SourceT=%9.6g  dt=%9.6g => dfrac=%12.9f\n",frac,R,SourceT,dt,-dfrac);
            int number_ref=1;
            while (fabs(dfrac)>1e-5 && number_ref<10)
            {
                // for subsequent refinement steps we use the change from the last step
                // to compute the drivative
                // before
                double Lastdt = dt;
                // after
                frac -=dfrac;
                SourceT = t1*(1.0-frac) + t2*frac;
                SourceX = X1*(1.0-frac) + X2*frac;
                R = (obs_pos-SourceX).norm();
                dt = SourceT + R / SpeedOfLight - obs_time;
                // how does dt change with frac
                dtdfrac = -(dt-Lastdt)/dfrac;
                dfrac = dt/dtdfrac; 
                // printf("frac=%12.9f  R=%9.6g  SourceT=%9.6g  dt=%9.6g => dfrac=%12.9f\n",frac,R,SourceT,dt,-dfrac);
                number_ref++;
            }
            // repeat with corrected frac
            frac -=dfrac;
            SourceT = t1*(1.0-frac) + t2*frac;
            SourceX = X1*(1.0-frac) + X2*frac;
            if (number_ref>=10)
            {
                printf("WARNING : no convergence in ChargedParticle::RetardedField\n");
                printf("Obs    : t=%9.6g  X=(%9.6g, %9.6g, %9.6g)\n",obs_time,obs_pos.x,obs_pos.y,obs_pos.z);
                printf("Source : t=%9.6g  X=(%9.6g, %9.6g, %9.6g)\n",SourceT,SourceX.x,SourceX.y,SourceX.z);
                R = (obs_pos-SourceX).norm();
                dt = SourceT + R / SpeedOfLight - obs_time;
                printf("frac=%12.9f  R=%9.6g  dt=%9.6g => dfrac=%12.9f\n",frac,R,dt,-dfrac);
            }
            SourceBeta = P[i1]*(1.0-frac) + P[i2]*frac;
            betagamma2 = SourceBeta.abs2nd();
            gamma2 = betagamma2 + 1.0;
            gamma = sqrt(gamma2);
            SourceBeta/=gamma;
            SourceBetaPrime = A[i1]*(1.0-frac) + A[i2]*frac;
            SourceBetaPrime/=gamma;
        }
        // printf("ChargedParticle::RetardedField - SourceX = (%9.6g, %9.6g, %9.6g)\n",SourceX.x,SourceX.y,SourceX.z);
        // compute the distance and direction to the observer
        Vector RVec = obs_pos - SourceX;
        // printf("ChargedParticle::RetardedField - RVec = (%9.6g, %9.6g, %9.6g)\n",RVec.x,RVec.y,RVec.z);
        double R = RVec.norm();
        Vector N = RVec;
        N.normalize();
        // compute the radiated field
        double scale = Charge*ElementaryCharge/(4.0*Pi*EpsNull);
        double bn3rd = pow(1.0 - dot(SourceBeta, N), 3.0);
        // velocity term
        EField = (N - SourceBeta) / (R*R*bn3rd*gamma2);
        // acceleration term
        EField += cross(N, cross(N - SourceBeta, SourceBetaPrime)) / (R*bn3rd*SpeedOfLight);
        EField *= scale;
        BField = cross(N,EField) / SpeedOfLight;
    };
    return ElMagField(EField, BField);
}

void ChargedParticle::integrateFieldTrace(
    Vector ObservationPoint,
    double t0,
    double dt,
    int nots,
    std::vector<ElMagField> *ObservationField)
{
    // cout << "ChargedParticle::integrateFieldTrace() NP=" << NP << endl;
    double tmax=t0+dt*nots;
    // if there is no trajectory of at least one step we do nothing
    if (NP<2) return;
    // loop index i over the time steps of the trajectory
    for (int i=0; i<NP-1; i++)
    {
	double ts1 = RetardedTime(i,ObservationPoint);
	double ts2 = RetardedTime(i+1,ObservationPoint);
	// cout << "time step " << ts1 << " ... " << ts2 << endl;
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
		    // cout << " idx=" << idx;
		};
		// handle the last (partially covered) bucket
		if (idx2<nots)
		{
		    double t_start = t0+idx2*dt;
		    ElMagField f_start = f1*((ts2-t_start)/dts) + f2*((t_start-ts1)/dts);
		    field = (f_start+f2)*0.5*((ts2-t_start)/dt);
		    ObservationField->at(idx2) += field;
		    // cout << " idx2=" << idx2;
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
		    // cout << " idx12=" << idx1;
		}
		else
		{
		    // handle the first (partially covered) bucket
		    double t_end = t0+(idx1+1)*dt;
		    ElMagField f_end = f1*((ts2-t_end)/dts) + f2*((t_end-ts1)/dts);
		    field = (f1+f_end)*0.5*((t_end-ts1)/dt);
		    ObservationField->at(idx1) += field;
		    // cout << " idx1=" << idx1;
		    // handle all fully covered buckets
		    for (int idx=idx1+1; idx<idx2; idx++)
		    {
			double t_center = t0 + (idx+0.5)*dt;
			field = f1*((ts2-t_center)/dts) + f2*((t_center-ts1)/dts);
			ObservationField->at(idx) += field;
			// cout << " idx=" << idx;
		    };
		    // handle the last (partially covered) bucket
		    if (idx2<nots)
		    {
			double t_start = t0+idx2*dt;
			ElMagField f_start = f1*((ts2-t_start)/dts) + f2*((t_start-ts1)/dts);
			field = (f_start+f2)*0.5*((ts2-t_start)/dt);
			ObservationField->at(idx2) += field;
			// cout << " idx2=" << idx2;
		    };
		}
	    };
	};
    };
    // cout << endl;
}

void ChargedParticle::integrateFieldTrace(
        Vector ObservationPoint,
        FieldTrace *trace)
{
    double t0 = trace->get_t0();
    double dt = trace->get_dt();
    int nots = trace->get_N();
    if (nots<2) return;
    // the FieldTrace spans the time from tmin to tmax
    double tmin=t0-dt*0.5;
    double tmax=t0+dt*((double)nots-0.5);
    // if there is no trajectory of at least one step we do nothing
    if (NP<2) return;
    // loop index step over the time steps of the trajectory
    for (int step=0; step<NP-1; step++)
    {
        // handle one trajectory time step ts1...ts2 with span dts at the observation point
	    double ts1 = RetardedTime(step,ObservationPoint);
	    double ts2 = RetardedTime(step+1,ObservationPoint);
	    double dts = ts2-ts1;
	    // field component to be added to the trace fields
	    ElMagField field;
	    if ((ts2>tmin) && (ts1<tmax))
	    // else - time step is completely outside observation range
	    {
	        // these are the fields observed at the ends of the trajectory segment
	        ElMagField f1 = RetardedField(step,ObservationPoint);
	        ElMagField f2 = RetardedField(step+1,ObservationPoint);
	        if (ts1<=tmin)
		    // if (ts2<=tmax) time step starts before but enters the observation range
		    // else time step is covering the entire observation range
		    // both cases are handled with the same code
	        {
		        // observation bucket in which the step end falls
		        int idx2 = round((ts2-t0)/dt);
		        if (idx2>nots) idx2=nots;
		        // handle all fully covered buckets
		        for (int idx=0; idx<idx2; idx++)
		        {
		            // center time of the bucket
		            double t_center = t0+idx*dt;
		            field = f1*((ts2-t_center)/dts) + f2*((t_center-ts1)/dts);
		            trace->add(idx,field);
		        };
		        // handle the last (partially covered) bucket
		        if (idx2<nots)
		        {
		            // start time of the bucket
		            double t_start = t0+(idx2-0.5)*dt;
		            ElMagField f_start = f1*((ts2-t_start)/dts) + f2*((t_start-ts1)/dts);
		            field = (f_start+f2)*0.5*((ts2-t_start)/dt);
		            trace->add(idx2,field);
		        };
	        }
	        else
		    // if (ts2<=tmax) time step is fully inside observation range
		    // else time step is leaving the observation range
		    // both cases are handled with the same code
	        {
		        // observation bucket in which the step start falls
		        int idx1 = round((ts1-t0)/dt);
		        // observation bucket in which the step end falls
		        int idx2 = round((ts2-t0)/dt);
		        if (idx2>nots) idx2=nots;
		        if (idx1==idx2)
		        // time step if fully within one bucket
		        {
		            field = (f1+f2)*0.5*(dts/dt);
		            trace->add(idx1,field);
		        }
		        else
		        {
		            // handle the first (partially covered) bucket
		            double t_end = t0+(idx1+0.5)*dt;
		            ElMagField f_end = f1*((ts2-t_end)/dts) + f2*((t_end-ts1)/dts);
		            field = (f1+f_end)*0.5*((t_end-ts1)/dt);
		            trace->add(idx1,field);
		            // handle all fully covered buckets
		            for (int idx=idx1+1; idx<idx2; idx++)
		            {
			            double t_center = t0+idx*dt;
			            field = f1*((ts2-t_center)/dts) + f2*((t_center-ts1)/dts);
			            trace->add(idx,field);
		            };
		            // handle the last (partially covered) bucket
		            if (idx2<nots)
		            {
			            double t_start = t0+(idx2-0.5)*dt;
			            ElMagField f_start = f1*((ts2-t_start)/dts) + f2*((t_start-ts1)/dts);
			            field = (f_start+f2)*0.5*((ts2-t_start)/dt);
			            trace->add(idx2,field);
		            };
		        };
	        };
	    };
    };
}

