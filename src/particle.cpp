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
#include "SDDS.h"

ChargedParticle::ChargedParticle()
{
    Charge = -1.0;
    Mass = 1.0;
    // no trajectory points
    NP = 0;
}

ChargedParticle::ChargedParticle(double charge, double mass)
{
    Charge = charge;
    Mass = mass;
    // no trajectory points
    NP = 0;
}

ChargedParticle::ChargedParticle(const ChargedParticle* Part)
{
    Charge = Part->Charge;
    Mass = Part->Mass;
    NP = Part->NP;
    if (NP > 0)
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

int ChargedParticle::getNP()
{
    return NP;
}

int ChargedParticle::getCharge()
{
    return Charge;
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

void ChargedParticle::TrackEuler(
    int Nstep,       // number of timesteps
    double tstep,    // time step size
    Vector X0,       // initial position
    Vector P0,       // initial momentum
    GeneralField* lattice)
{
    NP = Nstep + 1;
    Time.clear();
    X.clear();
    P.clear();
    A.clear();
    Time.push_back(0.0);
    X.push_back(X0);
    P.push_back(P0);
    double qm = Charge/Mass * InvRestMass ;    // charge over mass
    double betagamma2 = P0.abs2nd();
    double gamma = sqrt(betagamma2 + 1.0);
    Vector beta = P0 / gamma;
    ElMagField EB = lattice->Field(0.0, X0);
    Vector efield = EB.E();
    Vector bfield = EB.B();
    Vector force = cross(beta, bfield) + efield / SpeedOfLight;
    Vector dP_dt = force * qm;
    // store acceleration
    A.push_back(dP_dt);
    for (int i = 0; i < NP - 1; i++)
    {
	// compute derivatives of the particle motion
	// used to solve the equation of motion inside an undulator
	Vector p = P[i];
	betagamma2 = p.abs2nd();
	gamma = sqrt(betagamma2 + 1.0);
	beta = p / gamma;
	EB = lattice->Field(Time[i], X[i]);
	efield = EB.E();
	bfield = EB.B();
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

void ChargedParticle::TrackVay(
    int Nstep,
    double tstep,
    Vector X0,
    Vector P0,
    GeneralField* lattice)
{
    NP = Nstep + 1;
    Time.clear();
    Time.resize(NP);
    X.clear();
    X.resize(NP);
    P.clear();
    P.resize(NP);
    A.clear();
    A.resize(NP);
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
    ElMagField EB_h = lattice->Field(t_h, x_h);
    Vector E_h = EB_h.E();
    Vector B_h = EB_h.B();
    Vector dp_dt = (cross(beta_h, B_h) + E_h / SpeedOfLight) * qm;
    A[0] = dp_dt;
    // The leap-frog algorithm starts one half-step "before" the initial values
    // velocity p^(i+1) to be computed at the integer time step i+1
    // we store this initial velocity in p_i1 which is copied into p_i when we actually do the time step
    Vector p_i1 = p_h - dp_dt * 0.5 * tstep;
    double gamma_i1 = sqrt(p_i1.abs2nd() + 1.0);
    Vector beta_i1 = p_i1 / gamma_i1;
    for (int i = 0; i < NP; i++)
    {
	// velocity u^i at the integer time step i
	// has been computed in the last time step
	Vector p_i = p_i1;
	Vector beta_i = beta_i1;
	// compute the velocity change over the integer step
	EB_h = lattice->Field(t_h, x_h);
	E_h = EB_h.E();
	B_h = EB_h.B();
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
	Time[i] = t_h;
	X[i] = x_h;
	P[i] = p_h;
	A[i] = dp_dt;
	// compute the change in position and time
	// for the next step
	x_h += beta_i1 * SpeedOfLight * tstep;
	t_h += tstep;
    };
}

/*!
 * Implementation details:
 * ----------------------
 * 
 * The trajectory data is stored for the half-integer time steps of the algorithm.
 * 
 * This initialization routine sets the position and velocity of the first storage
 * point to the given start values and computes the velocity at the end of the
 * time step (half way to the next storage point). The given inital momentum is
 * treated as \f$\mathbf{u}^{i+1/2}\f$ of the algorithm. The computed momentum at the end of
 * the time step is stored in the ChargedParticle::VY_p_i1 variable for use in subsequent time steps.
 * (Note, that instead of \f$\mathbf{u}\f$ we are using \f$\mathbf{p}=\mathbf{u}/c\f$). In addition, we store the relativistic
 * factor belonging to this momentum ChargedParticle::VY_gamma_i1, the charge over mass ratio ChargedParticle::qm
 * and charge over mass divided by the double time step ChargedParticle::qmt2. These values as well
 * as the time step ChargedParticle::dt are constant during the tracking.
 */
void ChargedParticle::InitVay(
    double t0,			// time at the half-step point
	   Vector X0,		// position at the half-step point
	   Vector P0,		// momentum at the half-step point
	   double tstep,
	   GeneralField* field)
{
    // clear all old data
    Time.clear();
    X.clear();
    P.clear();
    A.clear();
    // preset the fixed values
    dt = tstep;
    qm = Charge/Mass * InvRestMass;
    qmt2 = qm * 0.5 * dt;
    // compute the fields at the center position of the step (which is given)
    ElMagField EB_h = field->Field(t0, X0);
    Vector E_h = EB_h.E();
    Vector B_h = EB_h.B();
    // perform the second half-step to compute the
    // velocity at the end of the time step
    Vector p_prime = P0 + E_h / SpeedOfLight * qmt2;
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
    Time.push_back(t0);
    X.push_back(X0);
    P.push_back(P0);
    double gamma_h = sqrt(P0.abs2nd() + 1.0);
    Vector beta_h = P0 / gamma_h;
    A.push_back((cross(beta_h, B_h) + E_h / SpeedOfLight) * qm);
    // after this step exactly one trajectory point exists
    NP = 1;
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
    NP += 1;
}

void ChargedParticle::Translate(Vector R)
{
    for (int i = 0; i < NP; i++)
    {
        X[i] = X[i] + R;
    };
}

void ChargedParticle::MirrorY(double MirrorY)
{
    Charge = -Charge;
    for (int i = 0; i < NP; i++)
    {
        X[i].y = 2.0 * MirrorY - X[i].y;
        P[i].y = -P[i].y;
        A[i].y = -A[i].y;
    };
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

ElMagField ChargedParticle::RetardedField(double time, Vector ObservationPoint)
{
    Vector EField = Vector(0.0, 0.0, 0.0);
    Vector BField = Vector(0.0, 0.0, 0.0);
    
    // return zero if there is no trajectory
    if (NP<1)
    {
	return ElMagField(EField, BField);
    };

    int i1 = 0;    // index of the first trajectory point
    Vector RVec = ObservationPoint - X[i1];
    double R = RVec.norm();
    double t1 = Time[i1] + R / SpeedOfLight;    // retarded observation time
    
    int i2 = NP - 1;    // index of the last trajectory point
    RVec = ObservationPoint - X[i2];
    R = RVec.norm();
    double t2 = Time[i2] + R / SpeedOfLight;    // retarded observation time
    
    // the field is different from zero only if the observation
    // time is within the possible retarded time interval
    if ((time >= t1) && (time <= t2))
    {
	// reduce the interval until the trajectory segment is found
	while (i2 - i1 > 1)
	{
	    int i = (i2 + i1) / 2;
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
	double frac = (time - t1) / (t2 - t1);
	Vector SourceX = X[i1] * (1.0 - frac) + X[i2] * frac;
	Vector SourceBeta = P[i1] * (1.0 - frac) + P[i2] * frac;
	double betagamma2 = SourceBeta.abs2nd();
	double gammasqr = betagamma2 + 1.0;
	double gamma = sqrt(gammasqr);
	SourceBeta /= gamma;
	Vector SourceBetaPrime = A[i1] * (1.0 - frac) + A[i2] * frac;
	SourceBetaPrime /= gamma;
	// now compute the field emitted from an interpolated source point
	RVec = ObservationPoint - SourceX;
	R = RVec.norm();
	Vector N = RVec;
	N.normalize();
	double scale = Charge*ElementaryCharge/(4.0*Pi*EpsNull);
	double bn3rd = pow(1.0 - dot(SourceBeta, N), 3.0);
	// velocity term
	EField += (N - SourceBeta) / (R * R * bn3rd) / (gammasqr);
	// acceleration term
	EField += cross(N, cross(N - SourceBeta, SourceBetaPrime)) / (R*bn3rd*SpeedOfLight);
	EField *= scale;
	BField = cross(N,EField)/SpeedOfLight;
    }
    
    return ElMagField(EField, BField);
}

int ChargedParticle::TimeDomainField(
    Vector ObservationPoint,
    std::vector<double> *ObservationTime,
    std::vector<ElMagField> *ObservationField)
{
    double scale=(Charge*ElementaryCharge/(4.0*Pi*EpsNull));
    // delete all possibly existing data
    ObservationTime->clear();
    ObservationField->clear();
    for( int i=0; i<NP; i++)
    {
	// loop over all source points stored in the trajectory
	Vector SourceX = X[i];
	Vector SourceBeta = P[i];
	Vector SourceBetaPrime = A[i];
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
	// compute the retarded time
	double rtime = Time[i] + R/SpeedOfLight;
	// now compute the radiated field
	double bn3rd = pow(1.0 - dot(SourceBeta, N), 3.0);
	// velocity term
	Vector EField = (N - SourceBeta) / (R*R*bn3rd*gamma2);
	// acceleration term
	EField += cross(N, cross(N - SourceBeta, SourceBetaPrime)) / (R*bn3rd*SpeedOfLight);
	EField *= scale;
	Vector BField = cross(N,EField) / SpeedOfLight;
	// store the data
	ObservationTime->push_back(rtime);
	ObservationField->push_back(ElMagField(EField, BField));
    }
    return NP;
}

int ChargedParticle::WriteSDDS(const char *filename)
{
    cout << "writing SDDS file " << filename << endl;
    SDDS_DATASET data;
    if (1 != SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,filename))
    {
	cout << "ChargedParticle::WriteSDDS - error initializing output\n";
	return 1;
    }
    if  (
	SDDS_DefineSimpleParameter(&data,"NumberTimeSteps","", SDDS_LONG)!=1 || 
	SDDS_DefineSimpleParameter(&data,"Charge","e", SDDS_DOUBLE)!=1 ||
	SDDS_DefineSimpleParameter(&data,"Mass","m_e", SDDS_DOUBLE)!=1
	)
    {
	cout << "ChargedParticle::WriteSDDS - error defining parameters\n";
	return 2;
    }
    if  (
	SDDS_DefineColumn(&data,"t\0","t\0","s\0","TimeInSeconds\0",NULL, SDDS_DOUBLE,0)   ==-1 || 
	SDDS_DefineColumn(&data,"x\0","x\0","m\0","DisplacementInX\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"y\0","y\0","m\0","DisplacementInY\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"z\0","z\0","m\0","DisplacementInZ\0",NULL, SDDS_DOUBLE,0) == -1 || 
	SDDS_DefineColumn(&data,"px\0","px\0",NULL,"Gammabetax\0",NULL, SDDS_DOUBLE,0)== -1 || 
	SDDS_DefineColumn(&data,"py\0","py\0",NULL,"Gammabetay\0",NULL,SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"pz\0","pz\0",NULL,"Gammabetaz\0",NULL,SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"gamma\0","gamma\0",NULL,"RelativisticFactor\0",NULL,SDDS_DOUBLE,0)==-1
	)
    {
	cout << "ChargedParticle::WriteSDDS - error defining data columns\n";
	return 3;
    }
    if (SDDS_WriteLayout(&data) != 1)
    {
	cout << "ChargedParticle::WriteSDDS - error writing layout\n";
	return 4;
    }
    // start a page with number of lines equal to the number of trajectory points
    cout << "SDDS start page" << endl;
    if (SDDS_StartPage(&data,(int32_t)NP) !=1 )
    {
	cout << "ChargedParticle::WriteSDDS - error starting page\n";
	return 5;
    }
    // write the single valued variables
    cout << "SDDS write parameters" << endl;
    if( SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
	"NumberTimeSteps",NP,
	"Charge",Charge,
	"Mass",Mass,
	NULL ) !=1
	)
    {
	cout << "ChargedParticle::WriteSDDS - error setting parameters\n";
	return 6;
    }
    // write the table of trajectory data
    cout << "SDDS writing " << NP << " trajectory points" << endl;
    for( int i=0; i<NP; i++)
    {
	if( SDDS_SetRowValues(&data,
	    SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,
	    "t",Time[i],
	    "x",X[i].x,
	    "y",X[i].y,
	    "z",X[i].z,
	    "px",P[i].x,
	    "py",P[i].y,
	    "pz",P[i].z,
	    "gamma",sqrt(1.0+P[i].abs2nd()),
	    NULL) != 1
	    )
	{
	    cout << "ChargedParticle::WriteSDDS - error writing data columns\n";
	    return 7;
	}
    }
    if( SDDS_WritePage(&data) != 1)
    {
	cout << "ChargedParticle::WriteSDDS - error writing page\n";
	return 8;
    }
    // finalize the file
    if (SDDS_Terminate(&data) !=1 )
    {
	cout << "ChargedParticle::WriteSDDS - error terminating data file\n";
	return 9;
    }	
    // no errors have occured if we made it 'til here
    cout << "writing SDDS done." << endl;
    return 0;
}
