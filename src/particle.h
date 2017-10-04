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

#pragma once

#include "fields.h"
#include "vector.h"

/*!
 * \class ChargedParticle
 * \brief Particle with electric charge.
 * @author Ulf Lehnert
 * @date 7.4.2017
 * 
 * This class holds the trajectory data of a single charged particle with given mass and charge.
 * The trajectory comprises time [s], position [m], momentum [] and acceleration [1/s] in
 * the laboratory frame.
 * 
 * For tracking this class provides a number of algorithms to perform a single time-step.
 * 
 * All charged particles are surrounded by electromagnetic fields and
 * accelerated particle produce radiation. The fields produced by the particle
 * at a given time and position in laboratory frame can be computed
 * taking into account the full retardation effect.
 * 
 * \todo Introduce particle ID
 * 
 */
class ChargedParticle
{
    
public:
    
    /*! Default constructor:<br>
     * Charge is set to electron charge, mass to electron mass.<br>
     * No trajectory is allocated
     */
    ChargedParticle();
    
    /*! Standard constructor:<br>
     * Charge is is given in units of elementary charge, mass im units of electron mass.<br>
     * No trajectory is allocated.
     */
    ChargedParticle(double charge, double mass);
    
    /*! Copy constructor:<br>
     * A given particle is identically reproduced, thereby creating copies of all trajectory data
     */
    ChargedParticle(const ChargedParticle* Part);
    
    /*! Destructor: free all memory of trajectory storage */
    ~ChargedParticle();
    
    /*! Return the number of points in the trajectory */
    int getNP();
    
    /*! Return the charge of the particle */
    int getCharge();
    
    /*! Return the time [s] for one point of the trajectory */
    double TrajTime(int step);
    
    /*! Return the position [m] for one point of the trajectory */
    Vector TrajPoint(int step);
    
    /*! Return the momentum \f$\beta*\gamma\f$ for one point of the trajectory */
    Vector TrajMomentum(int step);
    
    // track the particle through a given field
    // use the most simple Euler algorithm
    void TrackEuler(int Nstep,       // number of timesteps
		    double tstep,    // time step size
		    Vector X0,       // initial position
		    Vector P0,       // initial momentum
		    GeneralField* field);
    
    /*! @brief Track the particle through a given field.
     * 
     * The algorithm follows J.-L.Vay PHYSICS OF PLASMAS 15, 056701 (2008).
     * 
     * @param Nstep the number of tracking steps to be performed.
     * @param tstep the length of the time step.
     * @param X0 the start position of the particle
     * @param P0 the start momentum \f$\beta \gamma\f$ of the particle
     * @param field the field through which the particle is tracked
     */
    void TrackVay(int Nstep,       // number of timesteps
		  double tstep,    // time step size
		  Vector X0,       // initial position
		  Vector P0,       // initial momentum
		  GeneralField* field);
    
    /*! @brief Setup for tracking the particle using the Vay algorithm.
     * 
     * The algorithm follows J.-L.Vay PHYSICS OF PLASMAS 15, 056701 (2008).
     * 
     * @param t0 the start time
     * @param X0 the start position of the particle
     * @param P0 the start momentum \f$\beta \gamma\f$ of the particle
     * @param tstep the length of the time step.
     * @param field the field through which the particle is tracked
     */
    void InitVay(double t0,
		 Vector X0,
		 Vector P0,
		 double tstep,
		 GeneralField* field);

    /*! @brief Perform one tracking step using the Vay algorithm.
     * 
     * The algorithm follows J.-L.Vay PHYSICS OF PLASMAS 15, 056701 (2008).
     * 
     * @param field the field through which the particle is tracked
     */
    void StepVay(GeneralField* field);
    
    //! translate a given particle trajectory
    void Translate(Vector R);
    
    /*! mirroring a given particle trajectory
     *  on a plane y=MirrorY
     *  (this also reverses the charge)
     */
    void MirrorY(double MirrorY);
    
    /*! Particle phase space position at a given time.
     *  Interpolates between time steps after tracking.
     *  Can also be used when the tracking is only initalized to
     *  retreive the inital particle position at t=0
     */
    void CoordinatesAtTime(double time, Vector *position, Vector *momentum);
    
    /*! Electromagnetic field radiated by the particle
     * at a given observation point at time t in lab frame
     * retardation is properly accounted for
     * 
     * \param[in] time absolute time in s.
     * \param[in] ObservationPoint The position [m] of the observer.
     */
    ElMagField RetardedField(double time, Vector ObservationPoint);
    
    /*! Compute the electromagnetic field radiated by the particle
     * seen at the observation point. The field is given in time domain
     * starting at t0 with NOTS equidistant time steps of dt length.
     * The first time step starts at \f$t0\f$ and ends at \f$t0+dt\f$.
     * So the center (reference) time for the first sample is \f$t0+dt/2\f$.
     * The total timespan covererd ends at \f$t0+n*dt\f$.
     * 
     * This method first uses ChargedParticle::FieldTrace()
     * to obtain the non-interpolated data covering the requested range.
     * A step-wise linear funtion is used to interpolate from the non-equidistant
     * time steps stored in the particle trajectory, thereby, exactly preserving
     * the (first order) \f$\int E dt\f$ integral.
     * 
     * \param[in] ObservationPoint The position [m] of the observer.
     * \param[in] t0 Start time [s] of the trace (at the observer).
     * \param[in] dt Duration of one time segment [s] of the observation trace.
     * \param[in] nots Number of time segments.
     * \param[out] ObservationField Electromagnetic field integrals at the observation point.
     */
    void getTimeDomainField(
	Vector ObservationPoint,
	double t0,
	double dt,
	int nots,
	std::vector<ElMagField> *ObservationField);

    /*! Write all information in particular the trajectory data into an SDDS file.
     * 
     * returns values for error checks:
     *	 
     *	0  -  file successfully written \n
     *	1  -  error in SDDS_InitializeOutput \n
     *	2  -  error in SDDS_DefineSimpleParameter \n
     *	3  -  error in SDDS_DefineColumn \n
     *	4  -  error in SDDS_WriteLayout \n
     *	5  -  error in SDDS_StartPage \n
     *	6  -  error in SDDS_SetParameters \n
     *	7  -  error in SDDS_SetRowValues \n
     *	8  -  error in SDDS_WritePage \n
     *	9  -  error in SDDS_Terminate \n
     * 
     */
    int WriteTrajectorySDDS(
	const char *filename);
    
    /*! Write the time-domain field trace into an SDDS file.
     * 
     * All file writing methods should normally be implemented
     * for an observer. This one is an exception because the
     * single (indexed) field traces should not be public.
     * 
     * This method first computes the observation times and fields
     * for the given observation point.
     * These are written to an SDDS file with given name.
     * 
     * The file contains one table with 7 columns
     * - observation time [s]
     * - 3 componenets of the electric field [V/m]
     * - 3 componenets of the magnetic field [T]
     * 
     * @return values for error checks:
     *	
     * 10  -  no data \n
     *	0  -  successfully Written the file \n
     *	1  -  error in SDDS_InitializeOutput \n
     *	2  -  error in SDDS_DefineSimpleParameter \n
     *	3  -  error in SDDS_DefineColumn \n
     *	4  -  error in SDDS_WriteLayout \n
     *	5  -  error in SDDS_StartPage \n
     *	6  -  error in SDDS_SetParameters \n
     *	7  -  error in SDDS_SetRowValues \n
     *	8  -  error in SDDS_WritePage \n
     *	9  -  error in SDDS_Terminate \n
     * 
     */
    int WriteFieldTraceSDDS(
	Vector ObservationPoint,
	const char *filename);

private:
    
    //! number of trajectory points
    int NP;
    
    //! charge in units of ElementaryCharge
    double Charge;
    
    //! mass in unit of the electron rest mass
    double Mass;
    
    //! time in lab frame [s]
    std::vector<double> Time;
    
    //! position in lab frame [m]
    std::vector<Vector> X;
    
    //! dimensionless momentum in lab frame : \f$\beta\gamma = \frac{c p}{mc^2}\f$
    std::vector<Vector> P;
    
    //! acceleration [1/s] in lab frame \f$a = \frac{d}{dt}(\beta\gamma)\f$
    std::vector<Vector> A;
    
    //! time step for tracking - this will remain constant after being set at the start of tracking
    double dt;
    
    //! charge over mass ratio - this will remain constant during tracking
    double qm;

    //! charge over mass divided by the double time step - this will remain constant during tracking
    double qmt2;
    
    //! \f$u^{i+1}/c\f$ of the Vay algorithm
    Vector VY_p_i1;
    
    //! \f$\gamma^{i+1}/c\f$ of the Vay algorithm
    double VY_gamma_i1;
    
    /*! Compute the time when a signal emiited from the trajectory point
     * with given index is observed at the ObservationPoint.
     * 
     * watch out! - validity of the index is not checked
     */
    double RetardedTime(int index,
	Vector ObservationPoint);
	
    /*! Compute the field emiited from the trajectory point
     * with given index observed at the ObservationPoint.
     * 
     * watch out! - validity of the index is not checked
     * 
     * This method is private - no methods with direct indexing of
     * the trajectory points should be outside-visible. For the outside
     * world there exists a method with the same name but the
     * observation time as a parameter.
     */
    ElMagField RetardedField(int index,
	Vector ObservationPoint);

    /*! Compute the electromagnetic field radiated by a particle
     * seen from a given observation point. The field is given in time domain
     * at a number of time steps corresponding to the time steps of
     * the trajectory of the particle. The sample points at the observation
     * position, thus, are not predefined and not equi-distant. The samples
     * are delayed by the retardation corresponding to the observation distance.
     * Only those samples spanning (t1, t2) in time are delivered
     * including one sample before t1 and one sample after t2 if available.
     * It is not guarantied, that the whole span is covered in case of
     * time mismatches also 0 samples can be returned.
     * 
     * \return
     * Observation time and field are stored into given vectors.
     * The number of samples is returned.
     * 
     * \param[in] ObservationPoint Position [m] of the observer in space.
     * \param[in] t1 Start [s] of the range of interest in time.
     * \param[in] t2 End [s] of the range of interest in time.
     * \param[out] ObservationTime Sample time at the observation point.
     * \param[out] ObservationField Electromagnetic field samples at the observation point.
     */
    int FieldTrace(
	Vector ObservationPoint,
	double t1,
	double t2,
	std::vector<double> *ObservationTime,
	std::vector<ElMagField> *ObservationField);

    };
