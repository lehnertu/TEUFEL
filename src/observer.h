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

#include <complex>
#include "bunch.h"
#include "vector.h"
#include "particle.h"

/*
 * \class Observation
 * \brief Abstract class defining an observation of the electromagnetic
 * field created by an object (particle, bunch or beam).
 * @author Ulf Lehnert
 * @date 21.9.2017
 * 
 * template<class Obj> class Observation
 * {
 * }
 */

/*!
 * \class PointObserver
 * \brief Observer of emitted radiation at a single point in space.
 * @author Ulf Lehnert
 * @date 22.4.2017
 * 
 * This class handles the computation and storage of emitted electromagnetic
 * radiation from particles or ensembles of particles.
 */
class PointObserver
{

public:

    /*! Standard constructor:<br>
     * The observer is placed at a given point in space.
     */
    PointObserver(Vector position);

    /*! Compute the electromagnetic field radiated by a source
     * seen at the observation point. The field is given in time domain
     * starting at t0 with NOTS equidistant time steps of dt length.
     * A step-wise linear funtion is used to interpolate from the non-equidistant
     * time steps stored in the particle trajectory, thus, exactly preserving
     * the (first order) \f$\int E dt\f$ integral. The first time step
     * starts at \f$t0\f$ and ends at \f$t0+dt\f$.
     * The total timespan covererd ends at \f$t0+n*dt\f$.
     * 
     * This method is templated and spezialized for ChargedParticle()
     * Bunch() or Beam() being the field source.
     * 
     * The array PointObserver::TimeDomainField is cleared
     * and filled with the new values for the interpolated output field data.
     * 
     * \param[in] source The particle generating the field.
     * \param[in] t0 Start time of the trace.
     * \param[in] dt Time step of the observation trace.
     * \param[in] nots Number of time steps.
     */
    template <class Src>
    void ComputeTimeDomainField(
	Src *source,
	double t0,
	double dt,
	int nots);

    /*! Complex amplitude of the observation for a given frequency [Hz].
     * 
     * The 3 components of the electric field are analyzed.
     * 
     * This method requires that the a trace of observed field values has been
     * collected before and stored in PointObserver::TimeDomainField
     * Such a trace can be used for multiple
     * calls of the method in sequence - it is not altered.
     */
    void FrequencyObservation(
	double freq,
	std::complex<double> *Ex,
	std::complex<double> *Ey,
	std::complex<double> *Ez
    );

    /*! Write the time-domain field trace into an SDDS file.
     * 
     * This method requires that a trace of observed field values has been
     * collected before and stored in PointObserver::TimeDomainField.
     * 
     * The file contains one table with 7 columns
     * - observation time [s]
     * - 3 componenets of the electric field [V/m]
     * - 3 componenets of the magnetic field [T]
     * 
     * @return values for error checks:
     *	 
     *	0  -  successfully Written the file\n
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
    int WriteTimeDomainFieldSDDS(const char *filename);
    
    /*! Write a frequency spectrum into an SDDS file.
     * 
     * This method calls FrequencyObservation() and therefore requires
     * that the a trace of observed field values has been
     * collected before and stored in PointObserver::TimeDomainField
     * 
     * The file contains one table with 7 columns
     * - observation Frequency f [Hz]
     * - Amplitude Ax and phase Px of the x component of the electric field
     * - Amplitude Ay and phase Py of the y component of the electric field
     * - Amplitude Az and phase Pz of the z component of the electric field
     * 
     * @param filename name of the SDDS file to be created
     * @param fstart first frequency [Hz] to write into the file
     * @param fstop last frequency [Hz] to write into the file
     * @param fstep frequency step [Hz] of the list
     * 
     * @return values for error checks:
     *	 
     *	0  -  successfully Written the file\n
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
    int WriteSpectrumSDDS(
	const char *filename,
	double fstart,
	double fstop,
	double fstep);
 
private:

    //! the position of the observer
    Vector Pos;

    //! number of time steps in the interpolated trace
    int NOTS;

    //! the start of the field trace
    double t0_obs;

    //! the step of the field trace
    double dt_obs;

    //! the interpolated electromagnetic field
    std::vector<ElMagField> TimeDomainField;

};

/*!
 * \class ScreenObserver
 * \brief Observer of emitted radiation on a rectangular planar grid.
 * @author Ulf Lehnert
 * @date 17.9.2017
 * 
 * This class handles the computation and storage of emitted electromagnetic
 * radiation from particles or ensembles of particles.
 */
class ScreenObserver
{
    
public:
    
    /*! Standard constructor:<br>
     * The center of the screen is placed at a given origin in space.
     * From this center a grid is constructed spaning the given length
     * in two directions. (The two directions need not to be perpendicular
     * to each other.) The two span directions are subdivided by
     * nx/ny grid points (minimum 2). Only if both subdivisions are
     * odd numbers the origin also is one of the observation grid points.
     */
    ScreenObserver(
	Vector origin,
	Vector span_x,
	Vector span_y,
	unsigned int nx,
	unsigned int ny);
};

/*!
 * \class VolumeObserver
 * \brief Observer of emitted radiation on a rectangular volume grid.
 * @author Ulf Lehnert
 * @date 17.9.2017
 * 
 * This class handles the computation and storage of emitted electromagnetic
 * radiation from particles or ensembles of particles.
 */
class VolumeObserver
{
};