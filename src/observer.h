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
#include "fields.h"
#include "vector.h"

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
 * radiation from different sources (bunches and beams).
 */
template <class sourceT>
class PointObserver
{

public:

    /*! Standard constructor:<br>
     * Compute the electromagnetic field radiated by a given source.
     * 
     * The observer is placed at a given position in space.
     * The field is given in time domain
     * starting at t0 with nots equidistant time steps of dt length.
     * The first time step starts at \f$t0\f$ and ends at \f$t0+dt\f$.
     * The total timespan covererd ends at \f$t0+n*dt\f$.
     * 
     * This class is templated and spezialized for
     * Bunch() or Beam() being the field source.
     * 
     * \param[in] src The source generating the field.
     * \param[in] t0 Start time of the trace.
     * \param[in] dt Time step of the observation trace.
     * \param[in] nots Number of time steps.
     */
    PointObserver(
	sourceT *src,
	Vector position,
	double t0,
	double dt,
	int nots);

    /* The source has advanced one time step. Integrate the
     * fields emitted by the source during this latest step.
     */
    void integrate();

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

private:

    //! the field source
    sourceT *Source;
    
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
