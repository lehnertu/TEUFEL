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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

#include "beam.h"
#include "bunch.h"
#include "global.h"
#include "fields.h"
#include "observer.h"
#include "vector.h"

/*!
 * \class PointObserver
 * \brief Observer of emitted radiation time trace at a single point in space.
 * @author Ulf Lehnert
 * @date 22.4.2017
 * 
 * This class handles the computation and storage of emitted electromagnetic
 * radiation from different sources (particles, bunches and beams).
 * It provides a time-trace of electromagnetic field values.
 */
class PointObserver : public Observer
{

public:

    /*! Standard constructor:<br>
     * Compute the electromagnetic field radiated by a given source.
     * 
     * The observer is placed at a given position in space.
     * The field is recorded in time domain
     * starting at t0 with nots equidistant time steps of dt length.
     * The first time step starts at \f$t0\f$ and ends at \f$t0+dt\f$.
     * The total timespan covererd ends at \f$t0+n*dt\f$.
     * 
     * This class is templated and spezialized for
     * Bunch(), Beam() or Lattice() being the field source.
     * 
     * \param[in] filename The name of the generated output file.
     * \param[in] t0 Start time of the trace.
     * \param[in] dt Time step of the observation trace.
     * \param[in] nots Number of time steps.
     */
    PointObserver(
        std::string filename,
        Vector position,
        double t0,
        double dt,
        unsigned int nots);

    /*! Integrate the fields emitted by the source
     *  during all of its history, falling onto the time frame
     *  of observation.
     *  Should be called once after tracking all particles.
     * 
     * \param[in] src The source generating the field.
     *
     *  This method is defined for Beam(), Bunch() and Lattice() as field sources.
     */
    virtual void integrate(Beam *src);
    virtual void integrate(Bunch *src);
    virtual void integrate(Lattice *src);

    /*! Return the field value stored in one time slice
     */
    ElMagField getField(unsigned int idx);
    
    /*! Set one particular field value stored in one time slice with index it.
     *  This method throws an exception in case of an out-of-range index.
     *  This exception should not be caught as it represents an internal
     *  coding error (should never happen).
     */
    void setField(
        unsigned int it,
        ElMagField field);

    /*! The method gives the size of the buffer necessary to store
     *  the complete field information as a number of doubles (not bytes!).
     */
    virtual unsigned int getBufferSize();
    
    /*! Return all field values in a newly allocated buffer.
     *  Memory for the buffer is allocated by this method and must be freeed
     *  by the caller. Returns a pointer to the allocated memory.
     *  An exception is thrown if the alloaction of the buffer fails.
     */
    virtual double* getBuffer();

    /*! Set all field values from a given allocated buffer.
     *  The count value gives the size of the buffer as a number of doubles.
     *  An exception is thrown if it doesn't match the actual field size.
     */
    virtual void fromBuffer(double *buffer, unsigned int size);

    /*! @brief Write the time-domain field trace into an SDDS file.
     * 
     * The file name has been defined when creating the observer object.
     * 
     * This method requires that a trace of observed field values has been
     * collected before and stored in PointObserver::TimeDomainField.
     * 
     * The file contains one table with 7 columns
     * - observation time [s]
     * - 3 componenets of the electric field [V/m]
     * - 3 componenets of the magnetic field [T]
     * 
     * @throws IOexception
     */
    void WriteTimeDomainFieldSDDS();

    /*! Generate the output file(s) from this observation.
     *  The present code just calls  WriteTimeDomainFieldSDDS().
     */
    virtual void generateOutput();

private:

    //! file name for the final output
    std::string FileName;

    //! the position of the observer
    Vector Pos;

    //! number of time steps in the interpolated trace
    unsigned int NOTS;

    //! the start of the field trace
    double t0_obs;

    //! the step of the field trace
    double dt_obs;

    //! the interpolated electromagnetic field
    std::vector<ElMagField> TimeDomainField;

};

