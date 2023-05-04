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
#include "vector.h"

/*!
 * \class SourceScreen
 * \brief Use radiation fields recorded on a screen as an external field for tracking.
 * @author Ulf Lehnert
 * @date 14.12.2021
 * 
 * This class handles previously stored electromagnetic on a screen
 * as an external (lattice) field acting on the tracked particles.
 * All data is read from an *.h5 file generated in a previous TEUFEL run
 * or by external scripts (e.g. GaussianPulse.ipynb).
 *
 * It is aasumed that the radiation impinges on the screen essentially
 * from the normal direction with the source residing sufficiently far away.
 * Then all fields on the screen should be transversal, all normal
 * components are canceled. Given that, the normal and time derivatives
 * are related as \f$ \frac{\partial}{\partial n}V = -\frac{1}{c} \frac{\partial}{\partial t}V \f$
 * and the Kirchhoff theorem simplifies to 2 contributions, only.
 *
 * The screen is placed at a given position in space.
 * From there a grid of size (nx,ny) is created
 * extending along the dx and dy vectors.
 * These vectors define the size of one single grid cell.
 * The screen normal is perpendicular to both oriented such
 * that (x,y,-n) form a right-handed coordinate system.
 * I.e. if a screen is seen such that dx points left and dy points up
 * than the normal points at the observer.
 * 
 * The field is recorded in time domain
 * starting at t0 with nots equidistant time steps of dt length.
 * The first time step starts at \f$t0\f$ and ends at \f$t0+dt\f$.
 * The total timespan covererd ends at \f$t0+n*dt\f$.
 * The start time t0 is defined for the center of the grid.
 * The actual start time is cell dependent. It can vary by dtx/dty
 * (time shift per cell index) linearly depending on the position of the cell.
 * 
 * Usually, the radiation on a screen is recorded impinging from
 * the side where the normal vector is pointing to. A SourceScreen
 * propagates the fields further without change in the same direction
 * it was propagating when recorded, that is, in negative normal direction.
 */
class SourceScreen : public LocalizedField
{
    
public:
    
    /*! Standard constructor:<br>
     * Create a source of electromagnetic fields from a screen file.
     *
     * For a documentation of the file format see Screen::WriteTimeDomainFieldHDF5().
     * 
     * The geometry is fully read from the file attributes.
     * The origin of the screen can be shifted in space and time, however.
     * It is placed at the coordinaes given to the constructor (the parent
     * class LocalizedField handles that). In local coordinates the origin
     * of the screen is placed at the origin of the coordinate system and t=0.
     * The position and time information recorded in the file is discarded.
     *
     * \param[in] filename The name of the input file.
     * \param[in] position The new position of the screen.
     * \param[in] start_time The new start time of the screen (at center).
     */
    SourceScreen(
        std::string filename,
        Vector position,
        double time);
    
    //! We need a destructor to release the memory for the field traces
    ~SourceScreen();
    
    /*! The method gives the size of the buffer necessary to store
     *  the field information of either the traces or one of the derivatives
     *  as a number of doubles (not bytes!).
     */
    std::size_t getBufferSize() { return Nx*Ny*NOTS*6; };

    /*! Get the energy flow density vector of one trace
     *  integrated over time
     */
    Vector Poynting(unsigned int ix, unsigned int iy) { return Traces[ix][iy]->Poynting(); };

    /*! Compute the total electromagnetic energy flowing through the screen
     *  Energy flow vectors opposite the normal vector are counted positive.
     */
    double totalEnergy();

private:

	/*! Compute the electromagnetic fields in local coordinates:
     * This is the field in the coordinate system that was used in the file.
     * For field computation LocalizedField::Field() is used which performs the
     * shift of origin in space and time.
     */
    ElMagField LocalField(double t, Vector X) override;

    /*! Compute the center position of the grid cell
     *  from its indexes in local coordinates
     */
    Vector CellPosition(unsigned int ix, unsigned int iy);
    
    /*! Compute the trace start time of the grid cell
     *  from its indexes in local time
     */
    double CellTimeZero(unsigned int ix, unsigned int iy);

    /*! Write one of the fields into a scree-like HDF5 file.
     *  THIS IS FOR DEBUGGING PURPOSE ONLY
     *  @todo remove this
     *  A pointer to the field is given.
     */
    void WriteField(
        std::string filename,
        std::vector<std::vector<FieldTrace*>> Field
    );

    /*! Return all field values in a newly allocated buffer.
     *  Memory for the buffer is allocated by this method and must be freeed
     *  by the caller. Returns a pointer to the allocated memory.
     *  An exception is thrown if the alloaction of the buffer fails.
     *
     *  THIS IS FOR DEBUGGING PURPOSE ONLY
     *  @todo remove this
     */
    double* getBuffer(
        std::vector<std::vector<FieldTrace*>> Field
    );
    
private:
    
    //! file name for the final output
    std::string FileName;

    //! the x-direction vector of the grid
    Vector dX;
    
    //! the absolute value of the x grid spacing
    double dx;

    //! the unit vector of the x direction
    Vector ex;
        
    //! the y-direction vector of the grid
    Vector dY;

    //! the absolute value of the y grid spacing
    double dy;
    
    //! the unit vector of the y direction
    Vector ey;

    //! the normal direction vector of the screen
    Vector normal;
    
    //! the number of grid cells in x direction
    unsigned int Nx;
    
    //! the number of grid cells in x direction
    unsigned int Ny;
    
    //! number of time steps in the interpolated trace
    unsigned int NOTS;
    
    //! the time step of the field trace
    double dt_src;

    //! the start time change per index x
    double dtx_src;

    //! the start time change per index y
    double dty_src;
    
    /*! the electromagnetic field is stored in
     *  a 2-dimensional array of FieldTrace.
     *  The first (outermost) index is ix.
     *  The second index is iy.
     *  We store pointers to the traces for every grid cell.
     */ 
    std::vector<std::vector<FieldTrace*>> Traces;

    /*! the time derivative of the electromagnetic field is stored in
     *  a 2-dimensional array of FieldTrace.
     */ 
    std::vector<std::vector<FieldTrace*>> dt_Traces;

};
