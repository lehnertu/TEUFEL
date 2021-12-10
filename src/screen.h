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
 * \class ScreenObserver
 * \brief Observer of emitted radiation time trace on a gridded screen.
 * @author Ulf Lehnert
 * @date 18.10.2017
 * 
 * This class handles the computation and storage of emitted electromagnetic
 * radiation from different sources (lattice, bunches and beams).
 * For every point on a 2-dimensional screen time-traces of
 * the observed electromagnetic field are generated.
 */
class ScreenObserver : public Observer
{
    
public:
    
    /*! Standard constructor:<br>
     * Compute the electromagnetic field radiated by a given source.
     * 
     * The screen is placed at a given position in space.
     * From there a grid of size (nx,ny) is created
     * extending along the dx and dy vectors.
     * These vectors define the size of one single grid cell.
     * The screen normal is perpendicular to both oriented such
     * that (x,y,-n) form a right-handed coordinate system.
     * I.e. if a screen is seen such that dx points left and dy points up
     * than the normal points at the observer.
     * Only fields that are created from sources placed at
     * positive normal direction are registered by the screen
     * limiting the incidence to the front half plane.
     *
     * If the dx and dy ectors are not perpendicular to each other
     * dx will be modified to fulfill that condition.
     * If nxy are odd then the (nxy-1)/2 indexed grid cell (i.e. center=2 of 5)
     * (index running 0...nxy-1) will have its center exactly
     * at the screen center position.
     * If the index range is even than the first half of the array (i.e. 0...2 of 6) 
     * will be at negative displacements and the first cell
     * of the second half (i.e. center=3 of 6) will be placed exactly at the origin.
     * 
     * The field is recorded in time domain
     * starting at t0 with nots equidistant time steps of dt length.
     * The first time step starts at \f$t0\f$ and ends at \f$t0+dt\f$.
     * The total timespan covererd ends at \f$t0+n*dt\f$.
     * The start time t0 is defined for the center of the grid.
     * The actual start time is cell dependent. It can vary by dtx/dty
     * (time shift per cell index) linearly depending on the position of the cell.
     * The field is computed for the center point of the grid cells only
     * and assumed constant over the whole area of the cell for
     * power flow calculation.
     * 
     * \param[in] filename The name of the generated output file.
     * \param[in] position The center of the screen.
     * \param[in] dx The x direction/spacing of the grid.
     * \param[in] dy The y direction/spacing of the grid.
     * \param[in] nx The number of grid cells in x direction.
     * \param[in] ny The number of grid cells in y direction.
     * \param[in] t0 Start time of the trace.
     * \param[in] dt Time step of the observation trace.
     * \param[in] dtx Start time increment per index in x
     * \param[in] dty Start time increment per index in y
     * \param[in] nots Number of time steps.
     */
    ScreenObserver(
        std::string filename,
        Vector position,
        Vector dx,
        Vector dy,
        unsigned int nx,
        unsigned int ny,
        double t0,
        double dt,
        double dtx,
        double dty,
        unsigned int nots);
    
    //! We need a destructor to release the memory for the field traces
    ~ScreenObserver();
    
    //! Set the source of the fields to be recorded
    virtual void setSource(RadSource s);

    //! Report the source of the fields to be recorded
    virtual RadSource getSource();

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
     *  with index it from the grid cell with indices ix, iy.
     *  This method throws an exception in case of an out-of-range index.
     *  This exception should not be caught as it represents an internal
     *  coding error (should never happen).
     */
    ElMagField getField(
        unsigned int ix,
        unsigned int iy,
        unsigned int it);

    /*! Set one particular field value stored in one time slice
     *  with index it from the grid cell with indices ix, iy.
     *  This method throws an exception in case of an out-of-range index.
     *  This exception should not be caught as it represents an internal
     *  coding error (should never happen).
     */
    void setField(
        unsigned int ix,
        unsigned int iy,
        unsigned int it,
        ElMagField field);

    /*! The method gives the size of the buffer necessary to store
     *  the complete field information as a number of doubles (not bytes!).
     */
    virtual std::size_t getBufferSize() { return Nx*Ny*NOTS*6; };
    
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
    virtual void fromBuffer(double *buffer, std::size_t size);

    /*! @brief Write all the time-domain field traces into an HDF5 file.
     * 
     *  The file name was defined when creating the oberver object.
     * 
     *  Data is divided into 2 data sets
     *  At present 3 data sets are provided to maintain full compatibility with old code.
     * 
     *  "ObservationPosition" :
     *  - 2D array of cell center position (Vector) [m]
     *  - Attributes : Nx, Ny
     * 
     *  "ElMagField" :
     *  - 2D array [Nx,Ny] of field traces
     *    - [NOTS] array of ElMagField
     *      - 3 componenets of the electric field [V/m]
     *      - 3 componenets of the magnetic field [T]
     *  - Atributes : t0_obs, dt_obs, NOTS
     *
     *  new:
     *  all information is given in 2 data sets
     *
     *  1) "Screen" :
     *  - data: 4 vectors (3 doubles each) :
     *    - origin [m] : center position of the screen
     *    - normal [m] : backward normal of the screen
     *    - dx [m] : cell spacing
     *    - dy [m] : cell spacing
     *  - Attributes : 
     *    - Nx : number of cells
     *    - Ny : Number of cells
     *    - NOTS : number of time steps
     *    - t0 : start time [s] of the trace at the center cell
     *    - dt : time step [s] of the field traces
     *    - dtx : start time change [s/pixel] in x direction
     *    - dty : start time change [s/pixel] in y direction
     *   
     *  2) "ElMagField" :
     *  - data: 2D array [Nx,Ny] of field traces
     *    - [NOTS] array of ElMagField
     *      - 3 componenets of the electric field [V/m]
     *      - 3 componenets of the magnetic field [T]
     *  - Attributes :
     *    - none
     * 
     * @throws IOexception
     */
    void WriteTimeDomainFieldHDF5();

    /*! Generate the output file(s) from this observation.
     *  The present code just calls  WriteTimeDomainFieldHDF()
     */
    virtual void generateOutput();

private:
    
    /*! Compute the center position of the grid cell
     *  from its indexes
     */
    Vector CellPosition(unsigned int ix, unsigned int iy);
    
    /*! Compute the trace start time of the grid cell
     *  from its indexes
     */
    double CellTimeZero(unsigned int ix, unsigned int iy);
    
private:
    
    //! file name for the final output
    std::string FileName;

    //! the origin of the grid
    Vector Origin;
    
    //! the x-direction vector of the grid
    Vector dX;
    
    //! the y-direction vector of the grid
    Vector dY;

    //! the normal direction vector of the screen
    Vector normal;
    
    //! the number of grid cells in x direction
    unsigned int Nx;
    
    //! the number of grid cells in x direction
    unsigned int Ny;
    
    //! number of time steps in the interpolated trace
    unsigned int NOTS;
    
    //! the start time of the field trace
    double t0_obs;
    
    //! the time step of the field trace
    double dt_obs;

    //! the start time change per index x
    double dtx_obs;

    //! the start time change per index y
    double dty_obs;
    
    /*! the interpolated electromagnetic field is stored in
     *  a 2-dimensional array of FieldTrace.
     *  The first (outermost) index is ix.
     *  The second index is iy.
     *  We store pointers to the traces for every grid cell.
     */ 
    std::vector<std::vector<FieldTrace*>> Traces;

};
