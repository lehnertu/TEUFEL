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

#include <vector>
#include "global.h"
#include "fields.h"
#include "vector.h"

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
     * The field is recorded in time domain
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
	unsigned int nots);

    /*! The source has advanced one time step. Integrate the
     * fields emitted by the source during this latest step.
     */
    void integrate();

    /*! Return the field value stored in one time slice
     */
    ElMagField getField(int idx);
    
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
     * @throws IOexception
     */
    void WriteTimeDomainFieldSDDS(const char *filename);

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

/*!
 * \class ScreenObserver
 * \brief Observer of emitted radiation on a gridded screen.
 * @author Ulf Lehnert
 * @date 18.10.2017
 * 
 * This class handles the computation and storage of emitted electromagnetic
 * radiation from different sources (bunches and beams).
 */
template <class sourceT>
class ScreenObserver
{
    
public:
    
    /*! Standard constructor:<br>
     * Compute the electromagnetic field radiated by a given source.
     * 
     * The screen is placed at a given position in space.
     * From there a grid of size (nx,ny) is created
     * extending along the dx and dy vectors.
     * These vectors define the size of one single grid cell.
     * They need not neccessarily be orthogonal to each other.
     * The screen normal is perpendicular to both oriented such
     * that (x,y,n) form a right-handed coordinate system.
     * If nxy are odd then the (nxy-1)/2 indexed grid cell
     * (index running 0...nxy-1) will have its center exactly
     * at the screen center position.
     * 
     * The field is recorded in time domain
     * starting at t0 with nots equidistant time steps of dt length.
     * The first time step starts at \f$t0\f$ and ends at \f$t0+dt\f$.
     * The total timespan covererd ends at \f$t0+n*dt\f$.
     * The field is computed for the center point of the grid cells only
     * and assumed constant over the whole area of the cell for
     * power flow calculation.
     * 
     * This class is templated and spezialized for
     * Bunch() or Beam() being the field source.
     * 
     * \param[in] src The source generating the field.
     * \param[in] position The center of the screen.
     * \param[in] dx The x direction/spacing of the grid.
     * \param[in] dy The y direction/spacing of the grid.
     * \param[in] nx The number of grid cells in x direction.
     * \param[in] ny The number of grid cells in y direction.
     * \param[in] t0 Start time of the trace.
     * \param[in] dt Time step of the observation trace.
     * \param[in] nots Number of time steps.
     */
    ScreenObserver(
    	sourceT *src,
	    Vector position,
	    Vector dx,
	    Vector dy,
	    unsigned int nx,
	    unsigned int ny,
	    double t0,
	    double dt,
	    unsigned int nots);
    
    /*! The source has advanced one time step. Integrate the
     * fields emitted by the source during this latest step.
     */
    void integrate();
    
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
    unsigned int getCount();
    
    /*! Return all field values in a newly allocated buffer.
     *  Memory for the buffer is allocated by this method and must be freeed
     *  by the caller. Returns a pointer to the allocated memory.
     *  An exception is thrown if the alloaction of the buffer fails.
     */
    double* getBuffer();

    /*! Set all field values from a given allocated buffer.
     *  The count value gives the size of the buffer as a number of doubles.
     *  An exception is thrown if it doesn't match the actual field size.
     */
    void fromBuffer(double *buffer, unsigned int count);

    /*! Write all the time-domain field traces into an HDF5 file.
     *  Data is divided into 2 data sets
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
     * @throws IOexception
     */
    void WriteTimeDomainFieldHDF5(const char *filename);

private:
    
    /*! Compute the center position of the grid cell
     *  from its indeces
     */
    Vector CellPosition(int ix, int iy);
    
private:
    
    //! the field source
    sourceT *Source;
    
    //! the origin of the grid
    Vector O;
    
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
    
    /*! the interpolated electromagnetic field is stored in
     *  a 3-dimensional vector.
     *  The first (outermost) index is ix.
     *  The second index is iy.
     *  The third (innermost) index is the time index.
     */ 
    std::vector<std::vector<std::vector<ElMagField>>> TimeDomainField;

};
