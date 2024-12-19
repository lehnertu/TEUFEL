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
 * \class SnapshotObserver
 * \brief Observer of emitted radiation at a given point in time.
 * @author Ulf Lehnert
 * @date 22.4.2017
 * 
 * This class handles the computation and storage of emitted electromagnetic
 * radiation from different sources (particles, bunches and beams).
 * It provides the electromagnetic fields at a single moment in time
 * observed on a 2-dimensional grid in space.
 */
template <typename SourceT>
class SnapshotObserver : public Observer<SourceT>
{

public:

    /*! Standard constructor:<br>
     * Compute the electromagnetic field observed at a given point in time.
     * 
     * The field is computed on a grid in space.
     * The center of the grid is given with the position argument.
     * From there a grid of size (nx,ny) is created
     * extending along the dx and dy vectors.
     * These vectors define the size of one single grid cell.
     * They need not neccessarily be orthogonal to each other.
     * If nxy are odd then the (nxy-1)/2 indexed grid cell
     * (index running 0...nxy-1) will have its center exactly
     * at the given center position.
     * 
     * \param[in] filename The name of the generated output file.
     * \param[in] position The center of the grid.
     * \param[in] dx The x direction/spacing of the grid.
     * \param[in] dy The y direction/spacing of the grid.
     * \param[in] nx The number of grid cells in x direction.
     * \param[in] ny The number of grid cells in y direction.
     * \param[in] t  Time of observation.
     */
    SnapshotObserver(
        std::string filename,
        Vector position,
        Vector dx,
        Vector dy,
        unsigned int nx,
        unsigned int ny,
        double t);

    /*! Destructor: if necessary delete particleBuffer
     */
    ~SnapshotObserver() override;

    /*! Integrate the fields emitted by all particles of the source
     *  during all of their history, induced a the time and point of observation.
     *  Should be called after tracking all particles at least to the time of observation.
     * 
     * \param[in] src The source generating the field.
     *
     *  This method is defined for Beam(), Bunch() and Lattice() as field sources.
     */
    void integrate() override;
    
    /*! Get the number of particles contributing to the observation
     * 
     * \param[in] src The source generating the field.
     *
     *  This method is defined for Beam() and Bunch() as field sources.
     */
    int getNumParticles();

    /*! Store all particle positions interpolated to the time of observation
     *  \return The number of particles actually written to the buffer
     *  This method is only defined for Beam() as field source.
     *  There are 6 double values per particle (position and momentum=Beta*Gamma) stored into the buffer.
     */
    int storeParticlePositions();
    
    /*! Return the field value stored in one grid cell with indices ix, iy.
     *  This method throws an exception in case of an out-of-range index.
     *  This exception should not be caught as it represents an internal
     *  coding error (should never happen).
     */
    ElMagField getField(
        unsigned int ix,
        unsigned int iy);

    /*! Set the field value stored in one grid cell with indices ix, iy.
     *  This method throws an exception in case of an out-of-range index.
     *  This exception should not be caught as it represents an internal
     *  coding error (should never happen).
     */
    void setField(
        unsigned int ix,
        unsigned int iy,
        ElMagField field);

    /*! The method gives the size of the buffer necessary to store
     *  the complete field information as a number of doubles (not bytes!).
     */
    std::size_t getBufferSize() override { return Nx*Ny*6; };
    
    /*! Return all field values in a newly allocated buffer.
     *  Memory for the buffer is allocated by this method and must be freed
     *  by the caller. Returns a pointer to the allocated memory.
     *  An exception is thrown if the alloaction of the buffer fails.
     */
    double* getBuffer() override;

    /*! Set all field values from a given allocated buffer.
     *  The count value gives the size of the buffer as a number of doubles.
     *  An exception is thrown if it doesn't match the actual field size.
     */
    void fromBuffer(double *buffer, std::size_t size) override;

    /*! @brief Write all the field values into an HDF5 file.
     * 
     *  The file name was defined when creating the oberver object.
     * 
     *  Data is divided into 2 data sets
     * 
     *  "ObservationPosition" :
     *  - 2D array [Nx,Ny] of cell center position (Vector) [m]
     *  - Attributes : Nx, Ny
     * 
     *  "ElMagField" :
     *  - 2D array [Nx,Ny] of ElMagField
     *      - 3 componenets of the electric field [V/m]
     *      - 3 componenets of the magnetic field [T]
     *  - Atributes : t_obs
     * 
     *  If the source of observation consists of particles another data set
     *  containing the particle positions and momenta is generated
     *
     *  "ParticleCoordinates" :
     *  - 2D array [N,6] of double
     *      - 3 components of position [m]
     *      - 3 components of Beta*Gamma
     *  - Attributes : N_particles
     *
     * @throws IOexception
     */
    void WriteFieldHDF5();

    /*! Generate the output file(s) from this observation.
     *  The present code just calls  WriteFieldHDF5()
     */
    void generateOutput() override;

private:
    
    /*! Compute the center position of the grid cell
     *  from its indeces
     */
    Vector CellPosition(int ix, int iy);
    
private:
    
    //! file name for the final output
    std::string FileName;

    //! the origin of the grid
    Vector O;
    
    //! the x-direction vector of the grid
    Vector dX;
    
    //! the y-direction vector of the grid
    Vector dY;

    //! the number of grid cells in x direction
    unsigned int Nx;
    
    //! the number of grid cells in y direction
    unsigned int Ny;
    
    //! the observation time
    double t_obs;
    
    /*! the electromagnetic field is stored in a 2-dimensional vector.
     *  The first (outermost) index is ix.
     *  The second index is iy.
     */ 
    std::vector<std::vector<ElMagField>> FieldArray;

    //! buffer for particle coordinates
    double *particleBuffer;
    int particlesStored;
    
};

