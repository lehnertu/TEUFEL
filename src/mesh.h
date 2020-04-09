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

#include "beam.h"
#include "bunch.h"
#include "global.h"
#include "fields.h"
#include "fieldtrace.h"
#include "observer.h"
#include "vector.h"

/*! data type for tringle references */
struct tri_ref {
    int p1;
    int p2;
    int p3;
};

/*!
 * \class MeshedScreen
 * \brief Observer of emitted radiation time trace on a gridded screen.
 * @author Ulf Lehnert
 * @date 18.10.2017
 * 
 * This class handles the computation and storage of emitted electromagnetic
 * radiation from different sources (bunches and beams).
 * The electromagnetic fields are recorded on a triangulated mesh of arbitrary shape.
 *
 * For every triangle center point time-traces of
 * the observed electromagnetic field are generated.
 * The time traces for every observation point
 * have the same number and value of the time steps.
 * The start time t0 is individual for each trace.
 */
class MeshedScreen : public Observer
{
    
public:
    
    /*! Constructor for MeshedScreen objects from an HDF5 file:<br>
     *
     * The geometry and recording time information of the mesh
     * are read from the given file. The file format is identical to
     * that used in the TimeDomainTHz project (HDF5). When creating the
     * object all field informations that may be present in the file are disacarded.
     * These fields are recreated empty and can be filled by subsequently
     * calculating the received radiation. 
     *
     * For definition of the mesh 5 datasets are read from the file:
     * 
     *   "MeshCornerPoints"
     *      attributes: Ncp
     *      type: array Ncp*3 double
     *      content: point coordinates
     *
     *   "ObservationPosition"
     *      attributes: Np
     *      type: array Np*3 double
     *      content: point coordinates
     *
     *   "MeshTriangles"
     *      type: array Np*3 uint32
     *      content: for every triangle references to the 3 corner points
     *
     *   "ObservationTime"
     *      attributes: Nt, dt
     *      type: array Np double
     *      content: t0 for every trace
     *      It is no error if this dataset is missing - intialize to zero.
     *
     *   "ElMagField"
     *      type: array Np*Nt*6 double
     *      content: fields for every trace
     *      It is no error if this dataset is missing - intialize to zero.
     */
    MeshedScreen(std::string filename);
    
    /*! This method computes a number of internal data
     *  like cell areas, normals and performs some sanity checks.
     *  It should be called whenever a new screen object has been constructed.
     *
     *  A local coordinate system is set up for every mesh cell.
     *  The direction vectors (xi, eta, n) form a right-handed cartesic system.
     */
    void init();
    
    /*! Default destructor:<br>
     */
    ~MeshedScreen();

    /*! Get the position of one grid cell */
    Vector get_point(int ip) { return field_points[ip]; }

    /*! Get the start time of one trace */
    double get_t0(int ip) { return A[ip]->get_t0(); }

    /*! Get the cell-local coordinate system */
    Vector get_xi(int ip) { return xi[ip]; }
    Vector get_eta(int ip) { return eta[ip]; }
    Vector get_normal(int ip) { return normal[ip]; }

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
    
    /*! Integrate the fields emitted by the source
     *  during all of its history, falling onto the time frame
     *  of observation.
     *  Should be called once after tracking all particles.
     *
     *  This version is intended to be run on distributed systems.
     *  The computation is distributed (based on the transverse
     *  indices on the screen) over a number of calls (in parallel on
     *  a number of cores). The fields obtained from all of these
     *  cores have to be added up by an external routine.
     * 
     * \param[in] src The source generating the field.
     * \param[in] NumCores The number of core over which to distribute the computation.
     * \param[in] CoreId The index of the core [0...NumCores-1] for the present call.
     *
     *  This method is defined for Beam(), Bunch() and Lattice() as field sources.
     */
    virtual void integrate_mp(Beam *src, unsigned int NumCores, unsigned int CoreId);
    virtual void integrate_mp(Bunch *src, unsigned int NumCores, unsigned int CoreId);
    virtual void integrate_mp(Lattice *src, unsigned int NumCores, unsigned int CoreId);

    /*! The method gives the size of the buffer necessary to store
     *  the complete field information as a number of doubles (not bytes!).
     */
    virtual std::size_t getBufferSize() { return Np*Nt*6; };
    
    /*! Return all field values in a newly allocated buffer.
     *  Memory for the buffer is allocated by this method and must be freed
     *  by the caller. Returns a pointer to the allocated memory.
     *  An exception is thrown if the alloaction of the buffer fails.
     */
    virtual double* getBuffer();

    /*! Set all field values from a given allocated buffer.
     *  The size value gives the size of the buffer as a number of doubles.
     *  An exception is thrown if it doesn't match the actual field size.
     */
    virtual void fromBuffer(double *buffer, std::size_t size);

    /*! Write all data of the object into an HDF5 file.
     * 
     * @throws IOexception
     */
    void writeFile(std::string filename);

    /*! Write all data of the object into an HDF5 file.
     *  The file name was defined when creating the observer object.
     *  The existing file is overwritten.
     * 
     * @throws IOexception
     */
    void writeFile();

    /*! Compute the total electromagnetic energy flowing through the screen
     *  Energy flow vectors opposite the normal vector are counted positive.
     */
    double totalEnergy();

    /*! Write a report of the screen geometry and parameters
     *  including some summary field data
     *  onto an output stream
     */
    void writeReport(std::ostream *st);

    /*! Generate the output file(s) from this observation.
     *  The present code just calls  WriteTimeDomainFieldHDF()
     */
    virtual void generateOutput();

private:
    
    //! file name for the input and final output
    std::string FileName;

    /*! number of corner points of the grid cells */
    int Ncp;
    
    /*! number of the grid cells / field observation positions */
    int Np;

    /*! number of the time steps per trace */
    int Nt;
    
    /*! time-step value */
    double dt;

    /*! the coordinates of the mesh corner points */
    std::vector<Vector> triangle_points;
    
    /*! lists of references which points belong to the triangles */
    std::vector<tri_ref> triangles;
    
    /*! the coordinates of the mesh center points where fields are defined */
    std::vector<Vector> field_points;

    /*! the local coordinate system of the cells - computed by init() */
    std::vector<Vector> xi;
    std::vector<Vector> eta;
    std::vector<Vector> normal;

    /*! the cell areas - computed by init() */
    std::vector<double> area;
    
    /*! the total screen area - computed by init() */
    double total_area;

    /*! pointers to the field traces for every mesh cell */
    std::vector<FieldTrace*> A;

};
