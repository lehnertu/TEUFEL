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

#include "bunch.h"
#include "beam.h"
#include "global.h"
#include "observer.h"
#include "screen.h"

#include <iostream>
#include <math.h>
#include <omp.h>
#include "SDDS.h"
#include "hdf5.h"

ScreenObserver::ScreenObserver(
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
    unsigned int nots)
{
    FileName = filename,
    Origin = position;
    dX = dx;
    dY = dy;
    normal = cross(dX,dY)*-1.0;
    normal.normalize();
    Nx = nx;
    Ny = ny;
    t0_obs = t0;
    dt_obs = dt;
    dtx_obs = dtx;
    dty_obs = dty;
    NOTS = nots;
    if (teufel::rank==0)
    {
        std::cout << "screen observer : " << filename << std::endl;
        std::cout << "allocating " << (double)(6*Nx*Ny*NOTS) * sizeof(double) / 1e6 << " MB of memory per node";
        std::cout << std::endl << std::endl;
    }    

    // set the field sizes and fill the field with zeros
    Traces.resize(Nx);
    for (unsigned int ix = 0; ix < Nx; ix++) {
    	Traces[ix].resize(Ny);
    	for (unsigned int iy = 0; iy < Ny; iy++)
    	{
    	    FieldTrace* tr = new FieldTrace(CellTimeZero(ix,iy),dt_obs,NOTS);
    	    if (tr==0) throw(IOexception("MeshedScreen - error allocating memory."));
    	    Traces[ix][iy] = tr;
    	};
    }

    // set a default source
    source = BeamObservation;
}

ScreenObserver::~ScreenObserver()
{
    for (unsigned int ix = 0; ix < Nx; ix++)
    	for (unsigned int iy = 0; iy < Ny; iy++)
            delete Traces[ix][iy];
}

void ScreenObserver::setSource(RadSource s)
{
    source = s;
}

RadSource ScreenObserver::getSource()
{ 
    return source;
}

Vector ScreenObserver::CellPosition(unsigned int ix, unsigned int iy)
{
    // center index is Nx/2 Ny/2 (integer division!)
    return Origin + dX*(double)((int)ix-(int)Nx/2) + dY*(double)((int)iy-(int)Ny/2);
}

double ScreenObserver::CellTimeZero(unsigned int ix, unsigned int iy)
{
    // center index is Nx/2 Ny/2 (integer division!)
    return t0_obs + dtx_obs*(double)((int)ix-(int)Nx/2) + dty_obs*(double)((int)iy-(int)Ny/2);
}

void ScreenObserver::integrate(Beam *src)
{
    int counter = 0;
    double now = current_time();
    double print_time = now;
    #pragma omp parallel shared(counter)
    {
        for (unsigned int ix = 0; ix < Nx; ix++)
        {
            #pragma omp for
            for (unsigned int iy = 0; iy < Ny; iy++)
            {
                #pragma omp atomic
                counter++;
                // all timing information is already stored in the trace
                src->integrateFieldTrace(CellPosition(ix, iy), Traces[ix][iy]);
            };
            if (omp_get_thread_num() == 0)
            {
                now = current_time();
                if (now-print_time>60.0)
                {
                    print_time = now;
                    std::cout << "node " << teufel::rank << " : computed ";
                    std::cout << counter << " of " << Nx*Ny << " cells";
                    std::cout << " using " << omp_get_num_threads() << " threads." << std::endl;
                };
            };
        }
    }
}

void ScreenObserver::integrate(Bunch *src)
{
    int counter = 0;
    double now = current_time();
    double print_time = now;
    #pragma omp parallel shared(counter)
    {
        for (unsigned int ix = 0; ix < Nx; ix++)
        {
            #pragma omp for
            for (unsigned int iy = 0; iy < Ny; iy++)
            {
                #pragma omp atomic
                counter++;
                // all timing information is already stored in the trace
                src->integrateFieldTrace(CellPosition(ix, iy), Traces[ix][iy]);
            };
            if (omp_get_thread_num() == 0)
            {
                now = current_time();
                if (now-print_time>60.0)
                {
                    print_time = now;
                    std::cout << "node " << teufel::rank << " : computed ";
                    std::cout << counter << " of " << Nx*Ny << " cells";
                    std::cout << " using " << omp_get_num_threads() << " threads." << std::endl;
                };
            };
        };
    }
}

void ScreenObserver::integrate(Lattice *src)
{
    int counter = 0;
    double now = current_time();
    double print_time = now;
    #pragma omp parallel shared(counter)
    {
        for (unsigned int ix = 0; ix < Nx; ix++)
        {
            #pragma omp for
            for (unsigned int iy = 0; iy < Ny; iy++)
            {
                #pragma omp atomic
                counter++;
                for (unsigned int it = 0; it < NOTS; it++)
                {
                    Traces[ix][iy]->add(it,
                        src->Field(CellTimeZero(ix,iy)+(double)it*dt_obs, CellPosition(ix,iy)) );
                }
                if (omp_get_thread_num() == 0)
                {
                    now = current_time();
                    if (now-print_time>60.0)
                    {
                        print_time = now;
                        std::cout << "node " << teufel::rank << " : computed ";
                        std::cout << counter << " of " << Nx*Ny << " cells";
                        std::cout << " using " << omp_get_num_threads() << " threads." << std::endl;
                    };
                };
            };
        };
    };
}

ElMagField ScreenObserver::getField(
	unsigned int ix,
	unsigned int iy,
	unsigned int it)
{
    // default zero
    ElMagField field;
    if ((ix<Nx) && (iy<Ny) && (it<NOTS))
	    field = Traces[ix][iy]->get_field((std::size_t)it);
	else
	    throw(IOexception("ScreenObserver::getField() - requested index out of range."));
    return field;
}

void ScreenObserver::setField(
	unsigned int ix,
	unsigned int iy,
	unsigned int it,
	ElMagField field)
{
    if ((ix<Nx) && (iy<Ny) && (it<NOTS))
	   Traces[ix][iy]->set_field((std::size_t)it, field);
	else
	    throw(IOexception("ScreenObserver::setField() - requested index out of range."));
}

double* ScreenObserver::getBuffer()
{
    double* buffer = new double[getBufferSize()];
    if (buffer==0)
        throw(IOexception("ScreenObserver::getBuffer() - error allocating memory."));
    else
    {
        double* bp = buffer;
        for (unsigned int ix = 0; ix < Nx; ix++)
        	for (unsigned int iy = 0; iy < Ny; iy++)
        	    for (unsigned int it = 0; it < NOTS; it++)
        	    {
            		ElMagField field = getField(ix,iy,it);
            		*bp++ = field.E().x;
            		*bp++ = field.E().y;
            		*bp++ = field.E().z;
            		*bp++ = field.B().x;
            		*bp++ = field.B().y;
            		*bp++ = field.B().z;
        	    };
    };
    return buffer;
}

void ScreenObserver::fromBuffer(double *buffer, std::size_t size)
{
    if (size==Nx*Ny*NOTS*6)
    {
        double *bp = buffer;
        for (unsigned int ix = 0; ix < Nx; ix++)
        	for (unsigned int iy = 0; iy < Ny; iy++)
        	    for (unsigned int it = 0; it < NOTS; it++)
        	    {
                    Vector E, B;
                    E.x = *bp++;
                    E.y = *bp++;
                    E.z = *bp++;
                    B.x = *bp++;
                    B.y = *bp++;
                    B.z = *bp++;
            		setField(ix,iy,it,ElMagField(E,B));
        	    };
    } else {
	    throw(IOexception("ScreenObserver::fromBuffer() - buffer size mismatch."));
    }
}

void ScreenObserver::WriteTimeDomainFieldHDF5()
{
    hid_t atts, attr;
    herr_t status;
    cout << "writing HDF5 file " << FileName << endl;
    // Create a new file using the default properties.
    hid_t file = H5Fcreate (FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Fcreate()"));

// OLD

    // Create dataspace for the observation positions.
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hsize_t pdims[3];
    pdims[0] = Nx;
    pdims[1] = Ny;
    pdims[2] = 3;
    hid_t pspace = H5Screate_simple (3, pdims, NULL);
    if (pspace<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(pspace) - old"));
    // buffer the data
    double *pbuffer = new double[Nx*Ny*3];
    double *bp = pbuffer;
    for (unsigned int ix = 0; ix < Nx; ix++)
	for (unsigned int iy = 0; iy < Ny; iy++)
	{
	    Vector pos = CellPosition(ix, iy);
	    *bp++ = pos.x;
	    *bp++ = pos.y;
	    *bp++ = pos.z;
	};
    // Create the dataset creation property list
    hid_t pdcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (pdcpl<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pcreate(pdcpl) - old"));
    // Create the dataset.
    hid_t pdset = H5Dcreate(file,
        "ObservationPosition",		// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        pspace, H5P_DEFAULT,
        pdcpl, H5P_DEFAULT);
    if (pdset<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dcreate(pdset) - old"));
    // Write the data to the dataset
    status = H5Dwrite (pdset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			// mem space id
        pspace,
        H5P_DEFAULT,			// data transfer properties
        pbuffer);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dwrite(pdset)"));
    // attach scalar attribute Nx
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(Nx) - old"));
    attr = H5Acreate2(pdset, "Nx", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(Nx) - old"));
    status = H5Awrite(attr, H5T_NATIVE_INT, &Nx);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(Nx) - old"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(Nx) - old"));
    // attach scalar attribute Nx
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(Ny) - old"));
    attr = H5Acreate2(pdset, "Ny", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(Ny) - old"));
    status = H5Awrite(attr, H5T_NATIVE_INT, &Ny);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(Ny) - old"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(Ny) - old"));
    // Close and release resources.
    status = H5Pclose (pdcpl);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pclose(pdcpl) - old"));
    status = H5Dclose (pdset);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dclose(pdset) - old"));
    status = H5Sclose (pspace);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dclose(pspace) - old"));
    delete[] pbuffer;

// NEW

    // Create dataspace for the screen metadata.
    // data field 3 vectors
    hsize_t s_dims[2];
    s_dims[0] = 4;
    s_dims[1] = 3;
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hid_t s_space = H5Screate_simple (2, s_dims, NULL);
    if (s_space<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(s_space)"));
    // buffer the data
    double *s_buffer = new double[4*3];
    double *s_bp = s_buffer;
    *s_bp++ = Origin.x;
    *s_bp++ = Origin.y;
    *s_bp++ = Origin.z;
    *s_bp++ = normal.x;
    *s_bp++ = normal.y;
    *s_bp++ = normal.z;
    *s_bp++ = dX.x;
    *s_bp++ = dX.y;
    *s_bp++ = dX.z;
    *s_bp++ = dY.x;
    *s_bp++ = dY.y;
    *s_bp++ = dY.z;
    // Create the dataset creation property list
    hid_t s_dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (s_dcpl<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pcreate(s_dcpl)"));
    // Create the dataset.
    hid_t s_dset = H5Dcreate(file,
        "Screen",		        // dataset name
        H5T_NATIVE_DOUBLE,		// data type
        s_space, H5P_DEFAULT,
        s_dcpl, H5P_DEFAULT);
    if (s_dset<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dcreate(s_dset)"));
    // Write the data to the dataset
    status = H5Dwrite (s_dset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			    // mem space id
        s_space,
        H5P_DEFAULT,			// data transfer properties
        s_buffer);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dwrite(s_dset)"));
    // attach scalar attribute Nx
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(Nx)"));
    attr = H5Acreate2(s_dset, "Nx", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(Nx)"));
    status = H5Awrite(attr, H5T_NATIVE_INT, &Nx);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(Nx)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(Nx)"));
    // attach scalar attribute Nx
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(Ny)"));
    attr = H5Acreate2(s_dset, "Ny", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(Ny)"));
    status = H5Awrite(attr, H5T_NATIVE_INT, &Ny);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(Ny)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(Ny)"));
    // attach scalar attribute NOTS
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(NOTS)"));
    attr = H5Acreate2(s_dset, "NOTS", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(NOTS)"));
    status = H5Awrite(attr, H5T_NATIVE_INT, &NOTS);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(NOTS)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(NOTS)"));
    // attach scalar attribute t0
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(t0)"));
    attr = H5Acreate2(s_dset, "t0", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(t0)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &t0_obs);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(t0)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(t0)"));
    // attach scalar attribute dt
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(dt)"));
    attr = H5Acreate2(s_dset, "dt", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(dt)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dt_obs);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(dt)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(dt)"));
    // attach scalar attribute dtx
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(dtx)"));
    attr = H5Acreate2(s_dset, "dtx", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(dtx)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dtx_obs);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(dtx)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(dtx)"));
    // attach scalar attribute dty
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(dty)"));
    attr = H5Acreate2(s_dset, "dty", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(dty)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dty_obs);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(dty)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(dty)"));
    // Close and release resources.
    status = H5Pclose (s_dcpl);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pclose(s_dcpl)"));
    status = H5Dclose (s_dset);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dclose(s_dset)"));
    status = H5Sclose (s_space);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dclose(s_space)"));
    delete s_buffer;

// end NEW

    // Create dataspace for the field data.
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hsize_t dims[4];
    dims[0] = Nx;
    dims[1] = Ny;
    dims[2] = NOTS;
    dims[3] = 6;
    hid_t space = H5Screate_simple (4, dims, NULL);
    if (space<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(space)"));
    // buffer the data
    double *buffer = getBuffer();
    // Create the dataset creation property list
    hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pcreate(dcpl)"));
    // Create the dataset.
    hid_t dset = H5Dcreate (file,
        "ElMagField", 			// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        space, H5P_DEFAULT,
        dcpl, H5P_DEFAULT);
    if (dset<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dcreate(dset)"));
    // Write the data to the dataset
    status = H5Dwrite (dset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			// mem space id
        space,
        H5P_DEFAULT,			// data transfer properties
        buffer);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dwrite(dset)"));
    // attach scalar attributes
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(t0)"));
    attr = H5Acreate2(dset, "t0", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(t0)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &t0_obs);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(t0)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(t0)"));
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(dt)"));
    attr = H5Acreate2(dset, "dt", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(dt)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dt_obs);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(dt)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(dt)"));
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(NOTS)"));
    attr = H5Acreate2(dset, "NOTS", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(NOTS)"));
    status = H5Awrite(attr, H5T_NATIVE_INT, &NOTS);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(NOTS)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(NOTS)"));
    // Close and release resources.
    status = H5Pclose (dcpl);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pclose(dcpl)"));
    status = H5Dclose (dset);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dclose(dset)"));
    status = H5Sclose (space);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(space)"));
    delete buffer;

    status = H5Fclose (file);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Fclose()"));
    // no errors have occured if we made it 'til here
    cout << "writing HDF5 done." << endl;
    return;
}

void ScreenObserver::generateOutput()
{
    WriteTimeDomainFieldHDF5();
}
