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
    unsigned int nots)
{
    FileName = filename,
    O = position;
    dX = dx;
    dY = dy;
    normal = cross(dX,dY);
    normal.normalize();
    Nx = nx;
    Ny = ny;
    t0_obs = t0;
    dt_obs = dt;
    NOTS = nots;
    // set the field sizes and fill the field with zeros
    ElMagField field;
    TimeDomainField.resize(Nx);
    for (unsigned int ix = 0; ix < Nx; ix++) {
    	TimeDomainField[ix].resize(Ny);
    	for (unsigned int iy = 0; iy < Ny; iy++)
    	{
    	    TimeDomainField[ix][iy].clear();
    	    for (unsigned int i=0; i<NOTS; i++) TimeDomainField[ix][iy].push_back(field);
    	};
    }
}

Vector ScreenObserver::CellPosition(int ix, int iy)
{
    return O + dX*((double)ix-0.5*((double)Nx-1.0)) + dY*((double)iy-0.5*((double)Ny-1.0));
}

void ScreenObserver::integrate(Beam *src)
{
    for (unsigned int ix = 0; ix < Nx; ix++) {
    	for (unsigned int iy = 0; iy < Ny; iy++)
	    {
	        src->integrateFieldTrace(
	    	    CellPosition(ix, iy), t0_obs, dt_obs, NOTS, &TimeDomainField[ix][iy]);
	    };
    }
}

void ScreenObserver::integrate(Bunch *src)
{
    for (unsigned int ix = 0; ix < Nx; ix++) {
        for (unsigned int iy = 0; iy < Ny; iy++)
        {
            src->integrateFieldTrace(
                CellPosition(ix, iy), t0_obs, dt_obs, NOTS, &TimeDomainField[ix][iy]);
        };
    }
}

void ScreenObserver::integrate(Lattice *src)
{
    for (unsigned int ix = 0; ix < Nx; ix++)
        for (unsigned int iy = 0; iy < Ny; iy++)
            for (unsigned int it = 0; it < NOTS; it ++ )
            {
                TimeDomainField[ix][iy][it] += src->Field(t0_obs+it*dt_obs, CellPosition(ix, iy));
            }
}

void ScreenObserver::integrate_mp(Beam *src, unsigned int NumCores, unsigned int CoreId)
{
    for (unsigned int ix = 0; ix < Nx; ix++)
    	for (unsigned int iy = 0; iy < Ny; iy++)
    	    if (CoreId == (ix*Ny+iy) % NumCores)
    	    {
	            src->integrateFieldTrace(
	        	    CellPosition(ix, iy), t0_obs, dt_obs, NOTS, &TimeDomainField[ix][iy]);
	        };
}

void ScreenObserver::integrate_mp(Bunch *src, unsigned int NumCores, unsigned int CoreId)
{
    for (unsigned int ix = 0; ix < Nx; ix++)
        for (unsigned int iy = 0; iy < Ny; iy++)
    	    if (CoreId == (ix*Ny+iy) % NumCores)
            {
                src->integrateFieldTrace(
                    CellPosition(ix, iy), t0_obs, dt_obs, NOTS, &TimeDomainField[ix][iy]);
            };
}

void ScreenObserver::integrate_mp(Lattice *src, unsigned int NumCores, unsigned int CoreId)
{
    for (unsigned int ix = 0; ix < Nx; ix++)
        for (unsigned int iy = 0; iy < Ny; iy++)
    	    if (CoreId == (ix*Ny+iy) % NumCores) {
                for (unsigned int it = 0; it < NOTS; it ++ )
                {
                    TimeDomainField[ix][iy][it] += src->Field(t0_obs+it*dt_obs, CellPosition(ix, iy));
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
	    field = TimeDomainField[ix][iy][it];
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
	   TimeDomainField[ix][iy][it] = field;
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
    herr_t status;
    cout << "writing HDF5 file " << FileName << endl;
    // Create a new file using the default properties.
    hid_t file = H5Fcreate (FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Fcreate()"));
    
    // Create dataspace for the observation positions.
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hsize_t pdims[3];
    pdims[0] = Nx;
    pdims[1] = Ny;
    pdims[2] = 3;
    hid_t pspace = H5Screate_simple (3, pdims, NULL);
    if (pspace<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(pspace)"));
    // buffer the data
    double *buffer = new double[Nx*Ny*3];
    double *bp = buffer;
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
    if (pdcpl<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pcreate(pdcpl)"));
    // Create the dataset.
    hid_t pdset = H5Dcreate(file,
        "ObservationPosition",		// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        pspace, H5P_DEFAULT,
        pdcpl, H5P_DEFAULT);
    if (pdset<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dcreate(pdset)"));
    // Write the data to the dataset
    status = H5Dwrite (pdset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			// mem space id
        pspace,
        H5P_DEFAULT,			// data transfer properties
        buffer);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dwrite(pdset)"));
    // attach scalar attributes
    hid_t atts1  = H5Screate(H5S_SCALAR);
    if (atts1<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(atts1)"));
    hid_t attr1 = H5Acreate2(pdset, "Nx", H5T_NATIVE_INT, atts1, H5P_DEFAULT, H5P_DEFAULT);
    if (attr1<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(attr1)"));
    status = H5Awrite(attr1, H5T_NATIVE_INT, &Nx);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(attr1)"));
    status = H5Sclose (atts1);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(atts1)"));
    hid_t atts2  = H5Screate(H5S_SCALAR);
    if (atts2<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(atts2)"));
    hid_t attr2 = H5Acreate2(pdset, "Ny", H5T_NATIVE_INT, atts2, H5P_DEFAULT, H5P_DEFAULT);
    if (attr2<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(attr2)"));
    status = H5Awrite(attr2, H5T_NATIVE_INT, &Ny);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(attr2)"));
    status = H5Sclose (atts2);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(atts2)"));
    // Close and release resources.
    status = H5Pclose (pdcpl);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pclose(pdcpl)"));
    status = H5Dclose (pdset);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dclose(pdset)"));
    status = H5Sclose (pspace);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dclose(pspace)"));
    delete[] buffer;

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
    buffer = getBuffer();
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
    hid_t atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(t0)"));
    hid_t attr = H5Acreate2(dset, "t0", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
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
    delete[] buffer;

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
