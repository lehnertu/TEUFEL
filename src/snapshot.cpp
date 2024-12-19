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
#include "snapshot.h"

#include <iostream>
#include <math.h>
#include <omp.h>
#include "SDDS.h"
#include "hdf5.h"

template <typename SourceT>
SnapshotObserver<SourceT>::SnapshotObserver(
    std::string filename,
    Vector position,
    Vector dx,
    Vector dy,
    unsigned int nx,
    unsigned int ny,
    double t)
{
    FileName = filename,
    O = position;
    dX = dx;
    dY = dy;
    Nx = nx;
    Ny = ny;
    t_obs = t;

    // set the field sizes and fill the field with zeros
    ElMagField field;
    FieldArray.resize(Nx);
    for (unsigned int ix = 0; ix < Nx; ix++)
    	for (unsigned int iy = 0; iy < Ny; iy++)
    	    FieldArray[ix].push_back(field);
    // set a default source
    this->source_type = BeamObservation;
    // set a NULL pointer indicating no particle coordinates stored
    particleBuffer = 0;
    particlesStored = 0;
}

template <typename SourceT>
SnapshotObserver<SourceT>::~SnapshotObserver()
{
    if (particlesStored>0) delete particleBuffer;
}

template <>
void SnapshotObserver<Beam>::integrate()
{
    #pragma omp parallel for collapse(2)
    for (unsigned int ix = 0; ix < Nx; ix++) {
    	for (unsigned int iy = 0; iy < Ny; iy++) {
	        FieldArray[ix][iy] = source->RetardedField(t_obs, CellPosition(ix, iy));
	    }
	}
}

template <>
void SnapshotObserver<Bunch>::integrate()
{
    #pragma omp parallel for collapse(2)
    for (unsigned int ix = 0; ix < Nx; ix++) {
    	for (unsigned int iy = 0; iy < Ny; iy++) {
	        FieldArray[ix][iy] = source->RetardedField(t_obs, CellPosition(ix, iy));
	    }
	}
}

template <>
void SnapshotObserver<Lattice>::integrate()
{
    #pragma omp parallel for collapse(2)
    for (unsigned int ix = 0; ix < Nx; ix++) {
    	for (unsigned int iy = 0; iy < Ny; iy++) {
	        FieldArray[ix][iy] = source->Field(t_obs, CellPosition(ix, iy));
	    }
	}
}

template <>
int SnapshotObserver<Beam>::getNumParticles()
{
    if (source)
        return source->getNOP();
    else
        return 0;
}

template <>
int SnapshotObserver<Beam>::storeParticlePositions()
{
    Vector X, BG;
    int NOP = getNumParticles();
    particleBuffer = new double[6*NOP];
    if (particleBuffer==0)
    {
        throw(IOexception("SnapshotObserver::storeParticlePositions() - error allocating memory."));
        particlesStored = 0;
    }
    double *bp = particleBuffer;
    if (source)
    {
        int NOB = source->getNOB();
        for (int bunch_index=0; bunch_index<NOB; bunch_index++)
        {
            Bunch *b = source->getBunch(bunch_index);
            int b_NOP = b->getNOP();
            for (int i = 0; i < b_NOP; i++)
            {
                ChargedParticle* p = b->getParticle(i);
                p->CoordinatesAtTime(t_obs, &X, &BG);
                *bp++ = X.x;
                *bp++ = X.y;
                *bp++ = X.z;
                *bp++ = BG.x;
                *bp++ = BG.y;
                *bp++ = BG.z;
            }
        }
        particlesStored = NOP;
    }
    else
        particlesStored = 0;
    return particlesStored;
}

template <typename SourceT>
ElMagField SnapshotObserver<SourceT>::getField(
        unsigned int ix,
        unsigned int iy)
{
    // default zero
    ElMagField field;
    if ((ix<Nx) && (iy<Ny))
	    field = FieldArray[ix][iy];
	else
	    throw(IOexception("SnapshotObserver::getField() - requested index out of range."));
    return field;
}

template <typename SourceT>
void SnapshotObserver<SourceT>::setField(
        unsigned int ix,
        unsigned int iy,
        ElMagField field)
{
    if ((ix<Nx) && (iy<Ny))
	   FieldArray[ix][iy] = field;
	else
	    throw(IOexception("SnapshotObserver::setField() - requested index out of range."));
}

template <typename SourceT>
double* SnapshotObserver<SourceT>::getBuffer()
{
    double* buffer = new double[getBufferSize()];
    if (buffer==0)
        throw(IOexception("SnapshotObserver::getBuffer() - error allocating memory."));
    else
    {
        double* bp = buffer;
        for (unsigned int ix = 0; ix < Nx; ix++)
        	for (unsigned int iy = 0; iy < Ny; iy++)
       	    {
           		ElMagField field = getField(ix,iy);
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

template <typename SourceT>
void SnapshotObserver<SourceT>::fromBuffer(double *buffer, std::size_t size)
{
    if (size==Nx*Ny*6)
    {
        ElMagField field;
        double *bp = buffer;
        for (unsigned int ix = 0; ix < Nx; ix++)
        	for (unsigned int iy = 0; iy < Ny; iy++)
       	    {
                Vector E, B;
                E.x = *bp++;
                E.y = *bp++;
                E.z = *bp++;
                B.x = *bp++;
                B.y = *bp++;
                B.z = *bp++;
           		setField(ix,iy,ElMagField(E,B));
       	    };
    } else {
	    throw(IOexception("SnapshotObserver::fromBuffer() - buffer size mismatch."));
    }
}

template <typename SourceT>
void SnapshotObserver<SourceT>::WriteFieldHDF5()
{
    herr_t status;
    cout << "writing HDF5 file " << FileName << endl;
    // Create a new file using the default properties.
    hid_t file = H5Fcreate (FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Fcreate()"));
    
    // Create dataspace for the observation positions.
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hsize_t pdims[3];
    pdims[0] = Nx;
    pdims[1] = Ny;
    pdims[2] = 3;
    hid_t pspace = H5Screate_simple (3, pdims, NULL);
    if (pspace<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Screate(pspace)"));
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
    if (pdcpl<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Pcreate(pdcpl)"));
    // Create the dataset.
    hid_t pdset = H5Dcreate(file,
        "ObservationPosition",		// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        pspace, H5P_DEFAULT,
        pdcpl, H5P_DEFAULT);
    if (pdset<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dcreate(pdset)"));
    // Write the data to the dataset
    status = H5Dwrite (pdset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			// mem space id
        pspace,
        H5P_DEFAULT,			// data transfer properties
        buffer);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dwrite(pdset)"));
    // attach scalar attributes
    hid_t atts1  = H5Screate(H5S_SCALAR);
    if (atts1<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Screate(atts1)"));
    hid_t attr1 = H5Acreate2(pdset, "Nx", H5T_NATIVE_INT, atts1, H5P_DEFAULT, H5P_DEFAULT);
    if (attr1<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Acreate2(attr1)"));
    status = H5Awrite(attr1, H5T_NATIVE_INT, &Nx);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Awrite(attr1)"));
    status = H5Sclose (atts1);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Sclose(atts1)"));
    hid_t atts2  = H5Screate(H5S_SCALAR);
    if (atts2<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Screate(atts2)"));
    hid_t attr2 = H5Acreate2(pdset, "Ny", H5T_NATIVE_INT, atts2, H5P_DEFAULT, H5P_DEFAULT);
    if (attr2<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Acreate2(attr2)"));
    status = H5Awrite(attr2, H5T_NATIVE_INT, &Ny);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Awrite(attr2)"));
    status = H5Sclose (atts2);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Sclose(atts2)"));
    // Close and release resources.
    status = H5Pclose (pdcpl);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Pclose(pdcpl)"));
    status = H5Dclose (pdset);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dclose(pdset)"));
    status = H5Sclose (pspace);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dclose(pspace)"));
    delete[] buffer;

    // Create dataspace for the field data.
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hsize_t dims[3];
    dims[0] = Nx;
    dims[1] = Ny;
    dims[2] = 6;
    hid_t space = H5Screate_simple (3, dims, NULL);
    if (space<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Screate(space)"));
    // buffer the data
    buffer = getBuffer();
    // Create the dataset creation property list
    hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Pcreate(dcpl)"));
    // Create the dataset.
    hid_t dset = H5Dcreate (file,
        "ElMagField", 			// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        space, H5P_DEFAULT,
        dcpl, H5P_DEFAULT);
    if (dset<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dcreate(dset)"));
    // Write the data to the dataset
    status = H5Dwrite (dset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			// mem space id
        space,
        H5P_DEFAULT,			// data transfer properties
        buffer);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dwrite(dset)"));
    // attach scalar attributes
    hid_t atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Screate(t0)"));
    hid_t attr = H5Acreate2(dset, "t0", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Acreate2(t0)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &t_obs);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Awrite(t0)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Sclose(t0)"));
    // Close and release resources.
    status = H5Pclose (dcpl);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Pclose(dcpl)"));
    status = H5Dclose (dset);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dclose(dset)"));
    status = H5Sclose (space);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Sclose(space)"));
    delete[] buffer;

    // If available create a dataset with the particle coordinates
    if (particleBuffer and (particlesStored > 0))
    {
        if (DEBUGLEVEL>=2) std::cout << "SnaphotObserver saving coordinates of " << particlesStored << " particles - ";
        // Create dataspace for the field data.
        // Setting maximum size to NULL sets the maximum size to be the current size.
        hsize_t cdims[2];
        cdims[0] = particlesStored;
        cdims[1] = 6;
        space = H5Screate_simple (2, cdims, NULL);
        if (space<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Screate(space)"));
        // Create the dataset creation property list
        dcpl = H5Pcreate (H5P_DATASET_CREATE);
        if (dcpl<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Pcreate(dcpl)"));
        // Create the dataset.
        dset = H5Dcreate (file,
            "ParticleCoordinates",  // dataset name
            H5T_NATIVE_DOUBLE,		// data type
            space, H5P_DEFAULT,
            dcpl, H5P_DEFAULT);
        if (dset<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dcreate(dset)"));
        // Write the data to the dataset
        status = H5Dwrite (dset,
            H5T_NATIVE_DOUBLE, 		// mem type id
            H5S_ALL, 			    // mem space id
            space,
            H5P_DEFAULT,			// data transfer properties
            particleBuffer);
        if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dwrite(dset)"));
        // attach scalar attributes
        atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Screate(N_particles)"));
        attr = H5Acreate2(dset, "N_particles", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Acreate2(N_particles)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &particlesStored);
        if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Awrite(N_particles)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Sclose(N_particles)"));
        // Close and release resources.
        status = H5Pclose (dcpl);
        if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Pclose(dcpl)"));
        status = H5Dclose (dset);
        if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dclose(dset)"));
        status = H5Sclose (space);
        if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Sclose(space)"));
        if (DEBUGLEVEL>=2) std::cout << "done." << std::endl << std::flush;
    }
    
    status = H5Fclose (file);
    if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Fclose()"));
    // no errors have occured if we made it 'til here
    cout << "writing HDF5 done." << endl;
    return;
}

template <typename SourceT>
void SnapshotObserver<SourceT>::generateOutput()
{
    WriteFieldHDF5();
}

template <typename SourceT>
Vector SnapshotObserver<SourceT>::CellPosition(int ix, int iy)
{
    return O + dX*((double)ix-0.5*((double)Nx-1.0)) + dY*((double)iy-0.5*((double)Ny-1.0));
}

template class SnapshotObserver<Beam>;
template class SnapshotObserver<Bunch>;
template class SnapshotObserver<Lattice>;

