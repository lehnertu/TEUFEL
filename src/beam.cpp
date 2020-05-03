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

#include "beam.h"
#include "bunch.h"
#include "particle.h"
#include "global.h"
#include "hdf5.h"
#include <iostream>

Beam::Beam()
{
    NOB = 0;
    tracker = TRACKING_NONE;
}

Beam::~Beam()
{
    for(int i=0; i<NOB; i++)
        delete B[i];
}

void Beam::Add(Bunch *bunch)
{
    NOB++;
    B.push_back(bunch);
}

void Beam::clear()
{
    // as B is only a list of pointers we have to delete all bunches individually
    for(int i=0; i<NOB; i++)
        delete B[i];
    NOB = 0;
    B.clear();
}

Bunch* Beam::getBunch(int i)
{
    Bunch *b = 0;
    if (i>=0 && i<NOB) b=B[i];
    return b;
}

int Beam::getNOB()
{
    return NOB;
}

int Beam::getNOP()
{
    int nop = 0;
    for(int i=0; i<NOB; i++)
    {
        nop += B[i]->getNOP();
    }
    return nop;
}

double Beam::getTotalCharge()
{
    double charge = 0.0;
    for(int i=0; i<NOB; i++)
    {
        charge += B[i]->getTotalCharge();
    }
    return charge;
}

void Beam::setupTracking(GeneralField *field)
{
    switch (tracker)
    {
        case TRACKING_NONE:
            throw(IOexception("Beam::setupTracking - no tracking method provided."));
            break;
        case TRACKING_EULER:
            throw(IOexception("Beam::setupTracking - EULER tracking method not yet implemented."));
            break;
        case TRACKING_VAY:
            InitVay(field);
            break;
        default:
            throw(IOexception("Beam::setupTracking - unknown tracking method."));
    }
}

void Beam::doStep(GeneralField *field)
{
    switch (tracker)
    {
        case TRACKING_NONE:
            throw(IOexception("Beam::doStep - no tracking method provided."));
            break;
        case TRACKING_EULER:
            throw(IOexception("Beam::doStep - EULER tracking method not yet implemented."));
            break;
        case TRACKING_VAY:
            StepVay(field);
            break;
        default:
            throw(IOexception("Beam::doStep - unknown tracking method."));
    }
}

void Beam::InitVay(GeneralField *field)
{
    for(int i=0; i<NOB; i++)
    {
        B[i]->InitVay(dt, field);
    }
}

void Beam::StepVay(GeneralField *field)
{
    for(int i=0; i<NOB; i++)
    {
        B[i]->StepVay(field);
    }
}

//! @todo how to store several bunches hierarchically ?
int Beam::WriteWatchPointHDF5(std::string filename)
{
    int nop = getNOP();
    // buffer coordinates of all particles
    double *buffer = new double[nop*6];
    double *bp = buffer;
    // remaining buffer capacity
    int remaining = nop;
    for(int i=0; i<NOB; i++)
    {
        bp = B[i]->bufferCoordinates(bp, remaining);
        remaining = nop - (bp-buffer)/6;
    }
    if (remaining != 0)
        std::cout << "warning in Beam::WriteWatchPointHDF5() - buffer size mismatch" << std::endl;
    herr_t status;
    cout << "writing HDF5 file " << filename << endl;
    // Create a new file using the default properties.
    hid_t file = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0 )
    {
        std::cout << "Beam::WriteWatchPointHDF5 - error creating file" << std::endl;
        return -1;
    };
    // Create dataspace. Setting maximum size to NULL sets the maximum
    // size to be the current size.
    hsize_t dims[2];
    dims[0] = nop;
    dims[1] = 6;
    hid_t space = H5Screate_simple (2, dims, NULL);
    if (space < 0 )
    {
        std::cout << "Beam::WriteWatchPointHDF5 - error creating dataspace" << std::endl;
        return -1;
    };
    // Create the dataset creation property list
    hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl < 0 )
    {
        std::cout << "Beam::WriteWatchPointHDF5 - error creating property list" << std::endl;
        return -1;
    };
    // Create the dataset.
    hid_t dset = H5Dcreate (file,
        "electrons", 			// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        space, H5P_DEFAULT,
        dcpl, H5P_DEFAULT);
    if (dset < 0 )
    {
        std::cout << "Beam::WriteWatchPointHDF5 - error creating dataset" << std::endl;
        return -1;
    };
    // Write the data to the dataset
    status = H5Dwrite (dset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			// mem space id
        space,
        H5P_DEFAULT,			// data transfer properties
        buffer);
    if (status < 0 )
    {
        std::cout << "Beam::WriteWatchPointHDF5 - error writing dataset" << std::endl;
        return -1;
    }
    // Close and release resources.
    status = H5Pclose (dcpl);
    if (status < 0 )
    {
        std::cout << "Beam::WriteWatchPointHDF5 - error releasing property list" << std::endl;
        return -1;
    }
    status = H5Dclose (dset);
    if (status < 0 )
    {
        std::cout << "Beam::WriteWatchPointHDF5 - error releasing dataset" << std::endl;
        return -1;
    }
    status = H5Sclose (space);
    if (status < 0 )
    {
        std::cout << "Beam::WriteWatchPointHDF5 - error releasing dataspace" << std::endl;
        return -1;
    }
    status = H5Fclose (file);
    if (status < 0 )
    {
        std::cout << "Beam::WriteWatchPointHDF5 - error closing file" << std::endl;
        return -1;
    }
    // no errors have occured if we made it 'til here
    std::cout << "writing HDF5 done." << std::endl;
    return nop;
}

ElMagField Beam::RetardedField(double obs_time, Vector obs_pos)
{
    // initialize to zero
    ElMagField field;
    // integrate for all bunches
    for (int i=0; i<NOB; i++)
        field += B[i]->RetardedField(obs_time, obs_pos);
    return field;
}

void Beam::integrateFieldTrace(
    Vector ObservationPoint,
    double t0,
    double dt,
    int nots,
    std::vector<ElMagField> *ObservationField)
{
    // just do the summation over all te bunches
    for(int i=0; i<NOB; i++)
    {
        B[i]->integrateFieldTrace(ObservationPoint,t0,dt,nots,ObservationField);
    }
}

void Beam::integrateFieldTrace(
        Vector ObservationPoint,
        FieldTrace *trace)
{
    // just do the summation over all the bunches
    for(int i=0; i<NOB; i++)
        B[i]->integrateFieldTrace(ObservationPoint, trace);
}

