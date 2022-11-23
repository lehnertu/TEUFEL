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
#include "SDDS.h"
#include <math.h>
#include <iostream>

Beam::Beam()
{
    NOB = 0;
    NOTS = 0;
    dt = 0.0;
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
    NOTS = 0;
    B.clear();
}

void Beam::clearTrajectories()
{
    NOTS = 0;
    for(int i=0; i<NOB; i++)
        B[i]->clearTrajectories();
}

void Beam::preAllocate(int nTraj)
{
    for(int i=0; i<NOB; i++)
        B[i]->preAllocate(nTraj);
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
    NOTS++;
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
    NOTS = 0;
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

int Beam::getStepBufferSize()
{
    return getNOP()*10;
}

void Beam::bufferStep(double *buffer)
{
    double *buf = buffer;
    for(int i=0; i<NOB; i++)
    {
        // the pointer is advanced as the bunch is filled into the buffer
        buf = B[i]->bufferStep(buf);
    }
}

void Beam::setStepFromBuffer(double *buffer)
{
    double *buf = buffer;
    for(int i=0; i<NOB; i++)
    {
        // the pointer is advanced as the buffer is consumed
        buf = B[i]->setStepFromBuffer(buf);
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
    delete[] buffer;
    return nop;
}

int Beam::WriteWatchPointSDDS(std::string filename)
{
    cout << "writing SDDS file " << filename << endl;
    SDDS_DATASET data;
    if (1 != SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,filename.c_str()))
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error initializing output" << std::endl;
        return 1;
    }
    if  (SDDS_DefineSimpleParameter(&data,"NumberOfParticles","", SDDS_LONG)!=1)
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error defining parameters" << std::endl;
        return 2;
    }
    if  (
        SDDS_DefineColumn(&data,"t\0","t\0","s\0","Time\0",NULL, SDDS_DOUBLE,0)   ==-1 || 
        SDDS_DefineColumn(&data,"x\0","x\0","m\0","PositionX\0",NULL, SDDS_DOUBLE,0) == -1 ||
        SDDS_DefineColumn(&data,"y\0","y\0","m\0","PositionY\0",NULL, SDDS_DOUBLE,0) == -1 ||
        SDDS_DefineColumn(&data,"z\0","z\0","m\0","PositionZ\0",NULL, SDDS_DOUBLE,0) == -1 || 
        SDDS_DefineColumn(&data,"px\0","px\0",NULL,"BetaGammaX\0",NULL, SDDS_DOUBLE,0)== -1 || 
        SDDS_DefineColumn(&data,"py\0","py\0",NULL,"BetaGammaY\0",NULL,SDDS_DOUBLE,0) == -1 ||
        SDDS_DefineColumn(&data,"pz\0","pz\0",NULL,"BetaGammaZ\0",NULL,SDDS_DOUBLE,0) == -1 ||
        SDDS_DefineColumn(&data,"p\0","p\0",NULL,"BetaGamma\0",NULL,SDDS_DOUBLE,0) == -1 ||
        SDDS_DefineColumn(&data,"xp\0","xp\0",NULL,"Xprime\0",NULL, SDDS_DOUBLE,0)== -1 || 
        SDDS_DefineColumn(&data,"yp\0","yp\0",NULL,"Bprime\0",NULL,SDDS_DOUBLE,0) == -1 ||
        SDDS_DefineColumn(&data,"gamma\0","gamma\0",NULL,"RelativisticFactor\0",NULL,SDDS_DOUBLE,0)==-1
        )
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error defining data columns" << std::endl;
        return 3;
    }
    if (SDDS_WriteLayout(&data) != 1)
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error writing layout" << std::endl;
        return 4;
    }
    // start a page with number of lines equal to the number of particles
    int NOP = getNOP();
    // cout << "SDDS start page" << endl;
    if (SDDS_StartPage(&data,(int32_t)NOP) !=1 )
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error starting page" << std::endl;
        return 5;
    }
    // write the single valued variables
    // cout << "SDDS write parameters" << endl;
    if( SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
        "NumberOfParticles",NOP,
        NULL ) !=1
        )
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error setting parameters" << std::endl;
        return 6;
    }
    // write the table of particle data
    // cout << "SDDS writing " << NOP << " particles" << endl;
    int index=0;
    for(int b=0; b<NOB; b++)
    {
        Bunch *bunch = B[b];
        int nop = bunch->getNOP();
        for( int i=0; i<nop; i++)
        {
            ChargedParticle *P = bunch->getParticle(i);
            double t = P->getTime();
            Vector X = P->getPosition();
            Vector BG = P->getMomentum();
            if( SDDS_SetRowValues(&data,
                SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, index,
                "t",t,
                "x",X.x,
                "y",X.y,
                "z",X.z,
                "px",BG.x,
                "py",BG.y,
                "pz",BG.z,
                "p",BG.norm(),
                "xp",BG.x/BG.z,
                "yp",BG.y/BG.z,
                "gamma",sqrt(1.0+BG.abs2nd()),
                NULL) != 1
                )
            {
                std::cout << "Bunch::WriteWatchPointSDDS - error writing data columns" << std::endl;
                return 7;
            }
            index++;
        }
    }

    if( SDDS_WritePage(&data) != 1)
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error writing page" << std::endl;
        return 8;
    }
    // finalize the file
    if (SDDS_Terminate(&data) !=1 )
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error terminating data file" << std::endl;
        return 9;
    }	
    // no errors have occured if we made it 'til here
    // cout << "writing SDDS done." << endl;
    return 0;
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
        FieldTrace *trace)
{
    // just do the summation over all the bunches
    for(int i=0; i<NOB; i++)
        B[i]->integrateFieldTrace(ObservationPoint, trace);
}

