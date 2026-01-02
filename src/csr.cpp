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

#include "csr.h"
#include "global.h"

#include <iostream>
#include <math.h>
#include "hdf5.h"

CSR::CSR()
{
    is_initialized = false;
    numSlices = 0;
    createOutput = false;
}

CSR::CSR(
    const pugi::xml_node node,
    InputParser *parser )
{
    // no source definet yet
    is_initialized = false;
    pugi::xml_attribute att = node.attribute("N_slices");
    if (!att)
        throw(IOexception("InputParser::CSR - attribute N_slices not found."));
    numSlices = parser->parseInt(att);
    // define file output if requested
    pugi::xml_node lognode = node.child("log");
    if (lognode)
    {
        pugi::xml_attribute fn = lognode.attribute("file");
        if (!fn) throw(IOexception("InputParser::CSR - <log> filename for log not found."));
        FileName = fn.as_string();
        createOutput = true;
    } else {
        createOutput = false;
    }
    if (teufel::rank==0)
    {
        std::cout << "CSR interaction: " << numSlices << " slices" << std::endl;
    }
}

void CSR::init()
{
    if (teufel::rank==0)
    {
        std::cout << "CSR::init()" << std::endl;
    }
    if (is_initialized)
        throw(IOexception("error - CSR::init() called twice."));
    is_initialized = true;
}

CSR::~CSR()
{
// ChatGPT created
    for (Snapshot* s : history) {
        delete s;
    }
    history.clear();
// ChatGPT created
}

void CSR::update(Beam *beam, double tracking_time)
{
    int NOP = beam->getNOP();
    size_t bufsize = beam->getStepBufferSize();
    if (teufel::rank==0)
    {
        std::cout << "CSR::update() at tracking time " << tracking_time << " s";
        std::cout << "   NOP=" << NOP << " BUF=" << bufsize << std::endl;
    }
    
    // get the particle coordinates
    // the data obtained here correspond to the half-step positions
    // which are stored 
    double *buffer = new double[bufsize];
    beam->bufferStep(buffer);
    double *ptime = new double[NOP];
    Vector *position = new Vector[NOP];
    Vector *momentum = new Vector[NOP];
    Vector *accel = new Vector[NOP];
    double *bp = buffer;
    for(int i=0; i<NOP; i++)
    {
        ptime[i] = *bp++;
        position[i].x = *bp++;
        position[i].y = *bp++;
        position[i].z = *bp++;
        momentum[i].x = *bp++;
        momentum[i].y = *bp++;
        momentum[i].z = *bp++;
        accel[i].x = *bp++;
        accel[i].y = *bp++;
        accel[i].z = *bp++;
    };
    //! @todo all particles should have the same time stamp anyway - maybe better check
    double avg_time = 0;
    Vector avg_pos = VectorZero;
    Vector avg_mom = VectorZero;
    for(int i=0; i<NOP; i++)
    {
        avg_time += ptime[i];
        avg_pos += position[i];
        avg_mom += momentum[i];
    }
    avg_time /= NOP;
    avg_pos /= NOP;
    avg_mom /= NOP;
        
    // create a snapshot of the beam
    Snapshot* snap = new Snapshot{
        avg_time,
        avg_pos,
        avg_mom,
        {}      // slices (start empty)
    };

    /* Fill slices  - ChatGPT created
    snap->slices.reserve(numSlices);
    for (unsigned int i = 0; i < numSlices; ++i) {
        snap->slices.push_back(Slice{
            1.0,                    // charge
            Vector{0.0, 0.0, 0.0},  // position
            Vector{0.0, 0.0, 0.0},  // momentum
            Vector{0.0, 0.0, 0.0},  // accel
            0.1,                    // length
            0.01                    // radius
        });
    }
    */
    
    // Append to history
    history.push_back(snap);
    
    delete[] buffer;
    delete[] ptime;
    delete[] position;
    delete[] momentum;
    delete[] accel;
}

ElMagField CSR::Field(double t, Vector X)
{
    //! @todo returning zero fields for now
    return ElMagField();
}

void CSR::write_output()
{
    if (createOutput)
    {
        cout << "CSR : writing slice data to " << FileName << endl;

        herr_t status;
        // Create a new file using the default properties.
        hid_t file = H5Fcreate (FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file<0) throw(IOexception("CSR::write_output() - error in H5Fcreate()"));

        // --------- write snapshot data ----------

        // Create dataspace for the center coordinates (time, position, momentum).
        // Setting maximum size to NULL sets the maximum size to be the current size.
        int NOS = history.size();
        hsize_t snap_dims[2];
        snap_dims[0] = NOS;
        snap_dims[1] = 7;
        hid_t snap_space = H5Screate_simple (2, snap_dims, NULL);
        if (snap_space<0) throw(IOexception("CSR::write_output() - error in H5Screate(snap_space)"));

        // buffer the data
        double *snap_buffer = new double[NOS*7];
        double *bp = (double *)snap_buffer;
        for(int i=0; i<NOS; i++)
        {
            *bp++ = history[i]->tracking_time;
            *bp++ = history[i]->central_position.x;
            *bp++ = history[i]->central_position.y;
            *bp++ = history[i]->central_position.z;
            *bp++ = history[i]->central_momentum.x;
            *bp++ = history[i]->central_momentum.y;
            *bp++ = history[i]->central_momentum.z;
        }

        // Create the dataset creation property list
        hid_t snap_dcpl = H5Pcreate (H5P_DATASET_CREATE);
        if (snap_dcpl<0) throw(IOexception("CSR::write_output() - error in H5Pcreate(snap_dcpl)"));
        // Create the dataset
        hid_t snap_dset = H5Dcreate(file,
            "Snapshots",	     	    // dataset name
            H5T_NATIVE_DOUBLE,		// data type
            snap_space, H5P_DEFAULT,
            snap_dcpl, H5P_DEFAULT);
        if (snap_dset<0) throw(IOexception("CSR::write_output() - error in H5Dcreate(snap_dset)"));
        // Write the data to the dataset
        status = H5Dwrite (snap_dset,
            H5T_NATIVE_DOUBLE, 		// mem type id
            H5S_ALL, 			    // mem space id
            snap_space,
            H5P_DEFAULT,			// data transfer properties
            snap_buffer);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Dwrite(snap_dset)"));

        // attach scalar attributes
        hid_t atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("CSR::write_output() - error in H5Screate(N_steps)"));
        hid_t attr = H5Acreate2(snap_dset, "N_steps", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("CSR::write_output() - error in H5Acreate2(N_steps)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &NOS);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Awrite(N_steps)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Sclose(N_steps)"));

        // Close and release resources.
        status = H5Pclose (snap_dcpl);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Pclose(snap_dcpl)"));
        status = H5Dclose (snap_dset);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Dclose(snap_dset)"));
        status = H5Sclose (snap_space);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Dclose(snap_space)"));

        //! @todo write slice data

        status = H5Fclose(file);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Fclose()"));

        // no errors have occured if we made it 'til here
        cout << "writing HDF5 done." << endl;
    }
}

