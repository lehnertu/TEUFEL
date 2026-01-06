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

#include <cstddef>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <math.h>
#include "hdf5.h"
#include "particle.h"
#include "vector.h"

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

void CSR::init(Beam *beam)
{
    if (teufel::rank==0) std::cout << "CSR::init()" << std::endl;
    if (is_initialized)
        throw(IOexception("error - CSR::init() called twice."));
    source_beam = beam;
    // iterate over the beam to create a list of all particles
    // making up the field source
    NoP = source_beam->getNOP();
    size_t NoB = source_beam->getNOB();
    for (size_t ib=0; ib<NoB; ib++)
    {
        Bunch *B = source_beam->getBunch(ib);
        size_t B_NoP = B->getNOP();
        std::cout << "reading bunch No. " << ib << " with " << B_NoP << " particles." << std::endl;
        for (size_t ip=0; ip<B_NoP; ip++)
        {
            ChargedParticle *p = B->getParticle(ip);
            particles.push_back(p);
        };
    }
    if (teufel::rank==0) std::cout << "CSR::init() stored references to " << particles.size() << " particles." << std::endl;
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

void CSR::update(double tracking_time)
{
    if (teufel::rank==0)
        std::cout << "CSR::update() at tracking time " << tracking_time << " s" << std::endl;
    if ((int)NoP != source_beam->getNOP())
        throw(IOexception("error - CSR::update() mismatch of particle numbers."));
    
    // get the particle coordinates
    // the data obtained here correspond to the half-step positions
    //! @todo all particles should have the same time stamp anyway - maybe better check
    double avg_time = 0;
    Vector avg_pos = VectorZero;
    Vector avg_mom = VectorZero;
    Vector avg_acc = VectorZero;
    for(size_t i=0; i<NoP; i++)
    {
        avg_time += particles[i]->getTime();
        avg_pos += particles[i]->getPosition();
        avg_mom += particles[i]->getMomentum();
    }
    avg_time /= NoP;
    avg_pos /= NoP;
    avg_mom /= NoP;
        
    // create a snapshot of the beam
    Snapshot* snap = new Snapshot{
        avg_time,
        avg_pos,
        avg_mom,
        {}      // slices (start empty)
    };

    // compute particle distance from center reference plane
    Vector forward = avg_mom;
    forward.normalize();
    double *s = new double[NoP];
    for(size_t  i=0; i<NoP; i++)
        s[i] = dot(particles[i]->getPosition()-avg_pos, forward);
    // sort by longitudinal position
    size_t *sorting = new size_t[NoP];
    // fill sorting with indices 0..NOP-1
    std::iota(sorting, sorting + NoP, 0);
    // sort indices by comparing s[]
    std::sort(sorting, sorting + NoP,
          [&](size_t i, size_t j) { return s[i] < s[j]; });
    if (DEBUGLEVEL>=2)
    {
        std::cout << "   s[0]=" << s[sorting[0]] << " s[N]=" << s[sorting[NoP-1]] << std::endl;
    };
    
    // the first N_rem slices contain N_mod+1 particles, the rest N_mod
    size_t N_mod = NoP / numSlices;
    size_t N_rem = NoP % numSlices;
    size_t p_index = 0;

    // distribute the particles over the slices
    for (size_t sl=0; sl<numSlices; sl++)
    {
        // number of particles for this slice
        size_t n_sl = N_mod;
        if (sl<N_rem) n_sl++;
        // compute slice properties - weighted average by charge
        double total_charge = 0.0;
        double min_s = s[sorting[p_index]];
        double max_s = s[sorting[p_index]];
        avg_pos = VectorZero;
        avg_mom = VectorZero;
        avg_acc = VectorZero;
        size_t p_index_before_slice = p_index;
        for (size_t i_sl=0; i_sl<n_sl; i_sl++)
        {
            ChargedParticle *p = particles[sorting[p_index]];
            double p_charge = p->getCharge();
            total_charge += p_charge;
            double p_s = s[sorting[p_index]];
            if (min_s > p_s) min_s = p_s;
            if (max_s < p_s) max_s = p_s;
            avg_pos += p->getPosition() * p_charge;
            avg_mom += p->getMomentum() * p_charge;
            avg_acc += p->getAccel() * p_charge;
            p_index++;
        }
        avg_pos /= total_charge;
        avg_mom /= total_charge;
        avg_acc /= total_charge;
        // scan the same slice again to compute the rms radius
        p_index = p_index_before_slice;
        double rms_radius = 0.0;
        for (size_t i_sl=0; i_sl<n_sl; i_sl++)
        {
            ChargedParticle *p = particles[sorting[p_index]];
            double p_charge = p->getCharge();
            Vector rad = cross(p->getPosition()-avg_pos, forward);
            double r_sq = rad.abs2nd();
            rms_radius += p_charge * r_sq;
            p_index++;
        }
        rms_radius = sqrt(rms_radius/total_charge);
        // append the slice to the snapshot
        snap->slices.push_back(
            Slice{
                .s_min = min_s,
                .s_max = max_s,
                .charge = total_charge,
                .position = avg_pos,
                .momentum = avg_mom,
                .accel = avg_acc,
                .radius = rms_radius}
        );
    }
    if (DEBUGLEVEL>=2)
    {
        std::cout << "   particle index after slicing = " << p_index << std::endl;
    };

    // Append the snapshot to history
    history.push_back(snap);
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
        double *bp = snap_buffer;
        for(int i=0; i<NOS; i++)
        {
            Snapshot* snap = history[i];
            *bp++ = snap->tracking_time;
            *bp++ = snap->central_position.x;
            *bp++ = snap->central_position.y;
            *bp++ = snap->central_position.z;
            *bp++ = snap->central_momentum.x;
            *bp++ = snap->central_momentum.y;
            *bp++ = snap->central_momentum.z;
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

        // --------- write slice data ----------

        // Create dataspace for the slice properties.
        // Setting maximum size to NULL sets the maximum size to be the current size.
        size_t NoSteps = history.size();
        hsize_t slc_dims[3];
        slc_dims[0] = NoSteps;
        slc_dims[1] = numSlices;
        slc_dims[2] = 13; // the size of the Slice struct in doubles
        hid_t slc_space = H5Screate_simple (3, slc_dims, NULL);
        if (slc_space<0) throw(IOexception("CSR::write_output() - error in H5Screate(slc_space)"));

        // buffer the data
        double *slc_buffer = new double[NoSteps*numSlices*13];
        bp = slc_buffer;
        for (size_t i_st=0; i_st<NoSteps; i_st++)
        {
            Snapshot* snap = history[i_st];
            for (size_t i_slc=0; i_slc<numSlices; i_slc++)
            {
                Slice slc = snap->slices[i_slc];
                *bp++ = slc.s_min;
                *bp++ = slc.s_max;
                *bp++ = slc.charge;
                *bp++ = slc.position.x;
                *bp++ = slc.position.y;
                *bp++ = slc.position.z;
                *bp++ = slc.momentum.x;
                *bp++ = slc.momentum.y;
                *bp++ = slc.momentum.z;
                *bp++ = slc.accel.x;
                *bp++ = slc.accel.x;
                *bp++ = slc.accel.x;
                *bp++ = slc.radius;
            }
        }

        // Create the dataset creation property list
        hid_t slc_dcpl = H5Pcreate (H5P_DATASET_CREATE);
        if (slc_dcpl<0) throw(IOexception("CSR::write_output() - error in H5Pcreate(slc_dcpl)"));
        // Create the dataset
        hid_t slc_dset = H5Dcreate(file,
            "Slices",  	     	    // dataset name
            H5T_NATIVE_DOUBLE,		// data type
            slc_space, H5P_DEFAULT,
            slc_dcpl, H5P_DEFAULT);
        if (slc_dset<0) throw(IOexception("CSR::write_output() - error in H5Dcreate(slc_dset)"));

        // Write the data to the dataset
        status = H5Dwrite (slc_dset,
            H5T_NATIVE_DOUBLE, 		// mem type id
            H5S_ALL, 			    // mem space id
            slc_space,
            H5P_DEFAULT,			// data transfer properties
            slc_buffer);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Dwrite(slc_dset)"));

        // attach scalar attributes
        atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("CSR::write_output() - error in H5Screate(N_steps)"));
        attr = H5Acreate2(slc_dset, "N_steps", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("CSR::write_output() - error in H5Acreate2(N_steps)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &NoSteps);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Awrite(N_steps)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Sclose(N_steps)"));

        atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("CSR::write_output() - error in H5Screate(N_slices)"));
        attr = H5Acreate2(slc_dset, "N_slices", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("CSR::write_output() - error in H5Acreate2(N_slices)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &numSlices);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Awrite(N_slices)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Sclose(N_slices)"));

        // Close and release resources.
        status = H5Pclose (slc_dcpl);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Pclose(slc_dcpl)"));
        status = H5Dclose (slc_dset);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Dclose(slc_dset)"));
        status = H5Sclose (slc_space);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Dclose(slc_space)"));

        status = H5Fclose(file);
        if (status<0) throw(IOexception("CSR::write_output() - error in H5Fclose()"));
        
        delete[] snap_buffer;
        delete[] slc_buffer;

        // no errors have occured if we made it 'til here
        cout << "writing HDF5 done." << endl;
    }
}

