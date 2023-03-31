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

#include "fel1d.h"
#include "global.h"

#include <iostream>
#include <math.h>
#include "hdf5.h"

FEL_1D::FEL_1D(
    double time_step,
    int field_size )
{
    N_field = field_size;
    dt = time_step;
    N_steps = 0;
    prop = Vector(0.0, 0.0, 1.0);
    e_x = Vector(1.0, 0.0, 0.0);
    e_x = Vector(0.0, 1.0, 0.0);
    createOutput = false;
    setup();
}

FEL_1D::FEL_1D(
    double time_step,
    const pugi::xml_node node,
    InputParser *parser )
{
    dt = time_step;
    N_steps = 0;
    pugi::xml_attribute att = node.attribute("N");
    if (!att)
        throw(IOexception("InputParser::FEL1D - attribute N not found."));
    N_field = parser->parseInt(att);
    pugi::xml_node vec = node.child("position");
    if (!vec)
        throw(IOexception("InputParser::FEL1D - <position> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(vec.attribute("x"));
        y = parser->parseDouble(vec.attribute("y"));
        z = parser->parseDouble(vec.attribute("z"));
        head = Vector(x,y,z);
    }
    vec = node.child("prop");
    if (!vec)
        throw(IOexception("InputParser::FEL1D - <prop> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(vec.attribute("x"));
        y = parser->parseDouble(vec.attribute("y"));
        z = parser->parseDouble(vec.attribute("z"));
        prop = Vector(x,y,z);
    }
    prop.normalize();
    vec = node.child("pol");
    if (!vec)
        throw(IOexception("InputParser::FEL1D - <pol> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(vec.attribute("x"));
        y = parser->parseDouble(vec.attribute("y"));
        z = parser->parseDouble(vec.attribute("z"));
        e_x = Vector(x,y,z);
    }
    e_x -= prop*dot(e_x,prop);
    e_x.normalize();
    e_y = cross(prop, e_x);
    // define file output if requested
    pugi::xml_node lognode = node.child("log");
    if (lognode)
    {
        pugi::xml_attribute fn = lognode.attribute("file");
        if (!fn) throw(IOexception("InputParser::FEL1D - <log> filename for log not found."));
        FileName = fn.as_string();
        step_Output = 1;
        pugi::xml_attribute st = lognode.attribute("step");
        if (st) step_Output = st.as_int();
        createOutput = true;
    } else {
        createOutput = false;
    }
    setup();
    // seed the field if requested
    // this must be done after setup() so the fields have been created
    pugi::xml_node seednode = node.child("seed");
    if (seednode)
    {
        double E0 = parser->parseDouble(seednode.attribute("E0"));
        double lambda = parser->parseDouble(seednode.attribute("lambda"));
        double tau = parser->parseDouble(seednode.attribute("tau"));
        double tstart = parser->parseDouble(seednode.attribute("tstart"));
        seed(E0, lambda, tau, tstart);
    }
}

void FEL_1D::setup()
{
    // proper orthonormalization of (e_x, e_y, prop) is assumed
    dz = prop * (-dt*SpeedOfLight);
    field_E = std::vector<double>(N_field, 0.0);
    // the initialized field is the first stored array
    field_storage.clear();
    field_storage.push_back(field_E);
    N_steps = 0;
    // print properties
    if (teufel::rank==0)
    {
        std::cout << "FEL-1D interaction  N = " << N_field << ",   dz = " << dz.norm() << " m" << std::endl;
        std::cout << "  propagation = (" << prop.x << ", " << prop.y << ", " << prop.z << ")" << std::endl;
        std::cout << "  dz = (" << dz.x << ", " << dz.y << ", " << dz.z << ")" << std::endl;
        std::cout << "  e_E = (" << e_x.x << ", " << e_x.y << ", " << e_x.z << ")" << std::endl;
        std::cout << "  e_b = (" << e_y.x << ", " << e_y.y << ", " << e_y.z << ")" << std::endl;
        std::cout << "  logging to " << FileName << " every " << step_Output << " steps." << std::endl;
    }
}

void FEL_1D::seed(double E0, double lambda, double tau, double t_start)
{
    double omega = 2.0*Pi*SpeedOfLight/lambda;
    for (int i=0; i<N_field; i++)
    {
        double t = i*dt-t_start;
        field_E[i] = E0*cos(omega*t)*exp(-0.5*pow(t/tau,2));
    }
    // if we are seeding this one replaces the first storage array created by setup()
    field_storage.clear();
    field_storage.push_back(field_E);
    N_steps = 0;
    if (teufel::rank==0)
    {
        std::cout << "  seeding with E0 = " << E0 << "  lambda = " << lambda << "  tau = " << tau << "  t_start = " << t_start << std::endl;
    }
}

void FEL_1D::step(Beam *beam)
{
    //! @todo needs to be implemented

    N_steps++;
    
    // every number of steps we record the fields
    if (0 == N_steps % step_Output)
    {
        field_storage.push_back(field_E);
    }
}

void FEL_1D::write_output()
{
    //! @todo needs to be implemented
    if (createOutput)
    {
    
        cout << "FEL1D : writing HDF5 file " << FileName << endl;

        herr_t status;
        // Create a new file using the default properties.
        hid_t file = H5Fcreate (FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file<0) throw(IOexception("FEL1D::write_output - error in H5Fcreate()"));
        
        // Create dataspace for the observation positions.
        // Setting maximum size to NULL sets the maximum size to be the current size.
        int N_steps = field_storage.size();
        hsize_t pdims[2];
        pdims[0] = N_steps;
        pdims[1] = N_field;
        hid_t pspace = H5Screate_simple (2, pdims, NULL);
        if (pspace<0) throw(IOexception("FEL1D::write_output - error in H5Screate(pspace)"));
        // buffer the data
        double *buffer = new double[N_field*N_steps];
        double *bp = buffer;
        for (int ix=0; ix<N_steps; ix++)
        {
            std::vector<double> E = field_storage[ix];
	        for (int iy = 0; iy < N_field; iy++)
	            *bp++ = E[iy];
	    }
        // Create the dataset creation property list
        hid_t pdcpl = H5Pcreate (H5P_DATASET_CREATE);
        if (pdcpl<0) throw(IOexception("FEL1D::write_output - error in H5Pcreate(pdcpl)"));
        // Create the dataset.
        hid_t pdset = H5Dcreate(file,
            "InteractionField",		// dataset name
            H5T_NATIVE_DOUBLE,		// data type
            pspace, H5P_DEFAULT,
            pdcpl, H5P_DEFAULT);
        if (pdset<0) throw(IOexception("FEL1D::write_output - error in H5Dcreate(pdset)"));
        // Write the data to the dataset
        status = H5Dwrite (pdset,
            H5T_NATIVE_DOUBLE, 		// mem type id
            H5S_ALL, 			    // mem space id
            pspace,
            H5P_DEFAULT,			// data transfer properties
            buffer);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Dwrite(pdset)"));
        
        // attach scalar attributes
        hid_t atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("FEL1D::write_output - error in H5Screate(N_field)"));
        hid_t attr = H5Acreate2(pdset, "N_field", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("FEL1D::write_output - error in H5Acreate2(N_field)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &N_field);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Awrite(N_field)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Sclose(N_field)"));
        
        atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("FEL1D::write_output - error in H5Screate(N_steps)"));
        attr = H5Acreate2(pdset, "N_steps", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("FEL1D::write_output - error in H5Acreate2(N_steps)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &N_steps);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Awrite(N_steps)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Sclose(N_steps)"));

        atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("FEL1D::write_output - error in H5Screate(dt)"));
        attr = H5Acreate2(pdset, "dt", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("FEL1D::write_output - error in H5Acreate2(dt)"));
        status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dt);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Awrite(dt)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Sclose(dt)"));

        // Close and release resources.
        status = H5Pclose (pdcpl);
        if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Pclose(pdcpl)"));
        status = H5Dclose (pdset);
        if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dclose(pdset)"));
        status = H5Sclose (pspace);
        if (status<0) throw(IOexception("SnapshotObserver::WriteFieldHDF5 - error in H5Dclose(pspace)"));
        delete buffer;

        status = H5Fclose(file);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Fclose()"));
        // no errors have occured if we made it 'til here
        cout << "writing HDF5 done." << endl;
        return;
    }
}
    
FEL_1D::~FEL_1D()
{
}

ElMagField FEL_1D::Field(double t, Vector X)
//! @todo implementation missing
{
    return ElMagFieldZero;
}
