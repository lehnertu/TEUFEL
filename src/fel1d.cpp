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
    ds = time_step * SpeedOfLight;
    time = 0.0;
    N_steps = 0;
    step_Output = 1;
    prop = Vector(0.0, 0.0, 1.0);
    e_x = Vector(1.0, 0.0, 0.0);
    e_x = Vector(0.0, 1.0, 0.0);
    head = Vector(0.0, 0.0, 0.0);
    origin = head;
    createOutput = false;
    // assuming a 1mm waist with ZR=1m at 1m from the start
    setup(1.0, 1.0e-3, 1.0);
}

FEL_1D::FEL_1D(
    double time_step,
    const pugi::xml_node node,
    InputParser *parser )
{
    dt = time_step;
    ds = time_step * SpeedOfLight;
    time = 0.0;
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
        origin = head;
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
    // read the optical mode definition
    pugi::xml_node mode = node.child("mode");
    if (mode)
    {
        double zR = parser->parseDouble(mode.attribute("zR"));
        double w = parser->parseDouble(mode.attribute("w0"));
        double zw = parser->parseDouble(mode.attribute("zwaist"));
        setup(zR, w, zw);
    }
    else
        throw(IOexception("InputParser::FEL1D - <mode> not found."));
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

void FEL_1D::setup(double ZR, double w0, double z_w)
{
    // proper orthonormalization of (e_x, e_y, prop) is assumed
    dz = prop * (-dt*SpeedOfLight);
    field_E = std::vector<double>(N_field, 0.0);
    previous_E = field_E;
    // the initialized field is the first stored array (all zeros)
    field_storage.clear();
    field_storage.push_back(field_E);
    J_storage.clear();
    J_storage.push_back(field_E);
    N_steps = 0;
    time = 0.0;
    
    // the initialization of the envelope
    z_Rayleigh = ZR;
    w_0 = w0;
    z_waist = z_w;
    
    // print properties
    if (teufel::rank==0)
    {
        std::cout << "FEL-1D interaction  N=" << N_field << ",   dz=" << dz.norm() << " m" << std::endl;
        std::cout << "  propagation = (" << prop.x << ", " << prop.y << ", " << prop.z << ")" << std::endl;
        std::cout << "  dz = (" << dz.x << ", " << dz.y << ", " << dz.z << ")" << std::endl;
        std::cout << "  e_E = (" << e_x.x << ", " << e_x.y << ", " << e_x.z << ")" << std::endl;
        std::cout << "  e_B = (" << e_y.x << ", " << e_y.y << ", " << e_y.z << ")" << std::endl;
        std::cout << "  opt. mode zR=" << z_Rayleigh << " m,  w0=" << w_0 << " m  at s=" << z_waist << " m" << std::endl;
        if (createOutput)
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
    // the seed is a stable wave packet
    previous_E = field_E;
    // @todo make intentionally wrong, I want to see the field split into 2 counter-propagating fractions
    // previous_E.erase(previous_E.begin()); // pop the first
    // previous_E.push_back(0.0); // extend the end
    // if we are seeding this one replaces the first storage array created by setup()
    field_storage.clear();
    field_storage.push_back(field_E);
    J_storage.clear();
    J_storage.push_back(std::vector<double>(N_field, 0.0));
    N_steps = 0;
    if (teufel::rank==0)
    {
        std::cout << "  seeding with E0 = " << E0 << "  lambda = " << lambda << "  tau = " << tau << "  t_start = " << t_start << std::endl;
    }
}

void FEL_1D::init(Beam *beam)
{
    // sanity checks
    if (1e-9 < fabs(beam->getTimeStep()-dt)/dt)
        throw(IOexception("FEL1D::init() - time-step mismach to beam"));
        
    // the field after this step
    // it is shifted in time by one index
    std::vector<double> next_E = std::vector<double>(N_field, 0.0);
    
    // index 0 of next_E (the head) always remains zero (no field moving in from the front)
    // index 1 of next_E is index 0 of field_E
    next_E[1] = field_E[0];
    // @todo check if the derivatives really can be neglected
    // next_E[1] = field_E[1];
    // we loop over the indices of the current field
    for (int i=1; i<N_field-1; i++)
        next_E[i+1] = field_E[i-1] + field_E[i+1] - previous_E[i-1];
    
    // the time step has been computed, now we move forward
    N_steps++;
    time += dt;
    head += prop*ds;
    previous_E = field_E;
    field_E = next_E;
}

void FEL_1D::step(Beam *beam)
{
    // check if the beam is at the same time as the field
    int step_beam = beam->getNOTS();
    if (fabs(step_beam*dt-time)>0.25*dt)
    {
        std::cout << "FEL1D::step() - time mismatch" << std::endl;
        std::cout << "  step=" << N_steps << "   time=" << 1e9*time << " ns" << std::endl;
        std::cout << "  beam is at step=" << step_beam << "  t=" << 1e9*step_beam*dt << " ns" << std::endl;
        std::cout << "  first bunch particles average time t=" << 1e9*beam->getBunch(0)->avgTime() << " ns" << std::endl;
        throw(IOexception("FEL1D::step() - beam is at different time"));
    }

    // optical mode size
    double z = N_steps * ds;
    double w = w_0 * sqrt(1.0 + pow((z-z_waist)/z_Rayleigh,2));
    double A_opt = 0.5*Pi*w*w;
    
    // the transverse current density
    std::vector<double> J_x = std::vector<double>(N_field, 0.0);
    
    //! @todo get all particle coordinates
    //! @todo add the beam-induced fields
    int NOB = beam->getNOB();
    for(int ib=0; ib<NOB; ib++)
    {
        Bunch* B = beam->getBunch(ib);
        int NOP = B->getNOP();
        for(int ip=0; ip<NOP; ip++)
        {
            ChargedParticle* P = B->getParticle(ip);
            double t = P->getTime();
            //! @todo the particle time should be at the center of this time-step
            //! @todo actually the beam should be half a time step ahead of the field here
            if (0.25*dt < fabs(t-time))
            {
                std::cout << "FEL1D::step() - time mismatch" << std::endl;
                std::cout << "  step=" << N_steps << "   time =" << 1e9*time << " ns";
                std::cout << "   dt=" << 1e9*dt << " ns" << std::endl;
                std::cout << "  particle is at t=" << 1e9*t << " ns" << std::endl;
                throw(IOexception("FEL1D::step() - time mismatch"));
            }
            Vector X = P->getPosition();
            // find the interaction field (slice) index from the longitudinal particle position
            // small time deviations in time are ignored in the field computation (no interpolation)
            double delta_z = dot(head-X,prop);
            int index = rint(delta_z/ds);
            if ((index<0) || (index>=N_field))
            {
                std::cout << "FEL1D::step() - particle out off grid bounds" << std::endl;
                std::cout << "  step=" << N_steps << "  delta_z=" << delta_z << "   index=" << index << std::endl;
                throw(IOexception("FEL1D::step() - particle out off grid bounds"));
            }
            Vector BetaGamma = P->getMomentum();
            double bg2 = BetaGamma.abs2nd();
            double gamma = sqrt(bg2 + 1.0);
            Vector Beta_prime = P->getAccel()/gamma;
            double Q = P->getCharge();
            double j_x = dot(Beta_prime,e_x)*Q/A_opt/dt*MuNull*SpeedOfLight*SpeedOfLight;
            J_x[index] += j_x;
        }
    }
    
    // the field after this step
    // it is shifted in time by one index
    std::vector<double> next_E = std::vector<double>(N_field, 0.0);
    
    // index 0 of next_E (the head) always remains zero (no field moving in from the front)
    // index 1 of next_E is index 0 of field_E
    next_E[1] = field_E[0];
    // @todo check if the derivatives really can be neglected
    // next_E[1] = field_E[1];
    // we loop over the indices of the current field
    for (int i=1; i<N_field-1; i++)
        next_E[i+1] = field_E[i-1] + field_E[i+1] - previous_E[i-1];

    // the time step has been computed, now we move forward
    N_steps++;
    time += dt;
    head += prop*ds;
    previous_E = field_E;
    field_E = next_E;
    
    //! @todo we have to scale the fields according to the change in mode size
    
    // every number of steps we record the fields
    if (0 == N_steps % step_Output)
    {
        field_storage.push_back(field_E);
        J_storage.push_back(J_x);
    }
}

ElMagField FEL_1D::Field(double t, Vector X)
{
    // issue a warning if the query is too far in time (quarter step) from the current time step
    if (0.25*dt < fabs(t-time))
    {
        std::cout << "FEL1D::Field() - time mismatch" << std::endl;
        std::cout << "  step=" << N_steps << "   time=" << 1e9*time << " ns" << "   dt=" << 1e9*dt << " ns" << std::endl;
        std::cout << "  requested t=" << 1e9*t << " ns" << std::endl;
        throw(IOexception("FEL1D::Field() - time mismatch"));
    }

    // small time deviations in time are ignored in the field computation (no interpolation)
    double delta_z = dot(head-X,prop);
    int index = rint(delta_z/ds);
    if ((index<0) || (index>=N_field))
    {
        std::cout << "FEL1D::Field() - particle out off grid bounds" << std::endl;
        std::cout << "  step=" << N_steps << "  delta_z=" << delta_z << "   index=" << index << std::endl;
        throw(IOexception("FEL1D::Field() - particle out off grid bounds"));
    }
    
    //! todo the transverse shape of the field is not yet included
    
    // field is discretized to grid (not interpolated)
    Vector E = e_x * field_E[index];
    Vector B = e_y * field_E[index]/SpeedOfLight;
    
    return ElMagField(E,B);
}

void FEL_1D::write_output()
{
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
        
        // This declares a lambda, which can be called just like a function
        // some external variables are captured by reference
        auto set_int_attribute = [&pdset](std::string name, int* data) 
        { 
            hid_t atts  = H5Screate(H5S_SCALAR);
            if (atts<0) throw(IOexception(("FEL1D::write_output - error in H5Screate("+name+")").c_str()));
            hid_t attr = H5Acreate2(pdset, name.c_str(), H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
            if (attr<0) throw(IOexception(("FEL1D::write_output - error in H5Acreate2("+name+")").c_str()));
            herr_t status = H5Awrite(attr, H5T_NATIVE_INT, data);
            if (status<0) throw(IOexception(("FEL1D::write_output - error in H5Awrite("+name+")").c_str()));
            status = H5Sclose (atts);
            if (status<0) throw(IOexception(("FEL1D::write_output - error in H5Sclose("+name+")").c_str()));
        };
        auto set_double_attribute = [&pdset](std::string name, double* data) 
        { 
            hid_t atts  = H5Screate(H5S_SCALAR);
            if (atts<0) throw(IOexception(("FEL1D::write_output - error in H5Screate("+name+")").c_str()));
            hid_t attr = H5Acreate2(pdset, name.c_str(), H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
            if (attr<0) throw(IOexception(("FEL1D::write_output - error in H5Acreate2("+name+")").c_str()));
            herr_t status = H5Awrite(attr, H5T_NATIVE_DOUBLE, data);
            if (status<0) throw(IOexception(("FEL1D::write_output - error in H5Awrite("+name+")").c_str()));
            status = H5Sclose (atts);
            if (status<0) throw(IOexception(("FEL1D::write_output - error in H5Sclose("+name+")").c_str()));
        };

        // attach scalar attributes
        set_int_attribute("N_field", &N_field);
        set_int_attribute("N_steps", &N_steps);
        set_int_attribute("N_output", &step_Output);
        set_double_attribute("dt", &dt);
        set_double_attribute("w0", &w_0);
        set_double_attribute("z_Rayl", &z_Rayleigh);
        set_double_attribute("z_waist", &z_waist);
        set_double_attribute("origin.x", &origin.x);
        set_double_attribute("origin.y", &origin.y);
        set_double_attribute("origin.z", &origin.z);
        set_double_attribute("prop.x", &dz.x);
        set_double_attribute("prop.y", &dz.y);
        set_double_attribute("prop.z", &dz.z);
        
        // Close and release resources.
        status = H5Pclose (pdcpl);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Pclose(pdcpl)"));
        status = H5Dclose (pdset);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Dclose(pdset)"));
        status = H5Sclose (pspace);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Dclose(pspace)"));
        delete[] buffer;

        status = H5Fclose(file);
        if (status<0) throw(IOexception("FEL1D::write_output - error in H5Fclose()"));
        // no errors have occured if we made it 'til here
        cout << "writing HDF5 done." << endl;
    }
}
    
FEL_1D::~FEL_1D()
{
}

