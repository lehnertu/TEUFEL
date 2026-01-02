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

#include "csr2d.h"
#include "global.h"

#include <iostream>
#include <math.h>
#include "hdf5.h"

CSR_2D::CSR_2D()
{
    is_initialized = false;
    N_long = 0;
    N_trans = 0;
    e_normal = Vector(0.0, 1.0, 0.0);
    createOutput = false;
    step_Output = 1;
}

CSR_2D::CSR_2D(
    const pugi::xml_node node,
    InputParser *parser )
{
    // no source definet yet
    is_initialized = false;
    pugi::xml_attribute att = node.attribute("N_long");
    if (!att)
        throw(IOexception("InputParser::CSR_2D - attribute N_long not found."));
    N_long = parser->parseInt(att);
    att = node.attribute("N_trans");
    if (!att)
        throw(IOexception("InputParser::CSR_2D - attribute N_trans not found."));
    N_trans = parser->parseInt(att);
    pugi::xml_node vec = node.child("normal");
    if (!vec)
        throw(IOexception("InputParser::CSR_2D - <normal> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(vec.attribute("x"));
        y = parser->parseDouble(vec.attribute("y"));
        z = parser->parseDouble(vec.attribute("z"));
        e_normal = Vector(x,y,z);
        e_normal.normalize();
    }
    // define file output if requested
    pugi::xml_node lognode = node.child("log");
    if (lognode)
    {
        pugi::xml_attribute fn = lognode.attribute("file");
        if (!fn) throw(IOexception("InputParser::CSR_2D - <log> filename for log not found."));
        FileName = fn.as_string();
        step_Output = 1;
        pugi::xml_attribute st = lognode.attribute("step");
        if (st) step_Output = st.as_int();
        createOutput = true;
    } else {
        createOutput = false;
    }
    if (teufel::rank==0)
    {
        std::cout << "CSR-2D interaction" << std::endl;
        std::cout << "  " << N_long << " x " << N_trans << " grid,  ";
        std::cout << "normal = (" << e_normal.x << ", " << e_normal.y << ", " << e_normal.z << ")";
        std::cout << std::endl;
    }
}

void CSR_2D::init()
{
    if (is_initialized)
        throw(IOexception("error - CSR_2D::init() called twice."));
    if ( (N_long>1) and (N_trans>1))
    {
        interaction_field = new ElMagField[N_long*N_trans];
        if (teufel::rank==0)
        {
            std::cout << "CSR-2D::init() allocated " << N_long << " x " << N_trans << " grid." << std::endl;
        }
    }
    else
    {
        throw(IOexception("CSR_2D::init() - cannot allocate field map of zero size."));
    }
    is_initialized = true;
    step_counter = 0;
}

CSR_2D::~CSR_2D()
{
    // free field map memory
    if (is_initialized)
        delete[] interaction_field;
    // free output storage memory
    for (ElMagField *ptr : field_storage) delete ptr;
    for (Vector *ptr : position_storage) delete ptr;
    for (Vector *ptr : momentum_storage) delete ptr;
}

void CSR_2D::update(Beam *beam, double tracking_time)
{
    update_time = tracking_time;
    //! @todo remove debugging check
    if (DEBUGLEVEL>=2)
    {   
        if (std::isnan(update_time))
            throw std::runtime_error("CSR_2D::update(): update_time value is NaN!");
    }
    /*
    if (teufel::rank==0)
    {
        std::cout << "CSR-2D::update() at tracking time " << update_time << " s" << std::endl;
    }
    */
    // get the particle coordinates
    NOP = beam->getNOP();
    double *buffer = new double[NOP*6];
    beam->bufferCoordinates(buffer, NOP);
    Vector *position = new Vector[NOP];
    Vector *momentum = new Vector[NOP];
    double *bp = buffer;
    for(int i=0; i<NOP; i++)
    {
        position[i].x = *bp++;
        position[i].y = *bp++;
        position[i].z = *bp++;
        momentum[i].x = *bp++;
        momentum[i].y = *bp++;
        momentum[i].z = *bp++;
    };
    // determine orientation of the field map
    Vector e_long = VectorZero;
    for(int i=0; i<NOP; i++) e_long += momentum[i];
    e_long -= e_normal*dot(e_long,e_normal);
    e_long.normalize();
    Vector e_trans = cross(e_normal,e_long);
    e_trans.normalize();
    // determine the extensions of the needed field map
    double l_min = 1.0e300;
    double l_max = -1.0e300;
    double t_min = 1.0e300;
    double t_max = -1.0e300;
    double n_avg = 0.0;
    for(int i=0; i<NOP; i++)
    {
        double l = dot(position[i],e_long);
        if (l<l_min) l_min=l;
        if (l>l_max) l_max=l;
        double t = dot(position[i],e_trans);
        if (t<t_min) t_min=t;
        if (t>t_max) t_max=t;
        n_avg += dot(position[i],e_normal);
    }
    n_avg /= NOP;
    // the field map extends 10% beyond the range of particle coordinates
    double l_orig = 0.5*(l_max+l_min) - 0.6*(l_max-l_min);
    double t_orig = 0.5*(t_max+t_min) - 0.6*(t_max-t_min);
    origin = e_normal*n_avg + e_long*l_orig + e_trans*t_orig;
    // there are N grid nodes and N-1 cells in each direction
    // cell size
    double l_long = 1.2*(l_max-l_min)/(N_long-1);
    double l_trans = 1.2*(t_max-t_min)/(N_trans-1);
    // grid vector
    d_long = e_long * l_long;
    d_trans = e_trans * l_trans;
    
    // zero the interaction field
    for (int il=0; il<N_long; il++)
        for (int it=0; it<N_trans; it++)
            interaction_field[il*N_trans+it] = ElMagFieldZero;
    // compute fields for all particles of all bunches in the beam
    for(int i_bunch=0; i_bunch<beam->getNOB(); i_bunch++)
    {
        Bunch *b = beam->getBunch(i_bunch);
        for(int i_part=0; i_part<b->getNOP(); i_part++)
        {
            ChargedParticle* p = b->getParticle(i_part);
            // the current particle position
            Vector X = p->getPosition();
            // fractional cell index of particle position
            double f_long = dot(X-origin, e_long)/l_long;
            double index_l = std::floor(f_long);
            f_long -= index_l;
            double f_trans = dot(X-origin, e_trans)/l_trans;
            double index_t = std::floor(f_trans);
            f_trans -= index_t;
            // the intercation field shall be computed as if the particle was
            // at the center position of one grid cell
            // we cannot esily shift the whole trajectory so we
            // shift the grid nodes in opposite direction just for the computation
            Vector shift = (d_long*f_long+d_trans*f_trans) - (d_long*0.5+d_trans*0.5);
            //! @todo remove debugging check
            if (DEBUGLEVEL>=2)
            {   
                if (std::isnan(shift.x) || std::isnan(shift.y) || std::isnan(shift.z))
                    throw std::runtime_error("CSR_2D::update(): shift value is NaN!");
            }
            //! @todo OMP parallelize the loop over the grid
            for (int il=0; il<N_long; il ++)
                for (int it=0; it<N_trans; it ++)
                {
                    //! @todo use the shifted grid position - temporarily removed for debugging
                    // it does not change the outcome - only the artifacts appear at different positions
                    // Vector grid = origin + d_long*il + d_trans*it + shift;
                    Vector grid = origin + d_long*il + d_trans*it;
                    
                    //! @todo remove debugging check
                    if (DEBUGLEVEL>=2)
                    {   
                        if (std::isnan(grid.x) || std::isnan(grid.y) || std::isnan(grid.z))
                            throw std::runtime_error("CSR_2D::update(): grid value is NaN!");
                    }
                    ElMagField p_field = p->RetardedField(tracking_time, grid);
                    //! @todo remove debugging check
                    if (DEBUGLEVEL>=2)
                    {   
                        double *check_ptr = (double *)&p_field;
                        for (int check_i=0; check_i<6; check_i++)
                            if (std::isnan(*check_ptr++))
                            {
                                std::cout << "CSR_2D::update(): step=" << step_counter;
                                std::cout << " i_bunch=" << i_bunch << " i_part=" << i_part << std::endl;
                                throw std::runtime_error("CSR_2D::update(): p_field value is NaN!");
                            }
                    }
                    // the interaction field
                    interaction_field[il*N_trans+it] += p_field;
                };
        }
    }
    
    // when requested store field maps for file output
    if (0 == (step_counter % step_Output))
    {   
        // allocate new field map
        ElMagField *map = new ElMagField[N_long*N_trans];
        // copy current map data
        memcpy(map, interaction_field, N_long*N_trans*sizeof(ElMagField));
        // store pointer
        field_storage.push_back(map);
        // handle particle position
        Vector *pos = new Vector[NOP];
        memcpy(pos, position, NOP*sizeof(Vector));
        position_storage.push_back(pos);
        // handle particle momentum
        Vector *mom = new Vector[NOP];
        memcpy(mom, momentum, NOP*sizeof(Vector));
        momentum_storage.push_back(mom);
        // handle the field map geometry
        origin_storage.push_back(origin);
        dl_storage.push_back(d_long);
        dt_storage.push_back(d_trans);
        // the arrays remain allocated as part of the field_storage array
    }
    
    step_counter ++;
    delete[] buffer;
    delete[] position;
    delete[] momentum;
}

ElMagField CSR_2D::Field(double t, Vector X)
{
    //! @todo returning zero fields for now
    return ElMagField();
}

void CSR_2D::write_output()
{
    if (createOutput)
    {
        cout << "CSR_2D : writing interaction field data to " << FileName << endl;

        //! @todo write actual data

        herr_t status;
        // Create a new file using the default properties.
        hid_t file = H5Fcreate (FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file<0) throw(IOexception("CSR_2D::write_output() - error in H5Fcreate()"));

        // -------------------------------------------------------------------------------
        
        // Create dataspace for the field map positions (origin, d_long, d_trans).
        // Setting maximum size to NULL sets the maximum size to be the current size.
        int NOS = origin_storage.size();
        hsize_t g_dims[3];
        g_dims[0] = NOS;
        g_dims[1] = 3;
        g_dims[2] = 3;
        hid_t g_space = H5Screate_simple (3, g_dims, NULL);
        if (g_space<0) throw(IOexception("CSR_2D::write_output() - error in H5Screate(g_space)"));

        // buffer the data
        double *g_buffer = new double[NOS*3*sizeof(Vector)];
        Vector *bp = (Vector *)g_buffer;
        for(int i=0; i<NOS; i++)
        {
            *bp++ = origin_storage[i];
            *bp++ = dl_storage[i];
            *bp++ = dt_storage[i];
        }

        // Create the dataset creation property list
        hid_t g_dcpl = H5Pcreate (H5P_DATASET_CREATE);
        if (g_dcpl<0) throw(IOexception("CSR_2D::write_output() - error in H5Pcreate(g_dcpl)"));
        // Create the dataset.
        hid_t g_dset = H5Dcreate(file,
            "Geometry",	     	    // dataset name
            H5T_NATIVE_DOUBLE,		// data type
            g_space, H5P_DEFAULT,
            g_dcpl, H5P_DEFAULT);
        if (g_dset<0) throw(IOexception("CSR_2D::write_output() - error in H5Dcreate(g_dset)"));
        // Write the data to the dataset
        status = H5Dwrite (g_dset,
            H5T_NATIVE_DOUBLE, 		// mem type id
            H5S_ALL, 			    // mem space id
            g_space,
            H5P_DEFAULT,			// data transfer properties
            g_buffer);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Dwrite(g_dset)"));

        // attach scalar attributes
        hid_t atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("CSR_2D::write_output() - error in H5Screate(N_steps)"));
        hid_t attr = H5Acreate2(g_dset, "N_steps", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("CSR_2D::write_output() - error in H5Acreate2(N_steps)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &NOS);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Awrite(N_steps)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Sclose(N_steps)"));

        atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("CSR_2D::write_output() - error in H5Screate(N_long)"));
        attr = H5Acreate2(g_dset, "N_long", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("CSR_2D::write_output() - error in H5Acreate2(N_long)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &N_long);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Awrite(N_long)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Sclose(N_long)"));

        atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("CSR_2D::write_output() - error in H5Screate(N_trans)"));
        attr = H5Acreate2(g_dset, "N_trans", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("CSR_2D::write_output() - error in H5Acreate2(N_trans)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &N_trans);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Awrite(N_trans)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Sclose(N_trans)"));
        
        // Close and release resources.
        status = H5Pclose (g_dcpl);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Pclose(g_dcpl)"));
        status = H5Dclose (g_dset);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Dclose(g_dset)"));
        status = H5Sclose (g_space);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Dclose(g_space)"));

        // -------------------------------------------------------------------------------
        
        // Create dataspace for the particle positions.
        hsize_t p_dims[3];
        p_dims[0] = NOS;
        p_dims[1] = NOP;
        p_dims[2] = 6;
        hid_t p_space = H5Screate_simple (3, p_dims, NULL);
        if (p_space<0) throw(IOexception("CSR_2D::write_output() - error in H5Screate(p_space)"));

        // buffer the data
        double *p_buffer = new double[NOS*NOP*6];
        Vector *pp = (Vector *)p_buffer;
        for(int step=0; step<NOS; step++)
        {
            Vector *ps = position_storage[step];
            Vector *ms = momentum_storage[step];
            for(int part=0; part<NOP; part++)
            {
                *pp++ = ps[part];
                *pp++ = ms[part];
            }
        }

        // Create the dataset creation property list
        hid_t p_dcpl = H5Pcreate (H5P_DATASET_CREATE);
        if (p_dcpl<0) throw(IOexception("CSR_2D::write_output() - error in H5Pcreate(p_dcpl)"));
        // Create the dataset.
        hid_t p_dset = H5Dcreate(file,
            "Particles",    	    // dataset name
            H5T_NATIVE_DOUBLE,		// data type
            p_space, H5P_DEFAULT,
            p_dcpl, H5P_DEFAULT);
        if (p_dset<0) throw(IOexception("CSR_2D::write_output() - error in H5Dcreate(p_dset)"));
        // Write the data to the dataset
        status = H5Dwrite (p_dset,
            H5T_NATIVE_DOUBLE, 		// mem type id
            H5S_ALL, 			    // mem space id
            p_space,
            H5P_DEFAULT,			// data transfer properties
            p_buffer);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Dwrite(p_dset)"));

        // attach scalar attributes
        atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("CSR_2D::write_output() - error in H5Screate(NOP)"));
        attr = H5Acreate2(p_dset, "NOP", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("CSR_2D::write_output() - error in H5Acreate2(NOP)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &NOP);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Awrite(NOP)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Sclose(NOP)"));

        // Close and release resources.
        status = H5Pclose (p_dcpl);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Pclose(p_dcpl)"));
        status = H5Dclose (p_dset);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Dclose(p_dset)"));
        status = H5Sclose (p_space);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Dclose(p_space)"));

        // -------------------------------------------------------------------------------

        // Create dataspace for the field maps.
        hsize_t f_dims[4];
        f_dims[0] = NOS;
        f_dims[1] = N_long;
        f_dims[2] = N_trans;
        f_dims[3] = 6;
        hid_t f_space = H5Screate_simple (4, f_dims, NULL);
        if (f_space<0) throw(IOexception("CSR_2D::write_output() - error in H5Screate(f_space)"));

        // buffer the data
        double *f_buffer = new double[NOS*N_long*N_trans*6];
        ElMagField *fp = (ElMagField *)f_buffer;
        for(int step=0; step<NOS; step++)
        {
            ElMagField *map = field_storage[step];
            memcpy(fp, map, N_long*N_trans*sizeof(ElMagField));
            fp += N_long*N_trans;
        }
        
        // Create the dataset creation property list
        hid_t f_dcpl = H5Pcreate (H5P_DATASET_CREATE);
        if (f_dcpl<0) throw(IOexception("CSR_2D::write_output() - error in H5Pcreate(f_dcpl)"));
        // Create the dataset.
        hid_t f_dset = H5Dcreate(file,
            "Fields",       	    // dataset name
            H5T_NATIVE_DOUBLE,		// data type
            f_space, H5P_DEFAULT,
            f_dcpl, H5P_DEFAULT);
        if (f_dset<0) throw(IOexception("CSR_2D::write_output() - error in H5Dcreate(f_dset)"));
        // Write the data to the dataset
        status = H5Dwrite (f_dset,
            H5T_NATIVE_DOUBLE, 		// mem type id
            H5S_ALL, 			    // mem space id
            f_space,
            H5P_DEFAULT,			// data transfer properties
            f_buffer);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Dwrite(f_dset)"));

        // attach scalar attributes
        atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("CSR_2D::write_output() - error in H5Screate(NOS)"));
        attr = H5Acreate2(f_dset, "NOS", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("CSR_2D::write_output() - error in H5Acreate2(NOS)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &NOS);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Awrite(NOS)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Sclose(NOS)"));

        // attach scalar attributes
        atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("CSR_2D::write_output() - error in H5Screate(N_long)"));
        attr = H5Acreate2(f_dset, "N_long", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("CSR_2D::write_output() - error in H5Acreate2(N_long)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &N_long);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Awrite(N_long)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Sclose(N_long)"));

        // attach scalar attributes
        atts  = H5Screate(H5S_SCALAR);
        if (atts<0) throw(IOexception("CSR_2D::write_output() - error in H5Screate(N_trans)"));
        attr = H5Acreate2(f_dset, "N_trans", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
        if (attr<0) throw(IOexception("CSR_2D::write_output() - error in H5Acreate2(N_trans)"));
        status = H5Awrite(attr, H5T_NATIVE_INT, &N_trans);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Awrite(N_trans)"));
        status = H5Sclose (atts);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Sclose(N_trans)"));

        // Close and release resources.
        status = H5Pclose (f_dcpl);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Pclose(f_dcpl)"));
        status = H5Dclose (f_dset);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Dclose(f_dset)"));
        status = H5Sclose (f_space);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Dclose(f_space)"));

        // -------------------------------------------------------------------------------
        
        status = H5Fclose(file);
        if (status<0) throw(IOexception("CSR_2D::write_output() - error in H5Fclose()"));

        delete[] g_buffer;
        delete[] p_buffer;
        
        // no errors have occured if we made it 'til here
        cout << "writing HDF5 done." << endl;
    }
}

