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

#include "source_screen.h"
#include "global.h"

#include <iostream>
#include <math.h>
#include <omp.h>
#include "hdf5.h"

SourceScreen::SourceScreen(
        std::string filename,
        Vector position,
        double time) :
    LocalizedField(time,position)
{
    FileName = filename;
    if (teufel::rank==0) std::cout << "screen field source : " << filename << std::endl;
    herr_t status;
    hid_t dataset, attr;
    if (teufel::rank==0) std::cout << "reading HDF5 file " << filename << std::endl;
    hid_t file = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file<0) throw(IOexception("SourceScreen - error opening the file."));
    
    // read the screen geometry information
    dataset = H5Dopen2(file, "Screen", H5P_DEFAULT);
    if (dataset<0) throw(IOexception("SourceScreen - error opening dataset Screen"));
    attr = H5Aopen_name(dataset, "Nx");
    if (attr<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute Nx"));
    status =  H5Aread(attr, H5T_NATIVE_INT, &Nx);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute Nx"));
    status = H5Aclose(attr);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute Nx"));
    attr = H5Aopen_name(dataset, "Ny");
    if (attr<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute Ny"));
    status =  H5Aread(attr, H5T_NATIVE_INT, &Ny);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute Ny"));
    status = H5Aclose(attr);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute Ny"));
    attr = H5Aopen_name(dataset, "NOTS");
    if (attr<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute NOTS"));
    status =  H5Aread(attr, H5T_NATIVE_INT, &NOTS);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute NOTS"));
    status = H5Aclose(attr);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute NOTS"));

    // read the timing information
    attr = H5Aopen_name(dataset, "dt");
    if (attr<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute dt"));
    status =  H5Aread(attr, H5T_NATIVE_DOUBLE, &dt_src);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute dt"));
    status = H5Aclose(attr);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute dt"));
    attr = H5Aopen_name(dataset, "dtx");
    if (attr<0)
        {
            dtx_src = 0.0;
        }
    else
        {
            status =  H5Aread(attr, H5T_NATIVE_DOUBLE, &dtx_src);
            if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute dtx"));
            status = H5Aclose(attr);
            if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute dtx"));
        }
    attr = H5Aopen_name(dataset, "dty");
    if (attr<0)
        {
            dty_src = 0.0;
        }
    else
        {
            status =  H5Aread(attr, H5T_NATIVE_DOUBLE, &dty_src);
            if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute dty"));
            status = H5Aclose(attr);
            if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute dty"));
        }
        
    // read the vectors
    double* buf = new double[4*3];
    if (buf==0) throw(IOexception("SourceScreen - error allocating memory."));
    status = H5Dread (dataset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			    // mem space id
        H5S_ALL,
        H5P_DEFAULT,			// data transfer properties
        buf);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen"));
    // the origin from the file is ignored and replaced by the origin given to the constructor
    normal = Vector(buf[3],buf[4],buf[5]);
    dX = Vector(buf[6],buf[7],buf[8]);
    dY = Vector(buf[9],buf[10],buf[11]);
    // coax the grid vectors to form an orthogonal system
    // these values all correspond to the original orientation of the screen
    normal.normalize();
    dx = dX.norm();
    dX -= normal*dot(dX,normal);
    dX.normalize();
    ex = dX;
    dX *= dx;
    dy = dY.norm();
    dY -= normal*dot(dY,normal);
    dY -= dX*dot(dY,dX);
    dY.normalize();
    ey = dY;
    dY *= dy;
    delete buf;
    
    // done with the screen dataset
    status = H5Dclose(dataset);
    if (status<0) throw(IOexception("SourceScreen - error closing dataset Screen"));

    if (teufel::rank==0)
    {
        std::cout << "Nx=" << Nx << "  Ny=" << Ny << "  Nots=" << NOTS << std::endl;
        std::cout << "Origin = (" << origin.x << ", " << origin.y << ", " << origin.z << ") m" << std::endl;
        std::cout << "Normal = (" << normal.x << ", " << normal.y << ", " << normal.z << ") m" << std::endl;
        std::cout << "dX = (" << dX.x << ", " << dX.y << ", " << dX.z << ") m" << std::endl;
        std::cout << "dY = (" << dY.x << ", " << dY.y << ", " << dY.z << ") m" << std::endl;
        std::cout << "t0=" << t0 << "  dt=" << dt_src << "  dtx=" << dtx_src << "  dty=" << dty_src << std::endl;
        std::cout << "allocating " << (double)(3*6*Nx*Ny*NOTS) * sizeof(double) / 1e6 << " MB of memory per node" << std::endl;
        std::cout << std::endl;
    }    

    // read the field dataset
    ElMagField* field_buf = new ElMagField[getBufferSize()];
    if (field_buf==0) throw(IOexception("SourceScreen - error allocating memory."));
    dataset = H5Dopen2(file, "ElMagField", H5P_DEFAULT);
    
    // sanity check the field sizes
    hid_t dspace = H5Dget_space(dataset);
    int ndims = H5Sget_simple_extent_ndims(dspace);
    if (ndims != 4) throw(IOexception("SourceScreen - error reading dataset ElMagField - dimension != 4"));
    hsize_t dims[4];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    if (dims[0] != Nx) throw(IOexception("SourceScreen - error reading dataset ElMagField - size[0] != Nx"));
    if (dims[1] != Ny) throw(IOexception("SourceScreen - error reading dataset ElMagField - size[1] != Ny"));
    if (dims[2] != NOTS) throw(IOexception("SourceScreen - error reading dataset ElMagField - size[2] != NOTS"));
    if (dims[3] != 6) throw(IOexception("SourceScreen - error reading dataset ElMagField - size[3] != 6"));
    
    // buffer the data
    if (dataset<0) throw(IOexception("SourceScreen - error opening dataset ElMagField"));
    status = H5Dread (dataset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			    // mem space id
        H5S_ALL,
        H5P_DEFAULT,			// data transfer properties
        field_buf);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset ElMagField"));
    status = H5Dclose(dataset);
    if (status<0) throw(IOexception("SourceScreen - error closing dataset ElMagField"));
    
    // set the field sizes and
    // transfer the data into the internal data structures
    ElMagField *b = field_buf;
    Traces.resize(Nx);
    for (unsigned int ix = 0; ix < Nx; ix++) {
    	Traces[ix].resize(Ny);
    	for (unsigned int iy = 0; iy < Ny; iy++)
    	{
    	    // the field trace is taken from the file data
    	    FieldTrace* tr = new FieldTrace(CellTimeZero(ix,iy),dt_src,NOTS);
    	    if (tr==0) throw(IOexception("SourceScreen - error allocating memory."));
    	    tr -> set_buffer(b, NOTS);
    	    // the field data are expressed in screen-local coordinates
    	    tr -> transform(ex,ey,normal);
    	    // coax the fields to be purely transversal
    	    // (cancel the z component in the screen-local frame)
    	    tr -> cancel_long_fields(Vector(0.0,0.0,1.0));
    	    Traces[ix][iy] = tr;
    	    b+=NOTS;
    	};
    }
    
    // done with the field data
    delete field_buf;
        
    // done with the file
    status = H5Fclose (file);
    if (status<0) throw(IOexception("SourceScreen - error closing the file."));

    // compute the time derivatives
    dt_Traces.resize(Nx);
    for (unsigned int ix = 0; ix < Nx; ix++)
    {
    	dt_Traces[ix].resize(Ny);
    	for (unsigned int iy = 0; iy < Ny; iy++)
    	{
            // compute the time derivative (includes allocation of a trace)
    	    FieldTrace* dt_trace = Traces[ix][iy]->derivative();
    	    dt_Traces[ix][iy] = dt_trace;
	    }
    }
        
    if (teufel::rank==0)
    {
        std::cout << "total energy on source screen " << totalEnergy() << " J" << std::endl;
        std::cout << std::endl;
    }  
    
    // when we are done we write the field to files
    // for debugging purpose only
    // TODO: remove this
    WriteField("SourceScreenTraces_DEBUG.h5", Traces);
    WriteField("SourceScreenTimeDerivative_DEBUG.h5", dt_Traces);
  
}

SourceScreen::~SourceScreen()
{
    for (unsigned int ix = 0; ix < Nx; ix++)
    	for (unsigned int iy = 0; iy < Ny; iy++)
    	    {
                delete Traces[ix][iy];
                delete dt_Traces[ix][iy];
            };
}

Vector SourceScreen::CellPosition(unsigned int ix, unsigned int iy)
{
    // center index is Nx/2 Ny/2 (integer division!)
    // the origin in local coordinaes is always zero
    return dX*(double)((int)ix-(int)Nx/2) + dY*(double)((int)iy-(int)Ny/2);
}

double SourceScreen::CellTimeZero(unsigned int ix, unsigned int iy)
{
    // center index is Nx/2 Ny/2 (integer division!)
    // the center start time in local coordinates is zero
    return dtx_src*(double)((int)ix-(int)Nx/2) + dty_src*(double)((int)iy-(int)Ny/2);
}

ElMagField SourceScreen::LocalField(double t, Vector X)
{
    // integrate over the screen to collect all contributions
    // to the field at the requested point
    ElMagField total = ElMagFieldZero;
    for (unsigned int ix = 0; ix < Nx; ix++)
    	for (unsigned int iy = 0; iy < Ny; iy++)
        {
            Vector source_pos = CellPosition(ix,iy);
            FieldTrace* source_trace = Traces[ix][iy];
            Vector RVec = X - source_pos;
            double R = RVec.norm();
            // consider retardation
            double source_t = t - R/SpeedOfLight;
            if ((source_t>=source_trace->get_t0()) and (source_t<=source_trace->get_last_time()))
            {
                double R2 = R*R;
                double R3 = R2*R;
                ElMagField c1 = source_trace->get_field(source_t) * (-dot(RVec,normal)/R3);
                ElMagField c2 = dt_Traces[ix][iy]->get_field(source_t) * (1.0-dot(RVec,normal)/R) / (R*SpeedOfLight);
                total += c1+c2 ;
            };
        }
    return total*dx*dy/(4.0*Pi);
}

double SourceScreen::totalEnergy()
{
    double total = 0.0;
    for (unsigned int ix = 0; ix < Nx; ix++)
    	for (unsigned int iy = 0; iy < Ny; iy++)
            {
                Vector S = Traces[ix][iy]->Poynting();
                total -= dot(S,normal)*dX.norm()*dY.norm();
            }
    return total;
}

void SourceScreen::WriteField(
    std::string filename,
    std::vector<std::vector<FieldTrace*>> Field)
{
    hid_t atts, attr;
    herr_t status;
    cout << "writing HDF5 file " << filename << endl;
    // Create a new file using the default properties.
    hid_t file = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file<0) throw(IOexception("SourceScreen::WriteField - error in H5Fcreate()"));

    // Create dataspace for the screen metadata.
    // data field 3 vectors
    hsize_t s_dims[2];
    s_dims[0] = 4;
    s_dims[1] = 3;
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hid_t s_space = H5Screate_simple (2, s_dims, NULL);
    if (s_space<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(s_space)"));
    // buffer the data
    double *s_buffer = new double[4*3];
    double *s_bp = s_buffer;
    *s_bp++ = origin.x;
    *s_bp++ = origin.y;
    *s_bp++ = origin.z;
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
    if (s_dcpl<0) throw(IOexception("SourceScreen::WriteField - error in H5Pcreate(s_dcpl)"));
    // Create the dataset.
    hid_t s_dset = H5Dcreate(file,
        "Screen",		        // dataset name
        H5T_NATIVE_DOUBLE,		// data type
        s_space, H5P_DEFAULT,
        s_dcpl, H5P_DEFAULT);
    if (s_dset<0) throw(IOexception("SourceScreen::WriteField - error in H5Dcreate(s_dset)"));
    // Write the data to the dataset
    status = H5Dwrite (s_dset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			    // mem space id
        s_space,
        H5P_DEFAULT,			// data transfer properties
        s_buffer);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Dwrite(s_dset)"));
    // attach scalar attribute Nx
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(Nx)"));
    attr = H5Acreate2(s_dset, "Nx", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("SourceScreen::WriteField - error in H5Acreate2(Nx)"));
    status = H5Awrite(attr, H5T_NATIVE_INT, &Nx);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Awrite(Nx)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Sclose(Nx)"));
    // attach scalar attribute Nx
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(Ny)"));
    attr = H5Acreate2(s_dset, "Ny", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("SourceScreen::WriteField - error in H5Acreate2(Ny)"));
    status = H5Awrite(attr, H5T_NATIVE_INT, &Ny);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Awrite(Ny)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Sclose(Ny)"));
    // attach scalar attribute NOTS
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(NOTS)"));
    attr = H5Acreate2(s_dset, "NOTS", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("SourceScreen::WriteField - error in H5Acreate2(NOTS)"));
    status = H5Awrite(attr, H5T_NATIVE_INT, &NOTS);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Awrite(NOTS)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Sclose(NOTS)"));
    // attach scalar attribute t0
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(t0)"));
    attr = H5Acreate2(s_dset, "t0", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("SourceScreen::WriteField - error in H5Acreate2(t0)"));
    double t0_src = 0.0;
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &t0_src);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Awrite(t0)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Sclose(t0)"));
    // attach scalar attribute dt
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(dt)"));
    attr = H5Acreate2(s_dset, "dt", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("SourceScreen::WriteField - error in H5Acreate2(dt)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dt_src);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Awrite(dt)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Sclose(dt)"));
    // attach scalar attribute dtx
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(dtx)"));
    attr = H5Acreate2(s_dset, "dtx", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("SourceScreen::WriteField - error in H5Acreate2(dtx)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dtx_src);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Awrite(dtx)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Sclose(dtx)"));
    // attach scalar attribute dty
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(dty)"));
    attr = H5Acreate2(s_dset, "dty", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("SourceScreen::WriteField - error in H5Acreate2(dty)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dty_src);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Awrite(dty)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Sclose(dty)"));
    // Close and release resources.
    status = H5Pclose (s_dcpl);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Pclose(s_dcpl)"));
    status = H5Dclose (s_dset);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Dclose(s_dset)"));
    status = H5Sclose (s_space);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Dclose(s_space)"));
    delete s_buffer;

    // Create dataspace for the field data.
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hsize_t dims[4];
    dims[0] = Nx;
    dims[1] = Ny;
    dims[2] = NOTS;
    dims[3] = 6;
    hid_t space = H5Screate_simple (4, dims, NULL);
    if (space<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(space)"));
    // buffer the data
    double *buffer = getBuffer(Field);
    // Create the dataset creation property list
    hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(IOexception("SourceScreen::WriteField - error in H5Pcreate(dcpl)"));
    // Create the dataset.
    hid_t dset = H5Dcreate (file,
        "ElMagField", 			// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        space, H5P_DEFAULT,
        dcpl, H5P_DEFAULT);
    if (dset<0) throw(IOexception("SourceScreen::WriteField - error in H5Dcreate(dset)"));
    // Write the data to the dataset
    status = H5Dwrite (dset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			// mem space id
        space,
        H5P_DEFAULT,			// data transfer properties
        buffer);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Dwrite(dset)"));
    // attach scalar attributes
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(t0)"));
    attr = H5Acreate2(dset, "t0", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("SourceScreen::WriteField - error in H5Acreate2(t0)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &t0_src);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Awrite(t0)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Sclose(t0)"));
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(dt)"));
    attr = H5Acreate2(dset, "dt", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("SourceScreen::WriteField - error in H5Acreate2(dt)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dt_src);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Awrite(dt)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Sclose(dt)"));
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("SourceScreen::WriteField - error in H5Screate(NOTS)"));
    attr = H5Acreate2(dset, "NOTS", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("SourceScreen::WriteField - error in H5Acreate2(NOTS)"));
    status = H5Awrite(attr, H5T_NATIVE_INT, &NOTS);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Awrite(NOTS)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Sclose(NOTS)"));
    // Close and release resources.
    status = H5Pclose (dcpl);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Pclose(dcpl)"));
    status = H5Dclose (dset);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Dclose(dset)"));
    status = H5Sclose (space);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Sclose(space)"));
    delete buffer;

    status = H5Fclose (file);
    if (status<0) throw(IOexception("SourceScreen::WriteField - error in H5Fclose()"));
    // no errors have occured if we made it 'til here
    cout << "writing HDF5 done." << endl;
    return;
}

double* SourceScreen::getBuffer(
    std::vector<std::vector<FieldTrace*>> Field
)
{
    double* buffer = new double[getBufferSize()];
    if (buffer==0)
        throw(IOexception("SourceScreen::getBuffer() - error allocating memory."));
    else
    {
        double* bp = buffer;
        for (unsigned int ix = 0; ix < Nx; ix++)
        	for (unsigned int iy = 0; iy < Ny; iy++)
        	    for (unsigned int it = 0; it < NOTS; it++)
        	    {
        	        FieldTrace *trace = Field[ix][iy];
            		ElMagField field = trace->get_field((std::size_t)it);
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


