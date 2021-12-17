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
        double start_time) :
    ExternalField()
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
    attr = H5Aopen_name(dataset, "t0");
    if (attr<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute t0"));
    status =  H5Aread(attr, H5T_NATIVE_DOUBLE, &t0_src);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute t0"));
    status = H5Aclose(attr);
    if (status<0) throw(IOexception("SourceScreen - error reading dataset Screen attribute t0"));
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
    Origin = Vector(buf[0],buf[1],buf[2]);
    normal = Vector(buf[3],buf[4],buf[5]);
    dX = Vector(buf[6],buf[7],buf[8]);
    dY = Vector(buf[9],buf[10],buf[11]);
    delete buf;
    
    // done with the screen dataset
    status = H5Dclose(dataset);
    if (status<0) throw(IOexception("SourceScreen - error closing dataset Screen"));

    if (teufel::rank==0)
    {
        std::cout << "Nx=" << Nx << "  Ny=" << Ny << "  Nots=" << NOTS << std::endl;
        std::cout << "Origin = (" << Origin.x << ", " << Origin.y << ", " << Origin.z << ") m" << std::endl;
        std::cout << "Normal = (" << normal.x << ", " << normal.y << ", " << normal.z << ") m" << std::endl;
        std::cout << "dX = (" << dX.x << ", " << dX.y << ", " << dX.z << ") m" << std::endl;
        std::cout << "dY = (" << dY.x << ", " << dY.y << ", " << dY.z << ") m" << std::endl;
        std::cout << "t0=" << t0_src << "  dt=" << dt_src << "  dtx=" << dtx_src << "  dty=" << dty_src << std::endl;
        std::cout << "allocating " << (double)(3*6*Nx*Ny*NOTS) * sizeof(double) / 1e6 << " MB of memory per node" << std::endl;
    }    

    // read the field dataset
    // TODO: sanity check the field sizes
    ElMagField* field_buf = new ElMagField[getBufferSize()];
    if (field_buf==0) throw(IOexception("SourceScreen - error allocating memory."));
    dataset = H5Dopen2(file, "ElMagField", H5P_DEFAULT);
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
    // the derivatives are set to zero
    ElMagField *b = field_buf;
    Traces.resize(Nx);
    dt_Traces.resize(Nx);
    dn_Traces.resize(Nx);
    for (unsigned int ix = 0; ix < Nx; ix++) {
    	Traces[ix].resize(Ny);
    	dt_Traces[ix].resize(Ny);
    	dn_Traces[ix].resize(Ny);
    	for (unsigned int iy = 0; iy < Ny; iy++)
    	{
    	    FieldTrace* tr = new FieldTrace(CellTimeZero(ix,iy),dt_src,NOTS);
    	    if (tr==0) throw(IOexception("MeshedScreen - error allocating memory."));
    	    tr -> set_buffer(b, NOTS);
    	    b+=NOTS;
    	    Traces[ix][iy] = tr;
    	    tr = new FieldTrace(CellTimeZero(ix,iy),dt_src,NOTS);
    	    if (tr==0) throw(IOexception("MeshedScreen - error allocating memory."));
    	    dt_Traces[ix][iy] = tr;
    	    tr = new FieldTrace(CellTimeZero(ix,iy),dt_src,NOTS);
    	    if (tr==0) throw(IOexception("MeshedScreen - error allocating memory."));
    	    dn_Traces[ix][iy] = tr;
    	};
    }
    
    // done with the field data
    delete field_buf;
        
    // done with the file
    status = H5Fclose (file);
    if (status<0) throw(IOexception("MeshedScreen - error closing the file."));

    if (teufel::rank==0)
    {
        std::cout << "total energy on source screen " << totalEnergy() << " J" << std::endl;
        std::cout << std::endl;
    }    
}

SourceScreen::~SourceScreen()
{
    for (unsigned int ix = 0; ix < Nx; ix++)
    	for (unsigned int iy = 0; iy < Ny; iy++)
    	    {
                delete Traces[ix][iy];
                delete dt_Traces[ix][iy];
                delete dn_Traces[ix][iy];
            };
}

Vector SourceScreen::CellPosition(unsigned int ix, unsigned int iy)
{
    // center index is Nx/2 Ny/2 (integer division!)
    return Origin + dX*(double)((int)ix-(int)Nx/2) + dY*(double)((int)iy-(int)Ny/2);
}

double SourceScreen::CellTimeZero(unsigned int ix, unsigned int iy)
{
    // center index is Nx/2 Ny/2 (integer division!)
    return t0_src + dtx_src*(double)((int)ix-(int)Nx/2) + dty_src*(double)((int)iy-(int)Ny/2);
}

ElMagField SourceScreen::LocalField(double t, Vector X)
{
    return ElMagField(Vector(0.0,0.0,0.0),Vector(0.0,0.0,0.0));
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


