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
#include "point_observer.h"

#include <iostream>
#include <math.h>
#include "SDDS.h"
#include "hdf5.h"

PointObserver::PointObserver(
    std::string filename,
    Vector position,
    double t0,
    double dt,
    std::size_t nots)
{
    FileName = filename,
    Pos = position;
    trace = new FieldTrace(t0,dt,nots);
    // set a default source
    source = BeamObservation;
}

PointObserver::~PointObserver()
{
    delete trace;
}

void PointObserver::setSource(RadSource s)
{
    source = s;
}

RadSource PointObserver::getSource()
{ 
    return source;
}

void PointObserver::integrate(Beam *src)
{
    src->integrateFieldTrace(Pos, trace);
}

void PointObserver::integrate(Bunch *src)
{
    src->integrateFieldTrace(Pos, trace);
}

void PointObserver::integrate(Lattice *src)
{
    for (std::size_t it=0; it<trace->get_N(); it++)
        trace->add(it,src->Field(trace->get_time(it),Pos));
}

double* PointObserver::getBuffer()
{
    // buffer size of the trace is given in ElMagField units
    std::size_t Nb = trace->get_N();
    double* buffer = new double[6*Nb];
    if (buffer==0)
        throw(IOexception("PointObserver::getBuffer() - error allocating memory."));
    else
        trace->get_buffer((ElMagField *)buffer, Nb);
    return buffer;
}

void PointObserver::fromBuffer(double *buffer, std::size_t size)
{
    // buffer size of the trace is given in ElMagField units
    std::size_t Nb = trace->get_N();
    if (size==6*Nb)
    {
        ElMagField field;
        double *bp = buffer;
        for (std::size_t it=0; it<Nb; it++)
        {
            Vector E, B;
            E.x = *bp++;
            E.y = *bp++;
            E.z = *bp++;
            B.x = *bp++;
            B.y = *bp++;
            B.z = *bp++;
            trace->set_field(it,ElMagField(E,B));
        };
    } else {
        throw(IOexception("PointObserver::fromBuffer() - buffer size mismatch."));
    }
}

void PointObserver::WriteTimeDomainFieldSDDS()
{
    cout << "writing SDDS file " << FileName << endl;
    unsigned int NOTS = trace->get_N();
    SDDS_DATASET data;
    // if (1 != SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,filename))
    if (1 != SDDS_InitializeOutput(&data,SDDS_ASCII,1,NULL,NULL,FileName.c_str()))
    	throw(IOexception("PointObserver::WriteTimeDomainFieldSDDS - error in InitializeOutput()"));
    if  ( SDDS_DefineSimpleParameter(&data,"NumberTimeSteps","", SDDS_LONG) != 1 ||
    	SDDS_DefineSimpleParameter(&data,"t0","", SDDS_DOUBLE) != 1 ||
	    SDDS_DefineSimpleParameter(&data,"dt","", SDDS_DOUBLE) !=1
        )
        	throw(IOexception("PointObserver::WriteTimeDomainFieldSDDS - error in DefineSimpleParameter()"));
    if  (
    	SDDS_DefineColumn(&data,"t\0","t\0","s\0","time\0",NULL, SDDS_DOUBLE,0) == -1 ||
	    SDDS_DefineColumn(&data,"Ex\0","Ex\0","V/m\0","electric field\0",NULL, SDDS_DOUBLE,0) == -1 ||
	    SDDS_DefineColumn(&data,"Ey\0","Ey\0","V/m\0","electric field\0",NULL, SDDS_DOUBLE,0) == -1 ||
	    SDDS_DefineColumn(&data,"Ez\0","Ez\0","V/m\0","electric field\0",NULL, SDDS_DOUBLE,0) == -1 || 
	    SDDS_DefineColumn(&data,"Bx\0","Bx\0","T\0","magnetic field\0",NULL, SDDS_DOUBLE,0)== -1 || 
	    SDDS_DefineColumn(&data,"By\0","By\0","T\0","magnetic field\0",NULL,SDDS_DOUBLE,0) == -1 ||
	    SDDS_DefineColumn(&data,"Bz\0","Bz\0","T\0","magnetic field\0",NULL,SDDS_DOUBLE,0) == -1
        )
	        throw(IOexception("PointObserver::WriteTimeDomainFieldSDDS - error in DefineColumn()"));
    if (SDDS_WriteLayout(&data) != 1)
	    throw(IOexception("PointObserver::WriteTimeDomainFieldSDDS - error in WriteLayout()"));
    // start a page with number of lines equal to the number of trajectory points
    // cout << "SDDS start page" << endl;
    if (SDDS_StartPage(&data,(int32_t)NOTS) !=1 )
	    throw(IOexception("PointObserver::WriteTimeDomainFieldSDDS - error in StartPage()"));
    // write the single valued variables
    // cout << "SDDS write parameters" << endl;
    if  ( SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "NumberTimeSteps",NOTS, NULL ) != 1 || 
	    SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "t0",trace->get_t0(), NULL ) != 1 ||
	    SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "dt",trace->get_dt(), NULL ) != 1
        )
	        throw(IOexception("PointObserver::WriteTimeDomainFieldSDDS - error in SetParameters()"));
    // write the table of trajectory data
    // cout << "SDDS writing " << NOTS << " field values" << endl;
    for(unsigned int i=0; i<NOTS; i++)
    {
	    if (SDDS_SetRowValues(&data,
	        SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,
                "t",trace->get_t0()+(double)i*trace->get_dt(),
                "Ex",trace->get_field((std::size_t)i).E().x,
                "Ey",trace->get_field((std::size_t)i).E().y,
                "Ez",trace->get_field((std::size_t)i).E().z,
                "Bx",trace->get_field((std::size_t)i).B().x,
                "By",trace->get_field((std::size_t)i).B().y,
                "Bz",trace->get_field((std::size_t)i).B().z,
                NULL) != 1
	        )
            throw(IOexception("PointObserver::WriteTimeDomainFieldSDDS - error in SetRowValues()"));
    }
    if( SDDS_WritePage(&data) != 1)
    	throw(IOexception("PointObserver::WriteTimeDomainFieldSDDS - error in WritePage()"));
    // finalize the file
    if (SDDS_Terminate(&data) !=1 )
    	throw(IOexception("PointObserver::WriteTimeDomainFieldSDDS - error in Terminate()"));
    // no errors have occured if we made it 'til here
    // cout << "writing SDDS done." << endl;
    return;
}

void PointObserver::generateOutput()
{
    WriteTimeDomainFieldSDDS();
}

