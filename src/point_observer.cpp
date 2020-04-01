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
    unsigned int nots)
{
    FileName = filename,
    Pos = position;
    t0_obs = t0;
    dt_obs = dt;
    NOTS = nots;
    ElMagField field;
    // fill the field with zeros
    for (unsigned int i=0; i<NOTS; i++) TimeDomainField.push_back(field);
}

void PointObserver::integrate(Beam *src)
{
    src->integrateFieldTrace(
	Pos, t0_obs, dt_obs, NOTS, &TimeDomainField);
}

void PointObserver::integrate(Bunch *src)
{
    src->integrateFieldTrace(
        Pos, t0_obs, dt_obs, NOTS, &TimeDomainField);
}

void PointObserver::integrate(Lattice *src)
{
    for (unsigned int it = 0; it < NOTS; it ++ )
    {
        TimeDomainField[it] += src->Field(t0_obs+it*dt_obs, Pos);
    }
}

ElMagField PointObserver::getField(unsigned int idx)
{
    ElMagField field;
    if ( idx>=0 && idx<NOTS) field=TimeDomainField[idx];
    return field;
}

void PointObserver::setField(
        unsigned int it,
        ElMagField field)
{
    if (it<NOTS)
        TimeDomainField[it] = field;
    else
        throw(IOexception("PointObserver::setField() - requested index out of range."));
}

unsigned int PointObserver::getBufferSize()
{
    return NOTS*6;
}

double* PointObserver::getBuffer()
{
    double* buffer = new double[getBufferSize()];
    if (buffer==0)
        throw(IOexception("PointObserver::getBuffer() - error allocating memory."));
    else
    {
        double* bp = buffer;
        for (unsigned int it = 0; it < NOTS; it++)
        {
            ElMagField field = getField(it);
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

void PointObserver::fromBuffer(double *buffer, unsigned int size)
{
    if (size==NOTS*6)
    {
        ElMagField field;
        double *bp = buffer;
        for (unsigned int it = 0; it < NOTS; it++)
        {
            Vector E, B;
            E.x = *bp++;
            E.y = *bp++;
            E.z = *bp++;
            B.x = *bp++;
            B.y = *bp++;
            B.z = *bp++;
            setField(it,ElMagField(E,B));
        };
    } else {
        throw(IOexception("PointObserver::fromBuffer() - buffer size mismatch."));
    }
}

void PointObserver::WriteTimeDomainFieldSDDS()
{
    cout << "writing SDDS file " << FileName << endl;
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
	SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "t0",t0_obs, NULL ) != 1 ||
	SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "dt",dt_obs, NULL ) != 1
    )
	throw(IOexception("PointObserver::WriteTimeDomainFieldSDDS - error in SetParameters()"));
    // write the table of trajectory data
    // cout << "SDDS writing " << NOTS << " field values" << endl;
    for(unsigned int i=0; i<NOTS; i++)
    {
	if (SDDS_SetRowValues(&data,
	    SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,
	    "t",t0_obs+((double)i+0.5)*dt_obs,
			      "Ex",TimeDomainField[i].E().x,
			      "Ey",TimeDomainField[i].E().y,
			      "Ez",TimeDomainField[i].E().z,
			      "Bx",TimeDomainField[i].B().x,
			      "By",TimeDomainField[i].B().y,
			      "Bz",TimeDomainField[i].B().z,
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

