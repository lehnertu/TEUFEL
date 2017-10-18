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

#include "observer.h"
#include "bunch.h"
#include "beam.h"
#include "global.h"

#include <iostream>
#include <math.h>
#include "SDDS.h"

template <class sourceT>
PointObserver<sourceT>::PointObserver(
    sourceT *src,
    Vector position,
    double t0,
    double dt,
    int nots)
{
    Source = src;
    Pos = position;
    t0_obs = t0;
    dt_obs = dt;
    NOTS = nots;
    ElMagField field;
    // fill the field with zeros
    for (int i=0; i<NOTS; i++) TimeDomainField.push_back(field);
}

template <class sourceT>
void PointObserver<sourceT>::integrate()
{
    Source->integrateFieldTrace(
	Pos, t0_obs, dt_obs, NOTS, &TimeDomainField);
}

template <class sourceT>
ElMagField PointObserver<sourceT>::getField(int idx)
{
    ElMagField field;
    if ( idx>=0 && idx<NOTS) field=TimeDomainField[idx];
    return field;
}

template <class sourceT>
int PointObserver<sourceT>::WriteTimeDomainFieldSDDS(const char *filename)
{
    cout << "writing SDDS file " << filename << endl;
    SDDS_DATASET data;
    // if (1 != SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,filename))
    if (1 != SDDS_InitializeOutput(&data,SDDS_ASCII,1,NULL,NULL,filename))
    {
	cout << "WriteSDDS - error initializing output\n";
	return 1;
    }
    if  ( SDDS_DefineSimpleParameter(&data,"NumberTimeSteps","", SDDS_LONG) != 1 ||
	SDDS_DefineSimpleParameter(&data,"t0","", SDDS_DOUBLE) != 1 ||
	SDDS_DefineSimpleParameter(&data,"dt","", SDDS_DOUBLE) !=1
    )
    {
	cout << "WriteSDDS - error defining parameters\n";
	return 2;
    }
    if  (
	SDDS_DefineColumn(&data,"t\0","t\0","s\0","time\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"Ex\0","Ex\0","V/m\0","electric field\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"Ey\0","Ey\0","V/m\0","electric field\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"Ez\0","Ez\0","V/m\0","electric field\0",NULL, SDDS_DOUBLE,0) == -1 || 
	SDDS_DefineColumn(&data,"Bx\0","Bx\0","T\0","magnetic field\0",NULL, SDDS_DOUBLE,0)== -1 || 
	SDDS_DefineColumn(&data,"By\0","By\0","T\0","magnetic field\0",NULL,SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"Bz\0","Bz\0","T\0","magnetic field\0",NULL,SDDS_DOUBLE,0) == -1
    )
    {
	cout << "WriteSDDS - error defining data columns\n";
	return 3;
    }
    if (SDDS_WriteLayout(&data) != 1)
    {
	cout << "WriteSDDS - error writing layout\n";
	return 4;
    }
    // start a page with number of lines equal to the number of trajectory points
    // cout << "SDDS start page" << endl;
    if (SDDS_StartPage(&data,(int32_t)NOTS) !=1 )
    {
	cout << "WriteSDDS - error starting page\n";
	return 5;
    }
    // write the single valued variables
    // cout << "SDDS write parameters" << endl;
    if  ( SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "NumberTimeSteps",NOTS, NULL ) != 1 || 
	SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "t0",t0_obs, NULL ) != 1 ||
	SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "dt",dt_obs, NULL ) != 1
    )
    {
	cout << "ChargedParticle::WriteSDDS - error setting parameters\n";
	return 6;
    }
    // write the table of trajectory data
    // cout << "SDDS writing " << NOTS << " field values" << endl;
    for( int i=0; i<NOTS; i++)
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
	{
	    cout << "WriteSDDS - error writing data columns\n";
	    return 7;
	}
    }
    if( SDDS_WritePage(&data) != 1)
    {
	cout << "WriteSDDS - error writing page\n";
	return 8;
    }
    // finalize the file
    if (SDDS_Terminate(&data) !=1 )
    {
	cout << "WriteSDDS - error terminating data file\n";
	return 9;
    }	
    // no errors have occured if we made it 'til here
    // cout << "writing SDDS done." << endl;
    return 0;
}

// we have to instantiate the class for every possible source type
template class PointObserver<Bunch>;
template class PointObserver<Beam>;