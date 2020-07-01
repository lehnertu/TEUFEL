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

#include "logger.h"
#include "bunch.h"
#include "beam.h"
#include "global.h"

#include <iostream>
#include <math.h>
#include "SDDS.h"

template <class objectT>
TrackingLogger<objectT>::TrackingLogger(objectT *obj, const char *filename)
{
    Beam = obj;
    NOTS = 0;
    // no need to initalize the (static) fields
}

template <class objectT>
void TrackingLogger<objectT>::update()
{
    NOTS++;
    Time.push_back(Beam->getTime());
    Pos.push_back(Beam->avgPosition());
    BG.push_back(Beam->avgMomentum());
    PosRMS.push_back(Beam->rmsPosition());
    BGRMS.push_back(Vector(0.0,0.0,0.0));
    PosBG.push_back(Vector(0.0,0.0,0.0));
}

template <class objectT>
int TrackingLogger<objectT>::WriteBeamParametersSDDS(const char *filename)
{
    cout << "writing SDDS file " << filename << endl;
    SDDS_DATASET data;
    // if (1 != SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,filename))
    if (1 != SDDS_InitializeOutput(&data,SDDS_ASCII,1,NULL,NULL,filename))
    {
	cout << "WriteSDDS - error initializing output\n";
	return 1;
    }
    if  ( SDDS_DefineSimpleParameter(&data,"NumberTimeSteps","", SDDS_LONG) != 1 )
    {
	cout << "WriteSDDS - error defining parameters\n";
	return 2;
    }
    if  (
	SDDS_DefineColumn(&data,"t\0","t\0","s\0","time\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"x_av\0","x_av\0","m\0","average position\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"y_av\0","y_av\0","m\0","average position\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"z_av\0","z_av\0","m\0","average position\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"bgx_av\0","bgx_av\0","m\0","average momentum\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"bgy_av\0","bgy_av\0","m\0","average momentum\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"bgz_av\0","bgz_av\0","m\0","average momentum\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"x_rms\0","x_rms\0","m\0","r.m.s. position\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"y_rms\0","y_rms\0","m\0","r.m.s. position\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"z_rms\0","z_rms\0","m\0","r.m.s. position\0",NULL, SDDS_DOUBLE,0) == -1
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
    if  ( SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "NumberTimeSteps",NOTS, NULL ) != 1 )
    {
	cout << "WriteSDDS - error setting parameters\n";
	return 6;
    }
    // write the table of trajectory data
    // cout << "SDDS writing " << NOTS << " field values" << endl;
    for(unsigned int i=0; i<NOTS; i++)
    {
	if (SDDS_SetRowValues(&data,
	    SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,
	    "t",Time[i],
	    "x_av",Pos[i].x,
	    "y_av",Pos[i].y,
	    "z_av",Pos[i].z,
	    "bgx_av",BG[i].x,
	    "bgy_av",BG[i].y,
	    "bgz_av",BG[i].z,
	    "x_rms",PosRMS[i].x,
	    "y_rms",PosRMS[i].y,
	    "z_rms",PosRMS[i].z,
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

/*
template <class sourceT>
int PointObserver<sourceT>::WriteTimeDomainFieldSDDS(const char *filename)
{
    return 0;
}
*/

// we have to instantiate the class for every possible source type
template class TrackingLogger<Bunch>;
// template class TrackingLogger<Beam>;
