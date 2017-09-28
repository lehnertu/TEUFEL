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
#include "particle.h"
#include "global.h"

#include <iostream>
#include <math.h>
#include "SDDS.h"

PointObserver::PointObserver(Vector position)
{
    Pos = position;
}

template <class Src>
void PointObserver::ComputeTimeDomainField(
    Src *source,
    double t0,
    double dt,
    int nots)
{
    t0_obs = t0;
    dt_obs = dt;
    NOTS = nots;
    source->getTimeDomainField(
	Pos, t0_obs, dt_obs, NOTS, &TimeDomainField);
}

// we have to specialize the templated function for all possibles types
template <> void PointObserver::ComputeTimeDomainField<>(
    ChargedParticle *source, double t0, double dt, int nots);
template <> void PointObserver::ComputeTimeDomainField<>(
    Bunch *source, double t0, double dt, int nots);
/* not yet implemented
template void PointObserver::ComputeTimeDomainField<>(Beam *source, double t0, double dt, int nots);
*/

void PointObserver::FrequencyObservation(
    double freq,
    std::complex<double> *Ex,
    std::complex<double> *Ey,
    std::complex<double> *Ez)
{
    *Ex = {0.0, 0.0};
    *Ey = {0.0, 0.0};
    *Ez = {0.0, 0.0};
    double scale=0.5/sqrt(2.0*Pi);
    for (int i=1; i<NOTS; i++)
    {
	double t2 = t0_obs + i*dt_obs;
	double t1 = t2 - dt_obs;
	double cost1=cos(2.0*Pi*freq*t1);
	double cost2=cos(2.0*Pi*freq*t2);
	double sint1=sin(2.0*Pi*freq*t1);
	double sint2=sin(2.0*Pi*freq*t2);
	std::complex<double> exp_iwt1(cost1, sint1);
	std::complex<double> exp_iwt2(cost2, sint2);
	*Ex += exp_iwt1 * TimeDomainField[i-1].E().x * scale*dt_obs;
	*Ex += exp_iwt2 * TimeDomainField[i].E().x * scale*dt_obs;
	*Ey += exp_iwt1 * TimeDomainField[i-1].E().y * scale*dt_obs;
	*Ey += exp_iwt2 * TimeDomainField[i].E().y * scale*dt_obs;
	*Ez += exp_iwt1 * TimeDomainField[i-1].E().z * scale*dt_obs;
	*Ez += exp_iwt2 * TimeDomainField[i].E().z * scale*dt_obs;
    };
}

int PointObserver::WriteTimeDomainFieldSDDS(const char *filename)
{
    cout << "writing SDDS file " << filename << endl;
    SDDS_DATASET data;
    if (1 != SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,filename))
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
    cout << "SDDS start page" << endl;
    if (SDDS_StartPage(&data,(int32_t)NOTS) !=1 )
    {
	cout << "WriteSDDS - error starting page\n";
	return 5;
    }
    // write the single valued variables
    cout << "SDDS write parameters" << endl;
    if  ( SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "NumberTimeSteps",NOTS, NULL ) != 1 || 
	  SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "t0",t0_obs, NULL ) != 1 ||
	  SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "dt",dt_obs, NULL ) != 1
	)
    {
	cout << "ChargedParticle::WriteSDDS - error setting parameters\n";
	return 6;
    }
    // write the table of trajectory data
    cout << "SDDS writing " << NOTS << " field values" << endl;
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
    cout << "writing SDDS done." << endl;
    return 0;
}

int PointObserver::WriteSpectrumSDDS(
    const char *filename,
    double fstart,
    double fstop,
    double fstep)
{
    cout << "writing SDDS file " << filename << endl;
    SDDS_DATASET data;
    if (1 != SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,filename))
    {
	cout << "WriteSDDS - error initializing output\n";
	return 1;
    }
    if  (SDDS_DefineSimpleParameter(&data,"NumberTimeSteps","", SDDS_LONG)!=1)
    {
	cout << "WriteSDDS - error defining parameters\n";
	return 2;
    }
    if  (
	SDDS_DefineColumn(&data,"f\0","f\0","Hz\0","frequency\0",NULL, SDDS_DOUBLE,0)   ==-1 || 
	SDDS_DefineColumn(&data,"Ax\0","Ax\0","","amplitude\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"Px\0","Px\0","rad\0","phase\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"Ay\0","Ay\0","","amplitude\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"Py\0","Py\0","rad\0","phase\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"Az\0","Az\0","","amplitude\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"Pz\0","Pz\0","rad\0","phase\0",NULL, SDDS_DOUBLE,0) == -1
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
    cout << "SDDS start page" << endl;
    if (SDDS_StartPage(&data,(int32_t)NOTS) !=1 )
    {
	cout << "WriteSDDS - error starting page\n";
	return 5;
    }
    // write the single valued variables
    cout << "SDDS write parameters" << endl;
    if  (SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
	"NumberTimeSteps",NOTS,
	NULL ) !=1
    )
    {
	cout << "ChargedParticle::WriteSDDS - error setting parameters\n";
	return 6;
    }
    // compute and write the table of spectral data
    std::complex<double> Ex;
    std::complex<double> Ey;
    std::complex<double> Ez;
    int err=0;
    int i=0;	// row index in the table
    for (double f=fstart; f<=fstop; f+=fstep)
    {
	if (err==0)
	{
	    FrequencyObservation(f, &Ex, &Ey, &Ez);
	    if (SDDS_SetRowValues(&data,
		SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,
		"f",f,
		"Ax",std::abs(Ex),
		"Px",std::arg(Ex),
		"Ay",std::abs(Ey),
		"Py",std::arg(Ey),
		"Az",std::abs(Ez),
		"Pz",std::arg(Ez),
		NULL) != 1
		) err ++;
	}
	i++;
    }
    if (err==0)
    {
	cout << "SDDS written " << i << " spectral frequencies" << endl;
    } else {
	cout << "WriteSDDS - error writing data columns\n";
	return 7;
    }
    // finalize the file
    if( SDDS_WritePage(&data) != 1)
    {
	cout << "WriteSDDS - error writing page\n";
	return 8;
    }
    if (SDDS_Terminate(&data) !=1 )
    {
	cout << "WriteSDDS - error terminating data file\n";
	return 9;
    }	
    // no errors have occured if we made it 'til here
    cout << "writing SDDS done." << endl;
    return 0;
}
