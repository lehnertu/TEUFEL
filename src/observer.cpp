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
    NPT = 0;
}

int PointObserver::GetTimeDomainTrace(
    ChargedParticle *source)
{
    ObservationTime.clear();
    ObservationField.clear();
    NPT = source->TimeDomainField(Pos, &ObservationTime, &ObservationField);
    return NPT;
}

//! \todo untested code
void PointObserver::ComputeTimeDomainField(
    ChargedParticle *source,
    double t0,
    double dt,
    int nots)
{
    t0_int = t0;
    dt_int = dt;
    NOTS = nots;
    double t_max = t0+NOTS*dt;
    // the observed field during the time step (initialized zero)
    ElMagField field;
    // make the interpolated field an empty trace of the requested length
    InterpolatedField.clear();
    for (int i=0; i<NOTS; i++)
	InterpolatedField.push_back(field);
    
    GetTimeDomainTrace(source);
    // (source_t1, source_t2) is the current segment of the (source) trace
    // source_f1, source_f2 are the corresponding field values;
    double source_t1=0.0;
    double source_t2=0.0;
    ElMagField source_f1, source_f2;
    if (NPT>0)
    {
	source_t2=ObservationTime[0];
	source_f2=ObservationField[0];
    };
    // the index of the time/field point
    // where the current segment ends, that is the index ot ts2
    // must be smaller than (not equal) NPT
    int is2=0;

    // now we process the non-equidistant time trace segment by segment
    while ((is2+1<NPT) && (source_t2<t_max))
    {
	source_t1 = source_t2;
	source_f1 = source_f2;
	is2++;
	source_t2 = ObservationTime[is2];
	source_f2 = ObservationField[is2];
	// the indices of the steps in which the current segment ends lay
	int idx1 = floor((source_t1-t0)/dt);
	int idx2 = floor((source_t2-t0)/dt);
	if (idx1==idx2)
	{
	    // if the segment is fully contained in one time step
	    field = InterpolatedField[idx1];
	    field += (source_f1+source_f2)*(0.5*(source_t2-source_t1)/dt);
	    InterpolatedField[idx1] = field;
	} else {
	    // the segment is spanning more than one time step
	    // handle the first time step
	    field = InterpolatedField[idx1];
	    double dt_i1 = t0+(idx1+1)*dt - source_t1;
	    field += source_f1*(dt_i1/dt) + (source_f2-source_f1)*(0.5*dt_i1/(source_t2-source_t1));
	    InterpolatedField[idx1] = field;
	    // handle the intermediate time steps
	    for (int idx=idx1+1; idx<idx2; idx++)
	    {
		field = InterpolatedField[idx];
		double t_center = t0 + (idx+0.5)*dt;
		field += source_f1*((source_t2-t_center)/dt) + source_f2*((t_center-source_t1)/dt);
		InterpolatedField[idx] = field;
	    };
	    // handle the last time step
	    field = InterpolatedField[idx2];
	    double dt_i2 = source_t2 - (t0+idx2*dt);
	    field += source_f2*(dt_i2/dt) + (source_f1-source_f2)*(0.5*dt_i2/(source_t2-source_t1));
	    InterpolatedField[idx2] = field;
	}
    };
}

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
    for (int i=1; i<NPT; i++)
    {
	double cost1=cos(2.0*Pi*freq*ObservationTime[i-1]);
	double cost2=cos(2.0*Pi*freq*ObservationTime[i]);
	double sint1=sin(2.0*Pi*freq*ObservationTime[i-1]);
	double sint2=sin(2.0*Pi*freq*ObservationTime[i]);
	std::complex<double> exp_iwt1(cost1, sint1);
	std::complex<double> exp_iwt2(cost2, sint2);
	double dt=ObservationTime[i]-ObservationTime[i-1];
	*Ex += exp_iwt1 * ObservationField[i-1].E().x * scale*dt;
	*Ex += exp_iwt2 * ObservationField[i].E().x * scale*dt;
	*Ey += exp_iwt1 * ObservationField[i-1].E().y * scale*dt;
	*Ey += exp_iwt2 * ObservationField[i].E().y * scale*dt;
	*Ez += exp_iwt1 * ObservationField[i-1].E().z * scale*dt;
	*Ez += exp_iwt2 * ObservationField[i].E().z * scale*dt;
    };
}

int PointObserver::WriteTimeTraceSDDS(const char *filename)
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
	SDDS_DefineColumn(&data,"t\0","t\0","s\0","TimeInSeconds\0",NULL, SDDS_DOUBLE,0)   ==-1 || 
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
    if (SDDS_StartPage(&data,(int32_t)NPT) !=1 )
    {
	cout << "WriteSDDS - error starting page\n";
	return 5;
    }
    // write the single valued variables
    cout << "SDDS write parameters" << endl;
    if  (SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
	"NumberTimeSteps",NPT,
	NULL ) !=1
    )
    {
	cout << "ChargedParticle::WriteSDDS - error setting parameters\n";
	return 6;
    }
    // write the table of trajectory data
    cout << "SDDS writing " << NPT << " field values" << endl;
    for( int i=0; i<NPT; i++)
    {
	if (SDDS_SetRowValues(&data,
	    SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,
	    "t",ObservationTime[i],
	    "Ex",ObservationField[i].E().x,
	    "Ey",ObservationField[i].E().y,
	    "Ez",ObservationField[i].E().z,
	    "Bx",ObservationField[i].B().x,
	    "By",ObservationField[i].B().y,
	    "Bz",ObservationField[i].B().z,
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

int PointObserver::WriteTimeFieldSDDS(const char *filename)
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
    if (SDDS_StartPage(&data,(int32_t)NPT) !=1 )
    {
	cout << "WriteSDDS - error starting page\n";
	return 5;
    }
    // write the single valued variables
    cout << "SDDS write parameters" << endl;
    if  ( SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "NumberTimeSteps",NOTS, NULL ) != 1 || 
	  SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "t0",t0_int, NULL ) != 1 ||
	  SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "dt",dt_int, NULL ) != 1
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
	    "Ex",InterpolatedField[i].E().x,
	    "Ey",InterpolatedField[i].E().y,
	    "Ez",InterpolatedField[i].E().z,
	    "Bx",InterpolatedField[i].B().x,
	    "By",InterpolatedField[i].B().y,
	    "Bz",InterpolatedField[i].B().z,
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
    if (SDDS_StartPage(&data,(int32_t)NPT) !=1 )
    {
	cout << "WriteSDDS - error starting page\n";
	return 5;
    }
    // write the single valued variables
    cout << "SDDS write parameters" << endl;
    if  (SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
	"NumberTimeSteps",NPT,
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
