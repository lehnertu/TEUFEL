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

#include "beam.h"
#include "particle.h"
#include "global.h"
#include <random>
#include <iostream>
#include "SDDS.h"
#include "hdf5.h"


Distribution::Distribution(int dim, int nop)
{
    DIM = dim;
    NOP = nop;
    A = new double[DIM*NOP];
    // generate a random seed
    std::random_device rd;
    // pseudo-random generator of 32-bit numbers with a state size of 19937 bits
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i=0; i<DIM*NOP; i++)
	A[i] = dist(mt);
}

Distribution::~Distribution()
{
    delete[] A;
}

int Distribution::getDIM()
{
    return DIM;
}

int Distribution::getNOP()
{
    return NOP;
}
void Distribution::generateGaussian(double mean, double sigma, int dim)
{
    if (dim>=0 && dim<DIM)
    {
	// generate a random seed
	std::random_device rd;
	// pseudo-random generator of 32-bit numbers with a state size of 19937 bits
	std::mt19937 mt(rd());
	std::normal_distribution<double> dist(mean, sigma);
	for (int i=0; i<NOP; i++)
	    A[i*DIM+dim] = dist(mt);
    }
}

void Distribution::addCorrelation(int independent, int dependent, double factor)
{
    if (independent>=0 && independent<DIM &&
	dependent>=0 && dependent<DIM && dependent!=independent)
    {
	for (int i=0; i<NOP; i++)
	    A[i*DIM+dependent] += factor * A[i*DIM+independent];
    }
}

double Distribution::getCoordinate(int index, int dim)
{
    if (dim>=0 && dim<DIM && index>=0 && index<NOP)
	return A[index*DIM+dim];
    else
	return 0.0;
}

BunchedBeam::BunchedBeam()
{
    NOP = 0;
}

BunchedBeam::BunchedBeam(int N, double charge, double mass)
{
    NOP = N;
    for(int i=0; i<NOP; i++)
    {
	P.push_back(new ChargedParticle(charge,mass));
    }
}

BunchedBeam::BunchedBeam(BunchedBeam* origin)
{
    NOP = origin->getNOP();
    for(int i=0; i<NOP; i++)
    {
	P.push_back(new ChargedParticle(origin->getParticle(i)));
    }
}

BunchedBeam::~BunchedBeam()
{
}

void BunchedBeam::Add(ChargedParticle *part)
{
    NOP++;
    P.push_back(part);
}

int BunchedBeam::getNOP()
{
    return NOP;
}

double BunchedBeam::getTotalCharge()
{
    double charge = 0.0;
    for(int i=0; i<NOP; i++)
    {
	charge += getParticle(i)->getCharge();
    }
    return charge;
}

ChargedParticle* BunchedBeam::getParticle(int i)
{
    return P[i];
}

void BunchedBeam::InitVay(Distribution *dist,
		    double tstep,
		    GeneralField* field)
{
    dt = tstep;
    for(int i=0; i<NOP; i++)
    {
	Vector X0 = Vector(dist->getCoordinate(i,0),
			   dist->getCoordinate(i,1),
			   dist->getCoordinate(i,2));
	Vector P0 = Vector(dist->getCoordinate(i,3),
			   dist->getCoordinate(i,4),
			   dist->getCoordinate(i,5));
	double t0 = dist->getCoordinate(i,6);
	P[i]->InitVay(t0, X0, P0, tstep, field);
    }
}

void BunchedBeam::StepVay(GeneralField* field)
{
    for(int i=0; i<NOP; i++)
    {
	P[i]->StepVay(field);
    }
}

int BunchedBeam::WriteWatchPointSDDS(double time,
			       const char *filename)
{
    cout << "writing SDDS file " << filename << endl;
    SDDS_DATASET data;
    if (1 != SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,filename))
    {
	cout << "BunchedBeam::WriteWatchPointSDDS - error initializing output\n";
	return 1;
    }
    if  (SDDS_DefineSimpleParameter(&data,"NumberOfParticles","", SDDS_LONG)!=1)
    {
	cout << "BunchedBeam::WriteWatchPointSDDS - error defining parameters\n";
	return 2;
    }
    if  (
	SDDS_DefineColumn(&data,"t\0","t\0","s\0","Time\0",NULL, SDDS_DOUBLE,0)   ==-1 || 
	SDDS_DefineColumn(&data,"x\0","x\0","m\0","PositionX\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"y\0","y\0","m\0","PositionY\0",NULL, SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"z\0","z\0","m\0","PositionZ\0",NULL, SDDS_DOUBLE,0) == -1 || 
	SDDS_DefineColumn(&data,"px\0","px\0",NULL,"BetaGammaX\0",NULL, SDDS_DOUBLE,0)== -1 || 
	SDDS_DefineColumn(&data,"py\0","py\0",NULL,"BetaGammaY\0",NULL,SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"pz\0","pz\0",NULL,"BetaGammaZ\0",NULL,SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"p\0","p\0",NULL,"BetaGamma\0",NULL,SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"xp\0","xp\0",NULL,"Xprime\0",NULL, SDDS_DOUBLE,0)== -1 || 
	SDDS_DefineColumn(&data,"yp\0","yp\0",NULL,"Bprime\0",NULL,SDDS_DOUBLE,0) == -1 ||
	SDDS_DefineColumn(&data,"gamma\0","gamma\0",NULL,"RelativisticFactor\0",NULL,SDDS_DOUBLE,0)==-1
	)
    {
	cout << "BunchedBeam::WriteWatchPointSDDS - error defining data columns\n";
	return 3;
    }
    if (SDDS_WriteLayout(&data) != 1)
    {
	cout << "BunchedBeam::WriteWatchPointSDDS - error writing layout\n";
	return 4;
    }
    // start a page with number of lines equal to the number of trajectory points
    cout << "SDDS start page" << endl;
    if (SDDS_StartPage(&data,(int32_t)NOP) !=1 )
    {
	cout << "BunchedBeam::WriteWatchPointSDDS - error starting page\n";
	return 5;
    }
    // write the single valued variables
    cout << "SDDS write parameters" << endl;
    if( SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
	"NumberOfParticles",NOP,
	NULL ) !=1
	)
    {
	cout << "BunchedBeam::WriteWatchPointSDDS - error setting parameters\n";
	return 6;
    }
    // write the table of particle data
    cout << "SDDS writing " << NOP << " particles" << endl;
    for( int i=0; i<NOP; i++)
    {
	Vector X, BG;
	// query the particle for its coordinates at the given time
	P[i]->CoordinatesAtTime(time, &X, &BG);
	if( SDDS_SetRowValues(&data,
	    SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,
	    "t",time,
	    "x",X.x,
	    "y",X.y,
	    "z",X.z,
	    "px",BG.x,
	    "py",BG.y,
	    "pz",BG.z,
	    "p",BG.norm(),
	    "xp",BG.x/BG.z,
	    "yp",BG.y/BG.z,
	    "gamma",sqrt(1.0+BG.abs2nd()),
	    NULL) != 1
	    )
	{
	    cout << "BunchedBeam::WriteWatchPointSDDS - error writing data columns\n";
	    return 7;
	}
    }
    if( SDDS_WritePage(&data) != 1)
    {
	cout << "BunchedBeam::WriteWatchPointSDDS - error writing page\n";
	return 8;
    }
    // finalize the file
    if (SDDS_Terminate(&data) !=1 )
    {
	cout << "BunchedBeam::WriteWatchPointSDDS - error terminating data file\n";
	return 9;
    }	
    // no errors have occured if we made it 'til here
    cout << "writing SDDS done." << endl;
    return 0;
}

int BunchedBeam::WriteWatchPointHDF5(double time,
			const char *filename)
{
    herr_t status;
    cout << "writing HDF5 file " << filename << endl;
    // Create a new file using the default properties.
    hid_t file = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0 )
    {
	cout << "BunchedBeam::WriteWatchPointHDF5 - error crating file\n";
	return 1;
    };
    // Create dataspace. Setting maximum size to NULL sets the maximum
    // size to be the current size.
    hsize_t dims[2];
    dims[0] = NOP;
    dims[1] = 6;
    hid_t space = H5Screate_simple (2, dims, NULL);
    if (space < 0 )
    {
	cout << "BunchedBeam::WriteWatchPointHDF5 - error crating dataspace\n";
	return 2;
    };
    // buffer the data
    double *buffer = new double[NOP*6];
    double *bp = buffer;
    for( int i=0; i<NOP; i++)
    {
	Vector X, BG;
	// query the particle for its coordinates at the given time
	P[i]->CoordinatesAtTime(time, &X, &BG);
	*bp++ = X.x;
	*bp++ = X.y;
	*bp++ = X.z;
	*bp++ = BG.x;
	*bp++ = BG.y;
	*bp++ = BG.z;
    };
    // Create the dataset creation property list
    hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl < 0 )
    {
	cout << "BunchedBeam::WriteWatchPointHDF5 - error crating property list\n";
	return 3;
    };
    // Create the dataset.
    hid_t dset = H5Dcreate (file,
	"electrons", 			// dataset name
	H5T_NATIVE_DOUBLE,		// data type
	space, H5P_DEFAULT,
	dcpl, H5P_DEFAULT);
    if (dset < 0 )
    {
	cout << "BunchedBeam::WriteWatchPointHDF5 - error crating dataset\n";
	return 4;
    };
    // Write the data to the dataset
    status = H5Dwrite (dset,
	H5T_NATIVE_DOUBLE, 		// mem type id
	H5S_ALL, 			// mem space id
	space,
	H5P_DEFAULT,			// data transfer properties
	buffer);
    if (status < 0 )
    {
	cout << "BunchedBeam::WriteWatchPointHDF5 - error writing dataset\n";
	return 5;
    }	
    // Close and release resources.
    status = H5Pclose (dcpl);
    if (status < 0 )
    {
	cout << "BunchedBeam::WriteWatchPointHDF5 - error releasing property list\n";
	return 6;
    }	
    status = H5Dclose (dset);
    if (status < 0 )
    {
	cout << "BunchedBeam::WriteWatchPointHDF5 - error releasing dataset\n";
	return 7;
    }	
    status = H5Sclose (space);
    if (status < 0 )
    {
	cout << "BunchedBeam::WriteWatchPointHDF5 - error releasing dataspace\n";
	return 8;
    }	
    status = H5Fclose (file);
    if (status < 0 )
    {
	cout << "BunchedBeam::WriteWatchPointHDF5 - error closing file\n";
	return 9;
    }	
    // no errors have occured if we made it 'til here
    cout << "writing HDF5 done." << endl;
    return 0;
}