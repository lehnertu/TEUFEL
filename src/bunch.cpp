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
#include "particle.h"
#include "global.h"
#include <math.h>
#include <random>
#include <iostream>
#include <cstring>
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
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
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

void Distribution::scale( int dim, double factor)
{
    if (dim>=0 && dim<DIM)
    {
        for (int i=0; i<NOP; i++)
            A[i*DIM+dim] *= factor;
    }
}

void Distribution::generateGaussian(int dim, double sigma)
{
    if (dim>=0 && dim<DIM)
    {
        // generate a random seed
        std::random_device rd;
        // pseudo-random generator of 32-bit numbers with a state size of 19937 bits
        std::mt19937 mt(rd());
        std::normal_distribution<double> dist(0.0, sigma);
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

void Distribution::setCoordinate(int index, int dim, double value)
{
    if (dim>=0 && dim<DIM && index>=0 && index<NOP)
        A[index*DIM+dim] = value;
}

void Distribution::bufferData(double *buffer, int bufsize)
{
    if (bufsize == DIM*NOP)
    {
        std::memcpy(buffer,A,bufsize*sizeof(double));
    }
    else
    {
        throw(IOexception("Distribution::bufferData() - buffer size mismatch"));
    }
}

void Distribution::fromBuffer(double *buffer, int bufsize)
{
    if (bufsize == DIM*NOP)
    {
        std::memcpy(A,buffer,bufsize*sizeof(double));
    }
    else
    {
        throw(IOexception("Distribution::fromBuffer() - buffer size mismatch"));
    }
}

Distribution* Distribution::subDist(int index, int number)
{
    Distribution* sub = new Distribution(DIM, number);
    for (int i=0; i<number; i++)
        for (int k=0; k<DIM; k++)
        {
            sub->setCoordinate(i,k,getCoordinate(index+i,k));
        }
    return sub;
}

Bunch::Bunch()
{
    NOP = 0;
}

Bunch::Bunch(int N, double charge, double mass)
{
    NOP = N;
    for(int i=0; i<NOP; i++)
    {
        P.push_back(new ChargedParticle(charge,mass));
    }
}

Bunch::Bunch(Bunch* origin)
{
    NOP = origin->getNOP();
    for(int i=0; i<NOP; i++)
    {
        P.push_back(new ChargedParticle(origin->getParticle(i)));
    }
}

Bunch::Bunch(Distribution *dist, double reftime, Vector refpos, Vector refmom, double charge, double mass)
{
    NOP = dist->getNOP();
    for(int i=0; i<NOP; i++)
    {
        ChargedParticle *p = new ChargedParticle(charge,mass);
        double t0 = reftime + dist->getCoordinate(i,6);
        Vector X0 = refpos +
            Vector(dist->getCoordinate(i,0),
                   dist->getCoordinate(i,1),
                   dist->getCoordinate(i,2));
        Vector P0 = refmom +
            Vector(dist->getCoordinate(i,3),
                   dist->getCoordinate(i,4),
                   dist->getCoordinate(i,5));
        Vector A0 = Vector(0,0,0);
        p->initTrajectory(t0, X0, P0, A0);
        P.push_back(p);
    }
}

Bunch::Bunch(const char *filename)
{
    // set a safe default in case we return with read errors
    NOP = 0;
    std::cout << "reading particles from SDDS file " << filename << std::endl;
    // make a copy of the filename to make it changable
    char fn[100];
    strncpy (fn, filename, sizeof(fn)-1);
    SDDS_DATASET SDDS_dataset;
    // SDDS_InitializeInput - returns 1 on success, 0 on failure 
    if (SDDS_InitializeInput(&SDDS_dataset, fn ) != 1)
    {
        std::cout << "cannot read SDDS file " << filename;
        std::cout << " : creating empty bunch." << std::endl;
        return;
    };
    
    int page = SDDS_ReadPage(&SDDS_dataset);
    while (page >= 1)
    {
        page = SDDS_ReadPage(&SDDS_dataset);
    };
    
    char parName[] = "Particles";
    int32_t *nopData = SDDS_GetParameterAsLong(&SDDS_dataset, parName, NULL);
    if (nopData==NULL)
    {
        std::cout << "cannot read number of particles" << std::endl;
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        return;
    };

    // the file value is the total charge [C] unsigned value
    char chName[] = "Charge";
    double *chData = SDDS_GetParameterAsDouble(&SDDS_dataset, chName, NULL);
    if (chData==NULL)
    {
        std::cout << "cannot read charge" << std::endl;
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        return;
    };
    double charge = *chData/ElementaryCharge/(*nopData);

    char xName[] = "x";
    double *xData = SDDS_GetColumnInDoubles(&SDDS_dataset, xName);
    if (xData==NULL)
    {
        std::cout << "cannot read x column" << std::endl;
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        return;
    };

    char xpName[] = "xp";
    double *xpData = SDDS_GetColumnInDoubles(&SDDS_dataset, xpName);
    if (xpData==NULL)
    {
        std::cout << "cannot read xp column" << std::endl;
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        return;
    };

    char yName[] = "y";
    double *yData = SDDS_GetColumnInDoubles(&SDDS_dataset, yName);
    if (yData==NULL)
    {
        std::cout << "cannot read y column" << std::endl;
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        return;
    };

    char ypName[] = "yp";
    double *ypData = SDDS_GetColumnInDoubles(&SDDS_dataset, ypName);
    if (ypData==NULL)
    {
        std::cout << "cannot read yp column" << std::endl;
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        return;
    };

    char pName[] = "p";
    double *pData = SDDS_GetColumnInDoubles(&SDDS_dataset, pName);
    if (pData==NULL)
    {
        std::cout << "cannot read p column" << std::endl;
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        return;
    };

    char tName[] = "t";
    double *tData = SDDS_GetColumnInDoubles(&SDDS_dataset, tName);
    if (tData==NULL)
    {
        std::cout << "cannot read t column" << std::endl;
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        return;
    };

    NOP = *nopData;
    // t is the arrival time at a given z-plane
    // compute the mean arrival time of the bunch
    double sum = 0;
    for(int i=0; i<NOP; i++) sum += tData[i];
    double t0 = sum/NOP;
    // create the particles
    for(int i=0; i<NOP; i++)
    {
        // create negatively charged particles with electron charge-to-mass ratio
        ChargedParticle *p = new ChargedParticle(-charge,charge);
        // pData is the total momentum
        double betagamma = pData[i];
        double beta = sqrt( betagamma*betagamma / (betagamma*betagamma+1.0) );
        // xp and yp are the angles [rad] with respect to the axis (z)
        Vector direction = Vector(xpData[i], ypData[i], 1.0);
        direction.normalize();
        // relativistic velocity
        Vector v = direction * beta;
        // we translate the difference in arrival time into a position at t0
        Vector X0 = Vector(xData[i], yData[i], 0.0) + v*(t0-tData[i])*SpeedOfLight;
        Vector P0 = direction * betagamma;
        Vector A0 = Vector(0.0, 0.0, 0.0);
        p->initTrajectory(t0, X0, P0, A0);
        P.push_back(p);
    }
    
}

Bunch::~Bunch()
{
    for(int i=0; i<NOP; i++)
        delete P[i];
}

void Bunch::Add(ChargedParticle *part)
{
    NOP++;
    P.push_back(part);
}

void Bunch::clearTrajectories()
{
    for(int i=0; i<NOP; i++)
        P[i]->clearTrajectory();
}

void Bunch::preAllocate(int nTraj)
{
    for(int i=0; i<NOP; i++)
        P[i]->preAllocate(nTraj);
}

int Bunch::getNOP()
{
    return NOP;
}

double Bunch::getTotalCharge()
{
    double charge = 0.0;
    for(int i=0; i<NOP; i++)
    {
        charge += getParticle(i)->getCharge();
    }
    return charge;
}

ChargedParticle* Bunch::getParticle(int i)
{
    ChargedParticle *p = 0;
    if (i>=0 && i<NOP) p=P[i];
    return p;
}

void Bunch::replaceParticle(int i, ChargedParticle* part)
{
    if (i>=0 && i<NOP)
    {
        delete P[i];
        P[i] = part;
    }
}
    
void Bunch::InitVay(double tstep,
                    GeneralField *field)
{
    dt = tstep;
    for(int i=0; i<NOP; i++)
    {
        P[i]->InitVay(dt, field);
    }
}

void Bunch::StepVay(GeneralField *field)
{
    for(int i=0; i<NOP; i++)
    {
        P[i]->StepVay(field);
    }
}

double *Bunch::bufferStep(double *buffer)
{
    double *bp = buffer;
    int nop = NOP;
    for(int i=0; i<nop; i++)
    {
        double time = P[i]->getTime();
        Vector X = P[i]->getPosition();
        Vector BG = P[i]->getMomentum();
        Vector A = P[i]->getAccel();
        *bp++ = time;
        *bp++ = X.x;
        *bp++ = X.y;
        *bp++ = X.z;
        *bp++ = BG.x;
        *bp++ = BG.y;
        *bp++ = BG.z;
        *bp++ = A.x;
        *bp++ = A.y;
        *bp++ = A.z;
    }
    return bp;
}

double *Bunch::setStepFromBuffer(double *buffer)
{
    double *bp = buffer;
    int nop = NOP;
    double time;
    Vector X, BG, A;
    for(int i=0; i<nop; i++)
    {
        // extract values from buffer
        time = *bp++;
        X.x = *bp++;
        X.y = *bp++;
        X.z = *bp++;
        BG.x = *bp++;
        BG.y = *bp++;
        BG.z = *bp++;
        A.x = *bp++;
        A.y = *bp++;
        A.z = *bp++;
        // update the particle
        P[i]->setStep(time, X, BG, A);
    }
    return bp;
}

double Bunch::avgTime()
{
    double t = 0.0;
    for(int i=0; i<NOP; i++) t += P[i]->getTime();
    return t/NOP;
}

Vector Bunch::avgPosition()
{
    Vector pos = Vector(0.0, 0.0, 0.0);
    for(int i=0; i<NOP; i++) pos += P[i]->getPosition();
    return pos*(1.0/NOP);
}

Vector Bunch::rmsPosition()
{
    Vector sum = Vector(0.0, 0.0, 0.0);
    Vector mean = avgPosition();
    for(int i=0; i<NOP; i++)
    {
        Vector dev = P[i]->getPosition() - mean;
        sum += dev.square();
    };
    sum /= NOP;
    return sum.root();
}

Vector Bunch::avgMomentum()
{
    Vector mom = Vector(0.0, 0.0, 0.0);
    for(int i=0; i<NOP; i++)
        mom += P[i]->getMomentum();
    return mom*(1.0/NOP);
}

double* Bunch::bufferCoordinates(double *buffer, int size)
{
    double *bp = buffer;
    int nop = NOP;
    if (nop>size)
    {
        nop = size;
        std::cout << "warning in Bunch::bufferCoordinates() - buffer too small" << std::endl;
    }
    for(int i=0; i<nop; i++)
    {
        Vector X = P[i]->getPosition();
        Vector BG = P[i]->getMomentum();
        *bp++ = X.x;
        *bp++ = X.y;
        *bp++ = X.z;
        *bp++ = BG.x;
        *bp++ = BG.y;
        *bp++ = BG.z;
    }
    return bp;
}

int Bunch::WriteWatchPointSDDS(const char *filename)
{
    cout << "writing SDDS file " << filename << endl;
    SDDS_DATASET data;
    if (1 != SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,filename))
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error initializing output" << std::endl;
        return 1;
    }
    if  (SDDS_DefineSimpleParameter(&data,"NumberOfParticles","", SDDS_LONG)!=1)
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error defining parameters" << std::endl;
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
        std::cout << "Bunch::WriteWatchPointSDDS - error defining data columns" << std::endl;
        return 3;
    }
    if (SDDS_WriteLayout(&data) != 1)
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error writing layout" << std::endl;
        return 4;
    }
    // start a page with number of lines equal to the number of trajectory points
    // cout << "SDDS start page" << endl;
    if (SDDS_StartPage(&data,(int32_t)NOP) !=1 )
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error starting page" << std::endl;
        return 5;
    }
    // write the single valued variables
    // cout << "SDDS write parameters" << endl;
    if( SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
        "NumberOfParticles",NOP,
        NULL ) !=1
        )
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error setting parameters" << std::endl;
        return 6;
    }
    // write the table of particle data
    // cout << "SDDS writing " << NOP << " particles" << endl;
    for( int i=0; i<NOP; i++)
    {
        double t = P[i]->getTime();
        Vector X = P[i]->getPosition();
        Vector BG = P[i]->getMomentum();
        if( SDDS_SetRowValues(&data,
            SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,
            "t",t,
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
            std::cout << "Bunch::WriteWatchPointSDDS - error writing data columns" << std::endl;
            return 7;
        }
    }
    if( SDDS_WritePage(&data) != 1)
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error writing page" << std::endl;
        return 8;
    }
    // finalize the file
    if (SDDS_Terminate(&data) !=1 )
    {
        std::cout << "Bunch::WriteWatchPointSDDS - error terminating data file" << std::endl;
        return 9;
    }	
    // no errors have occured if we made it 'til here
    // cout << "writing SDDS done." << endl;
    return 0;
}

ElMagField Bunch::RetardedField(double obs_time, Vector obs_pos)
{
    // initialize to zero
    ElMagField field;
    // integrate for all particles
    for (int i=0; i<NOP; i++)
        field += P[i]->RetardedField(obs_time, obs_pos);
    return field;
}

void Bunch::integrateFieldTrace(
        Vector ObservationPoint,
        double t0,
        double dt,
        int nots,
        std::vector<ElMagField> *ObservationField)
{
    // cout << "Bunch::integrateFieldTrace() NOP=" << NOP << endl;
    // integrate for all particles
    for (int i=0; i<NOP; i++)
        P[i]->integrateFieldTrace(ObservationPoint, t0, dt, nots, ObservationField);
    // cout << endl;
}

void Bunch::integrateFieldTrace(
        Vector ObservationPoint,
        FieldTrace *trace)
{
    for (int i=0; i<NOP; i++)
        P[i]->integrateFieldTrace(ObservationPoint, trace);
}

