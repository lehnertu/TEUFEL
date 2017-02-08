/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/*
 $Log: not supported by cvs2svn $
 Revision 1.1  2006/03/15 17:10:22  shang
 first version, contains functions for physics apllications

*/

/* define structure for use with scanargs(), scanlist() */
#if !defined(PHYSICS_INCLUDED)
#define PHYSICS_INCLUDED 1

#include <assert.h>
#include "SDDS.h"
#include "mdb.h"
#undef epicsShareFuncOAGPHY

#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_OAGPHYLIB)
#define epicsShareFuncOAGPHY  __declspec(dllexport)
#else
#define epicsShareFuncOAGPHY
#endif
#else
#define epicsShareFuncOAGPHY
#endif

#define TWISS_DATA_TYPES 16
typedef struct {
  long beams; /* number of eletron beams = number of twiss input file pages.*/
  double **data; /* contains all the data for following types */
  double *betax, *alphax, *etax, *etaxp; 
  double *betay, *alphay, *etay, *etayp;
  double *ex0, *ey0, *Sdelta0, *pCentral;
  double *sigmax, *sigmaxp, *sigmay, *sigmayp;
} TWISS_PARAMETER;

epicsShareFuncOAGPHY extern long ReadTwissInput(char *inputfile,
                                                TWISS_PARAMETER *twiss,
                                                double coupling,
                                                double emitRatio,
                                                double period,
                                                long Nu);

epicsShareFuncOAGPHY extern long GetTwissValues(SDDS_DATASET *SDDSin,
                                                double *betax, double *alphax,
                                                double *etax, double *etaxp, 
                                                double *betay, double *alphay,
                                                double *etay, double *etayp, 
                                                double *ex0, double *ey0, 
                                                double *Sdelta0, double *pCentral, 
                                                double emitRatio, double coupling);
epicsShareFuncOAGPHY extern void FindPeak(double *E,double *spec,double *ep,double *sp,long n);
epicsShareFuncOAGPHY extern void ComputeBeamSize(double period, long Nu, double ex, 
                                                 double ey, double Sdelta0, 
                                                 double betax, double alphax, 
                                                 double etax, double etaxp,
                                                 double betay, double alphay, 
                                                 double etay, double etayp,
                                                 double *Sx, double *Sy, double *Sxp, double *Syp);

epicsShareFuncOAGPHY extern int Gauss_Convolve(double *E,double *spec,long *ns,double sigmaE, long *nstart, long *nend);
 
#endif

