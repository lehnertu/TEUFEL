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
 Revision 1.10  2006/10/19 17:55:40  soliday
 Updated to work with linux-x86_64.

 Revision 1.9  2005/06/13 21:57:13  shang
 added CalculatePhaseAndAmplitudeFromFreq() and adjustFrequencyHalfPlane

 Revision 1.8  2004/08/12 19:56:12  soliday
 Added Mike's changes

 Revision 1.7  2004/04/07 02:01:34  borland
 Added lowerFreqLimit and upperFreqLimit to PerformNAFF() prototype.

 Revision 1.6  2003/07/22 20:01:17  soliday
 Added support for Kylix.

 Revision 1.5  2002/08/14 15:40:15  soliday
 Added Open License

 Revision 1.4  2002/02/26 03:11:01  borland
 Added prototypes for NAFF functions (fftpackC.h) and
 1d parabolic optimization routine (mdb.h).  These are both due
 to removing the NAFF functions from sddsnaff.c

 Revision 1.3  1999/09/14 18:04:14  soliday
 Added export commands for WIN32 dll files.

 Revision 1.2  1995/09/05 21:15:04  saunders
 First test release of the SDDS1.5 package.

*/

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(FFTPACKC_INCLUDE)
#define FFTPACKC_INCLUDE 1

#undef epicsShareFuncFFTPACK
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_FFTPACK)
#define epicsShareFuncFFTPACK  __declspec(dllexport)
#else
#define epicsShareFuncFFTPACK
#endif
#else
#define epicsShareFuncFFTPACK
#endif

#define INVERSE_FFT 0x0001UL
#define MINUS_I_THETA 0x0002UL
/* If MINUS_I_THETA bit is set, then fftpack routines
   operate with a different convention than the normal
   one.  The Fourier sum is
   Sum[j=0,n-1]{ A(k)*exp[-i*j*k*2*pi/n] }
   This corresponds to a exp[+i ...] in the Fourier
   decomposition.
   However, by default, the signs of i in these expressions
   are reversed.  This is the more usual convention.
 */

epicsShareFuncFFTPACK void atexitFFTpack();
epicsShareFuncFFTPACK long realFFT(double *data, long n, unsigned long flags);
epicsShareFuncFFTPACK long complexFFT(double *data, long n, unsigned long flags);
epicsShareFuncFFTPACK long realFFT2(double *output, double *input, long n, unsigned long flags);
epicsShareFuncFFTPACK void FFTderivative(double *T, double *Y, long n_pts0,
    double **T_out, double **Y_out, long *n_out, long do_pad,
    long do_truncate, long zp_spectrum);
epicsShareFuncFFTPACK extern long power_of_2(long n);
epicsShareFuncFFTPACK long dp_pad_with_zeroes(double **t, double **f, long n);

#define NAFF_RMS_CHANGE_LIMIT    0x0001U
#define NAFF_FREQS_DESIRED       0x0002U
#define NAFF_MAX_FREQUENCIES     0x0004U
#define NAFF_FREQ_CYCLE_LIMIT    0x0008U
#define NAFF_FREQ_ACCURACY_LIMIT 0x0010U
#define NAFF_FREQ_FOUND          0x0100U
epicsShareFuncFFTPACK long CalculatePhaseAndAmplitudeFromFreq
  (double *hanning, long points, double NAFFdt, double frequency, double t0,
   double *phase, double *amplitude, double *significance, double *cosine, double *sine,
   unsigned long flags);

epicsShareFuncFFTPACK long PerformNAFF(double *frequency, double *amplitude, double *phase, 
				       double *significance, double t0, double dt,
				       double *data, long points, unsigned long flags,
				       double fracRMSChangeLimit, long maxFrequencies,
				       double freqCycleLimit, double fracFreqAccuracyLimit,
				       double lowerFreqLimit, double upperFreqLimit);

epicsShareFuncFFTPACK long simpleFFT(double *magnitude2, double *data, long points);
epicsShareFuncFFTPACK double adjustFrequencyHalfPlane(double frequency, 
                                                double phase0, double phase1, double dt);

#endif

#ifdef __cplusplus
}
#endif

