/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file   : mdb.h
 * purpose: definitions for general use, for mdblib, and for mdbmth.
 *
 * Michael Borland, 1988
 $Log: not supported by cvs2svn $
 Revision 1.131  2011/02/21 16:05:56  shang
 made the returnCode type of interp_short consistent with interp().

 Revision 1.130  2011/01/11 16:47:17  soliday
 The double_cmpdes function is now exported corretly with DLLs.

 Revision 1.129  2010/12/20 14:39:45  ywang25
  Added restartHaltonSequence and restartModHaltonSequence functions to reinitialize (optimized) Halton sequence.

 Revision 1.128  2010/10/19 14:29:09  ywang25
 Added the enforceVariableLimits function declaration for parallel swarm optimization.

 Revision 1.127  2010/06/23 18:48:29  shang
 added interp_short routine

 Revision 1.126  2010/06/11 14:25:05  borland
 Added prototype for index_min_max_long().

 Revision 1.125  2010/03/31 18:33:48  soliday
 Updated so that cpu_load and page_faults are exported properly for
 shared libraries.

 Revision 1.124  2010/02/04 23:42:53  soliday
 Updated so that functions that use complex variables can
 only be seen from c++

 Revision 1.123  2010/01/14 20:01:38  shang
 added modified halton sequence routines

 Revision 1.122  2009/12/18 21:12:14  shang
 moved definition of  czarray_2d(), free_czarray_2d(), and resize_czarray_2d() routines from elegant/track.h

 Revision 1.121  2009/12/07 18:32:00  soliday
 Added C99 support for Apple.

 Revision 1.120  2009/12/04 21:22:26  soliday
 The complex.h header is no longer included with C++ files because
 of a problem with sddspcas on linux.

 Revision 1.119  2009/12/03 17:03:38  soliday
 Fixed an issue on Solaris with the last change.

 Revision 1.118  2009/12/02 22:31:22  soliday
 Added complex number support for non C99 compilers.

 Revision 1.117  2009/12/01 21:38:01  soliday
 Added a workaround for C++

 Revision 1.116  2009/11/06 23:40:12  borland
 Changed for compatibility with C99 complex number extension.

 Revision 1.115  2009/09/23 21:46:10  borland
 Added prototype for mtimes().

 Revision 1.114  2009/06/21 18:12:38  borland
 Updated prototype for rk_odeint3_na().

 Revision 1.113  2009/04/27 17:11:48  soliday
 Added missing get_long1 declaration.

 Revision 1.112  2009/04/21 20:33:42  shang
 removed numerical recipe routines

 Revision 1.111  2009/04/21 19:05:42  shang
 removed dbeskv since it is from numerical receipes

 Revision 1.110  2009/04/21 15:27:25  shang
 added gy, k13, k23, polint, qromb, and trapzd

 Revision 1.109  2008/11/11 21:05:18  soliday
 Updated to fix an issue on vxWorks.

 Revision 1.108  2008/08/13 15:53:11  borland
 Added prototype for unweightedLinearFitSelect().

 Revision 1.107  2008/03/24 21:31:16  soliday
 Added strcmp_skip.

 Revision 1.106  2008/03/07 16:46:59  soliday
 Added the MakeMonthlyGenerationFilename definition.

 Revision 1.105  2007/07/13 15:03:19  soliday
 Renamed strcpy_s to strcpy_ss because Microsoft is already using this name.

 Revision 1.104  2007/04/17 17:29:43  soliday
 Fixed issue on WIN32

 Revision 1.103  2007/04/17 16:13:56  soliday
 Fixed isinf definition on Solaris 10 with gcc.

 Revision 1.102  2007/04/16 21:04:51  soliday
 Moved some procedures from mdbcommon to mdbmth.

 Revision 1.101  2007/04/03 19:40:09  soliday
 Fixed warning message on Solaris 10

 Revision 1.100  2007/02/02 18:36:01  borland
 Added prototype for strcpy_s() (Y. Wang).

 Revision 1.99  2007/02/01 21:28:54  soliday
 added a definition for int64_t for windows.

 Revision 1.98  2006/10/20 15:32:23  soliday
 Changed USE_TIMETAG.

 Revision 1.97  2006/10/20 15:22:21  soliday
 Updated some definitions used by SDDSepics programs.

 Revision 1.96  2006/10/19 17:55:40  soliday
 Updated to work with linux-x86_64.

 Revision 1.95  2006/09/12 14:10:42  borland
 Added macros for common integer powers of numbers.

 Revision 1.94  2006/08/31 15:06:56  soliday
 Updated to work with SDDS2

 Revision 1.93  2006/05/31 17:04:18  ywang25
 Updated for parallel random number generators.

 Revision 1.92  2006/02/08 23:02:24  soliday
 Added the random_oag, gauss_rn_oag and gauss_rn_lim_oag functions.

 Revision 1.91  2005/11/16 19:04:09  shang
 added wofz() -- complex error function for mdbmth

 Revision 1.90  2005/11/04 22:46:59  soliday
 Updated code to be compiled by a 64 bit processor.

 Revision 1.89  2005/04/07 19:33:21  borland
 Added WILDCARD_MATCH ability to match_string.  Used by sddsxref -wildMatch option.

 Revision 1.88  2005/03/03 17:22:19  soliday
 Updated to work on WIN32

 Revision 1.87  2005/02/02 16:07:27  soliday
 Moved a few routines from mdblib to mdbcommon

 Revision 1.86  2004/12/21 20:14:24  shang
 added routines for finding files between dates and sort_and_return_index() function.

 Revision 1.85  2004/12/17 20:34:51  soliday
 Updated declaration for rawread, tt_attach, tt_detach

 Revision 1.84  2004/12/03 17:42:28  soliday
 Put the sys/types.h and sys/stat.h includes back in because Linux was
 complaining without them.

 Revision 1.83  2004/12/02 23:02:01  soliday
 The section of code inside the _MATCH_STRING_ ifdef in both the match_string.h
 and mdb.h is now the same. So it will not matter which one is included first.

 Revision 1.82  2004/11/04 16:12:23  shang
 added functions for reading file links, gettting file state and checking if
 file is modified.

 Revision 1.81  2004/07/16 16:28:48  shang
 added countLimit argument to despikeData()

 Revision 1.80  2004/04/01 14:39:44  borland
 Added SIMPLEX_VERBOSE_LEVEL1 and SIMPLEX_VERBOSE_LEVEL2 flags for simplexMin()
 and simplexMinimization().

 Revision 1.79  2004/02/27 16:29:15  borland
 Added prototype for trapizoidIntegration1().

 Revision 1.78  2003/12/19 19:31:50  soliday
 Added strcmp_nh funciton.

 Revision 1.77  2003/11/07 16:47:45  borland
 Added prototype for solveQuadratic().

 Revision 1.76  2003/10/07 16:21:57  shang
 added modified bessel functions: dbeskv

 Revision 1.75  2003/08/28 22:38:38  soliday
 Added some definitions for vxWorks

 Revision 1.74  2003/08/28 15:50:14  soliday
 If MIN or MAX are already defined it will not redefine them.

 Revision 1.73  2003/07/23 16:21:20  soliday
 Added isinf for WIN32

 Revision 1.72  2003/07/22 21:09:45  soliday
 Removed reference to IEEE.

 Revision 1.71  2003/07/22 20:01:18  soliday
 Added support for Kylix.

 Revision 1.70  2003/07/09 19:16:51  soliday
 Fixed bessel function prototypes.

 Revision 1.69  2003/06/18 19:52:32  borland
 Removed the query() and query_e() macros, and replaced with queryn()
 and queryn_e() macros.  The former used gets() whereas the latter do
 not.

 Revision 1.68  2003/06/14 21:35:00  borland
 Added prototypes for double-precision bessel functions in mdbmth/dbessel.c

 Revision 1.67  2003/01/21 18:55:04  borland
 simplexMin() has another argument: a factor by which to multiply the
 parameter ranges at the end of a pass to get the initial step sizes for
 the next pass.  A value of 1.0 reproduces the old behavior.

 Revision 1.66  2003/01/16 20:07:05  soliday
 Added optimAbort

 Revision 1.65  2003/01/15 22:59:18  borland
 Added SIMPLEX_START_FROM_VERTEX1 flag for simplexMin().
 Added simplexDivisor argument for simplexMin(); default value is 3
 to reproduce old behavior.

 Revision 1.64  2003/01/08 22:42:26  borland
 Added prototype for binaryArraySearch().

 Revision 1.63  2003/01/08 19:33:18  borland
 Updated prototype for binaryIndexSearch().

 Revision 1.62  2002/10/28 17:10:47  shang
 added sorting function declarations

 Revision 1.61  2002/09/24 21:03:59  borland
 Added prototypes for randomSampleMin() and randomWalkMin().  Added
 target argument for grid_sample_opt() and grid_search_min().

 Revision 1.60  2002/09/09 19:33:30  soliday
 Added the missing TouchFile declaration.

 Revision 1.59  2002/08/15 16:52:08  soliday
 *** empty log message ***

 Revision 1.58  2002/08/14 15:40:15  soliday
 Added Open License

 Revision 1.57  2002/07/24 20:41:51  shang
 added search path related and file generation related functions

 Revision 1.56  2002/07/13 21:06:11  borland
 Added SIMPLEX_RANDOM_SIGNS macro.

 Revision 1.55  2002/06/26 16:28:05  soliday
 Fixed the round definition for negative numbers.

 Revision 1.54  2002/06/19 23:12:14  borland
 Fixed spelling of Savitzky-Golay routines.

 Revision 1.53  2002/06/18 13:43:40  borland
 Added prototype for trapazoidIntegration().

 Revision 1.52  2002/02/26 03:11:01  borland
 Added prototypes for NAFF functions (fftpackC.h) and
 1d parabolic optimization routine (mdb.h).  These are both due
 to removing the NAFF functions from sddsnaff.c

 Revision 1.51  2002/02/18 17:43:25  borland
 Added prototype for simplexMinAbort.

 Revision 1.50  2002/01/07 21:32:27  borland
 Added prototype for approximate_percentiles() function.

 Revision 1.49  2001/10/12 21:07:03  soliday
 Fixed function declaration for WIN32 and Linux.

 Revision 1.48  2001/09/28 19:51:45  shang
 add prototype of substituteTagValue()

 Revision 1.47  2001/07/31 20:47:25  borland
 Added prototype for OneDScanOptimize() and updated prototype for
 simplexMin().

 Revision 1.46  2001/07/13 14:51:19  soliday
 Added replaceString declaration.

 Revision 1.45  2001/05/31 03:18:27  borland
 Changes for simplexMin().

 Revision 1.44  2000/11/21 19:08:11  borland
 Updated prototype for powellMin().

 Revision 1.43  2000/11/06 17:59:15  borland
 Added prototype for powellMin()

 Revision 1.42  2000/11/04 17:49:03  borland
 Changed prototype for nextHaltonSequencePoint and added prototype for
 startHaltonSequence.

 Revision 1.41  2000/11/02 21:26:00  borland
 Revised prototypes for zeroIntHalve, zeroNewton, zeroInterp, and
 nextHaltonSequencePoint.

 Revision 1.40  2000/11/02 19:43:58  borland
 Added prototypes for nextHaltonSequencePoint and randomizeOrder.

 Revision 1.39  2000/10/31 21:25:27  soliday
 Fixed the declaration of computeMode so that it works with WIN32.

 Revision 1.38  2000/10/11 21:45:55  soliday
 Changed definition of isinf so that the sunmath library is no longer needed.

 Revision 1.37  2000/10/07 01:16:02  borland
 Added prototype for computeMedian function.

 Revision 1.36  2000/08/17 21:27:22  soliday
 computeCorrelations is now used by elegant so declaration was changed so that
 it can be called inside a DLL.

 Revision 1.35  2000/08/10 21:10:11  soliday
 Added definition for isinf on Solaris with gcc

 Revision 1.34  2000/08/09 21:59:58  borland
 Changed prototype for savitzky-golay smoother.

 Revision 1.33  2000/04/19 17:00:10  soliday
 Borland C no longer includes mdbtc.h

 Revision 1.32  2000/04/17 20:25:22  soliday
 Added option to define binaryInsert with Borland C.

 Revision 1.31  2000/04/17 19:27:11  soliday
 Removed binaryInsert prototype when compiling with Borland C.

 Revision 1.30  2000/04/13 16:10:27  soliday
 Changed WIN32 to _WIN32

 Revision 1.29  2000/04/11 16:19:24  soliday
 Modified prototypes to work with new mdbcommon library.

 Revision 1.28  2000/04/06 22:24:53  soliday
 Added support for Borland C.

 Revision 1.27  2000/03/27 20:25:59  borland
 Added prototype for random_4().

 Revision 1.26  2000/01/18 20:47:22  soliday
 Renamed compress to compressString to avoid a conflict with ZLIB.

 Revision 1.25  1999/09/14 18:04:58  soliday
 Added export commands for WIN32 dll files.

 Revision 1.24  1999/08/03 17:55:42  soliday
 Added keep_alloc_record declaration

 Revision 1.23  1999/07/29 21:23:23  borland
 Added prototype and macro for differential equation routine.

 Revision 1.22  1999/07/22 15:35:00  soliday
 Added macros for fopen modes.

 Revision 1.21  1999/07/09 14:24:42  soliday
 Borland added shiftedLinearCorrelationCoefficient

 Revision 1.20  1999/07/01 19:25:34  borland
 Added prototypes for Savitzky-Golay filters.

 Revision 1.19  1999/05/25 18:45:59  soliday
 Altered sleep macro for WIN32, also defined popen and pclose for WIN32

 Revision 1.18  1999/05/04 14:46:42  borland
 Added some WIN32-specific conditional compilation statements.

 Revision 1.17  1999/01/07 21:45:39  borland
 Modified prototypes for simplexMin and simplexMinimization.

 Revision 1.16  1998/08/26 14:49:41  borland
 Treatment of IEEE math function isinf is now uniform.  If on solaris
 and sunmath is missing, then modify mdb.h.

 Revision 1.15  1997/03/27 22:20:30  borland
 Added TimeToEpochText().

 Revision 1.14  1997/02/05 20:50:48  saunders
 Added 'extern "C" {}' and renamed some arguments in func prototypes.

 Revision 1.13  1997/02/03 21:21:09  borland
 Added prototype and definitions for renameRobust().

 Revision 1.12  1996/10/22 18:48:17  borland
 Added prototypes for poisson statistics significance level routines.

 Revision 1.11  1996/10/07 17:28:46  borland
 Changed prototype for despikeData to reflect long integer return value.

 Revision 1.10  1996/08/26 20:08:47  borland
 Added prototypes for new mdblib routines tokenIsInteger() and tokenIsNumber().

 Revision 1.9  1996/08/16 20:03:18  borland
 Added prototype for normSigLevel().

 * Revision 1.8  1996/03/28  04:58:43  borland
 * Changed time conversion routine prototypes (long integers are now short
 * integers).
 *
 * Revision 1.7  1996/03/19  23:59:51  borland
 * Added prototypes for new time conversion functions.
 *
 * Revision 1.6  1995/12/12  03:16:46  borland
 * Changed prototype for linearCorrelationCoefficient(); added prototype for
 * compute_percentiles().
 *
 * Revision 1.5  1995/12/02  02:15:15  borland
 * Added prototype for strslide().
 *
 * Revision 1.4  1995/11/13  16:08:35  borland
 * Added prototype for function mtime().
 *
 * Revision 1.3  1995/09/12  03:19:44  borland
 * Added prototypes for wild_match_ci, strchr_ci, strcmp_ci
 *
 * Revision 1.2  1995/09/05  21:15:14  saunders
 * First test release of the SDDS1.5 package.
 *
 */
#ifndef _MDB_
#define _MDB_ 1

#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#if defined(_WIN32) && !defined(_MINGW)
typedef __int16 int16_t;
typedef unsigned __int16 uint16_t;
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
#define PRId32 "ld"
#define SCNd32 "ld"
#define PRIu32 "lu"
#define SCNu32 "lu"
#if !defined(INT32_MAX)
#define INT32_MAX (2147483647)
#endif
#else
#if defined(vxWorks)
#define PRId32 "ld"
#define SCNd32 "ld"
#define PRIu32 "lu"
#define SCNu32 "lu"
#define INT32_MAX (2147483647)
#else
#include <inttypes.h>
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if defined(SUNOS4)
#define SUN_SPARC 1
#endif

#if defined(linux)
#define LINUX 1
#endif

  /*
#if !(defined(IEEE_MATH) && (defined(SUNOS4) || defined(SOLARIS) || defined(LINUX)))
#define isinf(x) (0)
#endif
  */

#if defined(SOLARIS)
#include <ieeefp.h>
#if defined(__GNUC__)
#define isinf(x) ((x==x) && !finite(x))
#else
#if (SOLARIS < 10)
#define isinf(x) ((x==x) && !finite(x))
#endif
#endif
#endif
#if defined(_WIN32) && !defined(_MINGW) && (_MSC_VER < 1800)
#define isnan(x) _isnan(x)
#define isinf(x) (0)
#endif
#if defined(vxWorks)
#include <private/mathP.h>
#define isinf(x) isInf(x)
#define isnan(x) isNan(x)
#endif

#include <string.h>
#include <stdio.h>

#define FOPEN_WRITE_MODE "wb"
#define FOPEN_READ_MODE  "rb"
#define FOPEN_READ_AND_WRITE_MODE "r+b"

#undef epicsShareFuncMDBLIB
#undef epicsShareFuncMDBMTH
#undef epicsShareFuncMDBCOMMON
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_MDBLIB)
#define epicsShareFuncMDBLIB  __declspec(dllexport)
#else
#define epicsShareFuncMDBLIB
#endif
#if defined(EXPORT_MDBMTH)
#define epicsShareFuncMDBMTH  __declspec(dllexport)
#else
#define epicsShareFuncMDBMTH
#endif
#if defined(EXPORT_MDBCOMMON)
#define epicsShareFuncMDBCOMMON  __declspec(dllexport)
#else
#define epicsShareFuncMDBCOMMON
#endif
#else
#define epicsShareFuncMDBLIB
#define epicsShareFuncMDBMTH
#define epicsShareFuncMDBCOMMON
#endif
  
#if defined(SUNOS4) && defined(GNU_C)
/* prototypes for functions not defined in stdio: */
extern int printf(const char *format_spec, ...);
extern int fprintf(FILE *file_ptr, const char *format_spec, ...);
/* int sprintf(char *str, const char *format_spec, ...); */
extern int scanf(const char *format_spec, ...);
extern int fscanf(FILE *file_ptr, const char *format_spec, ...);
extern int sscanf(char *str, const char *format_spec, ...);
extern int fputs(const char *string, FILE *file_ptr);
extern int puts(const char *string);
extern int fputc(char c, FILE *file_ptr);
extern int fclose(FILE *file_ptr);
extern int close(int file_descriptor);
extern void perror(char *s);
extern int fseek(FILE *file_ptr, long offset, int direction);
extern int fread(void *data, int size, int number, FILE *fp);
extern int fwrite(void *data, int size, int number, FILE *fp);
extern int fflush(FILE *file_ptr);
/* prototypes for functions not fully prototyped in math.h: */
extern double   acos(double x);
extern double   asin(double x);
extern double   atan(double x);
extern double   atan2(double y, double x);
extern double   ceil(double x);
extern double   cos(double x);
extern double   cosh(double x);
extern double   exp(double x);
extern double   fabs(double x);
extern double   floor(double x);
extern double   fmod(double x, double y);
extern double   frexp(double value, int *expo);
extern double   ldexp(double value, int expo);
extern double   log(double x);
extern double   log10(double x);
extern double   modf(double value, double *iptr);
extern double   pow(double x, double y);
extern double   sin(double x);
extern double   sinh(double x);
extern double   sqrt(double x);
extern double   tan(double x);
extern double   tanh(double x);
#endif

/* double-precision Bessel functions */
#define EPS 1.0e-16
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
#define PI 3.141592653589793
epicsShareFuncMDBMTH double dbesi0(double x);
epicsShareFuncMDBMTH double dbesi1(double x);
epicsShareFuncMDBMTH double dbesj0(double x);
epicsShareFuncMDBMTH double dbesj1(double x);
epicsShareFuncMDBMTH double dbesk0(double x);
epicsShareFuncMDBMTH double dbesk1(double x);
epicsShareFuncMDBMTH double dbesy0(double x);
epicsShareFuncMDBMTH double dbesy1(double x);
/*modified bessel function dbeskv */
epicsShareFuncMDBMTH double chebev(double a, double b, double c[], int m, double x);
epicsShareFuncMDBMTH void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi);


#include <stdlib.h>

epicsShareFuncMDBLIB long PackSuffixType(char *filename, char **unpackedName, unsigned long mode);
epicsShareFuncMDBLIB void *array_1d(long size_of_elem, long lower_index, long upper_index);
epicsShareFuncMDBLIB void **array_2d(long size_of_elem, long lower1, long upper1, long lower2,
                long upper2);
epicsShareFuncMDBLIB int free_array_1d(void *array, long size_of_elem, long lower_index,
                long upper_index);
epicsShareFuncMDBLIB int free_array_2d(void **array, long size_of_elem, long lower1, long upper1,
                long lower2, long upper2);
epicsShareFuncMDBLIB void **zarray_2d(long size, long n1, long n2);
epicsShareFuncMDBLIB void **resize_zarray_2d(long size, long old_n1, long old_n2,
            void **array, long n1, long n2);
epicsShareFuncMDBLIB int free_zarray_2d(void **array, long n1, long n2);
epicsShareFuncMDBLIB void **czarray_2d(long size, long n1, long n2);
epicsShareFuncMDBLIB int free_czarray_2d(void **array, long n1, long n2);
epicsShareFuncMDBLIB void **resize_czarray_2d(void **data, long size, long n1, long n2);
epicsShareFuncMDBLIB void zero_memory(void *memory, long n_bytes);
epicsShareFuncMDBLIB int tfree(void *ptr);
epicsShareFuncMDBLIB void keep_alloc_record(char *filename);

epicsShareFuncMDBLIB void fill_int_array(int *array, long n, int value);
epicsShareFuncMDBLIB void fill_short_array(short *array, long n, short value);
epicsShareFuncMDBLIB void fill_long_array(long *array, long n, long value);
epicsShareFuncMDBLIB void fill_float_array(float *array, long n, float value);
epicsShareFuncMDBLIB void fill_double_array(double *array, long n, double value);

epicsShareFuncMDBLIB void *tmalloc(unsigned long size_of_block);
epicsShareFuncMDBLIB void *trealloc(void *ptr, unsigned long size_of_block);

/* String-related macro definitions: */
#define chop_nl(m_s) ( ((m_s)[strlen(m_s)-1]=='\n') ? (m_s)[strlen(m_s)-1]=0 : 0)

#define queryn(s, t, n) ((*(t)=0),fputs(s,stdout),fgets(t,n,stdin),chop_nl(t))
#define queryn_e(s, t, n) ((*(t)=0),fputs(s,stderr),fgets(t,n,stdin),chop_l(t))

#define is_yes(c) ((c)=='y' || (c)=='Y')
#define is_no(c) ((c)=='n' || (c)=='N')

/*   -- Data-scanning routines: */
epicsShareFuncMDBLIB extern long   query_long(char *prompt, long default_value);
epicsShareFuncMDBLIB  extern int    query_int(char *prompt, int default_value);
epicsShareFuncMDBLIB  extern short  query_short(char *prompt, short default_value);
epicsShareFuncMDBLIB extern double query_double(char *prompt, double default_value);
epicsShareFuncMDBLIB extern float  query_float(char *prompt, float default_value);
epicsShareFuncMDBLIB extern int   get_double(double *target, char *source);
epicsShareFuncMDBLIB extern int   get_long(long *target, char *source);
epicsShareFuncMDBLIB extern int   get_long1(long *iptr, char *s);
epicsShareFuncMDBLIB extern int   get_short(short *target, char *source);
epicsShareFuncMDBLIB extern int   get_int(int *target, char *source);
epicsShareFuncMDBLIB extern int   get_float(float *target, char *source);
epicsShareFuncMDBLIB extern char  *get_token(char *source);
epicsShareFuncMDBLIB extern char  *get_token_buf(char *source, char *buffer, long buffer_length);
epicsShareFuncMDBLIB extern char  *get_token_t(char *source, char *token_delimiters);
epicsShareFuncMDBLIB extern char  *get_token_tq(char *source, char *token_start,
                           char *token_end, char *quote_start, char *quote_end);
epicsShareFuncMDBLIB long tokenIsInteger(char *token);
epicsShareFuncMDBLIB long tokenIsNumber(char *token);

/*   -- String routines: */
epicsShareFuncMDBLIB extern char *trim_spaces(char *s);
epicsShareFuncMDBLIB  extern char *replace_chars(char *string, char *from, char *to);
epicsShareFuncMDBLIB  extern char *rcdelete(char *string, char c_lower, char c_upper);
epicsShareFuncMDBLIB extern char *compressString(char *string, char *chars_to_compress);
epicsShareFuncMDBLIB extern char *delete_chars(char *s, char *t);
epicsShareFuncMDBLIB extern char *delete_bounding(char *string, char *chars_to_delete);
epicsShareFuncMDBLIB extern char *str_toupper(char *string);
epicsShareFuncMDBLIB extern char *str_tolower(char *string);
epicsShareFuncMDBLIB extern long is_blank(char *string);
epicsShareFuncMDBLIB extern char *str_in(char *string, char *sub_string);
epicsShareFuncMDBLIB extern char *strcpy_ss(char *dest, const char *src);

epicsShareFuncMDBLIB  extern char *str_inn(char *string, char *sub_string, long n_char_to_check);
epicsShareFuncMDBLIB  extern char *insert(char *place_to_insert, char *string_to_insert);
epicsShareFuncMDBLIB extern char *pad_with_spaces(char *s, int n_spaces);
epicsShareFuncMDBLIB extern char *cp_str(char **target, char *source);
epicsShareFuncMDBLIB  extern char *cpn_str(char **target, char *source, long n_characters);
epicsShareFuncMDBLIB extern long edit_string(char *text, char *edit);
epicsShareFuncMDBLIB extern void edit_strings(char **string, long strings, char *buffer, char *edit);
epicsShareFuncMDBLIB char *strslide(char *s, long distance);

/* ---search path routines-- */
epicsShareFuncMDBLIB extern void setSearchPath(char *input);
epicsShareFuncMDBLIB extern char *findFileInSearchPath(char *filename);

/* --file stat routines-- */
#include <sys/types.h>
#include <sys/stat.h>
epicsShareFuncMDBLIB extern char *dir_name (const char *path);
epicsShareFuncMDBLIB extern char *read_file_link(const char *filename);
epicsShareFuncMDBLIB extern const char *read_file_lastlink(const char *filename);
epicsShareFuncMDBLIB extern char *read_last_link_to_file(const char *filename);

epicsShareFuncMDBLIB extern long get_file_stat(const char *filename, const char *lastlink, struct stat *filestat);
epicsShareFuncMDBLIB extern long file_is_modified(const char *inputfile, char **final_file, struct stat *input_stat);

/* -- find files routines ---*/
#include <ctype.h>
#if !defined(_WIN32)
#include <dirent.h>
#endif
epicsShareFuncMDBCOMMON extern short make_four_digit_year (short year);
epicsShareFuncMDBCOMMON extern long is_leap_year (short year);
epicsShareFuncMDBCOMMON extern char **find_files_between_dates(char *directory, char *rootname, char *suffix, 
                            short startYear, short startMonth, short startDay, short startJDay,
                            short endYear, short endMonth, short endDay, short endJDay,
                            char *filter, char **extensionList, long extensions,
                            long tailsOnly, long *files, long increaseOrder);
void sort_files_by_start_time(char *directory, long isTail, char **fileList, long files, long increaseOrder);

epicsShareFuncMDBCOMMON extern char **ls_dir (char *path, char *matchstr, long tailsOnly, long *files);

#if !defined(__BORLANDC__) || defined(DefineBinaryInsert)
epicsShareFuncMDBLIB long binaryInsert(void **array, long members, void *newMember, 
             int (*compare)(const void *c1, const void *c2), int32_t *duplicate);
#endif
epicsShareFuncMDBLIB long binaryIndexSearch(void **array, long members, void *key, 
                       int (*compare)(const void *c1, const void *c2), long bracket);
epicsShareFuncMDBLIB long binaryArraySearch(void *array, size_t elemSize, long members, void *key, 
                                            int (*compare)(void *c1, void *c2), long bracket);

/* sort routines (previously sort.h) */
#if !defined(_MDBSORT_INCLUDED_)
#define _MDBSORT_INCLUDED_ 1 

/*following structs and function are moved from sddsxref.c for quick sorting. */
typedef struct {
  char *stringKey;
  double doubleKey;
  long rowIndex;
} KEYED_INDEX;

typedef struct {
  KEYED_INDEX **equivalent;
  long equivalents, nextIndex;
} KEYED_EQUIVALENT;

epicsShareFuncMDBLIB extern int CompareStringKeyedIndex(const void *ki1, const void *ki2);
epicsShareFuncMDBLIB extern int CompareDoubleKeyedIndex(const void *ki1, const void *ki2);
epicsShareFuncMDBLIB extern int CompareStringKeyedGroup(const void *kg1, const void *kg2);
epicsShareFuncMDBLIB extern int CompareDoubleKeyedGroup(const void *kg1, const void *kg2);
epicsShareFuncMDBLIB extern KEYED_EQUIVALENT **MakeSortedKeyGroups(long *keyGroups, long keyType, void *data, long points);
epicsShareFuncMDBLIB extern long FindMatchingKeyGroup(KEYED_EQUIVALENT **keyGroup, long keyGroups, long keyType,void *searchKeyData, long reuse);
epicsShareFuncMDBLIB extern long *sort_and_return_index(void *data, long type, long rows, long increaseOrder);

/* sort routines (previously sort.h) */
epicsShareFuncMDBLIB extern int double_cmpasc(const void *a, const void *b);
epicsShareFuncMDBLIB extern int double_cmpdes(const void *a, const void *b);
extern void double_copy(void *a, void *b);
extern int float_cmpasc(const void *a, const void *b);
extern int float_cmpdes(const void *a, const void *b);
extern void float_copy(void *a, void *b);
epicsShareFuncMDBLIB extern int long_cmpasc(const void *a, const void *b);
extern int long_cmpdes(const void *a, const void *b);
extern void long_copy(void *a, void *b);
epicsShareFuncMDBLIB extern int string_cmpasc(const void *a, const void *b);
extern int string_cmpdes(const void *a, const void *b);
epicsShareFuncMDBLIB extern void string_copy(void *a, void *b);
epicsShareFuncMDBLIB extern int row_compare(const void *a, const void *b);
extern void row_copy(void *a, void *b);
epicsShareFuncMDBLIB extern void set_up_row_sort(int sort_by_column, size_t n_columns,
    size_t element_size,
    int (*compare)(const void *a, const void *b));
epicsShareFuncMDBLIB extern int unique(void *base, size_t n_items, size_t size,
    int (*compare)(const void *a, const void *b),
    void (*copy)(void *a, void *b));
#endif

/* string array matching (previously match_string.h): */
#if !defined(_MATCH_STRING_)
epicsShareFuncMDBLIB extern long match_string(char *string, char **option_list, long n_options,
                        long match_mode_flags);

epicsShareFuncMDBLIB extern int strncmp_case_insensitive(char *s1, char *s2, long n);
epicsShareFuncMDBLIB extern int strcmp_case_insensitive(char *s1, char *s2);
#if defined(_WIN32)
#if defined(__BORLANDC__)
#define strcasecmp(s, t) stricmp(s, t)
#define strncasecmp(s, t, n) strnicmp(s, t, n)
#else
#define strcasecmp(s, t) _stricmp(s, t)
#define strncasecmp(s, t, n) _strnicmp(s, t, n)
#endif
#endif
#define _MATCH_STRING_ 1


#define DCL_STYLE_MATCH 0
#define UNIQUE_MATCH DCL_STYLE_MATCH
#define CASE_SENSITIVE 1
#define MATCH_WHOLE_STRING 2
#define RETURN_FIRST_MATCH 8
#define EXACT_MATCH (CASE_SENSITIVE|MATCH_WHOLE_STRING|RETURN_FIRST_MATCH)
#define WILDCARD_MATCH 16

#endif

epicsShareFuncMDBLIB extern char *clean_filename(char *filename);
epicsShareFuncMDBLIB extern long fexists(const char *filename);
#define RENAME_OVERWRITE 0x0001UL
extern long renameRobust(char *oldName, char *newName, unsigned long flags);

extern char *exp_notation(double x, long n1, long n2);
extern void add_to_headers(char **header, long n_headers, char **item,
    long min_width, long format_index);
long replaceFile(char *file, char *replacement);
epicsShareFuncMDBLIB long replaceFileAndBackUp(char *file, char *replacement);
epicsShareFuncMDBLIB void add_to_standard_headers(char *name_header, char *unit_header,
    char *printf_string, char *new_name, char *new_unit, char *new_format,
    long min_width);
extern long format_length(char *format_specifier);
extern char *sbinary(char *s, int len, long n);
epicsShareFuncMDBLIB extern long bitsSet(unsigned long data);
epicsShareFuncMDBLIB extern void interpret_escapes(char *s);
epicsShareFuncMDBLIB extern int replace_string(char *target, char *source, char *orig, char *newOne);
epicsShareFuncMDBLIB int replace_stringn(char *t, char *s, char *orig, char *repl, long count_limit);
epicsShareFuncMDBLIB void interpret_escaped_quotes(char *s);
epicsShareFuncMDBLIB extern int replaceString(char *t, char *s, char *orig, 
			 char *repl, long count_limit, long here);

extern char **wild_list(int *n_names_ret, int **origin_ret, char **item_list,
        int num_items);
epicsShareFuncMDBLIB int wild_match(char *string, char *tmplate);
epicsShareFuncMDBLIB int wild_match_ci(char *string, char *tmplate);
char *strchr_ci(char *s, char c);
epicsShareFuncMDBLIB int strcmp_ci(const char *s, const char *t);
epicsShareFuncMDBLIB char *expand_ranges(char *tmplate);
epicsShareFuncMDBLIB int has_wildcards(char *tmplate);
epicsShareFuncMDBLIB char *unescape_wildcards(char *tmplate);
epicsShareFuncMDBLIB int strcmp_nh(const char *s, const char *t);
epicsShareFuncMDBLIB int strcmp_skip(const char *s1, const char *s2, const char *skip);

/*   -- Routines for flagging and aborting on errors: */
epicsShareFuncMDBLIB extern void bomb(char *error_message, char *usage_message);
epicsShareFuncMDBLIB extern long bombre(char *error_message, char *usage_message, long return_value);
extern long err_mess(long status, char *routine_name, char *message);
extern long err_mess_sys(long status, char *routine_name, char *message);
extern void fatal_err(long error_code, char *error_message);

/*   -- IO routines: */
extern char *ffgets(char *target, long target_length, FILE *file_pointer);
epicsShareFuncMDBLIB extern FILE *fopen_e(char *file_name, char *open_mode, long error_mode);
#define FOPEN_EXIT_ON_ERROR 0
#define FOPEN_RETURN_ON_ERROR 1
#define FOPEN_INFORM_OF_OPEN  2
#define FOPEN_SAVE_IF_EXISTS  4
epicsShareFuncMDBLIB extern void backspace(long n);

/*   -- Run-time statistics: */
epicsShareFuncMDBLIB extern void init_stats(void);
epicsShareFuncMDBLIB extern void report_stats(FILE *fp, char *label);
epicsShareFuncMDBLIB extern double delapsed_time();
epicsShareFuncMDBLIB extern long memory_count(void);
epicsShareFuncMDBLIB extern long cpu_time(void);
epicsShareFuncMDBLIB extern long page_faults(void);

/*   -- Miscellaneous routines: */
extern long log_usage(char *log_file_name, char *program_name);
epicsShareFuncMDBLIB extern char *tmpname(char *target);
epicsShareFuncMDBLIB char *mtime(void);
epicsShareFuncMDBLIB char *mtimes(void);
short IsLeapYear(short year);
short JulianDayFromMonthDay(short month, short day, short year, short *julianDay);
short MonthDayFromJulianDay(short julianDay, short year, short *month, short *day);
epicsShareFuncMDBLIB short TimeEpochToBreakdown(short *year, short *jDay, short *month, short *day, double *hour, double epochTime);
epicsShareFuncMDBLIB short TimeEpochToText(char *text, double epochTime);
epicsShareFuncMDBLIB short TimeBreakdownToEpoch(short year, short jDay, short month, short day, double hour, double *epochTime);
epicsShareFuncMDBLIB int makedir (char *newdir);

/********************* end of mdblib routines *****************************/

/************************ mathematical stuff *******************************/
#ifndef _MDB_MTH_
#define _MDB_MTH_ 1

#define FABS(x) fabs(x)
#define SIGN(x) ((x)<0?-1:((x)>0?1:0))
#define IS_NEGATIVE(x) (SIGN(x)==-1)
#define IS_POSITIVE(x) (SIGN(x)==1)

/* mdbmth routines */
epicsShareFuncMDBMTH long gaussianQuadrature(double (*fn)(), double a, double b, long n, double err, double *result);
epicsShareFuncMDBCOMMON extern int fixcount(char *filename, long n_points);
epicsShareFuncMDBMTH extern long factorial(long n);
epicsShareFuncMDBMTH extern double dfactorial(long n);
epicsShareFuncMDBMTH extern double g_int(double (*function)(), double lower_limit,
                    double upper_limit, long n_panels, double error_limit);
epicsShareFuncMDBMTH extern double simpson(double (*function)(), double lower_limit,
                    double upper_limit, long n_panels);
epicsShareFuncMDBMTH long trapazoidIntegration(double *x, double *y, long n, double *integral);
epicsShareFuncMDBMTH long trapazoidIntegration1(double *x, double *y, long n, double *integral);
epicsShareFuncMDBMTH int GillMillerIntegration(double *integral, double *error, double *f, double *x, long n);
epicsShareFuncMDBMTH extern double ipow(double base, long power);
epicsShareFuncMDBMTH extern double amod(double x, double m);
epicsShareFuncMDBMTH extern double zeroInterp(double (*function)(), double value,
                                              double x_initial, double x_final,
                                              double x_step, double effective_zero);
epicsShareFuncMDBMTH extern double zeroIntHalve(double (*function)(), double value,
                                                double x_initial, double x_final,
                                                double x_step, double effective_zero);
epicsShareFuncMDBMTH extern double zeroNewton(double (*function)(), 
                                              double value, double x_initial, double dx_deriv,
                                              long n_passes, double effective_zero);
#if defined(__BORLANDC__)
epicsShareFuncMDBMTH extern double poly2(double *coef, long n_coefs, double x);
#else
epicsShareFuncMDBMTH extern double poly(double *coef, long n_coefs, double x);
#endif
epicsShareFuncMDBMTH extern double dpoly(double *coef, long n_coefs, double x);
epicsShareFuncMDBMTH extern double polyp(double *coef, long *power, long n_coefs, double x);
epicsShareFuncMDBMTH extern double dpolyp(double *coef, long *power, long n_coefs, double x);
epicsShareFuncMDBMTH extern int solveQuadratic(double a, double b, double c, double *solution);
epicsShareFuncMDBMTH extern double K_cei(double k);
epicsShareFuncMDBMTH extern double E_cei(double k);
epicsShareFuncMDBMTH extern double dK_cei(double k);
epicsShareFuncMDBMTH extern double dE_cei(double k);
epicsShareFuncMDBMTH extern double celliptic(double kc, double p, double a, double b);
epicsShareFuncMDBMTH extern float drand(long dummy);
epicsShareFuncMDBMTH extern double rdrand(double lower_limit, double upper_limit);
epicsShareFuncMDBMTH extern void r_theta_rand(double *r, double *theta, double r_min,
                           double r_max);
epicsShareFuncMDBMTH extern double random_1(long iseed);
epicsShareFuncMDBMTH extern double random_2(long iseed);
epicsShareFuncMDBMTH extern double random_3(long iseed);
epicsShareFuncMDBMTH extern double random_4(long iseed);
epicsShareFuncMDBMTH extern double urandom_gauss(long iseed);
epicsShareFuncMDBMTH extern double gauss_rn(long iseed, double (*urandom)(long iseed1));
epicsShareFuncMDBMTH extern double gauss_rn_lim(double mean, double sigma, double limit_in_sigmas, double (*urandom)(long iseed));

epicsShareFuncMDBMTH extern double random_oag(long iseed, long increment);
epicsShareFuncMDBMTH extern double gauss_rn_oag(long iseed, long increment, double (*urandom)(long iseed1, long increment));
epicsShareFuncMDBMTH extern double gauss_rn_lim_oag(double mean, double sigma, double limit_in_sigmas, long increment, double (*urandom)(long iseed, long increment));

epicsShareFuncMDBMTH extern long randomizeOrder(char *ptr, long size, long length, long iseed, double (*urandom)(long iseed1));
epicsShareFuncMDBMTH extern double nextHaltonSequencePoint(long ID);
epicsShareFuncMDBMTH extern int32_t startHaltonSequence(int32_t *radix, double value);
epicsShareFuncMDBMTH extern int32_t restartHaltonSequence(long ID, double value);
epicsShareFuncMDBMTH extern double nextModHaltonSequencePoint(long ID);
epicsShareFuncMDBMTH extern int32_t startModHaltonSequence(int32_t *radix, double value);
epicsShareFuncMDBMTH extern int32_t restartModHaltonSequence(long ID, double tiny);
epicsShareFuncMDBMTH extern long convertSequenceToGaussianDistribution(double *data, long points, double limit);
epicsShareFuncMDBMTH extern double KS_Qfunction(double lambda);
extern double twoVariableKStest(double *d1, long n1, double *d2, long n2, double *MaxCDFerror);
epicsShareFuncMDBMTH extern double linearCorrelationCoefficient(double *data1, double *data2, 
                                                                short *accept1, short *accept2, 
                                                                long rows, long *count);
epicsShareFuncMDBMTH extern double linearCorrelationSignificance(double r, long rows);
epicsShareFuncMDBMTH extern double shiftedLinearCorrelationCoefficient(double *data1, double *data2, 
                                                                       short *accept1, short *accept2,
                                                                       long rows, long *count, 
                                                                       long shift);
epicsShareFuncMDBMTH extern double betaInc(double x, double a, double b);
epicsShareFuncMDBMTH extern double betaComp(double a, double b);
epicsShareFuncMDBMTH extern double gammaP(double a, double x);
epicsShareFuncMDBMTH extern double gammaQ(double a, double x);
epicsShareFuncMDBMTH extern double tTailSigLevel(double t0, long nu, long tails);
epicsShareFuncMDBMTH extern double FSigLevel(double var1, double var2, long nu1, long nu2);
epicsShareFuncMDBMTH extern double rSigLevel(double r0, long nu);
epicsShareFuncMDBMTH double ChiSqrSigLevel(double ChiSquared0, long nu);
epicsShareFuncMDBMTH double normSigLevel(double z0, long tails);
epicsShareFuncMDBMTH double poissonSigLevel(long n, double n0);

epicsShareFuncMDBMTH extern long is_prime(long number);
epicsShareFuncMDBMTH extern long smallest_factor(long number);
epicsShareFuncMDBMTH extern long next_prime_factor(long *number);
epicsShareFuncMDBMTH long largest_prime_factor(long number);
epicsShareFuncMDBMTH extern void copy_sp_array(float  **c, float  *o, long n);
epicsShareFuncMDBMTH extern void copy_dp_array(double **c, double *o, long n);
epicsShareFuncMDBMTH extern double bessel_Jn(double z, long n);
epicsShareFuncMDBMTH extern double bessel_Yn(double z, long n);
epicsShareFuncMDBMTH extern double Ai(double x), Bi(double x), Aip(double x),
        Bip(double x);
epicsShareFuncMDBMTH extern int wofz(double *xi, double *yi, double *u, double *v, long *flag);


epicsShareFuncMDBCOMMON long lsfn(double *xd, double *yd, double *sy,
    long nd, long nf, double *coef, double *s_coef,
    double *chi, double *diff);
epicsShareFuncMDBCOMMON long lsfp(double *xd, double *yd, double *sy,
    long n_pts, long n_terms, long *power, double *coef, double *s_coef,
    double *chi, double *diff);
epicsShareFuncMDBCOMMON long lsfg(double *xd, double *yd, double *sy,
    long n_pts, long n_terms, int32_t *order,
    double *coef, double *s_coef, double *chi, double *diff,
    double (*fn)(double x, long ord));

/*functions for generation file names, moved from SDDSepics.c, May 8, 2002 */
/*i.e. the declarations of functions in logfile_generation.c */
epicsShareFuncMDBMTH extern double computeYearStartTime(double StartTime);
epicsShareFuncMDBMTH extern void getTimeBreakdown(double *Time, double *Day, 
                                                     double *Hour, double *JulianDay,
                                                     double *Year, double *Month, 
                                                     char **TimeStamp);
epicsShareFuncMDBMTH extern void makeTimeBreakdown(double Time, double *ptrTime,
                                                      double *ptrDay, double *ptrHour,
                                                      double *ptrJulianDay, double *ptrYear,
                                                      double *ptrMonth, char **ptrTimeStamp);
epicsShareFuncMDBMTH extern char *makeTimeStamp(double Time);
epicsShareFuncMDBMTH extern double getTimeInSecs(void);
epicsShareFuncMDBMTH extern double getHourOfDay(void);
epicsShareFuncMDBMTH extern char *getHourMinuteSecond ();
epicsShareFuncMDBMTH extern void checkGenerationFileLocks(char *match_date);
epicsShareFuncMDBMTH extern char *MakeGenerationFilename(char *rootname, long digits, 
                                                            char *delimiter, char *lastFile);
epicsShareFuncMDBMTH extern char *MakeDailyGenerationFilename(char *rootname, long digits,
                                                                 char *delimiter, long timetag);
epicsShareFuncMDBMTH extern char *MakeMonthlyGenerationFilename(char *rootname, long digits,
                                                                char *delimiter, long timetag);
epicsShareFuncMDBMTH extern char *MakeSCRDailyTimeGenerationFilename(char *rootname);
#define DEFAULT_GENERATIONS_DIGITS 4
#define USE_TIMETAG   0x0010U
epicsShareFuncMDBMTH extern void usleepSystemIndependent(long usec);
epicsShareFuncMDBMTH extern double k13(double z);
epicsShareFuncMDBMTH extern double k23(double z);
epicsShareFuncMDBMTH extern double gy(long n, double y);
epicsShareFuncMDBMTH extern double qromb(double (*func)(), long maxe,  double a, double b,  double eps);

#define DIFFEQ_EXIT_COND_FAILED -4
#define DIFFEQ_ZERO_STEPSIZE -3
#define DIFFEQ_CANT_TAKE_STEP -2
#define DIFFEQ_OUTSIDE_INTERVAL -1
#define DIFFEQ_XI_GT_XF 0
#define DIFFEQ_SOLVED 1
#define DIFFEQ_SOLVED_ALREADY 1
#define DIFFEQ_ZERO_FOUND 2
#define DIFFEQ_END_OF_INTERVAL 3
epicsShareFuncMDBMTH char *diffeq_result_description(long return_code);

epicsShareFuncMDBMTH extern long rk_odeint(
    double *y_i, void (*derivs)(double *yp, double *y, double x),
    long n_eq, double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec, double (*exfn)(double *yp, double *y, double x),
    double exit_accuracy, long n_to_skip,
    void (*store_data)(double *yp, double *y, double x, double exval)  );
epicsShareFuncMDBMTH extern long rk_odeint1(
    double *y_i, void (*derivs)(double *yp, double *y, double x),
    long n_eq, double *accuracy, long *accmode,
    double *tiny, long *misses, double *x0, double xf, double x_accuracy,
    double h_start, double h_max, double *h_rec );
epicsShareFuncMDBMTH extern long rk_odeint2(
    double *y_i, void (*derivs)(double *qp, double *q, double t), long n_eq,
    double *accuracy, long *accmode,
    double *tiny, long *misses, double *x0, double xf, double x_accuracy,
    double h_start, double h_max, double *h_rec, double exit_value,
    long i_exit_value, double exit_accuracy, long n_to_skip );
epicsShareFuncMDBMTH extern long rk_odeint3(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec, double (*exfn)(double *yp, double *y, double x),
    double exit_accuracy);
epicsShareFuncMDBMTH extern long rk_odeint4(
    double *y_i, void (*derivs)(double *qp, double *q, double t), long n_eq,
    double *accuracy, long *accmode,
    double *tiny, long *misses, double *x0, double xf, double x_accuracy,
    double h_start, double h_max, double *h_rec, double exit_value,
    long i_exit_value, double exit_accuracy, long n_to_skip,
    void (*store_data)(double *qp, double *q, double t, double exfn) );
epicsShareFuncMDBMTH extern long rk_odeint_na(
    double *y_i, void (*derivs)(double *yp, double *y, double x),
    long n_eq, double *null1, long *null2, double *null3, long *null4,
    double *x0, double xf, double dummy1,
    double h, double dummy2, double *dummy3 );
epicsShareFuncMDBMTH extern long rk_odeint3_na(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec, double (*exfn)(double *yp, double *y, double x),
    double exit_accuracy,  void (*stochastic)(double *y, double x, double h));
epicsShareFuncMDBMTH extern long bs_odeint(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec,
    double (*exfn)(double *yp, double *y, double x), double exit_accuracy, long n_to_skip,
    void (*store_data)(double *qp, double *q, double t, double exfn_value)  );
epicsShareFuncMDBMTH extern long bs_odeint1(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec );
epicsShareFuncMDBMTH extern long bs_odeint2(
    double *y_i, void (*derivs)(double *qp, double *q, double t), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec,
    double exit_value, long i_exit_value, double exit_accuracy, long n_to_skip );
epicsShareFuncMDBMTH extern long bs_odeint3(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec, double (*exfn)(double *yp, double *y, double x),
    double exit_accuracy);
epicsShareFuncMDBMTH extern long bs_odeint4(
    double *y_i, void (*derivs)(double *qp, double *q, double t), long n_eq,
    double *accuracy, long *accmode,
    double *tiny, long *misses, double *x0, double xf, double x_accuracy,
    double h_start, double h_max, double *h_rec, double exit_value,
    long i_exit_value, double exit_accuracy, long n_to_skip,
    void (*store_data)(double *qp, double *q, double t, double exfn) );
epicsShareFuncMDBMTH extern long bss_odeint(
    float *y_i, void (*derivs)(float *yp, float *y, float x), long n_eq,
    float *accuracy, long *accmode, float *tiny, long *misses,
    float *x0, float xf, float x_accuracy, float h_start, float h_max,
    float *h_rec, float (*exfn)(float *yp, float *y, float x),
    float exit_accuracy, long n_to_skip,
    void (*store_data)(float *yp, float *y, float x, float exval)  );
epicsShareFuncMDBMTH extern long bss_odeint1(
    float *y_i, void (*derivs)(float *yp, float *y, float x), long n_eq,
    float *accuracy, long *accmode, float *tiny, long *misses, float *x0,
    float xf, float x_accuracy, float h_start, float h_max, float *h_rec );
epicsShareFuncMDBMTH extern long bss_odeint3(
    float *y_i, void (*derivs)(float *yp, float *y, float x), long n_eq,
    float *accuracy, long *accmode, float *tiny, long *misses,
    float *x0, float xf, float x_accuracy, float h_start, float h_max,
    float *h_rec, float (*exfn)(float *yp, float *y, float x),
    float exit_accuracy);
epicsShareFuncMDBMTH extern long rks_odeint(
    float *y_i, void (*derivs)(float *dydx, float *y, float x),
    long n_eq, float *accuracy, long *accmode, float *tiny, long *misses,
    float *x0, float xf, float x_accuracy, float h_start, float h_max,
    float *h_rec, float (*exit_func)(float *dydx, float *y, float x),
    float exit_accuracy, long n_to_skip,
    void (*store_data)(float *dydx, float *y, float x, float exval)  );
epicsShareFuncMDBMTH extern long rks_odeint1(
    float *y_i, void (*derivs)(float *yp, float *y, float x), long n_eq,
    float *accuracy, long *accmode, float *tiny, long *misses, float *x0,
    float xf, float x_accuracy, float h_start, float h_max, float *h_rec);
epicsShareFuncMDBMTH extern long rks_odeint2(
    float *y_i, void (*derivs)(float *yp, float *y, float x), long n_eq,
    float *accuracy, long *accmode, float *tiny, long *misses, float *x0,
    float xf, float x_accuracy, float h_start, float h_max, float *h_rec,
    float exit_value, long i_exit_value, float exit_accuracy, long n_to_skip);
epicsShareFuncMDBMTH extern long rks_odeint_na(
    float *y_i, void (*derivs)(float *yp, float *y, float x),
    long n_eq, float *null1, long *null2, float *null3, long *null4,
    float *x0, float xf, float dummy1,
    float h, float dummy2, float *dummy3 );
epicsShareFuncMDBMTH extern long stiff_odeint1(
    double *y_i, void (*derivs)(double *yp, double *y, double x),
    long n_eq, double *accuracy, long *accmode,
    double *tiny, long *misses, double *x0, double xf, double x_accuracy,
    double h_start, double h_max, double *h_rec );
epicsShareFuncMDBMTH void mmid(double *y, double *dydx, long nvar, double xs, double htot,
    long nstep, double *yout, void (*derivs)(double *dydxa, double *ya, double xa)
    );
epicsShareFuncMDBMTH void mmid2(double *y, double *dydx, long nvar, double xs, double htot,
    long nstep, double *yout, void (*derivs)(double *dydxa, double *ya, double xa)
    );
epicsShareFuncMDBMTH extern long mmid_odeint3_na(
    double *y_i, void (*derivs)(double *yp, double *y, double x), long n_eq,
    double *accuracy, long *accmode, double *tiny, long *misses,
    double *x0, double xf, double x_accuracy, double h_start, double h_max,
    double *h_rec, double (*exfn)(double *yp, double *y, double x),
    double exit_accuracy);

epicsShareFuncMDBMTH void smoothData(double *data, long rows, long smoothPoints, long smoothPasses);
epicsShareFuncMDBMTH long despikeData(double *data, long rows, long neighbors, long passes, long averageOf,
    double threshold, long countLimit);
void SavitzkyGolayCoefficients(double *coef, long maxCoefs,
                              long order, long nLeft, long nRight,
                              long derivativeOrder, long wrapAround);
epicsShareFuncMDBCOMMON long SavitzkyGolaySmooth(double *data, long rows,
                                              long order, long nLeft, 
                                              long nRight, long derivativeOrder);
epicsShareFuncMDBMTH void TouchFile(char *filename);

#define SavitzyGolaySmooth(data, rows, order, nLeft, nRight, derivativeOrder) \
SavitzkyGolaySmooth(data, rows, order, nLeft, nRight, derivativeOrder)

epicsShareFuncMDBMTH extern long optimAbort(long abort);
epicsShareFuncMDBMTH extern double minc(double (*fn)(double *param), double *x, double *dx,
    double *dx_lim, double *xlo, double *xhi, long np, long ns_max,
    long p_flag);

epicsShareFuncMDBMTH void set_argument_offset(double offset);
epicsShareFuncMDBMTH void set_argument_scale(double scale);
epicsShareFuncMDBMTH double dtcheby(double x, long n);
epicsShareFuncMDBMTH double tcheby(double x, long n);
epicsShareFuncMDBMTH double ipower(double x, long n);
epicsShareFuncMDBMTH double dipower(double x, long n);
epicsShareFuncMDBMTH double eval_sum(double (*fn)(double x, long ord), double *coef, int32_t *order, long n_coefs, double x0);
epicsShareFuncMDBMTH long powellMin(double *yReturn, double *xGuess, double *dxGuess, double *xLowerLimit,
                                    double *xUpperLimit, long dims, double target, double tolerance,
                                    double (*func)(double *x, long *invalid), 
                                    void (*report)(double ymin, double *xmin, long pass, long evals, long dims),
                                    long maxPasses, long maxEvaluations, long linMinIterations);
#define SIMPLEX_NO_1D_SCANS        0x0001U
#define SIMPLEX_RANDOM_SIGNS       0x0002U
#define SIMPLEX_START_FROM_VERTEX1 0x0004U
#define SIMPLEX_VERBOSE_LEVEL1     0x0008U
#define SIMPLEX_VERBOSE_LEVEL2     0x0010U
epicsShareFuncMDBMTH long simplexMinAbort(long abort);
epicsShareFuncMDBMTH long simplexMin(double *yReturn, double *xGuess, double *dxGuess, double *xLowerLimit,
                double *xUpperLimit, short *disable,
                long dimensions, double target, double tolerance, double (*func)(double *x, long *invalid), 
                void (*report)(double ymin, double *xmin, long pass, long n_evals, long n_dim),
                long maxEvaluations, long maxPasses, long maxDivisions, double divisorFactor, 
                double passRangeFactor, unsigned long flags);
epicsShareFuncMDBMTH long simplexMinimization(double **simplexVector, double *fValue, double *coordLowerLimit,
                         double *coordUpperLimit, short *disable, long dimensions, long activeDimensions,
                         double target, double tolerance, long tolerance_mode,
                         double (*function)(double *x, long *invalid), long maxEvaluations, 
                         long *evaluations, unsigned long flags);
epicsShareFuncMDBMTH void enforceVariableLimits(double *x, double *xlo, double *xhi, long n);
#define ONEDSCANOPTIMIZE_REFRESH 0x0001U
epicsShareFuncMDBMTH long OneDScanOptimize (double *yReturn, double *xGuess,double *dxGuess,
                 double *xLowerLimit, double *xUpperLimit,short *disable,
                 long dimensions, 
                 double target,              /* will return if any value is <= this */
                 double tolerance,           /* <0 means fractional, >0 means absolute */
                 double (*func)(double *x, long *invalid), 
                 void (*report)(double ymin, double *xmin, long pass, long evals, long dims),
                 long maxSteps,
                 long maxDivsions, long maxRepeats,
                 unsigned long flags);                                            

epicsShareFuncMDBMTH long OneDParabolicOptimization
(double *yReturn, double *xGuess, double dx,
 double xLower, double xUpper, 
 double (*func)(double x, long *invalid),
 long maxCycles, double dxLimit, double tolerance,
 long maximize);

epicsShareFuncMDBMTH long grid_search_min(double *best_result, double *best_x, double *lower, double *upper,
    double *step, long n_dimen, double target, double (*func)(double *x, long *invalid));
epicsShareFuncMDBMTH long grid_sample_min(double *best_result, double *best_x, double *lower, double *upper,
    double *step, long n_dimen, double target, double (*func)(double *x, long *invalid), double sample_fraction,
    double (*random_f)(long iseed));
epicsShareFuncMDBMTH long randomSampleMin(double *best_result, double *best_x,
                                          double *lower, double *upper, long n_dimen,
                                          double target,
                                          double (*func)(double *x, long *invalid), long nSamples,
					  double (*random_f)(long iseed));
epicsShareFuncMDBMTH long randomWalkMin(double *best_result, double *best_x,
                                          double *lower, double *upper, double *range, long n_dimen,
                                          double target,
					  double (*func)(double *x, long *invalid), long nSamples,
					  double (*random_f)(long iseed));

epicsShareFuncMDBMTH long advance_values(double *value, long *value_index, double *initial, double *step, long n_values,
    long *counter, long *max_count, long n_indices);
epicsShareFuncMDBMTH long advance_counter(long *counter, long *max_count, long n_indices);

epicsShareFuncMDBMTH extern long compute_average(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long compute_middle(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long compute_median(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long compute_percentile(double *value, double *x, long n, double percentile);
epicsShareFuncMDBMTH extern long compute_percentiles(double *value, double *percent, long values, double *x, long n);
epicsShareFuncMDBMTH extern long approximate_percentiles(double *value, double *percent, long values, double *x, long n, long bins);

epicsShareFuncMDBMTH extern long find_average(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long find_middle(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long find_median(double *value, double *x, long n);
epicsShareFuncMDBMTH extern long find_percentile(double *value, double *x, long n, double percentile);
epicsShareFuncMDBMTH extern long find_median_of_row(double *value, double **x, long index, long n);

epicsShareFuncMDBMTH extern long make_histogram(double *hist, long n_bins, double lo, double hi, double *data,
    long n_pts, long new_start);
epicsShareFuncMDBMTH extern long make_histogram_weighted(double *hist, long n_bins, double lo, double hi, double *data,
    long n_pts, long new_start, double *weight);
epicsShareFuncMDBMTH long computeMode(double *result, double *data, long pts, double binSize, long bins);

epicsShareFuncMDBMTH long findCrossingPoint(long start, double *data, long points, double level, long direction,
                       double *interpData, double *result);
epicsShareFuncMDBMTH long findTopBaseLevels(double *top, double *base, double *data, long points,
                       long bins, double sigmasRequired);

epicsShareFuncMDBMTH extern double standardDeviation(double *x, long n);
epicsShareFuncMDBMTH long unweightedLinearFit(double *xData, double *yData, long nData, double *slope, double *intercept, double *variance);
epicsShareFuncMDBMTH long unweightedLinearFitSelect(double *xData, double *yData, short *select, long nData, double *slope, double *intercept, double *variance);
epicsShareFuncMDBMTH extern double rmsValue(double *y, long n);
epicsShareFuncMDBMTH extern double arithmeticAverage(double *y, long n);
epicsShareFuncMDBMTH extern double meanAbsoluteDeviation(double *y, long n);
epicsShareFuncMDBMTH extern long computeMoments(double *mean, double *rms, double *standardDev,
          double *meanAbsoluteDev, double *x, long n);
epicsShareFuncMDBMTH extern long computeCorrelations(double *C11, double *C12, double *C22, double *x, double *y, long n);
epicsShareFuncMDBMTH extern long computeWeightedMoments(double *mean, double *rms, double *standardDev,
          double *meanAbsoluteDev, double *x, double *w, long n);
extern long accumulateMoments(double *mean, double *rms, double *standardDev,
          double *x, long n, long reset);
extern long accumulateWeightedMoments(double *mean, double *rms, double *standardDev,
          double *x, double *w, long n, long reset);
epicsShareFuncMDBMTH extern double weightedAverage(double *y, double *w, long n);
epicsShareFuncMDBMTH extern double weightedRMS(double *y, double *w, long n);
epicsShareFuncMDBMTH extern double weightedMAD(double *y, double *w, long n);
epicsShareFuncMDBMTH extern double weightedStDev(double *y, double *w, long n);

epicsShareFuncMDBMTH extern int find_min_max(double *min, double *max, double *list, long n);
epicsShareFuncMDBMTH extern int index_min_max(long *imin, long *imax, double *list, long n);
epicsShareFuncMDBMTH extern int index_min_max_long(long *imin, long *imax, long *list, long n);
extern int assign_min_max(double *min, double *max, double val);
extern int find_min_max_2d(double *min, double *max, double **value,
            long n1, long n2);
extern int find_min_max_2d_float(float *min, float *max, float **value,
    long n1, long n2);
extern int find_min(double *min, double *loc, double *c1, double *c2, long n);
extern int find_max(double *max, double *loc, double *c1, double *c2, long n);
epicsShareFuncMDBMTH extern double max_in_array(double *array, long n);
epicsShareFuncMDBMTH extern double min_in_array(double *array, long n);
epicsShareFuncMDBMTH extern void median_filter(double *x, double *m, long n, long w);

/*interpolate functions from interp.c */
typedef struct {
  double value;
  unsigned long flags;
} OUTRANGE_CONTROL;
#define OUTRANGE_VALUE       0x00000001
#define OUTRANGE_SKIP        0x00000002
#define OUTRANGE_SATURATE    0x00000004
#define OUTRANGE_EXTRAPOLATE 0x00000008
#define OUTRANGE_ABORT       0x00000010
#define OUTRANGE_WARN        0x00000020
#define OUTRANGE_WRAP        0x00000040
epicsShareFuncMDBMTH extern double interpolate(double *f, double *x, long n, double xo, 
                                               OUTRANGE_CONTROL *belowRange, 
                                               OUTRANGE_CONTROL *aboveRange, 
                                               long order, unsigned long *returnCode, long M);
epicsShareFuncMDBMTH extern double interp(double *y, double *x, long n, double x0, long warn, long order, long *returnCode);
int interpolate_minimum(double *fmin, double *zmin, double *value, double z_lo,
    double z_hi, long n);
epicsShareFuncMDBMTH double LagrangeInterp(double *x, double *f, long order, double x0, long *returnCode);


epicsShareFuncMDBLIB extern void substituteTagValue(char *input, long buflen, 
                        char **macroTag, char **macroValue, long macros); 

epicsShareFuncMDBMTH short interp_short(short *f, double *x, long n, double xo, long warnings,
                                        short order, unsigned long *returnCode, long *next_start_pos);
#define iceil(x) ((int)ceil(x))
#define round(x) ( x < 0.0 ? ((int)((x)-.5)) : ((int)((x)+.5)) )

#ifndef MIN
#define MIN(x,y) ( ((x)>(y)) ? (y) : (x))
#endif
#ifndef MAX
#define MAX(x,y) ( ((x)<(y)) ? (y) : (x))
#endif

#define SWAP_LONG(x, y) {long tmp_swap_long; tmp_swap_long=(x); (x)=(y); (y)=tmp_swap_long; }
#define SWAP_INT(x, y) {int tmp_swap_int; tmp_swap_int=(x); (x)=(y); (y)=tmp_swap_int; }
#define SWAP_SHORT(x, y) {short tmp_swap_short; tmp_swap_short=(x); (x)=(y); (y)=tmp_swap_short; }
#define SWAP_DOUBLE(x, y) {double tmp_swap_double; tmp_swap_double=(x); (x)=(y); (y)=tmp_swap_double; }
#define SWAP_FLOAT(x, y) {float tmp_swap_float; tmp_swap_float=(x); (x)=(y); (y)=tmp_swap_float; }
#define SWAP_PTR(x, y) {void *tmp_swap_ptr; tmp_swap_ptr=(x); (x)=(y); (y)=tmp_swap_ptr; }

#define INTERPOLATE(y1, y2, x1, x2, x0) (((y2)-(y1))/((x2)-(x1))*((x0)-(x1)) + (y1))

#define sqr(x)  ipow(x, 2)
#define pow2(x) sqr(x)
#define pow3(x) ipow(x, 3)
#define pow4(x) ipow(x, 4)
#define pow5(x) ipow(x, 5)

#define SQR(x)  ((x)*(x))
#define POW2(x) ((x)*(x))
#define POW3(x) (POW2(x)*(x))
#define POW4(x) (POW2(POW2(x)))
#define POW5(x) (POW4(x)*(x))

#include "constants.h"

#endif  /* _MDB_MTH_ */


#ifdef __cplusplus
#include <complex>
epicsShareFuncMDBMTH std::complex <double> complexErf(std::complex <double> z, long *flag);
epicsShareFuncMDBMTH std::complex <double> cexpi(double p);
epicsShareFuncMDBMTH std::complex <double> cipowr(std::complex <double> a, int n);
#endif

epicsShareFuncMDBMTH void complex_multiply(
                      double *r0, double *i0,    /* result */
                      double  r1, double  i1,
                      double  r2, double  i2
                      );
epicsShareFuncMDBMTH void complex_divide(
                    double *r0, double *i0,    /* result */
                    double  r1, double  i1,
                    double  r2, double  i2,
                    double threshold
                    );


#if defined(_WIN32)
#include <windows.h>
#define sleep(sec) Sleep(sec * 1000)
#if !defined(_MINGW)
#define popen(x,y) _popen(x,y)
#define pclose(x) _pclose(x)
#endif
#endif

/* machine-specific include file: */
#ifdef VAX_VMS
#include "mdbvax.h"
#endif
#ifdef SUNOS4
#include "mdbsunos4.h"
#endif
#if defined(__TURBOC__) && !defined(__BORLANDC__)
#include "mdbtc.h"
#endif

#ifdef __cplusplus
}
#endif 

#endif /* _MDB_ */

