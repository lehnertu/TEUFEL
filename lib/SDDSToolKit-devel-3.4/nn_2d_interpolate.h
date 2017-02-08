#ifndef _NN_2D_INTERPOLATE_H
#define _NN_2D_INTERPOLATE_H

#include <limits.h>
#include <float.h>
#include <math.h>
#include <errno.h>
#include "nn.h"
#include "nan.h"
#include "minell.h"

#if !defined(NN_SERIAL)
#define NMAX 4096
#endif
#define STRBUFSIZE 256

#define NIMAX 2048
#define BUFSIZE 10240
#define NALLOCATED_START 1024

#define NN 0
#define CSA 1

#ifdef __cplusplus
extern "C" {
#endif

#undef epicsShareFuncMDBLIB
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_MDBLIB)
#define epicsShareFuncMDBLIB  __declspec(dllexport)
#else
#define epicsShareFuncMDBLIB
#endif
#else
#define epicsShareFuncMDBLIB
#endif



#if !defined(_POINT_STRUCT)
#define _POINT_STRUCT
typedef struct {
    double x;
    double y;
    double z;
} point;
#endif

#if !defined(_SPECS_STRUCT)
#define _SPECS_STRUCT
typedef struct {
    int generate_points, thin, nointerp, linear, invariant, square, range;
    int nx, ny, nxd, nyd, nppc;
    double rmax, wmin, zoom, xmin, xmax, ymin, ymax, dx, dy, k; 
    int npoints, method;
} specs;
#endif

void generate_output_points(specs *spec, int *nout, point **pout);
/*nn routine*/
epicsShareFuncMDBLIB void do_nn_2d_interpolate(specs *spec, int *nin, point **pin, int *nout, point **pout);
epicsShareFuncMDBLIB specs *specs_create(void);

#endif
