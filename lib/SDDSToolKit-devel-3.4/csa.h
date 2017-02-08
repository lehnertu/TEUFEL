/******************************************************************************
 *
 * File:           csa.h
 *
 * Created:        16/10/2002
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        A header for csa library (2D data approximation with
 *                 bivariate C1 cubic spline)
 *
 * Revisions:      None
 *
 *****************************************************************************/

#if !defined(_CSA_H)
#define _CSA_H

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

extern int csa_verbose;
extern char* csa_version;

struct csa;
typedef struct csa csa;

csa* csa_create();
void csa_destroy(csa* a);
void csa_addpoints(csa* a, int n, point points[]);
void csa_addstd(csa* a, int n, double variance[]);
void csa_calculatespline(csa* a);
void csa_approximatepoint(csa* a, point* p);
void csa_approximatepoints(csa* a, int n, point* points);

void csa_setnpmin(csa* a, int npmin);
void csa_setnpmax(csa* a, int npmax);
void csa_setk(csa* a, int k);
void csa_setnppc(csa* a, int nppc);

int do_csa_2d_interpolate(specs *spec, int nin, point *pin, int *nout, point **pout, double *std);

#endif
