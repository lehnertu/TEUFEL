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
 Revision 1.13  2006/10/19 21:49:24  soliday
 Fixed issue with tsetFlags

 Revision 1.12  2006/05/22 22:48:01  jiaox
 Added gray keyword to -shade option

 Revision 1.11  2005/04/26 22:21:40  shang
 added xlabelScale and ylabelScale to plot contour functions

 Revision 1.10  2005/03/07 22:48:49  shang
 added colorSymbol and colorUnits arguments to go_shade_grid()

 Revision 1.9  2004/09/14 19:05:10  soliday
 Added tsetFlags to go_plot_contours and go_shade_grid

 Revision 1.8  2002/12/04 17:55:59  soliday
 Modifed the do_plot_contours arguments

 Revision 1.7  2002/08/14 15:40:14  soliday
 Added Open License

 Revision 1.6  2002/01/10 13:45:26  borland
 Added arguments to support thickness option on sddscontour.

 Revision 1.5  2001/08/29 19:13:47  soliday
 Added the ability to use a dynamic dx and dy.

 Revision 1.4  2001/06/04 20:21:46  soliday
 Added the -layout option.

 Revision 1.3  2001/01/16 20:32:14  norume
 Make prototypes match functions.

 Revision 1.2  1995/09/05 21:15:00  saunders
 First test release of the SDDS1.5 package.

*/

#include <stdio.h>

/* prototypes for code in file contour6.c */
void swap_xy(char **label, double *xmin, double *ymin, double *dx, double *dy,
    int *nx, int *ny, double ***data_value);
#if defined(GNU_C)
extern char *getenv(char *env_var);
#endif

/* prototypes for code in file contour_mask.c */
int go_mask_contour_data(char *maskfile, double **data_value, long nx, 
    long ny, double xmin, double ymin, double dx, double dy);
int contour_mask(double **data, int nxd, int nyd, double xmind, double ymind, 
    double dxd, double dyd, double **mask, int nxm, int nym, 
    double xminm, double yminm, double dxm, double dym);

/* prototypes for code in file contour_read.c */
int contour_read(
    char *file,
    char **label,
    char **contour_quantity,
    double *xmin, double *ymin, double *dx, double *dy,
    int *nx, int *ny,
    double ***data,
    int *n_contours,
    double **contour_level
    );

/* prototypes for code in file contour_sample2.c */
double interpolate_2d(double f00, double f10, double f01, double f11, 
    double x00, double y00, double dx, double dy, double x, double y);
double weighted_sum(
    double f00, double f10, double f01, double f11, double x00, double y00,
    double dx, double dy, double x, double y
    );
int go_sample_line(char *sample_file, char *x_variable, char *y_variable, 
    char *topline, char *title, char *contour_quantity, char *inputfile, 
    double **data_value, double xmin_sc, double ymin_sc, double dx, 
    double dy, long nx, long ny, double slope, double intercept, 
    double xmin_samp, double xmax_samp, long n_samples, 
    long abscissa_code);
double interpolate_value(double x, double y, long ix, long iy, 
    double **data, double xmin, double ymin, double dx, double dy, 
    long nx, long ny);

/* prototypes for code in file contour_window.c */
int contour_window(double ***new_data, int *new_nx, int *new_ny, 
    double *new_xmin, double *new_xmax, double *new_ymin, 
    double *new_ymax, double **old_data, int old_nx, int old_ny, 
    double old_xmin, double old_ymin, double dx, double dy);

/* prototypes for code in file contour_write.c */
int contour_write(
    char *file,
    char **label,
    char *contour_quantity,
    double xmin, double ymin, double dx, double dy,
    int nx, int ny,
    double **data,   
    int n_contours,
    double *contour_level
    );

/* prototypes for code in file draw_contours.c */
void find_intersection(double x, double y, double f1, double f2, double cval,
    double xx[4], double yy[4], long *nn);
void draw_contours(double **fxy, double xmin, double xmax, double ymin, double ymax,
                   long nx, long ny, double *cval, long nc, long thickness);
void find_intersection(double x, double y, double f1, double f2, double cval,
    double xx[4], double yy[4], long *nn);

/* prototypes for code in file path_integral2.c */
int go_do_path_integral(char *pathfile, char *listfile);

/* prototypes for code in file plot_contours4.c */
void label_contour(
    double **data, long nx, long ny, double value, char *string, double xmin,
    double dx, double ymin, double dy, double error_limit, long label_mode);
void label_contours(
    double **data, long nx, long ny, double *contour_level, long n_contours,
    long label_interval, long label_offset, double xmin, double dx, double ymin,
    double dy);
void plot_shapes(char *shapes, long do_swap_xy, int linetype);
typedef struct {
  char *filename;
  char *xColumn, *yColumn;
  double **xData, **yData;
  long *nPoints, nPages;
  double scale;
  short plotSymbols, fill; /* if non-zero, lineType is interpreted as symbol type */
  long lineType, thickness;
} SHAPE_DATA;
void go_plot_contours(
    char *device, char *title, char *xvar, char *yvar, char *topline,
    double **data, double xmin, double xmax, double ymin, double ymax,
    double dx, double dy, long nx, long ny,
    double *contour_level, long n_contours, long contour_label_interval,
    long contour_label_offset, long layout[2], long ix, long iy,
    char *shapes, int *pen, long flags, long pause_interval,
    SHAPE_DATA *shape, long nshapes, unsigned long long tsetFlags, double xlabelScale, double ylabelScale, short no_setup,
    long thickness);
void label_contours(
    double **data,
    long nx,
    long ny,
    double *contour_level,
    long n_contours,
    long label_interval,
    long label_offset,
    double xmin,
    double dx,
    double ymin,
    double dy
    );
void label_contour(
    double **data,
    long nx,
    long ny,
    double value,
    char *string,
    double xmin,
    double dx,
    double ymin,
    double dy,
    double error_limit,
    long label_mode  
    );
void plot_shapes(char *shapes, long do_swap_xy, int linetype);

/* prototypes for code in file print_contour.c */
int go_dump_contour_data(char *dump_file, char *x_variable, char *y_variable, 
    char *topline, char *title, char *contour_quantity, char *inputfile,
    double **data_value, double xmin_sc, double ymin_sc, double dx, 
    double dy, long nx, long ny);

void go_shade_grid(
    char *device, char *title, char *xvar, char *yvar, char *topline,
    double **data, double xmin, double xmax, double ymin, double ymax,
    double *xintervals, double *yintervals, long nx, long ny,
    double minlev, double maxlev, long nlev,
    double hue0, double hue1, long layout[2], long ix, long iy,
    char *shapes, int *pen, long flags, long pause_interval,
    long thickness, unsigned long long tsetFlags, char *colorName, char *colorUnits, double xlabelScale, double ylabelScale,long gray);

void shade_grid(double **fxy, double xmin, double ymin, double dx, double dy, 
		double *xintervals, double *yintervals, long nx, long ny, 
                double *min, double *max, double hue0, double hue1, long nlev, long reverse, long flags);


double **fft_interpolation_index1(double **data, long nx, long ny, long nx_mult, long lowpass, long flags);
double **fft_interpolation_index2(double **data, long nx, long ny, long ny_mult, long lowpass, long flags);
#define CONTOUR_FLOOR 1
#define CONTOUR_CEILING 2
#define CONTOUR_ANTIRIPPLE 4

# if __WORDSIZE == 64
#define TICKSET_GIVEN          0x000000001UL
#define TICKSET_XGRID          0x000000002UL
#define TICKSET_YGRID          0x000000004UL
#define TICKSET_XLINETYPE      0x000000008UL
#define TICKSET_YLINETYPE      0x000000010UL
#define TICKSET_XFRACTION      0x000000020UL
#define TICKSET_YFRACTION      0x000000040UL
#define TICKSET_XDIVISIONS     0x000000080UL
#define TICKSET_YDIVISIONS     0x000000100UL
#define TICKSET_XSPACING       0x000000200UL
#define TICKSET_YSPACING       0x000000400UL
#define TICKSET_XLOGARITHMIC   0x000000800UL
#define TICKSET_YLOGARITHMIC   0x000001000UL
#define TICKSET_LINETYPE       0x000002000UL
#define TICKSET_FRACTION       0x000004000UL
#define TICKSET_XMODULUS       0x000008000UL
#define TICKSET_YMODULUS       0x000010000UL
#define TICKSET_XFACTOR        0x000020000UL
#define TICKSET_YFACTOR        0x000040000UL
#define TICKSET_XTIME          0x000080000UL
#define TICKSET_YTIME          0x000100000UL
#define TICKSET_XNONEXPLABELS  0x000200000UL
#define TICKSET_YNONEXPLABELS  0x000400000UL
#define TICKSET_XOFFSET        0x000800000UL
#define TICKSET_YOFFSET        0x001000000UL
#define TICKSET_XINVERT        0x002000000UL
#define TICKSET_YINVERT        0x004000000UL
#define TICKSET_XSCALECHAR     0x008000000UL
#define TICKSET_YSCALECHAR     0x010000000UL
#define TICKSET_XTHICKNESS     0x020000000UL
#define TICKSET_YTHICKNESS     0x040000000UL
#define TICKSET_THICKNESS      0x080000000UL
# else
#define TICKSET_GIVEN          0x000000001ULL
#define TICKSET_XGRID          0x000000002ULL
#define TICKSET_YGRID          0x000000004ULL
#define TICKSET_XLINETYPE      0x000000008ULL
#define TICKSET_YLINETYPE      0x000000010ULL
#define TICKSET_XFRACTION      0x000000020ULL
#define TICKSET_YFRACTION      0x000000040ULL
#define TICKSET_XDIVISIONS     0x000000080ULL
#define TICKSET_YDIVISIONS     0x000000100ULL
#define TICKSET_XSPACING       0x000000200ULL
#define TICKSET_YSPACING       0x000000400ULL
#define TICKSET_XLOGARITHMIC   0x000000800ULL
#define TICKSET_YLOGARITHMIC   0x000001000ULL
#define TICKSET_LINETYPE       0x000002000ULL
#define TICKSET_FRACTION       0x000004000ULL
#define TICKSET_XMODULUS       0x000008000ULL
#define TICKSET_YMODULUS       0x000010000ULL
#define TICKSET_XFACTOR        0x000020000ULL
#define TICKSET_YFACTOR        0x000040000ULL
#define TICKSET_XTIME          0x000080000ULL
#define TICKSET_YTIME          0x000100000ULL
#define TICKSET_XNONEXPLABELS  0x000200000ULL
#define TICKSET_YNONEXPLABELS  0x000400000ULL
#define TICKSET_XOFFSET        0x000800000ULL
#define TICKSET_YOFFSET        0x001000000ULL
#define TICKSET_XINVERT        0x002000000ULL
#define TICKSET_YINVERT        0x004000000ULL
#define TICKSET_XSCALECHAR     0x008000000ULL
#define TICKSET_YSCALECHAR     0x010000000ULL
#define TICKSET_XTHICKNESS     0x020000000ULL
#define TICKSET_YTHICKNESS     0x040000000ULL
#define TICKSET_THICKNESS      0x080000000ULL
# endif




