/************************************************************************/
/*                                                                      */
/*  TEUFEL - THz Emission from Undulators and Free-Electron Lasers      */
/*                                                                      */
/*  written by  U.Lehnert                               12/2016         */
/*                                                                      */
/************************************************************************/

#include "cavity.h"
#include "global.h"

#include <math.h>

PillboxCavity::PillboxCavity( double freq,
                              double grad,
                              double length,
                              double rcell,
                              int    nlong )
  {
    Freq = freq;
    Grad = grad;
    Length = length;
    R = rcell;
    R2 = R*R;
    Nlong = nlong;
    kz = Nlong*Pi/Length;
    omega = 2.0*Pi*freq;
  };

Vector PillboxCavity::ElementLocalEField(double t, Vector X)
  {
    Vector E = Vector(0.0, 0.0, 0.0);
    if ((X.z>=0.0) && (X.z<=Length))
      {
        double r2 = X.x*X.x + X.y*X.y;
        if (r2<R2)
          {
            // double r = sqrt(r2);
            double Ez = (1-r2/R2)*cos(kz*X.z);
            double Er = 0.5*(1-r2/(2.0*R2))*kz*sin(kz*X.z);     // actually Er/r
            E = Vector(X.x*Er, X.y*Er, Ez) * Grad * cos(omega*t);
          }
      };
    return(E);
  }

Vector PillboxCavity::ElementLocalBField(double t, Vector X)
  {
    return Vector(0.0 ,0.0 ,0.0);
  };
