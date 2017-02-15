/************************************************************************/
/*                                                                      */
/*  TEUFEL - THz Emission from Undulators and Free-Electron Lasers      */
/*                                                                      */
/*  written by  U.Lehnert & Vipul Joshi                 12/2016         */
/*                                                                      */
/************************************************************************/

#include "homogeneousmagnet.h"
#include "global.h"
#include <math.h>

HomogeneousMagnet::HomogeneousMagnet(Vector B)
{
  B0=B;
}

Vector HomogeneousMagnet::getB0()
{
  return B0;
}

Vector HomogeneousMagnet::ElementLocalEField(double t, Vector X)
{
  Vector p=Vector(0.0,0.0,0.0);
  return (p);
}

Vector HomogeneousMagnet::ElementLocalBField(double t, Vector X)
{
  return B0;
}


