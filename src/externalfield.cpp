/************************************************************************/
/*                                                                      */
/*  TEUFEL - THz Emission from Undulators and Free-Electron Lasers      */
/*                                                                      */
/*  written by  U.Lehnert                               12/2016         */
/*                                                                      */
/************************************************************************/

#include "externalfield.h"

ExternalField::ExternalField()
{
}

ExternalField::~ExternalField()
{
}

Vector ExternalField::EField(double t, Vector X)
{
  return ElementLocalEField(t, X);
}

Vector ExternalField::BField(double t, Vector X)
{
  return ElementLocalBField(t, X);
}


// ********** Lattice **********

Lattice::Lattice()
{
}

Lattice::~Lattice()
{
}

void Lattice::addElement(ExternalField *element)
{
  // elements->append(element);
  elements.push_back(element);
}

Vector Lattice::EField(double t, Vector X)
{
  Vector E = Vector(0.0, 0.0, 0.0);
  for (unsigned int i=0; i<elements.size(); i++)
  {
    E += elements.at(i)->EField(t, X);
  };
  return(E);
}

Vector Lattice::BField(double t, Vector X)
{
  Vector B = Vector(0.0, 0.0, 0.0);
  for (unsigned int i=0; i<elements.size(); i++)
  {
    B += elements.at(i)->BField(t, X);
  };
  return(B);
}
