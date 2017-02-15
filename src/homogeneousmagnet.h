#ifndef HOMOGENEOUSMAGNET_H
#define HOMOGENEOUSMAGNET_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "externalfield.h"
#include "vector.h"

using namespace std;

//a class for a homogenous magnetic field 
//the field would extend everywhere

class HomogeneousMagnet : public ExternalField
{

  public:

    // constructor
    HomogeneousMagnet(Vector B);        // field vector B[T]

    Vector getB0();			// report the field vector
    
  private:

    Vector ElementLocalEField(double t, Vector X);

    Vector ElementLocalBField(double t, Vector X);

  private:
      
    Vector B0;				// the constant field vector

};

#endif
