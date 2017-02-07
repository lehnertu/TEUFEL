/************************************************************************/
/*                                                                      */
/*  TEUFEL - THz Emission from Undulators and Free-Electron Lasers      */
/*                                                                      */
/*  written by  U.Lehnert                               12/2016         */
/*                                                                      */
/************************************************************************/

#ifndef VECTOR_H
#define VECTOR_H

using namespace std;

class Vector
{

  public:

    Vector();
    Vector(double x0, double y0, double z0);

    double x;
    double y;
    double z;

    // sum and difference of two vectors
    Vector operator+ (Vector v2);
    Vector& operator+= (Vector v2);
    Vector operator- (Vector v2);

    // negative vector
    Vector operator- ();

    // multiplication of a vector with a factor
    Vector operator* (double faktor);
    Vector& operator*= (double faktor);
    Vector operator/ (double faktor);
    Vector& operator/= (double faktor);

    // absolute value of a vector
    double norm();

    // squared absolute value of a vector
    double abs2nd();

    // make a vector unit length
    void normalize();

  private:

};

// dot product of two vectors
double dot(Vector a, Vector b);

// cross product of two vectors
Vector cross(Vector a, Vector b);

#endif