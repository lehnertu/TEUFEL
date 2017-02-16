/*=========================================================================
 * 
 *  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers
 * 
 *  Copyright (c) 2017 U. Lehnert
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * =========================================================================*/

#include "vector.h"

#include <math.h>

Vector::Vector()
{
}

Vector::Vector(double x0, double y0, double z0)
{
  x=x0;
  y=y0;
  z=z0;
}

Vector Vector::operator+ (Vector v2)
{
  Vector temp;
  temp.x = x + v2.x;
  temp.y = y + v2.y;
  temp.z = z + v2.z;
  return (temp);
}

Vector& Vector::operator+= (Vector v2)
{
  x+=v2.x;
  y+=v2.y;
  z+=v2.z;
  return (*this);
}

Vector Vector::operator- (Vector v2)
{
  Vector temp;
  temp.x = x - v2.x;
  temp.y = y - v2.y;
  temp.z = z - v2.z;
  return (temp);
}

Vector Vector::operator- ()
{
  Vector temp;
  temp.x = -x;
  temp.y = -y;
  temp.z = -z;
  return (temp);
}

Vector Vector::operator* (double factor)
{
  Vector temp;
  temp.x = x*factor;
  temp.y = y*factor;
  temp.z = z*factor;
  return (temp);
}

Vector& Vector::operator*= (double factor)
{
  x*=factor;
  y*=factor;
  z*=factor;
  return (*this);
}

Vector Vector::operator/ (double factor)
{
  Vector temp;
  temp.x = x/factor;
  temp.y = y/factor;
  temp.z = z/factor;
  return (temp);
}

Vector& Vector::operator/= (double factor)
{
  x/=factor;
  y/=factor;
  z/=factor;
  return (*this);
}

double Vector::norm()
{
  return(sqrt(x*x+y*y+z*z));
}

double Vector::abs2nd()
{
  return(x*x+y*y+z*z);
}

void Vector::normalize()
{
  double abs = norm();
  x /= abs;
  y /= abs;
  z /= abs;
}

double dot(Vector a, Vector b)
{
  double temp;
  temp = a.x * b.x;
  temp += a.y * b.y;
  temp += a.z * b.z;
  return (temp);
}

Vector cross(Vector a, Vector b)
{
  Vector r;
  r.x = a.y*b.z - a.z*b.y;
  r.y = a.z*b.x - a.x*b.z;
  r.z = a.x*b.y - a.y*b.x;
  return (r);
}

