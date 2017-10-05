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

#include "fields.h"

// ********** ElMagField **********

ElMagField::ElMagField()
{
    Zero();
}

ElMagField::ElMagField(Vector E, Vector B)
{
    vecE=E;
    vecB=B;
}

void ElMagField::Zero()
{
    vecE=Vector(0.0,0.0,0.0);
    vecB=Vector(0.0,0.0,0.0);
}

Vector ElMagField::E() { return vecE; }

Vector ElMagField::B() { return vecB; }

ElMagField ElMagField::operator+ (ElMagField other)
{
    ElMagField temp;
    temp.vecE = vecE + other.vecE;
    temp.vecB = vecB + other.vecB;
    return (temp);
}

ElMagField& ElMagField::operator+= (ElMagField other)
{
    vecE += other.vecE;
    vecB += other.vecB;
    return (*this);
}

ElMagField ElMagField::operator- (ElMagField other)
{
    ElMagField temp;
    temp.vecE = vecE - other.vecE;
    temp.vecB = vecB - other.vecB;
    return (temp);
}

ElMagField ElMagField::operator* (double factor)
{
    ElMagField temp;
    temp.vecE = vecE*factor;
    temp.vecB = vecB*factor;
    return (temp);
}

// ********** ElMagObs **********

ElMagObs::ElMagObs()
{
    t = 0.0;
    vecE = Vector(0.0,0.0,0.0);
    vecB = Vector(0.0,0.0,0.0);
}

ElMagObs::ElMagObs(double time, Vector E, Vector B)
{
    t = time;
    vecE = E;
    vecB = B;
}

double ElMagObs::Time()
{
    return t;
}

Vector ElMagObs::E()
{ 
    return vecE;
}

Vector ElMagObs::B()
{ 
    return vecB;
}

ElMagField ElMagObs::EB()
{ 
    return ElMagField(vecE,vecB);
}

// ********** GeneralField **********

// ********** HomogeneousField **********

HomogeneousField::HomogeneousField()
{
    EB = ElMagField(Vector(0.0,0.0,0.0),Vector(0.0,0.0,0.0));
}

HomogeneousField::HomogeneousField(Vector E, Vector B)
{
    EB = ElMagField(E,B);
}

ElMagField HomogeneousField::Field(double t, Vector X)
{
    return EB;
}

// ********** ExternalField **********

ExternalField::ExternalField()
{
    origin = Vector(0.0,0.0,0.0);
}

ExternalField::ExternalField(Vector pos)
{
    origin = pos;
}

ElMagField ExternalField::Field(double t, Vector X)
{
    return LocalField(t,X-origin);
}

ElMagField ExternalField::LocalField(double t, Vector X)
{
    return ElMagField(Vector(0.0,0.0,0.0),Vector(0.0,0.0,0.0));
}

// ********** Lattice **********

Lattice::Lattice()
{
}

Lattice::~Lattice()
{
    // destroy all objects belonging to the lattice
    while (!elements.empty())
    {
	delete elements.back();
	elements.pop_back();
    }
}

void Lattice::addElement(GeneralField *element)
{
    elements.push_back(element);
}

int Lattice::count() { return elements.size(); }

ElMagField Lattice::Field(double t, Vector X)
{
    // this initalizes with zero
    ElMagField total = ElMagField();
    // sum up all elementary fields
    for (unsigned int i=0; i<elements.size(); i++)
    {
	total += elements.at(i)->Field(t, X);
    };
    return total;
}
