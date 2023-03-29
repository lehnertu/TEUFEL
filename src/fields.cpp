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
#include "global.h"

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
    vecE=VectorZero;
    vecB=VectorZero;
}

Vector ElMagField::E() { return vecE; }

Vector ElMagField::B() { return vecB; }

Vector ElMagField::Poynting() { return cross(vecE, vecB) / MuNull; }

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

ElMagField& ElMagField::operator-= (ElMagField other)
{
    vecE -= other.vecE;
    vecB -= other.vecB;
    return (*this);
}

ElMagField ElMagField::operator* (double factor)
{
    ElMagField temp;
    temp.vecE = vecE*factor;
    temp.vecB = vecB*factor;
    return (temp);
}

ElMagField ElMagField::operator/ (double factor)
{
    ElMagField temp;
    temp.vecE = vecE/factor;
    temp.vecB = vecB/factor;
    return (temp);
}

ElMagField& ElMagField::operator*= (double factor)
{
    vecE *= factor;
    vecB *= factor;
    return (*this);
}

ElMagField& ElMagField::operator/= (double factor)
{
    vecE /= factor;
    vecB /= factor;
    return (*this);
}

// ********** ElMagObs **********

ElMagObs::ElMagObs()
{
    t = 0.0;
    vecE = VectorZero;
    vecB = VectorZero;
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
    EB = ElMagFieldZero;
}

HomogeneousField::HomogeneousField(Vector E, Vector B)
{
    EB = ElMagField(E,B);
}

ElMagField HomogeneousField::Field(double t, Vector X)
{
    return EB;
}

// ********** LocalizedField **********

LocalizedField::LocalizedField()
{
    t0 = 0.0;
    origin = VectorZero;
}

LocalizedField::LocalizedField(double time, Vector pos)
{
    t0 = time;
    origin = pos;
    t0 = time;
}

ElMagField LocalizedField::Field(double t, Vector X)
{
    return LocalField(t-t0,X-origin);
}

ElMagField LocalizedField::LocalField(double t, Vector X)
{
    return ElMagFieldZero;
}

// ********** InteractionField **********

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
