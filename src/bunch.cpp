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

#include "bunch.h"
#include "particle.h"
#include "global.h"
#include <random>

Distribution::Distribution(int dim, int nop)
{
    DIM = dim;
    NOP = nop;
    A = new double[DIM*NOP];
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i=0; i<DIM*NOP; i++)
	A[i] = dist(mt);
}

Distribution::~Distribution()
{
    delete[] A;
}

Bunch::Bunch()
{
    NOP = 0;
}

Bunch::Bunch(Bunch* origin)
{
    NOP = origin->getNOP();
    for(int i=0; i<NOP; i++)
    {
	P.push_back(new ChargedParticle(origin->getParticle(i)));
    }
}

Bunch::~Bunch()
{
}

void Bunch::Add(ChargedParticle *part)
{
    NOP++;
    P.push_back(part);
}

int Bunch::getNOP()
{
    return NOP;
}

double Bunch::getTotalCharge()
{
    double charge = 0.0;
    for(int i=0; i<NOP; i++)
    {
	charge += getParticle(i)->getCharge();
    }
    return charge;
}

ChargedParticle* Bunch::getParticle(int i)
{
    return P[i];
}
