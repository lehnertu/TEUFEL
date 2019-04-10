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

#include "undulator.h"
#include "global.h"

#include <math.h>

PlanarUndulator::PlanarUndulator() :
    ExternalField()
{
    Setup(0.0, 1.0, 1);
}

PlanarUndulator::PlanarUndulator(Vector pos) :
    ExternalField(pos)
{
    Setup(0.0, 1.0, 1);
}

PlanarUndulator::PlanarUndulator(const pugi::xml_node node, InputParser *parser) :
    ExternalField()
{
    parser->parseCalcChildren(node);
    pugi::xml_node position = node.child("position");
    if (!position)
        throw(IOexception("InputParser::PlanarUndulator - undulator <position> not found."));
    else
    {
        double x, y, z;
        x = parser->parseValue(position.attribute("x"));
        y = parser->parseValue(position.attribute("y"));
        z = parser->parseValue(position.attribute("z"));
        origin = Vector(x,y,z);
    }
    pugi::xml_node field = node.child("field");
    if (!field)
        throw(IOexception("InputParser::PlanarUndulator - undulator <field> not found."));
    else
    {
        double B = parser->parseValue(field.attribute("B"));
        double period = parser->parseValue(field.attribute("period"));
        int N = field.attribute("N").as_int(0);
        Setup(B, period, N);
    }
}

void PlanarUndulator::Setup(
    double B,
    double lambda,
    int N
    )
{
    BPeak = B;
    LambdaU = lambda;
    NPeriods = N;
    Krms = LambdaU * SpeedOfLight * BPeak / (2.0 * Pi * mecsquared) / sqrt(2.0);
    if (teufel::rank==0)
        std::cout << "planar undulator  N = " << NPeriods << ",  lambda = " << LambdaU << ",  K(rms) = " << Krms << std::endl;
    ky = 2.0 * Pi / LambdaU;
    kz = 2.0 * Pi / LambdaU;
}

double PlanarUndulator::GetBPeak()
{
    return BPeak;
}

double PlanarUndulator::GetLambdaU()
{
    return LambdaU;
}

int PlanarUndulator::GetNPeriods()
{
    return NPeriods;
}

double PlanarUndulator::GetKpeak()
{
    return Krms * sqrt(2.0);
}

double PlanarUndulator::GetKrms()
{
    return Krms;
}

ElMagField PlanarUndulator::LocalField(double t, Vector X)
{
    Vector E = Vector(0.0, 0.0, 0.0);
    Vector B = Vector(0.0, 0.0, 0.0);
    // outside the range (z1, z2) the field is zero
    double z1 = -(NPeriods+1)*LambdaU/2.0;
    double z2 = (NPeriods+1)*LambdaU/2.0;
    // ramp in over one LambdaU
    if (z1 <= X.z && X.z < z1+LambdaU)
    {
	B.y = ((X.z - z1) / LambdaU) * BPeak * sin(kz * X.z) * cosh(kz * X.y);
	B.z = ((X.z - z1) / LambdaU) * BPeak * cos(kz * X.z) * sinh(kz * X.y);
    }
    // central part of the field
    if (z1+LambdaU <= X.z && X.z < z2-LambdaU)
    {
        B.y = BPeak * sin(kz * X.z) * cosh(kz * X.y);
        B.z = BPeak * cos(kz * X.z) * sinh(kz * X.y);
    }
    // ramp out over one LambdaU
    if (z2-LambdaU <= X.z && X.z <= z2)
    {
	B.y = ((z2 - X.z) / LambdaU) * BPeak * sin(kz * X.z) * cosh(kz * X.y);
	B.z = ((z2 - X.z) / LambdaU) * BPeak * cos(kz * X.z) * sinh(kz * X.y);
    }
    return ElMagField(E, B);
}
