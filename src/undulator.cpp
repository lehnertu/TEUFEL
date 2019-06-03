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
    // non-zero part of the field
    if (z1 <= X.z && X.z < z2)
    {
        B.x = 0.0;
        B.y = BPeak * sin(kz * X.z) * cosh(kz * X.y);
        B.z = BPeak * cos(kz * X.z) * sinh(kz * X.y);
        // ramp in over one LambdaU
        if (X.z < z1+LambdaU) B *= (X.z - z1) / LambdaU;
        // ramp out over one LambdaU
        if (z2-LambdaU <= X.z) B *= (z2 - X.z) / LambdaU;
    };
    return ElMagField(E, B);
}

TransverseGradientUndulator::TransverseGradientUndulator() :
    ExternalField()
{
    Setup(0.0, 1.0, 0.0, 1.0, 1);
}

TransverseGradientUndulator::TransverseGradientUndulator(Vector pos) :
    ExternalField(pos)
{
    Setup(0.0, 1.0, 0.0, 1.0, 1);
}

TransverseGradientUndulator::TransverseGradientUndulator(const pugi::xml_node node, InputParser *parser) :
    ExternalField()
{
    parser->parseCalcChildren(node);
    pugi::xml_node position = node.child("position");
    if (!position)
        throw(IOexception("InputParser::TransverseGradientUndulator - undulator <position> not found."));
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
        throw(IOexception("InputParser::TransverseGradientUndulator - undulator <field> not found."));
    else
    {
        double B = parser->parseValue(field.attribute("B"));
        double kx = parser->parseValue(field.attribute("kx"));
        double mod = 0.0;
        pugi::xml_attribute modatt = field.attribute("modulation");
        if (modatt) mod=parser->parseValue(modatt);
        double period = parser->parseValue(field.attribute("period"));
        int N = field.attribute("N").as_int(0);
        Setup(B, kx, mod, period, N);
    }
}

void TransverseGradientUndulator::Setup(
    double B,
    double grad,
    double modu,
    double lambda,
    int N
    )
{
    BPeak = B;
    kx = grad;
    Modulation = modu;
    LambdaU = lambda;
    NPeriods = N;
    Krms = LambdaU * SpeedOfLight * BPeak / (2.0 * Pi * mecsquared) / sqrt(2.0);
    if (teufel::rank==0)
    {
        std::cout << "transverse gradient undulator  N = " << NPeriods << ",  lambda = " << LambdaU << " m" << std::endl;
        std::cout << "  K(rms) = " << Krms << ",  Bpeak = " << BPeak << " T,  grad = " << BPeak*kx << " T/m" << std::endl;
    }
    kz = 2.0 * Pi / LambdaU;
    ky = sqrt(kz*kz+kx*kx);
}

double TransverseGradientUndulator::GetBPeak()
{
    return BPeak;
}

double TransverseGradientUndulator::GetLambdaU()
{
    return LambdaU;
}

int TransverseGradientUndulator::GetNPeriods()
{
    return NPeriods;
}

double TransverseGradientUndulator::GetKpeak()
{
    return Krms * sqrt(2.0);
}

double TransverseGradientUndulator::GetKrms()
{
    return Krms;
}

ElMagField TransverseGradientUndulator::LocalField(double t, Vector X)
{
    Vector E = Vector(0.0, 0.0, 0.0);
    Vector B = Vector(0.0, 0.0, 0.0);
    // outside the range (z1, z2) the field is zero
    double z1 = -(NPeriods+1)*LambdaU/2.0;
    double z2 = (NPeriods+1)*LambdaU/2.0;
    // non-zero part of the field
    if (z1 <= X.z && X.z < z2)
    {
        // planar undulator field
        B.x = 0.0;
        B.y = BPeak * (1.0+Modulation*sin(kz*X.z)) * sin(kz * X.z) * cosh(kz * X.y);
        B.z = BPeak * (1.0+Modulation*sin(kz*X.z)) * cos(kz * X.z) * sinh(kz * X.y);
        // add gradient field
        B.x += BPeak * (1.0+Modulation*sin(kz*X.z)) * kx/ky * cos(kx * X.x) * sinh(ky * X.y) * sin(kz * X.z);
        B.y += BPeak * (1.0+Modulation*sin(kz*X.z)) * sin(kx * X.x) * cosh(ky * X.y) * sin(kz * X.z);
        B.z += BPeak * (1.0+Modulation*sin(kz*X.z)) * kz/ky * sin(kx * X.x) * sinh(ky * X.y) * cos(kz * X.z);
        // ramp in over one LambdaU
        if (X.z < z1+LambdaU) B *= (X.z - z1) / LambdaU;
        // ramp out over one LambdaU
        if (z2-LambdaU <= X.z) B *= (z2 - X.z) / LambdaU;
    };
    return ElMagField(E, B);
}
