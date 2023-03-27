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

#include "wave.h"
#include "global.h"

#include <math.h>

GaussianWave::GaussianWave() :
    LocalizedField()
{
    Setup(1.0, complex<double>(0.0, 0.0), 1.0);
}

GaussianWave::GaussianWave(double time, Vector pos) :
    LocalizedField(time,pos)
{
    Setup(1.0, complex<double>(0.0, 0.0), 1.0);
}

GaussianWave::GaussianWave(const pugi::xml_node node, InputParser *parser) :
    LocalizedField()
{
    // parser->parseCalcChildren(node);
    pugi::xml_node position = node.child("position");
    if (!position)
        throw(IOexception("InputParser::GaussianWave - wave <position> not found."));
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
        throw(IOexception("InputParser::GaussianWave - wave <field> not found."));
    else
    {
        double Ar = parser->parseValue(field.attribute("ReA"));
        double Ai = parser->parseValue(field.attribute("ImA"));
        double l = parser->parseValue(field.attribute("lambda"));
        double r = parser->parseValue(field.attribute("rayleigh"));
        Setup(l, complex<double>(Ar,Ai), r);
    }
}

void GaussianWave::Setup(
        double lda,
        complex<double> amplitude,
        double range
    )
{
    lambda = lda;
    k = 2.0*Pi/lambda;
    omega = k*SpeedOfLight;
    A = amplitude;
    zR = range;
    w0 = sqrt(lambda*zR/Pi);
    double I0 = norm(A)/2.0*EpsNull*SpeedOfLight;
    double Ptot = Pi/2.0 * I0*w0*w0;
    if (teufel::rank==0)
        std::cout << "gaussian wave  E0 = " << sqrt(norm(A)) << " V/m   w0 = " << w0 << " m   Ptot = " << Ptot << " W" << std::endl;
}

ElMagField GaussianWave::LocalField(double t, Vector X)
{
    // cylindrical coordinates
    double z = X.z;
    double r = Vector(X.x,X.y,0.0).norm();
    // mode size
    double w = w0 * sqrt(1.0+pow(z/zR,2.0));
	// curvature radius
	double R = z + zR*zR/z;
	// field density fraction
	// the phase factor is missing because it is contained in the complex amplitude
	complex<double> dens = w0/w * exp(complex<double>(-pow(r/w,2.0), omega*t -k*z -Pi*r*r/lambda/R));
	// std::cout << "z=" << z << " m    r=" << r << " m   ";
	// std::cout << "fac=(" << dens.real() << ", " << dens.imag() << ")" << std::endl;
	Vector E = Vector( (A*dens).real(), 0.0, 0.0 );
	Vector B = Vector(0.0, (A*dens).real()/SpeedOfLight, 0.0);
    return ElMagField(E, B);
}



GaussianWavePacket::GaussianWavePacket() :
    LocalizedField()
{
    Setup(1.0, complex<double>(0.0, 0.0), 1.0, 1.0);
}

GaussianWavePacket::GaussianWavePacket(double time, Vector pos) :
    LocalizedField(time,pos)
{
    Setup(1.0, complex<double>(0.0, 0.0), 1.0, 1.0);
}

GaussianWavePacket::GaussianWavePacket(const pugi::xml_node node, InputParser *parser) :
    LocalizedField()
{
    // parser->parseCalcChildren(node);
    pugi::xml_node position = node.child("position");
    if (!position)
        throw(IOexception("InputParser::GaussianWave - wave <position> not found."));
    else
    {
        double x, y, z, t;
        x = parser->parseValue(position.attribute("x"));
        y = parser->parseValue(position.attribute("y"));
        z = parser->parseValue(position.attribute("z"));
        t = parser->parseValue(position.attribute("t"));
        origin = Vector(x,y,z);
        t0 = t;
    }
    pugi::xml_node field = node.child("field");
    if (!field)
        throw(IOexception("InputParser::GaussianWave - wave <field> not found."));
    else
    {
        double Ar = parser->parseValue(field.attribute("ReA"));
        double Ai = parser->parseValue(field.attribute("ImA"));
        double l = parser->parseValue(field.attribute("lambda"));
        double r = parser->parseValue(field.attribute("rayleigh"));
        double tau = parser->parseValue(field.attribute("tau"));
        Setup(l, complex<double>(Ar,Ai), r, tau);
    }
}

void GaussianWavePacket::Setup(
        double lda,
        complex<double> amplitude,
        double range,
        double duration
    )
{
    // origin is already set
    lambda = lda;
    k = 2.0*Pi/lambda;
    omega = k*SpeedOfLight;
    A = amplitude;
    zR = range;
    w0 = sqrt(lambda*zR/Pi);
    tau = duration;
    double I0 = norm(A)/2.0*EpsNull*SpeedOfLight;
    double Ptot = Pi/2.0 * I0*w0*w0;
    if (teufel::rank==0) {
        std::cout << "gaussian wave packet  E_peak = " << sqrt(norm(A)) << " V/m   w0 = " << w0 << " m   P_peak = " << Ptot << " W" << std::endl;
        std::cout << "t0 = " << t0 << " s   tau = " << tau << " s" << std::endl;
    }
}

ElMagField GaussianWavePacket::LocalField(double t, Vector X)
{
    // cylindrical coordinates
    double z = X.z;
    double r = Vector(X.x,X.y,0.0).norm();
    // mode size
    double w = w0 * sqrt(1.0+pow(z/zR,2.0));
	// curvature radius
	double R = z + zR*zR/z;
	// field density fraction
	// the phase factor is missing because it is contained in the complex amplitude
	complex<double> dens = w0/w *
	    exp(complex<double>(-pow(r/w,2.0)-pow((t-z/SpeedOfLight)/tau,2.0), omega*t -k*z -Pi*r*r/lambda/R));
	// std::cout << "z=" << z << " m    r=" << r << " m   ";
	// std::cout << "fac=(" << dens.real() << ", " << dens.imag() << ")" << std::endl;
	Vector E = Vector( (A*dens).real(), 0.0, 0.0 );
	Vector B = Vector(0.0, (A*dens).real()/SpeedOfLight, 0.0);
    return ElMagField(E, B);
}

