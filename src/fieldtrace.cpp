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

#include <math.h>
#include "fields.h"
#include "global.h"
#include "fieldtrace.h"

FieldTrace::FieldTrace(double t0_p, double dt_p, std::size_t N_p)
{
    t0 = t0_p;
    dt = dt_p;
    N = N_p;
    trace = std::vector<ElMagField>(N,ElMagFieldZero);
}

FieldTrace::FieldTrace(FieldTrace *original)
{
    t0 = original->t0;
    dt = original->dt;
    N = original->N;
    trace = original->trace;
}

void FieldTrace::zero()
{
    for (std::size_t i=0; i<N; i++)
        trace[i] = ElMagFieldZero;
}

FieldTrace::~FieldTrace()
{
    // delete [] trace;
}

FieldTrace & FieldTrace::operator=(const FieldTrace &t)
{
    // check for "self assignment" and do nothing in that case
    if (this != &t)
    {
        // check for a size mismatch and throw an exception in that case
        if (N != t.N)
        {
            throw(IOexception("FieldTrace assignment - size mismatch."));
        }
        // copy all data from the other field trace
        else
        {
            t0 = t.t0;
            dt = t.dt;
            // for (int i=0; i<N; i++) trace[i]=t.trace[i];
            trace = t.trace;
        }
    }
    return *this;    
}

FieldTrace FieldTrace::operator* (double factor)
{
    FieldTrace temp = *this;
    for (std::size_t i=0; i<temp.N; i++)
        temp.trace[i] *= factor;
    return(temp);
}

FieldTrace FieldTrace::operator/ (double factor)
{
    FieldTrace temp = *this;
    for (std::size_t i=0; i<temp.N; i++)
        temp.trace[i] /= factor;
    return(temp);
}

FieldTrace FieldTrace::operator+ (FieldTrace other)
{
    if (N != other.N) throw(IOexception("FieldTrace::operator+ - size mismatch."));
    if (t0 != other.t0) throw(IOexception("FieldTrace::operator+ - start time mismatch."));
    if (dt != other.dt) throw(IOexception("FieldTrace::operator+ - time step mismatch."));
    FieldTrace temp = *this;
    for (std::size_t i=0; i<N; i++)
        temp.trace[i] = trace[i] + other.trace[i];
    return (temp);
}

FieldTrace& FieldTrace::operator+= (FieldTrace other)
{
    if (N != other.N) throw(IOexception("FieldTrace::operator+= - size mismatch."));
    if (t0 != other.t0) throw(IOexception("FieldTrace::operator+= - start time mismatch."));
    if (dt != other.dt) throw(IOexception("FieldTrace::operator+= - time step mismatch."));
    for (std::size_t i=0; i<N; i++)
        trace[i] += other.trace[i];
    return (*this);
}

FieldTrace FieldTrace::operator- (FieldTrace other)
{
    if (N != other.N) throw(IOexception("FieldTrace::operator- - size mismatch."));
    if (t0 != other.t0) throw(IOexception("FieldTrace::operator- - start time mismatch."));
    if (dt != other.dt) throw(IOexception("FieldTrace::operator- - time step mismatch."));
    FieldTrace temp = *this;
    for (std::size_t i=0; i<N; i++)
        temp.trace[i] = trace[i] - other.trace[i];
    return (temp);
}

void FieldTrace::get_buffer(ElMagField *buffer, std::size_t Nb)
{
    if (Nb != N)
        throw(IOexception("FieldTrace::get_buffer - size mismatch."));
    ElMagField *buf = buffer;
    for (std::size_t it=0; it<N; it++) *buf++ = get_field(it);
}

void FieldTrace::set_buffer(ElMagField *buffer, std::size_t Nb)
{
    if (Nb != N)
        throw(IOexception("FieldTrace::set - size mismatch."));
    ElMagField *buf = buffer;
    for (std::size_t it=0; it<N; it++) set_field(it, *buf++);
}

double FieldTrace::get_time(std::size_t index)
{
    if (index<0 || index>=N)
        throw(IOexception("FieldTrace::get_time - index out of range."));
    return t0+index*dt;
}

ElMagField FieldTrace::get_field(std::size_t index)
{
    if (index<0 || index>=N)
        throw(IOexception("FieldTrace::get_field - index out of range."));
    return trace[index];
}

void FieldTrace::set_field(std::size_t index, ElMagField f)
{
    if (index>=N)
        throw(IOexception("FieldTrace::set - index out of range."));
    trace[index] = f;
}

ElMagField FieldTrace::get_field(double time)
{
    double ref = (time-t0)/dt;
    int index = (int)floor(ref);
    double frac = ref-(double)index;
    ElMagField f1=ElMagFieldZero;
    if ((index>=0) && (index<(int)N)) f1=trace[(std::size_t)index];
    ElMagField f2=ElMagFieldZero;
    if ((index+1>=0) && (index+1<(int)N)) f2=trace[(std::size_t)(index+1)];
    return(f1*(1.0-frac)+f2*frac);
}

void FieldTrace::add(std::size_t index, ElMagField f)
{
    if (index>=N)
        throw(IOexception("FieldTrace::add - index out of range."));
    trace[index] += f;
};

Vector FieldTrace::Poynting()
{
    Vector sum(0.0,0.0,0.0);
    for (std::size_t i=0; i<N; i++)
        sum += trace[i].Poynting();
    return(sum*dt);
}

FieldTrace* FieldTrace::derivative()
{
    // create a new field trace as a copy
    FieldTrace* temp = new FieldTrace(*this);
    if (temp==0) throw(IOexception("FieldTrace::derivative() - error allocating memory."));
    // compute the new data
    temp->set_field(0, (trace[1]-trace[0])/dt );
    for (std::size_t i=1; i<N-1; i++)
        temp->set_field(i, (trace[i+1]-trace[i-1])/(2.0*dt) );
    temp->set_field(N-1, (trace[N-1]-trace[N-2])/dt );
    return(temp);
}

void FieldTrace::retard(double delta_t, FieldTrace *target)
{
    for (std::size_t i=0; i<target->get_N(); i++)
    {
        double t = target->get_time(i) - delta_t;
        target->trace[i] = get_field(t);
    };
    if (target->Poynting().norm()==0.0)
    {
        if (DEBUGLEVEL>=2)
        { 
            std::cout << "retard from t0=" << t0 << " N=" << N << " P=" << Poynting().norm() << std::endl;
            std::cout << "      by delta=" << delta_t << "s" << std::endl;
            std::cout << "retard   to t0=" << target->get_t0() << " N=" << target->get_N();
            std::cout << " P=" << target->Poynting().norm() << std::endl;
        }
        throw(IOexception("FieldTrace::retard - zero field, possibly time mismatch."));
    };
}

void FieldTrace::transform(Vector ex, Vector ey, Vector ez)
{
    for (std::size_t i=0; i<N; i++)
    {
        Vector E = trace[i].E();
        Vector B = trace[i].B();
        Vector E_new = Vector(dot(E,ex),dot(E,ey),dot(E,ez));
        Vector B_new = Vector(dot(B,ex),dot(B,ey),dot(B,ez));
        trace[i] = ElMagField(E_new, B_new);
    }
}

