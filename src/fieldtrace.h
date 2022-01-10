/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Copyright (c) 2017 U. Lehnert

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

#pragma once

#include <vector>
#include "vector.h"
#include "fields.h"

/*!
 * \class FieldTrace
 * \brief Type for time traces of electromagnetic fields.
 * @author Ulf Lehnert
 * @date 18.11.2019
 * 
 * This class holds a number of the electromagnetic field samples
 * covering a certain time window. Time of the first sample is t0.
 * All other samples are spaced by dt from each other.
 */
class FieldTrace
{

public:

    /*! Default constructor:<br>
     *  Allocate memory for N samples
     *  and initialize all components to zero.
     */
    FieldTrace(double t0_p, double dt_p, std::size_t N_p);

    /*! Copy constructor:<br>
     */
    FieldTrace(FieldTrace *original);

    /*! Erase the fields of a trace */
    void zero();
    
    /*! Default destructor:<br>
     *  Free memory
     */
    ~FieldTrace();

    /*! Assignment operator which copies all data of another field trace.
     *  An exception is thrown if the sizes don't match.
     */
    FieldTrace & operator=(const FieldTrace &t);

    /*! Get the number of doubles (not bytes!) for a buffer
     *  holding the trace data
     */
    // int getBufferSize() { return 6*N; };

    /*! Copy all entries of a field trace into a buffer.
     *  The memory has to be allocated by the caller beforehand.
     *  The number of elements (ElMagField) in the buffer needs to be provided for checking.
     */
    void get_buffer(ElMagField *buffer, std::size_t Nb);

    /*! Set all entries of a field trace from a buffer.
     *  The number of elements (ElMagField) in the buffer needs to be provided for checking.
     */
    void set_buffer(ElMagField *buffer, std::size_t Nb);

    /*! Multiplication of the field with a real factor */
    FieldTrace operator* (double factor);
    
    /*! Division of the field by a real factor */
    FieldTrace operator/ (double factor);
    
    /*! Sum of the fields of two traces.
     *  An exception is thrown in case of mismatch of the time or size definitions.
     */
    FieldTrace operator+ (FieldTrace other);

    /*! Accumulating sum of the fields of two traces.
     *  An exception is thrown in case of mismatch of the time or size definitions.
     */
    FieldTrace& operator+= (FieldTrace other);

    /*! Difference of the fields of two traces.
     *  An exception is thrown in case of mismatch of the time or size definitions.
     */
    FieldTrace operator- (FieldTrace other);

    /*! Get/set the start time of the trace */
    double get_t0() { return t0; };
    void set_t0(double t) { t0=t; };

    /*! Get the sampling definitions
     *  These are not expected to be changed, so, no setter routines are provided */
    double get_dt() { return dt; };
    std::size_t get_N() { return N; };
    
    /*! Get the time of one point of the trace */
    double get_time(std::size_t index);
    
    /*! Get the time of the last point of the trace */
    double get_last_time();
    
    /*! Get the center time of the trace */
    double get_center_time();
    
    /*! Get the field at one point of the trace */
    ElMagField get_field(std::size_t index);
    
    /*! Set one particular field value stored in one time slice with index it.
     *  This method throws an exception in case of an out-of-range index.
     *  This exception should not be caught as it represents an internal
     *  coding error (should never happen).
     */
    void set_field(std::size_t index, ElMagField field);

    /*! Get the field at one point in time.
        This interpolates linearly between samples.
        Outside the covered trace length zero is returned.
     */
    ElMagField get_field(double time);
    
    /*! Add some field to the value recorded for a particular time step */
    void add(std::size_t index, ElMagField f);
    
    /*! Poynting vector - time-integrated energy flow density */
    Vector Poynting();
    
    /*! compute a new trace containing time derivative of fields */
    FieldTrace* derivative();
    
    /*! Compute a retarded trace:
     *  @param delta_t - the amount of retardation in seconds
     *  @param target - the new trace (fields will be overwritten)
     *  The target trace timing will be retained.
     *  This method throws an exception if the returned field trace has zero radiation intensity.
     *  This is usually an indication for the time frame of the retarded trace
     *  not overlaping with the retarded source fields.
     */
    void retard(double delta_t, FieldTrace *target);

    /*! Transform the fields into another coordinate system.
     *  The vectors ex, ey and ez are supposed to form a right-handed cartesian system.
     */
    void transform(Vector ex, Vector ey, Vector ez);
    
private:

    std::size_t N;
    double t0;
    double dt;
    std::vector<ElMagField> trace;

};

