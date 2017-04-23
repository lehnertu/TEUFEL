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

#pragma once

#include <complex>
#include "vector.h"
#include "particle.h"

/*!
 * \class PointObserver
 * \brief Observer of emitted radiation at a single point in space.
 * @author Ulf Lehnert
 * @date 22.4.2017
 * 
 * This class handles the computation and storage of emitted electromagnetic
 * radiation from particles or ensembles of particles.
 */
class PointObserver
{
    
public:
    
    /*! Standard constructor:<br>
     * The observer is placed at a given point in space.
     */
    PointObserver(Vector position);
    
    /*! Collect the time-domain field data of a single particle.
     *  This includes both "static" and "radiative" fields.
     * 
     *  The spacing of the observation times is non-equidistant.
     *  It follows from the retardation of the tracking time steps of the
     *  observed particle.
     * 
     *  @return the number of observations in time
     */
    int GetTimeDomainTrace(
	ChargedParticle *source
    );
    
    /*! Complex amplitude of the observation for a given frequency [Hz].
     * 
     * The 3 components of the electric field are analyzed.
     * 
     * This method requires that the a trace of observed field values has been
     * collected before and stored in PointObserver::ObservationTime
     * and PointObserver::ObservationField. Such a trace can be used for multiple
     * calls of the method in sequence - it is not altered.
     */
    void FrequencyObservation(
	double freq,
	std::complex<double> *Ex,
	std::complex<double> *Ey,
	std::complex<double> *Ez
    );
    
    /*! Write a time-domain field trace into an SDDS file.
     * 
     * The file contains one table with 7 columns
     * - observation time [s]
     * - 3 componenets of the electric field [V/m]
     * - 3 componenets of the magnetic field [T]
     * 
     * @return values for error checks:
     *	 
     *	0  -  successfully Written the file\n
     *	1  -  error in SDDS_InitializeOutput \n
     *	2  -  error in SDDS_DefineSimpleParameter \n
     *	3  -  error in SDDS_DefineColumn \n
     *	4  -  error in SDDS_WriteLayout \n
     *	5  -  error in SDDS_StartPage \n
     *	6  -  error in SDDS_SetParameters \n
     *	7  -  error in SDDS_SetRowValues \n
     *	8  -  error in SDDS_WritePage \n
     *	9  -  error in SDDS_Terminate \n
     * 
     */
    int WriteTimeTraceSDDS(const char *filename);
    
    /*! Write a frequency spectrum into an SDDS file.
     * 
     * The file contains one table with 7 columns
     * - observation Frequency f [Hz]
     * - Amplitude Ax and phase Px of the x component of the electric field
     * - Amplitude Ay and phase Py of the y component of the electric field
     * - Amplitude Az and phase Pz of the z component of the electric field
     * 
     * @param filename name of the SDDS file to be created
     * @param fstart first frequency [Hz] to write into the file
     * @param fstop last frequency [Hz] to write into the file
     * @param fstep frequency step [Hz] of the list
     * 
     * @return values for error checks:
     *	 
     *	0  -  successfully Written the file\n
     *	1  -  error in SDDS_InitializeOutput \n
     *	2  -  error in SDDS_DefineSimpleParameter \n
     *	3  -  error in SDDS_DefineColumn \n
     *	4  -  error in SDDS_WriteLayout \n
     *	5  -  error in SDDS_StartPage \n
     *	6  -  error in SDDS_SetParameters \n
     *	7  -  error in SDDS_SetRowValues \n
     *	8  -  error in SDDS_WritePage \n
     *	9  -  error in SDDS_Terminate \n
     * 
     */
    int WriteSpectrumSDDS(
	const char *filename,
	double fstart,
	double fstop,
	double fstep);
 
private:
    
    //! the position of the observer
    Vector Pos;
    
    //! number of time steps in the observed trace
    int NOTS;
    
    //! the times at which observed fields are reported
    std::vector<double> ObservationTime;
    
    //! the observed electromagnetic field
    std::vector<ElMagField> ObservationField;
    
};
