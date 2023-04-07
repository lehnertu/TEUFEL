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

#include <vector>
#include "vector.h"

/*!
 * \class ParameterLogger
 * \brief Store beam parameters while tracking.
 * @author Ulf Lehnert
 * @date 19.10.2017
 * 
 * This class handles the computation and storage of beam parameters
 * accessible during the tracking.
 * @todo specialization for beam still missing, only for bunch implemented
 * @todo could (optionally) log lattice fields at bunch center
 */
template <class objectT>
class ParameterLogger
{

public:

    /*! Standard constructor:<br>
     * Computation and storage of beam parameters of a given beam object.
     * 
     * This class is templated and spezialized for
     * Bunch() or Beam() being observed.
     * 
     * \param[in] obj The object (bunch or beam) to be observed.
     * \param[in] filename The name of the file to write.
     * \param[in] the logging frequency in tracking steps
     */
    ParameterLogger(objectT *obj, const char *filename, int step);

    /*! Request the logger to record the bunching factor at a given frequency.
     *  The recorded value is a complex number including the phase information.
     *
     * \param[in] freq The frequency [Hz] at which to compute the bunching factor
     */
    void record_bunching(double freq);
    
    /*! Check whether the given tracking step should be logged */
    bool log_requested(int step) { return 0 == (step % stepsize); };
    
    /*! The source has advanced one time step.
     *  Compute and store the quantities of interest.
     */
    void update();

    /*! Write the collected data into an SDDS file.
     * 
     * The file contains one table with following columns
     * - observation time [s]
     * - mean beam momentum (beta*gamma)
     * - average particle position
     * - average particle momentum
     * - gamma - beam energy
     * - RMS beam size
     * - RMS beam momentum
     * - bunching factor at certain frequency
     * - correlation x-bgx y-bgy z-bgz
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
    int WriteBeamParametersSDDS();

private:

    //! the observed beam object
    objectT *Beam;
    
    //! the file name
    const char* fn;

    //! the number of stored datasets
    unsigned int NOTS;
    
    //! the logging frequency in tracking steps
    unsigned int stepsize;
    
    //! whether to include a bunching factor
    bool include_bunching;
    
    //! the frequency for the bunching factor
    double bunching_freq;
    
    //! the observation time
    std::vector<double> Time;

    //! the observation position
    std::vector<Vector> Pos;

    //! the observation beam energy
    std::vector<double> Gamma;
    
    //! the observation momentum
    std::vector<Vector> BG;
    
    //! the observation position spread
    std::vector<Vector> PosRMS;
    
    //! the observation momentum spread
    std::vector<Vector> BGRMS;

    //! the observation position-momentum correlation
    std::vector<Vector> PosBG;
    
    //! magnitude of observed bunching factors
    std::vector<double> BF;
    
    //! relative energy spread
    std::vector<double> Delta;
    
};

/*!
 * \class TrajectoryLogger
 * \brief Store particle coordinates while tracking.
 * @author Ulf Lehnert
 * @date 19.10.2017
 * 
 * This class handles the storage of particle trajectories.
 * @todo specialization for beam still missing, only for bunch implemented
 * @todo could (optionally) log lattice fields at bunch center
 */
template <class objectT>
class TrajectoryLogger
{

public:

    /*! Standard constructor:<br>
     * Storage of particle coordinates while tracking.
     * 
     * This class is templated and spezialized for
     * Bunch() or Beam() being observed.
     * 
     * \param[in] obj The object (bunch or beam) to be observed.
     * \param[in] filename The name of the file to write.
     * \param[in] the logging frequency in tracking steps
     */
    TrajectoryLogger(objectT *obj, const char *filename, int step);

    /*! Check whether the given tracking step should be logged */
    bool log_requested(int step) { return 0 == (step % stepsize); };
    
    /*! The source has advanced one time step.
     *  Compute and store the quantities of interest.
     */
    void update();

    /*! Write the collected data into an HDF5 file.
     * 
     */
    int WriteBeamParametersSDDS();

private:

    //! the observed beam object
    objectT *Beam;
    
    //! the file name
    const char* fn;

    //! the number of stored datasets
    unsigned int NOTS;
    
    //! the logging frequency in tracking steps
    unsigned int stepsize;
    
};
