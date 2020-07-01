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
 * \class TrackingLogger
 * \brief Store beam parameters while tracking.
 * @author Ulf Lehnert
 * @date 19.10.2017
 * 
 * This class handles the computation and storage of beam parameters
 * accessible during the tracking.
 * @todo specialization for beam still missing, only for bunch implemented
 */
template <class objectT>
class TrackingLogger
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
     */
    TrackingLogger(objectT *obj, const char *filename);

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
     * - RMS beam size
     * - RMS beam momentum
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
    int WriteBeamParametersSDDS(const char *filename);

private:

    //! the observed beam object
    objectT *Beam;

    //! the number of stored datasets
    unsigned int NOTS;
    
    //! the observation time
    std::vector<double> Time;

    //! the observation position
    std::vector<Vector> Pos;

    //! the observation momentum
    std::vector<Vector> BG;
    
    //! the observation position spread
    std::vector<Vector> PosRMS;
    
    //! the observation momentum spread
    std::vector<Vector> BGRMS;

    //! the observation position-momentum correlation
    std::vector<Vector> PosBG;
    
};
