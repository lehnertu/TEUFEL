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

#include "pugixml.hpp"
#include "vector.h"
#include "fields.h"
#include "beam.h"
#include "parser.h"
#include <vector>

/*!
 * \class FEL_1D
 * \brief An electromagnetic field interacting with the beam
 * @author Ulf Lehnert
 * @date 27.3.2023
 * 
 * This is an electromagnetic field described on a 1D grid.
 * The grid propagates with the wave at the speed of light.
 * The transverse extent of the field is defined by a prescribed optical mode definition.
 * One polarization direction perpendicular to the direction of propagation is described.
 * 
 * The index 0 grid cell is placed at the given position in space and time.
 * The grid extends spatially opposite to the direction of propagation, so,
 * all grid cells subsequentially pass through the origin point.
 * The cell size must be the time-step divided by the speed of light for
 * the diffraction algorithm to work properly.
 */
class FEL_1D : public InteractionField
{

public:

    /*! The default contructor just calls the default constructor of the base class
     *  and initalizes all variables with sane values. This will not yet produce any fields.
     */
    FEL_1D(
        double time_step,
        int field_size
        );

    /*! This constructor takes the information from an XML node
     *  describing all field properties. It will throw exceptions
     *  if necessary information is missing or cannot be interpreted.
     * 
     *  A reference to the input parser must be provided as it is
     *  necessary to run the input through the calculator.
     */
    FEL_1D(
        double time_step,
        const pugi::xml_node node,
        InputParser *parser
        );

    /*! All derived classes from GeneralField must provide a destructor */
    virtual ~FEL_1D();

    /*! Advance the electromagnetic field one time step.
     *  The particles of the given beam drive the interaction.
     *  This method must be called in a leap-frog sequence interleaved with
     *  the tracking steps of the beam.
     */
    virtual void step(Beam *beam);
    
    /*! Write logging data to file if requested (do nothing otherwise).
     *  This will be called after the tracking is finished.
     * 
     * This must be implemented for all derived interactions.
     */
    virtual void write_output();

    /*!
     * The electromagnetic field at a given time and point in space.
     * 
     * The coordinates [m] and the time [s] refer to the laboratory (rest) frame.
     * The field is returned as a tuple of electric field [V/m] and
     * magnetic field [T] vectors.
     */
    virtual ElMagField Field(double t, Vector X);

private:

    /*! Method called by all constructors for initialization and printout
     *
     *  Here we define the envelope of the electromagnetic field.
     *  The transverse properties are assumed to be constant along the stored field.
     *  This is justified by the assumption that the propagates field is short
     *  compared to the propagation length. The transverse profile does, however,
     *  vary along the propagation (all cells synchronuously).
     *
     *  The field assumes a Gaussian optical mode. The central wavelength
     *  is not known, so we declare the mode by specifying the Rayleigh range
     *  and the position and size of the optical waist (this implicitly defines a wavelength).
     *
     *  \param ZR the Rayleigh length [m]
     *  \param w0 the waist size [m]
     *  \param z_wa [m] the distance of the waist from the origin in propagation direction
     */
    void setup(double zR, double w0, double z_w);

    /*! Seeding the field with a Gaussian wave packet:
     *  \param E0 peak amplitude of the electric field [V/m]
     *  \param lambda central wavelength of the radiation [m]
     *  \param tau duration of the pulse [s]
     *  \param t_start the distance from the field head to the pulse center [s]
     */
    void seed(double E0, double lambda, double tau, double t_start);
    
    //! number of the field grid cells
    int  N_field;
    
    //! time step [s]
    double  dt;
    
    //! the time of the current simulation step [s]
    double  time;
    
    //! the number of steps that the field has been propagated
    int N_steps;
    
    //! the grid vector [m] (opposite to propagation direction)
    Vector  dz;
    
    //! the current position [m] of the first cell
    Vector  head;
    
    //! the position [m] of the first cell at start time
    Vector  origin;
    
    //! The propagation direction (unit vector) of the electromagnetic wave
    Vector prop;

    //! The propagation step size [m] of the electromagnetic wave
    double ds;
    
    /*! The polarization direction of the electric field.
     *  This vector is forced to be perpendicular to the propagation direction.
     */
    Vector e_x;
    
    /*! The polarization direction of the magnetic field.
     *  This vector is generated perpendicular to the propagation direction and the electric field.
     */
    Vector e_y;
    
    //! the Rayleight length [m] of the optical mode
    double z_Rayleigh;
    
    //! the waist size [m] of the optical mode
    double w_0;
    
    //! the distance of the waist from the origin
    double z_waist;
    
    //! the electric field [V/m] array (current values)
    std::vector<double> field_E;

    //! the electric field [V/m] array from the last time step
    std::vector<double> previous_E;

    //! storage of electric fields over the computed time-steps
    std::vector<std::vector<double>> field_storage;
    
    //! whether to create a field output file
    bool createOutput;
    
    //! how often to create a field output
    int step_Output;

    //! the output file name
    std::string FileName;

};

