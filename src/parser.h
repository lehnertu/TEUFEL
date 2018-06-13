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
#include "beam.h"
#include "muParser.h"

using namespace std;

/*!
 *    \class InputParser
 *    \brief Parse the XML input file.
 * 
 *    @author Ulf Lehnert
 *    @date 14.9.2017
 *    
 *    All input about the computation case is read from an XML input file.
 *    This class provides methods for parsing that file and
 *    creation of all objects needed to run the simulation.
 */
class InputParser
{

public:
    
    /*! The constructor of the simulation gets the top node of
     *  an open XML file. This will be the root node from which
     *  all necessary information will be parsed.
     */
    InputParser(const pugi::xml_node node);
    
    //! Destructor only used to free the memory
    ~InputParser();
    
    /*! The children of the given node are checked, whether they are <calc> nodes.
     *  In that case the provided statements are fed to the calculator.
     * 
     *  Possible combinations of attributes are:
     * 
     *  var="name" eq="..."
     *  A variable with given name is initialized with the given equation.
     * 
     *  print="comment" eq="..."
     *  The result of the given equation is printed along with the comment.
     *  
     *  All defined variables are persistent for the duration of the run
     *  and can be used in further equations.
     */
    void parseCalc(const pugi::xml_node node);
    
    /*! Parse the input file and create all described lattice elements.
     *  Each one is added to the given lattice object.
     *  If there is no child <lattice> under the root node
     *  the whole program is aborted.
     * 
     *  @return The number of elements added to the lattice.
     */
    int parseLattice(Lattice *lattice);

    /*! Parse the input file and create the beam.
     *  If there is no child <beam> under the root node
     *  the whole program is aborted.
     * 
     *  Every entry in the beam section creates a single bunch.
     *  Also for each <particle> entry an individual bunch is created.
     * 
     *  @return The number of bunches is returned.
     */
    int parseBeam(Beam *beam);
    
    /* Parse the input file and create the defined observers.
     * Will not cause problems, if no observers are defined (tracking only).
     */
    void parseObservations();

private:
    
    //! the root node of the input file
    pugi::xml_node root;
    
    //! the math equation parser
    mu::Parser *calc;
};

