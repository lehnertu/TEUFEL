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

#include "pugixml.hpp"
#include "beam.h"
#include "muParser.h"
#include "beam.h"
#include "fields.h"
#include "observer.h"

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
 * 
 *    The class incorporates a calculator based on the [muParser](http://beltoforion.de/article.php?a=muparser) library.
 *    All floating point input values are evaluated by this calculator.
 * 
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
    
    /*! The provided node is a <calc> node.
     *  The provided statements are fed to the calculator.
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
     * 
     *  A number of constants is pre-defined for use in calculations
     *  "_c", "_e", "_mec2", "_eps0", "_mu0"
     */
    void parseCalc(const pugi::xml_node node);
    
    /*! The children of the given node are checked, whether they are <calc> nodes.
     *  In that case the provided statements are fed to the calculator
     *  by calling parseCalc().
     */
    void parseCalcChildren(const pugi::xml_node node);

     /*! Convert a node attribute into a double value using the calculator.
     *  In case of an error 0.0 is returned.
     */
    double parseValue(const pugi::xml_attribute att);
    
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
    
    /*! Parse the input file and setup for tracking particles.
     *  If there is no child <tracking> under the root node
     *  the whole program is aborted.
     * 
     *  The parameters for tracking are stored in the beam object.
     * 
     *  Watch points where beam snapshots are stored during tracking
     *  are defined in this section. The definitions are appended to the
     *  list provided.
     * 
     *  @return The number of watch points
     */
    int parseTracking(Beam *beam, std::vector<watch_t> *watches);
    
    /*! Parse the input file and create the defined observers.
     *  All observers found are appended to the given list.
     *  Will not cause problems, if no observers are defined (tracking only)
     *  but the main nodes <observer> should still be present.
     * 
     *  For construction of the observer the observed beam object must be given.
     * 
     *  @return The number of observers is returned.
     */
    int parseObservers(std::vector<Observer*> *listObservers, Beam* beam);

private:
    
    //! the root node of the input file
    pugi::xml_node root;
    
    //! the math equation parser
    mu::Parser *calc;
};

