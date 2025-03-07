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

#include <iostream>
#include <iomanip>
#include <math.h>

#include "dipole.h"
#include "fel1d.h"
#include "fields.h"
#include "global.h"
#include "mesh.h"
#include "parser.h"
#include "point_observer.h"
#include "screen.h"
#include "snapshot.h"
#include "source_screen.h"
#include "undulator.h"
#include "vector.h"
#include "wave.h"

InputParser::InputParser(const pugi::xml_node node)
{
    // this is the <teufel> node
    root = node;
    calc = new mu::Parser();
    // we pre-define a number of constants in the calculator
    calc->DefineConst("_c", (double)SpeedOfLight);
    calc->DefineConst("_e", (double)ElementaryCharge);
    calc->DefineConst("_h", (double)PlanckH);
    calc->DefineConst("_mec2", (double)mecsquared);
    calc->DefineConst("_eps0", (double)EpsNull);
    calc->DefineConst("_pi", (double)Pi);
    // evaluate all <calc> nodes present as children of root
    parseCalcChildren(root);
}

InputParser::~InputParser()
{
    delete calc;
}

void InputParser::parseCalc(const pugi::xml_node node)
{
    std::string prt = node.attribute("print").value();
    std::string var = node.attribute("var").value();
    std::string eq = node.attribute("eq").value();
    // if this is a print statement
    if (prt.length() > 0)
    {
        try
        {
            calc->SetExpr(eq);
            if (teufel::rank==0)
                std::cout << "calc : " << prt << calc->Eval() << std::endl;
        }
        catch (mu::Parser::exception_type &e)
        {
            if (teufel::rank==0)
            {
                std::cout << "InputParser::parseCalc : can't evaluate " << eq << std::endl;
                std::cout << e.GetMsg() << endl;
            }
        }
    };
    // if this is a variable definition
    if (var.length() > 0)
    {
        try
        {
            calc->SetExpr(eq);
            calc->DefineConst(var, calc->Eval());
        }
        catch (mu::Parser::exception_type &e)
        {
            if (teufel::rank==0)
            {
                std::cout << "InputParser::parseCalc : can't evaluate " << eq << std::endl;
                std::cout << e.GetMsg() << endl;
            }
        }
    };
}

void InputParser::parseCalcChildren(const pugi::xml_node node)
{
    for (pugi::xml_node_iterator it = node.begin(); it != node.end(); ++it)
    {
        pugi::xml_node element = *it;
        std::string type = element.name();
        if ( type == "calc") parseCalc(element);
    }
}

double InputParser::parseDouble(const pugi::xml_attribute attr)
{
    double retval = 0.0;
    try
    {
        calc->SetExpr(attr.value());
        retval = calc->Eval();
    }
    catch (mu::Parser::exception_type &e)
    {
        if (teufel::rank==0)
        {
            std::cout << "InputParser::parseDouble : can't evaluate " << attr.name() << " from " << attr.value() << std::endl;
            std::cout << e.GetMsg() << endl;
        }
    }
    return retval;
}

//! @todo use parseInt for all integer inputs
double InputParser::parseInt(const pugi::xml_attribute attr)
{
    int retval = 0;
    try
    {
        calc->SetExpr(attr.value());
        retval = rint(calc->Eval());
    }
    catch (mu::Parser::exception_type &e)
    {
        if (teufel::rank==0)
        {
            std::cout << "InputParser::parseInt : can't evaluate " << attr.name() << " from " << attr.value() << std::endl;
            std::cout << e.GetMsg() << endl;
        }
    }
    return retval;
}

int InputParser::parseLattice(Lattice *lattice)
{
    int count = 0;
    std::string type;
    std::string name;
    pugi::xml_node latticenode = root.child("lattice");
    if (!latticenode)
        throw(IOexception("InputParser::parseLattice(Lattice - section <lattice> not found in input."));
    // loop over all children of the lattice node
    for (pugi::xml_node_iterator it = latticenode.begin(); it != latticenode.end(); ++it)
    {
        pugi::xml_node element = *it;
        type = element.name();
        if (type == "calc")
            parseCalc(element);
        else if (type == "undulator")
        {
            parseCalcChildren(element);
            type = element.attribute("type").value();
            name = element.attribute("name").value();
            if (type == "planar")
            {
                if (teufel::rank==0) std::cout << name << "::PlanarUndulator" << std::endl;
                // the undulator object parses its own input
                // we provide a reference to the parser
                PlanarUndulator* Undu = new PlanarUndulator(element, this);
                lattice->addElement(Undu);
            }
            else if (type == "transverse gradient")
            {
                if (teufel::rank==0) std::cout << name << "::TransverseGradientUndulator" << std::endl;
                // the undulator object parses its own input
                // we provide a reference to the parser
                TransverseGradientUndulator* Undu = new TransverseGradientUndulator(element, this);
                lattice->addElement(Undu);
            }
            else
                throw(IOexception("InputParser::parseLattice(Lattice - unknown undulator type."));
            count++;
        }
        else if (type == "wave")
        {
            parseCalcChildren(element);
            type = element.attribute("type").value();
            name = element.attribute("name").value();
            if (type == "gaussian")
            {
                if (teufel::rank==0) std::cout << name << "::GaussianWave" << std::endl;
                // the wave object parses its own input
                // we provide a reference to the parser
                GaussianWave *wave = new GaussianWave(element, this);
                lattice->addElement(wave);
            }
            else if (type == "packet")
            {
                if (teufel::rank==0) std::cout << name << "::GaussianWavePacket" << std::endl;
                // the wave packet object parses its own input
                // we provide a reference to the parser
                GaussianWavePacket *wave = new GaussianWavePacket(element, this);
                lattice->addElement(wave);
            }
            else
                throw(IOexception("InputParser::parseLattice(Lattice - unknown wave type."));
            count++;
        }
        else if (type == "screen")
        {
            parseCalcChildren(element);
            type = element.attribute("type").value();
            name = element.attribute("name").value();
            if (teufel::rank==0) std::cout << name << "::SourceScreen" << std::endl;
            pugi::xml_attribute fn = element.attribute("file");
            if (!fn)
                throw(IOexception("InputParser::parseLattice - <screen> attribute file not found."));                
            double x, y, z, t0;
            pugi::xml_node posnode = element.child("position");
            if (!posnode)
                throw(IOexception("InputParser::parseLattice - <screen> <position> not found."));
            x = parseDouble(posnode.attribute("x"));
            y = parseDouble(posnode.attribute("y"));
            z = parseDouble(posnode.attribute("z"));
            t0 = parseDouble(posnode.attribute("t"));
            Vector pos = Vector(x, y, z);
            SourceScreen *screen = new SourceScreen(
                fn.as_string(), pos, t0);
            lattice->addElement(screen);
            count++;
        }
        else if (type == "dipole")
        {
            parseCalcChildren(element);
            type = element.attribute("type").value();
            name = element.attribute("name").value();
            if (type == "hard edge")
            {
                if (teufel::rank==0) std::cout << name << "::HardEdgeDipole" << std::endl;
                // the dipole object parses its own input
                // we provide a reference to the parser
                HardEdgeDipole* dipole = new HardEdgeDipole(element, this);
                lattice->addElement(dipole);
            }
            else if (type == "soft edge")
            {
                if (teufel::rank==0) std::cout << name << "::SoftEdgeDipole" << std::endl;
                // the dipole object parses its own input
                // we provide a reference to the parser
                SoftEdgeDipole* dipole = new SoftEdgeDipole(element, this);
                lattice->addElement(dipole);
            }
            else
                throw(IOexception("InputParser::parseLattice(Lattice - unknown dipole type."));
            count++;
        }
        else if (type == "background")
        {
            parseCalcChildren(element);
            double x, y, z;
            pugi::xml_node Enode = element.child("E");
            Vector E;
            if (Enode)
            {
                x = parseDouble(Enode.attribute("x"));
                y = parseDouble(Enode.attribute("y"));
                z = parseDouble(Enode.attribute("z"));
                E = Vector(x, y, z);
            }
            pugi::xml_node Bnode = element.child("B");
            Vector B;
            if (Bnode)
            {
                x = parseDouble(Bnode.attribute("x"));
                y = parseDouble(Bnode.attribute("y"));
                z = parseDouble(Bnode.attribute("z"));
                B = Vector(x, y, z);
            }
            HomogeneousField *background = new HomogeneousField(E,B);
            lattice->addElement(background);
            count++;
        }
        else
            throw(IOexception("InputParser::parseLattice(Lattice - unknown lattice element."));
    }
    return count;
}

int InputParser::parseBeam(Beam *beam, std::vector<Logger*> *logs, ProbeInfo *probe)
{
    int count = 0;
    pugi::xml_node beamnode = root.child("beam");
    if (!beamnode)
    {
        throw(IOexception("InputParser::parseBeam - section <beam> not found in input."));
    }
    else
    {
        for (pugi::xml_node_iterator it = beamnode.begin(); it != beamnode.end(); ++it)
        {
            pugi::xml_node entry = *it;
            std::string type = entry.name();
            if (type == "calc")
                parseCalc(entry);
            else if (type == "particle")
            {
                double gamma = parseDouble(entry.attribute("gamma"));
                double betagamma = sqrt(gamma * gamma - 1.0);
                double charge = parseDouble(entry.attribute("charge"))/ElementaryCharge;
                double cmr = parseDouble(entry.attribute("cmr"));
                parseCalcChildren(entry);
                double t0;
                pugi::xml_node timenode = entry.child("time");
                if (!timenode)
                    t0 = 0.0;
                else
                    t0 = parseDouble(timenode.attribute("t0"));
                pugi::xml_node posnode = entry.child("position");
                if (!posnode)
                    throw(IOexception("InputParser::parseBeam - <particle> <position> not found."));
                double x, y, z;
                x = parseDouble(posnode.attribute("x"));
                y = parseDouble(posnode.attribute("y"));
                z = parseDouble(posnode.attribute("z"));
                Vector pos = Vector(x, y, z);
                pugi::xml_node dirnode = entry.child("momentum");
                if (!dirnode)
                    throw(IOexception("InputParser::parseBeam - <particle> <momentum> not found."));
                x = parseDouble(dirnode.attribute("x"));
                y = parseDouble(dirnode.attribute("y"));
                z = parseDouble(dirnode.attribute("z"));
                Vector mom = Vector(x, y, z);
                mom.normalize();
                mom *= betagamma;
                Vector acc = Vector(0.0, 0.0, 0.0);
                // now we have all information - create a particle
                ChargedParticle *p = new ChargedParticle(charge,charge/cmr);
                p->initTrajectory(t0, pos, mom, acc);
                // create a bunch with just this particle
                Bunch *single = new Bunch();
                single->Add(p);
                // create a logger for this bunch if defined in the input file
                pugi::xml_node lognode = entry.child("log");
                if (lognode)
                {
                    pugi::xml_attribute fn = lognode.attribute("file");
                    if (!fn) throw(IOexception("InputParser::parseBeam - <particle> filename for log not found."));
                    int step = 1;
                    pugi::xml_attribute st = lognode.attribute("step");
                    if (st) step = st.as_int();
                    ParameterLogger<Bunch>* logger = new ParameterLogger<Bunch>(single, fn.as_string(), step);
                    logs->push_back(logger);
                }
                // add this bunch to the beam
                beam->Add(single);
                count++;
            }
            else if (type == "bunch")
            {
                double gamma = parseDouble(entry.attribute("gamma"));
                double betagamma = sqrt(gamma * gamma - 1.0);
                double charge = parseDouble(entry.attribute("charge"))/ElementaryCharge;
                double cmr = parseDouble(entry.attribute("cmr"));
                int NoP = entry.attribute("n").as_int(1);
                parseCalcChildren(entry);
                double t0 = 0.0;
                double x, y, z;
                double xrms = 0.0;
                double yrms = 0.0;
                double zrms = 0.0;
                double zft = 0.0;
                // parse start time
                pugi::xml_node timenode = entry.child("time");
                if (timenode)
                {
                    parseCalcChildren(timenode);
                    t0 = parseDouble(timenode.attribute("t0"));
                }
                // parse position
                pugi::xml_node posnode = entry.child("position");
                if (!posnode)
                    throw(IOexception("InputParser::parseBeam - <bunch> <position> not found."));
                parseCalcChildren(posnode);
                // parse mean position
                x = parseDouble(posnode.attribute("x"));
                y = parseDouble(posnode.attribute("y"));
                z = parseDouble(posnode.attribute("z"));
                Vector pos = Vector(x, y, z);
                // parse position distribution
                pugi::xml_attribute xrmsatt = posnode.attribute("xrms");
                if (xrmsatt) xrms=parseDouble(xrmsatt);
                pugi::xml_attribute yrmsatt = posnode.attribute("yrms");
                if (yrmsatt) yrms=parseDouble(yrmsatt);
                pugi::xml_attribute zrmsatt = posnode.attribute("zrms");
                if (zrmsatt) zrms=parseDouble(zrmsatt);
                pugi::xml_attribute zftatt = posnode.attribute("zft");
                if (zftatt) zft=parseDouble(zftatt);
                // parse momentum
                double xprms = 0.0;
                double yprms = 0.0;
                double delta = 0.0;
                pugi::xml_node momnode = entry.child("momentum");
                if (!momnode)
                    throw(IOexception("InputParser::parseBeam - <bunch> <momentum> not found."));
                parseCalcChildren(momnode);
                // parse mean direction of momentum
                x = parseDouble(momnode.attribute("x"));
                y = parseDouble(momnode.attribute("y"));
                z = parseDouble(momnode.attribute("z"));
                Vector dir = Vector(x, y, z);
                dir.normalize();
                // parse momentum distribution
                pugi::xml_attribute xpatt = momnode.attribute("xrms");
                if (xpatt) xprms=parseDouble(xpatt);
                pugi::xml_attribute ypatt = momnode.attribute("yrms");
                if (ypatt) yprms=parseDouble(ypatt);
                pugi::xml_attribute delatt = momnode.attribute("delta");
                if (delatt) delta=parseDouble(delatt);
                // create the bunch particle distribution
                Distribution *dist = new Distribution(6, NoP);
                dist->generateGaussian(0, xrms);
                dist->generateGaussian(1, yrms);
                if (zft != 0.0)
                    dist->scale(2, zft);
                else
                    {
                        if (zrms != 0.0) dist->generateGaussian(2, zrms);
                        else dist->scale(2, 0.0);
                    };
                dist->generateGaussian(3, betagamma*xprms);
                dist->generateGaussian(4, betagamma*yprms);
                dist->generateGaussian(5, betagamma*delta);
                // parse and add correlations
                pugi::xml_node corrnode = entry.child("correlations");
                if (corrnode)
                {
                    parseCalcChildren(corrnode);
                    // all children are handled in sequence
                    for (pugi::xml_node_iterator cit = corrnode.begin(); cit != corrnode.end(); ++cit)
                    {
                        pugi::xml_node cn = *cit;
                        std::string ctype = cn.name();
                        if (ctype == "linear")
                        {
                            int ind = 0;
                            int dep = 0;
                            double fac = 0.0;
                            pugi::xml_attribute indatt = cn.attribute("indep");
                            if (indatt)
                                ind = int(parseDouble(indatt));
                            else
                                throw(IOexception("InputParser::parseBeam - correlation indep not found."));
                            pugi::xml_attribute depatt = cn.attribute("dep");
                            if (depatt)
                                dep = int(parseDouble(depatt));
                            else
                                throw(IOexception("InputParser::parseBeam - correlation dep not found."));
                            pugi::xml_attribute facatt = cn.attribute("factor");
                            if (facatt)
                                fac = parseDouble(facatt);
                            else
                                throw(IOexception("InputParser::parseBeam - correlation factor not found."));
                            // std::cout << "correlation " << fac << " " << ind << " -> " << dep << " added" << endl;
                            dist->addCorrelation(ind, dep, fac);
                        } else {
                            throw(IOexception("InputParser::parseBeam - unknown type of correlation."));
                        }
                    };
                };
                // now we can create particles according to the coordinate distribution
                Bunch *bunch = new Bunch(dist, t0, pos, dir*betagamma, charge/(double)NoP, charge/cmr/(double)NoP);
                // create a logger for this bunch if defined in the input file
                pugi::xml_node lognode = entry.child("log");
                if (lognode)
                {
                    pugi::xml_attribute fn = lognode.attribute("file");
                    if (!fn) throw(IOexception("InputParser::parseBeam - <bunch> filename for log not found."));
                    int step = 1;
                    //! @todo if step is not given it becomes 0 and leads to divide-by-zero errors
                    pugi::xml_attribute st = lognode.attribute("step");
                    if (st) step = st.as_int();
                    ParameterLogger<Bunch>* logger = new ParameterLogger<Bunch>(bunch, fn.as_string(), step);
                    pugi::xml_attribute bf = lognode.attribute("bunching-freq");
                    if (bf)
                    {
                        double freq = parseDouble(bf);
                        logger->record_bunching(freq);
                    };
                    logs->push_back(logger);
                };
                // create a particle trajectory logger for this bunch if defined in the input file
                pugi::xml_node trajnode = entry.child("trajectories");
                if (trajnode)
                {
                    pugi::xml_attribute fn = trajnode.attribute("file");
                    if (!fn) throw(IOexception("InputParser::parseBeam - <bunch> filename for trajectories not found."));
                    int step = 1;
                    pugi::xml_attribute st = trajnode.attribute("step");
                    if (st) step = st.as_int();
                    int particles = 1000000;
                    pugi::xml_attribute np = trajnode.attribute("particles");
                    if (np) particles = np.as_int();
                    TrajectoryLogger<Bunch>* traj_logger = new TrajectoryLogger<Bunch>(bunch, fn.as_string(), step, particles);
                    logs->push_back(traj_logger);
                };
                // add the bunch to the beam
                beam->Add(bunch);
                delete dist;
                count++;
            }
            else if (type == "SDDS")
            {
                pugi::xml_attribute sddsname = entry.attribute("file");
                if (!sddsname)
                    throw(IOexception("InputParser::parseBeam - <SDDS> attribute file not found."));   
                // parse propagation  direction of the input beam
                Vector dir = Vector(0.0, 0.0, 1.0);
                pugi::xml_node dirnode = entry.child("dir");
                if (dirnode)
                {
                    parseCalcChildren(dirnode);
                    double x = parseDouble(dirnode.attribute("x"));
                    double y = parseDouble(dirnode.attribute("y"));
                    double z = parseDouble(dirnode.attribute("z"));
                    dir = Vector(x, y, z);
                    dir.normalize();
                };
                // now we can create the bunch from the SDDS input file
                Bunch *bunch = new Bunch(sddsname.as_string(), dir);
                // parse and add correlations
                pugi::xml_node corrnode = entry.child("correlations");
                if (corrnode)
                {
                    parseCalcChildren(corrnode);
                    // all children are handled in sequence
                    for (pugi::xml_node_iterator cit = corrnode.begin(); cit != corrnode.end(); ++cit)
                    {
                        pugi::xml_node cn = *cit;
                        std::string ctype = cn.name();
                        if (ctype == "linear")
                        {
                            int ind = 0;
                            int dep = 0;
                            double fac = 0.0;
                            pugi::xml_attribute indatt = cn.attribute("indep");
                            if (indatt)
                                ind = int(parseDouble(indatt));
                            else
                                throw(IOexception("InputParser::parseBeam <SDDS> - correlation indep not found."));
                            pugi::xml_attribute depatt = cn.attribute("dep");
                            if (depatt)
                                dep = int(parseDouble(depatt));
                            else
                                throw(IOexception("InputParser::parseBeam <SDDS> - correlation dep not found."));
                            pugi::xml_attribute facatt = cn.attribute("factor");
                            if (facatt)
                                fac = parseDouble(facatt);
                            else
                                throw(IOexception("InputParser::parseBeam <SDDS> - correlation factor not found."));
                            // std::cout << "correlation " << fac << " " << ind << " -> " << dep << " added" << endl;
                            bunch->addCorrelation(ind, dep, fac);
                        } else {
                            throw(IOexception("InputParser::parseBeam <SDDS> - unknown type of correlation."));
                        }
                    };
                };
                // create a logger for this bunch if defined in the input file
                pugi::xml_node lognode = entry.child("log");
                if (lognode)
                {
                    pugi::xml_attribute fn = lognode.attribute("file");
                    if (!fn) throw(IOexception("InputParser::parseBeam - <bunch> filename for log not found."));
                    int step = 1;
                    pugi::xml_attribute st = lognode.attribute("step");
                    if (st) step = st.as_int();
                    ParameterLogger<Bunch>* logger = new ParameterLogger<Bunch>(bunch, fn.as_string(), step);
                    pugi::xml_attribute bf = lognode.attribute("bunching-freq");
                    if (bf)
                    {
                        double freq = parseDouble(bf);
                        logger->record_bunching(freq);
                    };
                    logs->push_back(logger);
                };
                // add the bunch to the beam
                beam->Add(bunch);
                count++;
            }
            else if (type == "probe")
            {
                probe->requested = true;
                pugi::xml_attribute filename = entry.attribute("file");
                if (!filename)
                    throw(IOexception("InputParser::parseBeam - <probe> attribute file not found."));
                probe->filename = filename.as_string();
                probe->step = 1;
                pugi::xml_attribute st = entry.attribute("step");
                if (st) probe->step = st.as_int();
                probe->number = 1;
                pugi::xml_attribute nm = entry.attribute("number");
                if (nm) probe->number = nm.as_int();
            }
            else
                throw(IOexception("InputParser::parseBeam - unknown type of beam entry."));
        };
    };
    return count;
}

void InputParser::parseTracking(
    Beam *beam,
    std::vector<watch_t> *watches,
    std::list<InteractionField*> *interactions
    )
{
    pugi::xml_node track = root.child("tracking");
    if (!track)
    {
        throw(IOexception("InputParser::parseTracking - section <tracking> not found in input."));
    }
    else
    {
        std::string method = track.attribute("method").value();
        if (method == "Vay")
        {
            beam->setTrackingMethod(TRACKING_VAY);
        }
        else
            throw(IOexception("InputParser::parseTracking - unknown tracking method."));
        pugi::xml_attribute timestep = track.attribute("delta_t");
        if (!timestep)
            throw(IOexception("InputParser::parseTracking - <tracking> attribute delta_t not found."));
        beam->setTimeStep(parseDouble(timestep));
        pugi::xml_attribute nst = track.attribute("n");
        if (!nst)
            throw(IOexception("InputParser::parseTracking - <tracking> attribute n not found."));
        beam->setNOTS(parseInt(nst));
        // parse all children of the <tracking> node
        for (pugi::xml_node_iterator it = track.begin(); it != track.end(); ++it)
        {
            pugi::xml_node child = *it;
            std::string type = child.name();
            if (type == "calc")
                parseCalc(child);
            else if (type == "fel1d")
            {
                parseCalcChildren(child);
                // the interaction object parses its own input, we provide a reference to the parser
                // the timestep attribute has already been evaluated - thats necessary
                FEL_1D* fel_interaction = new FEL_1D(parseDouble(timestep), child, this);
                interactions->push_back(fel_interaction);
            }
            else if (type == "watch")
            {
                pugi::xml_attribute step = child.attribute("step");
                if (!step)
                    throw(IOexception("InputParser::parseTracking - <watch> attribute step not found."));
                int n_step = int(parseInt(step));
                pugi::xml_attribute fn = child.attribute("file");
                if (!fn)
                    throw(IOexception("InputParser::parseTracking - <watch> attribute file not found."));
                watch_t w = { n_step, fn.as_string(), hdf5 };
                watches->push_back(w);
            }
            else if (type == "watch-sdds")
            {
                pugi::xml_attribute step = child.attribute("step");
                if (!step)
                    throw(IOexception("InputParser::parseTracking - <watch-sdds> attribute step not found."));
                int n_step = int(parseInt(step));
                pugi::xml_attribute fn = child.attribute("file");
                if (!fn)
                    throw(IOexception("InputParser::parseTracking - <watch-sdds> attribute file not found."));
                watch_t w = { n_step, fn.as_string(), sdds };
                watches->push_back(w);
            }
            else
                throw(IOexception("InputParser::parseTracking - unknown child under <tracking>."));
        }
    }
}

void InputParser::parseObservers(std::vector<BaseObserver*> *listObservers)
{
    pugi::xml_node obsroot = root.child("observer");
    if (!obsroot)
    {
        throw(IOexception("InputParser::parseObservers - section <observer> not found in input."));
    }
    else
    {
        for (pugi::xml_node_iterator it = obsroot.begin(); it != obsroot.end(); ++it)
        {
            pugi::xml_node obs = *it;
            std::string type = obs.name();
            if (type == "calc")
                parseCalc(obs);
            else if (type == "snapshot")
            {
                double x, y, z;
                pugi::xml_attribute fn = obs.attribute("file");
                if (!fn)
                    throw(IOexception("InputParser::parseObservers - <snapshot> attribute file not found."));                
                parseCalcChildren(obs);
                pugi::xml_node tnode = obs.child("time");
                if (!tnode)
                    throw(IOexception("InputParser::parseObservers - <snapshot> <time> not found."));
                double t0 = parseDouble(tnode.attribute("t0"));
                pugi::xml_node posnode = obs.child("position");
                if (!posnode)
                    throw(IOexception("InputParser::parseObservers - <snapshot> <position> not found."));
                x = parseDouble(posnode.attribute("x"));
                y = parseDouble(posnode.attribute("y"));
                z = parseDouble(posnode.attribute("z"));
                Vector pos = Vector(x, y, z);
                pugi::xml_node hnode = obs.child("hpitch");
                if (!hnode)
                    throw(IOexception("InputParser::parseObservers - <snapshot> <hpitch> not found."));
                x = parseDouble(hnode.attribute("x"));
                y = parseDouble(hnode.attribute("y"));
                z = parseDouble(hnode.attribute("z"));
                Vector h = Vector(x, y, z);
                pugi::xml_attribute nh = hnode.attribute("n");
                if (!nh)
                    throw(IOexception("InputParser::parseObservers - <snapshot> <hpitch> attribute n not found."));
                pugi::xml_node vnode = obs.child("vpitch");
                if (!vnode)
                    throw(IOexception("InputParser::parseObservers - <snapshot> <vpitch> not found."));
                x = parseDouble(vnode.attribute("x"));
                y = parseDouble(vnode.attribute("y"));
                z = parseDouble(vnode.attribute("z"));
                Vector v = Vector(x, y, z);
                pugi::xml_attribute nv = vnode.attribute("n");
                if (!nv)
                    throw(IOexception("InputParser::parseObservers - <snapshot> <vpitch> attribute n not found."));
                pugi::xml_attribute src = obs.attribute("source");
                BaseObserver *snapObs;
                if (src)
                {
                    if (std::string(src.value()) == "beam")
                    {
                        snapObs = new SnapshotObserver<Beam>(
                            fn.as_string(),
                            pos, h, v,
                            nh.as_int(),
                            nv.as_int(),
                            t0 );
                        snapObs->setSourceType(BeamObservation);
                    }
                    else if (std::string(src.value()) == "lattice")
                    {
                        snapObs = new SnapshotObserver<Lattice>(
                            fn.as_string(),
                            pos, h, v,
                            nh.as_int(),
                            nv.as_int(),
                            t0 );
                        snapObs->setSourceType(LatticeObservation);
                    }
                    else throw(IOexception("InputParser::parseObservers - <snapshot> illegal source attribute."));
                }
                else
                {
                    snapObs = new SnapshotObserver<Beam>(
                        fn.as_string(),
                        pos, h, v,
                        nh.as_int(),
                        nv.as_int(),
                        t0 );
                    snapObs->setSourceType(BeamObservation);
                }
                listObservers->push_back(snapObs);
            }
            else if (type == "screen")
            {
                double x, y, z;
                double t0, dt, dtx, dty;
                pugi::xml_attribute fn = obs.attribute("file");
                if (!fn)
                    throw(IOexception("InputParser::parseObservers - <screen> attribute file not found."));                
                parseCalcChildren(obs);
                pugi::xml_node posnode = obs.child("position");
                if (!posnode)
                    throw(IOexception("InputParser::parseObservers - <screen> <position> not found."));
                x = parseDouble(posnode.attribute("x"));
                y = parseDouble(posnode.attribute("y"));
                z = parseDouble(posnode.attribute("z"));
                Vector pos = Vector(x, y, z);
                pugi::xml_node hnode = obs.child("hpitch");
                if (!hnode)
                    throw(IOexception("InputParser::parseObservers - <screen> <hpitch> not found."));
                x = parseDouble(hnode.attribute("x"));
                y = parseDouble(hnode.attribute("y"));
                z = parseDouble(hnode.attribute("z"));
                pugi::xml_attribute att_dtx = hnode.attribute("dt");
                if (!att_dtx)
                    dtx = 0.0;
                else
                    dtx = parseDouble(att_dtx);
                Vector dx = Vector(x, y, z);
                pugi::xml_attribute nh = hnode.attribute("n");
                if (!nh)
                    throw(IOexception("InputParser::parseObservers - <screen> <hpitch> attribute n not found."));
                pugi::xml_node vnode = obs.child("vpitch");
                if (!vnode)
                    throw(IOexception("InputParser::parseObservers - <screen> <vpitch> not found."));
                x = parseDouble(vnode.attribute("x"));
                y = parseDouble(vnode.attribute("y"));
                z = parseDouble(vnode.attribute("z"));
                pugi::xml_attribute att_dty = vnode.attribute("dt");
                if (!att_dty)
                    dty = 0.0;
                else
                    dty = parseDouble(att_dty);
                Vector dy = Vector(x, y, z);
                pugi::xml_attribute nv = vnode.attribute("n");
                if (!nv)
                    throw(IOexception("InputParser::parseObservers - <screen> <vpitch> attribute n not found."));
                pugi::xml_node tnode = obs.child("time");
                if (!tnode)
                    throw(IOexception("InputParser::parseObservers - <screen> <time> not found."));
                t0 = parseDouble(tnode.attribute("t0"));
                dt = parseDouble(tnode.attribute("dt"));
                pugi::xml_attribute nt = tnode.attribute("n");
                if (!nt)
                    throw(IOexception("InputParser::parseObservers - <screen> <time> attribute n not found."));
                pugi::xml_attribute src = obs.attribute("source");
                BaseObserver *screenObs;
                if (src)
                {
                    if (std::string(src.value()) == "beam")
                    {
                        screenObs = new ScreenObserver<Beam>(
                            fn.as_string(),
                            pos, dx, dy,
                            nh.as_int(),
                            nv.as_int(),
                            t0, dt, dtx, dty, nt.as_int() );
                        screenObs->setSourceType(BeamObservation);
                    }
                    else if (std::string(src.value()) == "lattice")
                    {
                        screenObs = new ScreenObserver<Lattice>(
                            fn.as_string(),
                            pos, dx, dy,
                            nh.as_int(),
                            nv.as_int(),
                            t0, dt, dtx, dty, nt.as_int() );
                        screenObs->setSourceType(LatticeObservation);
                    }
                    else throw(IOexception("InputParser::parseObservers - <screen> illegal source attribute."));
                }
                else
                {
                    screenObs = new ScreenObserver<Beam>(
                        fn.as_string(),
                        pos, dx, dy,
                        nh.as_int(),
                        nv.as_int(),
                        t0, dt, dtx, dty, nt.as_int() );
                    screenObs->setSourceType(BeamObservation);
                }
                listObservers->push_back(screenObs);
            }
            else if (type == "mesh")
            {
                pugi::xml_attribute fn = obs.attribute("file");
                if (!fn)
                    throw(IOexception("InputParser::parseObservers - <mesh> attribute file not found."));                
                parseCalcChildren(obs);
                pugi::xml_attribute src = obs.attribute("source");
                if (src)
                {
                    if (std::string(src.value()) == "beam")
                    {
                        MeshedScreen<Beam> *obs = new MeshedScreen<Beam>(fn.as_string());
                        obs->init();
                        obs->zero();
                        obs->setSourceType(BeamObservation);
                        listObservers->push_back(obs);
                        if (teufel::rank==0) obs->writeReport(&cout); 
                    }
                    else if (std::string(src.value()) == "lattice")
                    {
                        MeshedScreen<Lattice> *obs = new MeshedScreen<Lattice>(fn.as_string());
                        obs->init();
                        obs->zero();
                        obs->setSourceType(LatticeObservation);
                        listObservers->push_back(obs);
                        if (teufel::rank==0) obs->writeReport(&cout);
                    }
                    else throw(IOexception("InputParser::parseObservers - <mesh> illegal source attribute."));
                }
                else
                {
                    MeshedScreen<Beam> *obs = new MeshedScreen<Beam>(fn.as_string());
                    obs->init();
                    obs->zero();
                    obs->setSourceType(BeamObservation);
                    listObservers->push_back(obs);
                    if (teufel::rank==0) obs->writeReport(&cout);
                }
            }
            else if (type == "point")
            {
                double x, y, z;
                double t0, dt;
                pugi::xml_attribute fn = obs.attribute("file");
                if (!fn)
                    throw(IOexception("InputParser::parseObservers - <point> attribute file not found."));
                parseCalcChildren(obs);
                pugi::xml_node posnode = obs.child("position");
                if (!posnode)
                    throw(IOexception("InputParser::parseObservers - <point> <position> not found."));
                x = parseDouble(posnode.attribute("x"));
                y = parseDouble(posnode.attribute("y"));
                z = parseDouble(posnode.attribute("z"));
                Vector pos = Vector(x, y, z);
                pugi::xml_node tnode = obs.child("time");
                if (!tnode)
                    throw(IOexception("InputParser::parseObservers - <point> <time> not found."));
                t0 = parseDouble(tnode.attribute("t0"));
                dt = parseDouble(tnode.attribute("dt"));
                pugi::xml_attribute nt = tnode.attribute("n");
                if (!nt)
                    throw(IOexception("InputParser::parseObservers - <point> <time> attribute n not found."));
                pugi::xml_attribute src = obs.attribute("source");
                if (src)
                {
                    if (std::string(src.value()) == "beam")
                    {
                        PointObserver<Beam> *pointObs = new PointObserver<Beam>(
                            fn.as_string(),
                            pos,
                            t0, dt, nt.as_int() );
                        pointObs->setSourceType(BeamObservation);
                        listObservers->push_back(pointObs);
                    }
                    else if (std::string(src.value()) == "lattice")
                    {
                        PointObserver<Lattice> *pointObs = new PointObserver<Lattice>(
                            fn.as_string(),
                            pos,
                            t0, dt, nt.as_int() );
                        pointObs->setSourceType(LatticeObservation);
                        listObservers->push_back(pointObs);
                    }
                    else throw(IOexception("InputParser::parseObservers - <mesh> illegal source attribute."));
                }
                else
                {
                    PointObserver<Beam> *pointObs = new PointObserver<Beam>(
                        fn.as_string(),
                        pos,
                        t0, dt, nt.as_int() );
                    pointObs->setSourceType(BeamObservation);
                    listObservers->push_back(pointObs);
                }
            }
            else
                throw(IOexception("InputParser::parseObservers - unknown observer type."));
        };
    }
}
