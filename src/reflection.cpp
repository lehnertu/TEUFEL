/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    executable for radiation propagation

  Copyright (c) 2020 U. Lehnert

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

/*!
    \brief TEUFEL general (screen-to-screen) radiation reflection tool.
    
    @author Ulf Lehnert
    @date 21.11.2019
    @file reflection.cpp
    
    The input beam is read from a HDF5 file defining the field traces on a meshed screen.
    The surface of the screen is considered a perfect electric conductor reflecting
    the incident fields. The reflected beam is propagated to an output screen.
    The geometry of the output screen is read from a second HDF5 file.
    From this second file only the geometry information is used.
    the computed output fields will be written to a third file.

*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "global.h"
#include "fields.h"
#include "fieldtrace.h"
#include "mesh.h"

#include <Eigen/Dense>

int main(int argc, char* argv[])
{
    // time-keeping variables used in all long-duration loops
    timespec start_time, stop_time, print_time, now_time;
    double elapsed;
    // loop counter for parallel loops
    volatile int counter;

    std::cout << std::endl;
    std::cout << "Propagate Reflected Radiation from Screen to Screen" << std::endl << std::endl;
    
    if (argc<4)
    {
        std::cout << "Usage: Propagate sourcefile targetfile outfile" << std::endl;
        return 1;
    }
    std::string sourcefile(argv[1]);
    std::string infile(argv[2]);
    std::string outfile(argv[3]);

    // **************************************
    // load the input field from a file
    // **************************************

    std::cout << std::endl << "=== Source Screen ===" << std::endl;
    MeshedScreen *source = new MeshedScreen(sourcefile);
    source->init();
    // print report
    source->writeReport(&cout);

    Vector avg_normal = Vector(0.0,0.0,0.0);
    Vector avg_xi = Vector(0.0,0.0,0.0);
    Vector avg_eta = Vector(0.0,0.0,0.0);
    for (int ip=0; ip<source->get_Np(); ip++)
    {
        avg_normal += source->get_normal(ip);
        avg_xi += source->get_xi(ip);
        avg_eta += source->get_eta(ip);
    };
    
    // smalles box in space and time containing all field points
    double bbx1 = 1e30;
    double bbx2 = -1e30;
    double bby1 = 1e30;
    double bby2 = -1e30;
    double bbz1 = 1e30;
    double bbz2 = -1e30;
    double bbt1 = 1e30;
    double bbt2 = -1e30;
    for (int ip=0; ip<source->get_Np(); ip++)
    {
        Vector p = source->get_point(ip);
        if (p.x<bbx1) bbx1=p.x;
        if (p.x>bbx2) bbx2=p.x;
        if (p.y<bby1) bby1=p.y;
        if (p.y>bby2) bby2=p.y;
        if (p.z<bbz1) bbz1=p.z;
        if (p.z>bbz2) bbz2=p.z;
        double t0 = source->get_t0(ip);
        if (t0<bbt1) bbt1=t0;
        double tmax = t0 + source->get_Nt() * source->get_dt();
        if (tmax>bbt2) bbt2=tmax;
    };
    std::cout << "bounding box:" << std::endl;
    std::cout << "x : (" << bbx1 << ", " << bbx2 << ") m" << std::endl;
    std::cout << "y : (" << bby1 << ", " << bby2 << ") m" << std::endl;
    std::cout << "z : (" << bbz1 << ", " << bbz2 << ") m" << std::endl;
    std::cout << "t : (" << 1e9*bbt1 << ", " << 1e9*bbt2 << ") ns" << std::endl;

    avg_normal /= source->get_Np();
    avg_xi /= source->get_Np();
    avg_eta /= source->get_Np();
    std::cout << "n =   (" << avg_normal.x << ", " << avg_normal.y << ", " << avg_normal.z << ")" << std::endl;
    std::cout << "xi =  (" << avg_xi.x << ", " << avg_xi.y << ", " << avg_xi.z << ")" << std::endl;
    std::cout << "eta = (" << avg_eta.x << ", " << avg_eta.y << ", " << avg_eta.z << ")" << std::endl;
    Vector S(0.0, 0.0, 0.0);
    for (int ip=0; ip<source->get_Np(); ip++)
        S += source->Poynting(ip)*source->get_area(ip);
    std::cout << "total energy flow vector : (" << S.x << ", " << S.y << ", " << S.z << ") J" << std::endl;


    // **************************************
    // compute the surface reflected fields
    // **************************************
    // Just above the surface of a perfect electric conductor only normal electric
    // and tangential magnetic fields can exist.
    // The sum of incident and reflected fields has to satisfy this condition.
    // Therefore, tangential electric and normal magnetic field components
    // are reversed upon reflection. Only the reflected fields are later propagated.

    std::cout << std::endl << "computing the reflected fields ..." << std::endl;
    // record the start time
    clock_gettime(CLOCK_REALTIME, &start_time);

    for (int index=0; index<source->get_Np(); index++)
    {
        // get one field trace
        FieldTrace trace = source->get_trace(index);
        // compute reflection
        Vector normal = source->get_normal(index);
        normal.normalize();
        for (std::size_t it=0; it<trace.get_N(); it++)
        {
            ElMagField EB = trace.get_field(it);
            Vector E = EB.E();
            Vector En = normal * dot(E, normal);
            Vector Et = E - En;
            Vector B = EB.B();
            Vector Bn = normal * dot(B, normal);
            Vector Bt = B - Bn;
            trace.set_field(it, ElMagField(En-Et,Bt-Bn));
        };
        // write it back to the source
        source->set_trace(index, trace);
    };

    // record the finish time
    clock_gettime(CLOCK_REALTIME, &stop_time);
    elapsed = stop_time.tv_sec-start_time.tv_sec +
        1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    std::cout << "time elapsed during computation : " << elapsed << " s" << std::endl;

    Vector Sr(0.0, 0.0, 0.0);
    for (int ip=0; ip<source->get_Np(); ip++)
        Sr += source->Poynting(ip)*source->get_area(ip);
    std::cout << "total energy flow vector : (" << Sr.x << ", " << Sr.y << ", " << Sr.z << ") J" << std::endl;

    // **************************************
    // compute the derivatives of the fields
    // **************************************

    std::cout << std::endl << "computing the time derivatives of the fields ..." << std::endl;
    // record the start time
    clock_gettime(CLOCK_REALTIME, &start_time);

    // compute a full set of time derivatives of the fields
    std::vector<FieldTrace*> source_dA_dt = std::vector<FieldTrace*>(source->get_Np());
    for (int ip=0; ip<source->get_Np(); ip++)
    {
        FieldTrace source_trace = source->get_trace(ip);
        FieldTrace *deriv = source_trace.derivative();
        source_dA_dt[ip] = deriv;
    }
    
    // record the finish time
    clock_gettime(CLOCK_REALTIME, &stop_time);
    elapsed = stop_time.tv_sec-start_time.tv_sec +
        1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    std::cout << "time elapsed during computation : " << elapsed << " s" << std::endl;

    // prepare a full set of normal derivatives of the fields
    // these will be overwritten with data while we compute the transverse derivatives
    std::vector<FieldTrace*> source_dA_dn = std::vector<FieldTrace*>(source->get_Np());
    for (int ip=0; ip<source->get_Np(); ip++)
    {
        FieldTrace *deriv = new FieldTrace(source->get_trace(ip));
        source_dA_dn[ip] = deriv;
    };
    
    std::cout << std::endl << "computing the spatial derivatives of the fields ..." << std::endl;
    // record the start time
    clock_gettime(CLOCK_REALTIME, &start_time);
    print_time = start_time;
    
    double cSquared = SpeedOfLight*SpeedOfLight;
    std::size_t NOTS = source->get_Nt();
    
    counter = 0;
    // parallel domain
    // optional parameter :  num_threads(32)
    // all variables declared outside this block are shared, e.g. source
    // all variables declared inside this block are private
    #pragma omp parallel shared(counter)
    {
        #pragma omp single
        {
            std::cout << "computing on " << omp_get_num_threads() << " parallel threads" << std::endl;
        }
        #pragma omp for
        for (int index=0; index<source->get_Np(); index++)
        {    
            // determine the neighbourhood of this cell
            int nbh[24];
            int count = source->get_Neighbourhood(index,nbh);
            // std::cout << "neighbourhood size = " << count << std::endl;
            // std::cout << "[";
            // for (int i=0; i<count; i++) std::cout << nbh[i] << " ";
            // std::cout << "]" << std::endl;
            
            // get the coordinates of the neighbours in the cell-local coordinate system
            // assemble these into a matrix A
            Vector local_ref = source->get_point(index);
            Vector local_xi = source->get_xi(index);
            Vector local_eta = source->get_eta(index);
            Vector local_n = source->get_normal(index);
            Eigen::VectorXd xi(count), eta(count), xi2(count), eta2(count), xieta(count), one(count);
            for (int i=0; i<count; i++)
            {
                Vector local_coo = source->get_point(nbh[i]) - local_ref;
                local_coo.transform(local_xi,local_eta,local_n);
                one(i) = 1.0;
                xi(i) = local_coo.x;
                eta(i) = local_coo.y;
                xi2(i) = local_coo.x * local_coo.x;
                eta2(i) = local_coo.y * local_coo.y;
                xieta(i) = local_coo.x * local_coo.y;
            };
            Eigen::MatrixXd A(count,6);
            A.col(0) << one;
            A.col(1) << xi;
            A.col(2) << eta;
            A.col(3) << xi2;
            A.col(4) << eta2;
            A.col(5) << xieta;
            
            // get the field traces of the neighbourhood in the local coordinate system
            FieldTrace* local_trace = new FieldTrace(source->get_trace(index));
            FieldTrace* nbh_trace[24];
            for (int i=0; i<count; i++)
            {
                nbh_trace[i] = new FieldTrace(source->get_trace(nbh[i]));
                nbh_trace[i]->transform(local_xi,local_eta,local_n);
            };
            
            // get the time derivative of the center field in local coordinates
            FieldTrace local_dA_dt = source_dA_dt[index];
            local_dA_dt.transform(local_xi,local_eta,local_n);
            
            // this is the trace of normal derivatives we have to compute
            FieldTrace *local_dA_dn = source_dA_dn[index];
            
            for (std::size_t it=0; it<NOTS; it++)
            {
                // we solve the equation A x = E, A x = B in a least-squares sense
                Eigen::VectorXd Exi(count);
                Eigen::VectorXd Eeta(count);
                Eigen::VectorXd En(count);
                Eigen::VectorXd Bxi(count);
                Eigen::VectorXd Beta(count);
                Eigen::VectorXd Bn(count);
                Eigen::VectorXd fit(6);
                double t = local_trace->get_time(it);
                for (int i=0; i<count; i++)
                {
                    // we do not query the fields for the same time index but for the same time
                    ElMagField F = nbh_trace[i]->get_field(t);
                    Exi(i) = F.E().x;
                    Eeta(i) = F.E().y;
                    En(i) = F.E().z;
                    Bxi(i) = F.B().x;
                    Beta(i) = F.B().y;
                    Bn(i) = F.B().z;
                };
                // the solution delivers the spatial derivatives of the fields
                fit = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Exi);
                double dx_Ex = fit(1);
                // double dy_Ex = fit(2);
                fit = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Eeta);
                // double dx_Ey = fit(1);
                double dy_Ey = fit(2);
                fit = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(En);
                double dx_Ez = fit(1);
                double dy_Ez = fit(2);
                fit = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Bxi);
                double dx_Bx = fit(1);
                // double dy_Bx = fit(2);
                fit = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Beta);
                // double dx_By = fit(1);
                double dy_By = fit(2);
                fit = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Bn);
                double dx_Bz = fit(1);
                double dy_Bz = fit(2);

                // obtain the time-derivatives of the fiels
                ElMagField dt_F = local_dA_dt.get_field(it);
                double dt_Ex = dt_F.E().x;
                double dt_Ey = dt_F.E().y;
                double dt_Bx = dt_F.B().x;
                double dt_By = dt_F.B().y;
                
                // use Maxwells equations to compute the normal derivatives of the fields
                double dz_Ex = dx_Ez - dt_By;
                double dz_Ey = dy_Ez + dt_Bx;
                double dz_Ez = -dx_Ex - dy_Ey;
                double dz_Bx = dx_Bz + dt_Ey/cSquared;
                double dz_By = dy_Bz - dt_Ex/cSquared;
                double dz_Bz = -dx_Bx - dy_By;

                // we set the normal derivatives by back-transforming the computed values
                // into the global coordinate system
                Vector E = local_xi*dz_Ex + local_eta*dz_Ey + local_n*dz_Ez;
                Vector B = local_xi*dz_Bx + local_eta*dz_By + local_n*dz_Bz;
                ElMagField f = ElMagField(E,B);
                local_dA_dn->set_field(it,f);
            };
            
            // clean up the traces
            for (int i=0; i<count; i++) delete nbh_trace[i];
            delete local_trace;

            // all threads increment the counter
            #pragma omp atomic
            counter++;
            // only the master thread keeps the time
            if (omp_get_thread_num()==0)
            {
                clock_gettime(CLOCK_REALTIME, &now_time);
                elapsed = now_time.tv_sec-print_time.tv_sec +
                    1e-9*(now_time.tv_nsec-print_time.tv_nsec);
                if (elapsed>30.0)
                {
                    print_time = now_time;
                    elapsed = now_time.tv_sec-start_time.tv_sec +
                    1e-9*(now_time.tv_nsec-start_time.tv_nsec);
                    std::cout << "completed " << std::setprecision(4) << 100.0*(double)counter/(double)(source->get_Np()) << "%";
                    std::cout << "   after " << elapsed << " seconds" << std::endl;
                };
            };
        };
    };
    // end parallel domain
    
    // record the finish time
    clock_gettime(CLOCK_REALTIME, &stop_time);
    elapsed = stop_time.tv_sec-start_time.tv_sec +
        1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    std::cout << "time elapsed during computation : " << elapsed << " s" << std::endl;
    
    // **************************************
    // load the target geometry from a file
    // **************************************

    std::cout << std::endl << "=== Target Screen ===" << std::endl;
    MeshedScreen *target = new MeshedScreen(infile);
    target->init();
    target->zero();
    // print report
    target->writeReport(&cout);
    
    // smalles box in space and time containing all field points
    bbx1 = 1e30;
    bbx2 = -1e30;
    bby1 = 1e30;
    bby2 = -1e30;
    bbz1 = 1e30;
    bbz2 = -1e30;
    bbt1 = 1e30;
    bbt2 = -1e30;
    for (int ip=0; ip<target->get_Np(); ip++)
    {
        Vector p = target->get_point(ip);
        if (p.x<bbx1) bbx1=p.x;
        if (p.x>bbx2) bbx2=p.x;
        if (p.y<bby1) bby1=p.y;
        if (p.y>bby2) bby2=p.y;
        if (p.z<bbz1) bbz1=p.z;
        if (p.z>bbz2) bbz2=p.z;
        double t0 = target->get_t0(ip);
        if (t0<bbt1) bbt1=t0;
        double tmax = t0 + target->get_Nt() * target->get_dt();
        if (tmax>bbt2) bbt2=tmax;
    };
    std::cout << "bounding box:" << std::endl;
    std::cout << "x : (" << bbx1 << ", " << bbx2 << ") m" << std::endl;
    std::cout << "y : (" << bby1 << ", " << bby2 << ") m" << std::endl;
    std::cout << "z : (" << bbz1 << ", " << bbz2 << ") m" << std::endl;
    std::cout << "t : (" << 1e9*bbt1 << ", " << 1e9*bbt2 << ") ns" << std::endl;

    avg_normal = Vector(0.0,0.0,0.0);
    avg_xi = Vector(0.0,0.0,0.0);
    avg_eta = Vector(0.0,0.0,0.0);
    for (int ip=0; ip<target->get_Np(); ip++)
    {
        avg_normal += target->get_normal(ip);
        avg_xi += target->get_xi(ip);
        avg_eta += target->get_eta(ip);
    };
    avg_normal /= target->get_Np();
    avg_xi /= target->get_Np();
    avg_eta /= target->get_Np();
    std::cout << "n =   (" << avg_normal.x << ", " << avg_normal.y << ", " << avg_normal.z << ")" << std::endl;
    std::cout << "xi =  (" << avg_xi.x << ", " << avg_xi.y << ", " << avg_xi.z << ")" << std::endl;
    std::cout << "eta = (" << avg_eta.x << ", " << avg_eta.y << ", " << avg_eta.z << ")" << std::endl;

    double min_t0 = 1.0e30;
    double max_t0 = -1.0e30;
    for (int i=0; i<target->get_Np(); i++)
    {
        double t0 = target->get_t0(i);
        if (t0>max_t0) max_t0=t0;
        if (t0<min_t0) min_t0=t0;
    }
    std::cout << "target start timing range (" << min_t0 << ", " << max_t0 << ")" << std::endl << std::endl;

    // **************************************
    // propagate the fields to the target
    // **************************************
    
    std::cout << std::endl << "propagating the fields ..." << std::endl;
    // record the start time
    clock_gettime(CLOCK_REALTIME, &start_time);
    print_time = start_time;

    // counter for the number of target cells computed
    counter = 0;
    // counter for the number of zero-field warnings encountered
    volatile int warnings_counter = 0;
    // parallel domain
    // optional parameter :  num_threads(32)
    // all variables declared outside this block are shared, e.g. source
    // all variables declared inside this block are private
    #pragma omp parallel shared(counter, warnings_counter)
    {
        #pragma omp single
        {
            std::cout << "computing on " << omp_get_num_threads() << " parallel threads" << std::endl;
        }
        #pragma omp for
        for (int index=0; index<target->get_Np(); index++)
        {
            Vector target_pos = target->get_point(index);
            // this is the result sum (copy of the target trace)
            FieldTrace propagated_trace = target->get_trace(index);
            propagated_trace.zero();
            // this is the component to be added to the result trace
            FieldTrace component = propagated_trace;
            
            for (int isource=0; isource<source->get_Np(); isource++)
            {
                Vector source_pos = source->get_point(isource);
                double dA = source->get_area(isource)/(4.0*Pi);
                Vector RVec = target_pos - source_pos;
                double R = RVec.norm();
                double R2 = R*R;
                double R3 = R2*R;
                Vector Normal = source->get_normal(isource);
                FieldTrace t1 = source->get_trace(isource);
                try
                {
                    t1.retard(R/SpeedOfLight, &component);
                    propagated_trace += component * (-dot(RVec,Normal)/R3) * dA;
                }
                catch(IOexception const & exc)
                {
                    #pragma omp atomic
                    warnings_counter++;
                    if (DEBUGLEVEL>=2)
                    {
                        std::cout << "FieldTrace 1st term zero propagating from i(source)=" << isource;
                        std::cout << " to i(target)=" << index << std::endl;
                        std::cout << "    x = (" << source_pos.x << ", " << source_pos.y << ", " << source_pos.z << ")";
                        std::cout << "    t = " << t1.get_t0() << " ... " << t1.get_t0()+t1.get_N()*t1.get_dt() << std::endl;
                        std::cout << "    x = (" << target_pos.x << ", " << target_pos.y << ", " << target_pos.z << ")";
                        std::cout << "    t = " << component.get_t0();
                        std::cout << " ... " << component.get_t0()+component.get_N()*component.get_dt() << std::endl;
                    }
                }
                try
                {
                    FieldTrace t2 = *source_dA_dt[isource];
                    t2.retard(R/SpeedOfLight, &component);
                    propagated_trace += component * (-dot(RVec,Normal)/(R2*SpeedOfLight)) * dA;
                }
                catch(IOexception const & exc)
                {
                    #pragma omp atomic
                    warnings_counter++;
                    if (DEBUGLEVEL>=2) std::cout << "FieldTrace 2nd term zero propagating from i(source)=" << isource;
                    if (DEBUGLEVEL>=2) std::cout << " to i(target)=" << index << std::endl;
                }
                try
                {
                    FieldTrace t3 = *source_dA_dn[isource];
                    t3.retard(R/SpeedOfLight, &component);
                    // TODO: why is this term positive - should be negative
                    propagated_trace += component * (1.0/R) * dA;
                }
                catch(IOexception const & exc)
                {
                    #pragma omp atomic
                    warnings_counter++;
                    if (DEBUGLEVEL>=2) std::cout << "FieldTrace 3rd term zero propagating from i(source)=" << isource;
                    if (DEBUGLEVEL>=2) std::cout << " to i(target)=" << index << std::endl;
                }
            };
            
            target->set_trace(index, propagated_trace);
            // target->writeTraceReport(&std::cout, index);
            
            // all threads increment the counter
            #pragma omp atomic
            counter++;
            // only the master thread keeps the time
            if (omp_get_thread_num()==0)
            {
                clock_gettime(CLOCK_REALTIME, &now_time);
                elapsed = now_time.tv_sec-print_time.tv_sec +
                    1e-9*(now_time.tv_nsec-print_time.tv_nsec);
                if (elapsed>30.0)
                {
                    print_time = now_time;
                    elapsed = now_time.tv_sec-start_time.tv_sec +
                    1e-9*(now_time.tv_nsec-start_time.tv_nsec);
                    std::cout << "completed " << std::setprecision(4) << 100.0*(double)counter/(double)(target->get_Np()) << "%";
                    std::cout << "   after " << elapsed << " seconds";
                    if (warnings_counter>0) std::cout << "   received " << warnings_counter << " 'zero field' warnings";
                    std::cout << std::endl;
                };
            };
        };
    };
    // end parallel domain

    // record the finish time
    clock_gettime(CLOCK_REALTIME, &stop_time);
    elapsed = stop_time.tv_sec-start_time.tv_sec +
        1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    std::cout << "time elapsed during computation : " << elapsed << " s" << std::endl;
    
    // print report
    std::cout << std::setprecision(6);
    target->writeReport(&cout);
    // write the target data to file
    target->writeFile(outfile);
    
    // release all allocated memory
    for (int ip=0; ip<source->get_Np(); ip++) delete source_dA_dt[ip];
    for (int ip=0; ip<source->get_Np(); ip++) delete source_dA_dn[ip];
    delete source;
    delete target;
    
    return 0;
}
