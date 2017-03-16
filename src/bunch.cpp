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

#include "bunch.h"
#include "particle.h"
#include "global.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <tuple>
#include <externalfield.h>
#include "omp.h"
#include "SDDS.h"

Bunch::Bunch()
{
	//Create a single charged particle with default charge, mass values.
	setNOP(1);
	b=new ChargedParticle[1];
	//Initial values has to be set at the initialization
	b[0].setInitialMomentum(Vector(0.0,0,0.0));
	b[0].setInitialPosition(Vector(0.0,0.00,0.0));
	b[0].setInitialTime(0.0);	
	b[0].setParticleID(0);
	b[0].setCharge(-1);
	b[0].setMass(1);

}

Bunch::Bunch(const char *filename, int NP, int charge, int mass)
{
	//set the number of particles equal to NP
	setNOP(NP);
	b=new ChargedParticle[NP];
	Charge=charge;
	Mass=mass;
	//check if the file has correct number of rows and columns
	//Rows must be equal to number of particles
	//columns must be equal to 7 to completely describe phase space
	// (x,px,y,py,z,pz,time)
	LoadBeamProfile(filename, b);

};


Bunch::Bunch(const Bunch *bunch)
{
	int NP = bunch->NOP;
	setNOP(NP);
	int NOTS = bunch->NT;
	b=new ChargedParticle[NP];
	InitializeTrajectory(NOTS);
	for (int i=0;i<NOP;i++)
	{
//#pragma omp parallel for
		for (int k =0;k<NOTS;k++)
		{
			b[i].setTrajPoint(k,(bunch->b[i]).TrajPoint(k));
			b[i].setTrajMomentum(k,(bunch->b[i]).TrajMomentum(k));
			b[i].setTrajTime(k,(bunch->b[i]).TrajTime(k));
			b[i].setTrajAcceleration(k,(bunch->b[i]).TrajAcceleration(k));
			b[i].setParticleID(bunch->b[i].getParticleID());
			b[i].setCharge(bunch->b[i].getCharge());
			b[i].setMass(bunch->b[i].getMass());
		}
	}

}
// set the number of particles the bunch will have
void Bunch:: setNOP(int NP)
{
	NOP=NP;

};


int Bunch::getNOP()
{
	return NOP;
};

int Bunch:: getNOTS()
{
	return NT;
}

int Bunch:: getCharge()
{
	return Charge;
}

int Bunch:: getMass()
{
	return Mass;
}
tuple<Vector,Vector> Bunch::MutualField(int stepnumber, int ParticleID, double t)
{
	Vector Ef = Vector(0.0,0.0,0.0);
	Vector Bf = Vector(0.0,0.0,0.0);
	int j=ParticleID;
	Vector Robs=b[j].TrajPoint(stepnumber);
	tuple<Vector,Vector>FF[NOP];
#pragma omp parallel for shared(j,stepnumber,t,Robs)
	for (int i=0; i<NOP;i++)
	{
		// for i == j; i.e. field due to particle i on its position is zero
		//implemented in the InteractionField routine of particles
		FF[i]= b[i].InteractionField(j,stepnumber,t, Robs);
	}
	// tuple of electric and magnetic fields have been obtained
	//add all of them together
	
	double Ex,Ey,Ez,Bx,By,Bz;
	Ex=0.0;
	Ey=0.0;
	Ez=0.0;
	Bx=0.0;
	By=0.0;
	Bz=0.0;
#pragma omp parallel for reduction(+:Ex,Ey,Ez)
	for(int i=0;i<NOP;i++)
	{
		Ex+=get<0>(FF[i]).x;
		Ey+=get<0>(FF[i]).y;
		Ez+=get<0>(FF[i]).z;
	}
#pragma omp parallel for reduction(+:Bx,By,Bz)
	for(int i=0;i<NOP;i++)
	{
		Bx+=get<1>(FF[i]).x;
		By+=get<1>(FF[i]).y;
		Bz+=get<1>(FF[i]).z;
	}

	Ef=Vector(Ex,Ey,Ez);
	Bf=Vector(Bx,By,Bz);

	return make_tuple(Ef,Bf);
}


//evaluate the radiation field profile
//sums up the electric and magnetic fields at any screen

// gives only the retarded field values
tuple<Vector,Vector> Bunch::RadiationField(Vector Robs, double t)
{
	Vector Ef = Vector(0.0,0.0,0.0);
	Vector Bf = Vector(0.0,0.0,0.0);
	tuple<Vector,Vector>FF[NOP];
#pragma omp parallel for shared(t,Robs)
	for (int i=0; i<NOP;i++)
	{
		// for i == j; i.e. field due to particle i on its position is zero
		//implemented in the InteractionField routine of particles
		FF[i]= b[i].RetardedEField(t, Robs);
	}
	// tuple of electric and magnetic fields have been obtained
	//add all of them together
	
	double Ex,Ey,Ez,Bx,By,Bz;
	Ex=0.0;
	Ey=0.0;
	Ez=0.0;
	Bx=0.0;
	By=0.0;
	Bz=0.0;
//#pragma omp parallel for reduction(+:Ex,Ey,Ez)
	for(int i=0;i<NOP;i++)
	{
		Ex+=get<0>(FF[i]).x;
		Ey+=get<0>(FF[i]).y;
		Ez+=get<0>(FF[i]).z;
	}
//#pragma omp parallel for reduction(+:Bx,By,Bz)
	for(int i=0;i<NOP;i++)
	{
		Bx+=get<1>(FF[i]).x;
		By+=get<1>(FF[i]).y;
		Bz+=get<1>(FF[i]).z;
	}

	Ef=Vector(Ex,Ey,Ez);
	Bf=Vector(Bx,By,Bz);

	return make_tuple(Ef,Bf);
}

// Track the particles in an external field 
//field is generated by the lattice
//the interaction fields are added with the lattice fields

void Bunch::Track_Euler(int NOTS, double tstep, Lattice *field)
{

	//allocate memory to every particle
	//to store trajectory details
	InitializeTrajectory(NOTS);
	NT = NOTS;
	TIMESTEP=tstep;
	double start=omp_get_wtime();
	for (int i=0;i<NOTS;i++)
	{
		
#pragma omp parallel for 
		for(int k=0;k<NOP;k++)
		{
			Vector X0=b[k].TrajPoint(i);
			Vector P0=b[k].TrajMomentum(i);
			double t0=b[k].TrajTime(i);
			double qm=b[k].getCharge()*InvRestMass/b[k].getMass();
			Vector p=b[k].TrajMomentum(i);
			double betagamma2=p.abs2nd();
			double gamma=sqrt(betagamma2+1.0);
			Vector beta = p/gamma;
			
			tuple<Vector,Vector>MField=MutualField(i,k,t0);
    			Vector efield = field->EField(t0, X0) + get<0>(MField);
    			Vector bfield = field->BField(t0,X0)+ get<1>(MField);
			//Vector efield = field->EField(t0, X0) ;
    			//Vector bfield = field->BField(t0,X0);
    			Vector force = cross(beta, bfield) + efield/SpeedOfLight;
   			Vector dX_dt = beta * SpeedOfLight;
   			Vector dP_dt = force * qm;
    			// store acceleration
   			b[k].setTrajAcceleration(i,dP_dt);
   			// integrator step
			b[k].setTrajPoint(i+1,X0+dX_dt*tstep);
			b[k].setTrajTime(i+1,t0+tstep);
			b[k].setTrajMomentum(i+1,P0+dP_dt*tstep);
		}
		
	}
	double end=omp_get_wtime();
	cout<<"\033[1;31m Tracking Completed in: \033[0m"<<(end-start)<<"\033[1;31m seconds\033[0m\n"<<endl;

}

void Bunch::Track_Vay(int NOTS, double tstep, Lattice *field)
{
	InitializeTrajectory(NOTS);
	NT = NOTS;
	TIMESTEP=tstep;
	double qm=b[0].getCharge()*InvRestMass/b[0].getMass();
	double t2=0.5*tstep;
	double qmt2=qm*t2;
	double t_h[NOP];
	Vector x_h[NOP];
	Vector p_h[NOP];
	double gamma_h[NOP];
	Vector beta_h[NOP];
	Vector E_h[NOP];
	Vector B_h[NOP];
	Vector dp_dt[NOP];
	Vector A[NOP];
	Vector p_i1[NOP];
	double gamma_i1[NOP];
	Vector beta_i1[NOP];
	tuple<Vector,Vector> F[NOP];
	double start=omp_get_wtime();
//first of all setup the initial information for half steps
//for all the particles	
#pragma parallel for private(t_h,x_h,p_h,gamma_h,beta_h,F,E_h,B_h,dp_dt,p_i1,gamma_i1,beta_i1,k,qmt2,qm,tstep)
	for(int k=0;k<NOP;k++)
		{
			t_h[k]=b[k].getInitialTime();
			x_h[k]=b[k].getInitialPosition();
			p_h[k]=b[k].getInitialMomentum();
			gamma_h[k]=sqrt(p_h[k].abs2nd()+1.0);
			beta_h[k]=p_h[k]/gamma_h[k];
			F[k]=MutualField(0, k, t_h[k]);
			E_h[k]=field->EField(t_h[k],x_h[k])+get<0>(F[k]);
			B_h[k]=field->BField(t_h[k],x_h[k])+get<1>(F[k]);
			//E_h[k] = field->EField(t_h[k],x_h[k]) ;
    			//B_h[k]=field->BField(t_h[k],x_h[k]);
			dp_dt[k]=(cross(beta_h[k],B_h[k])+E_h[k]/SpeedOfLight)*qm;
			p_i1[k]=p_h[k]-dp_dt[k]*0.5*tstep;
			gamma_i1[k]=sqrt(p_i1[k].abs2nd()+1.0);
			beta_i1[k]=p_i1[k]/gamma_i1[k];	
			b[k].setTrajAcceleration(0,dp_dt[k]);
			
		};

	for (int i=0;i<NOTS;i++)
	{
		
#pragma parallel for private(t2,i,t_h,x_h,p_h,gamma_h,beta_h,F,E_h,B_h,dp_dt,p_i1,gamma_i1,beta_i1,k,qmt2,qm,tstep)
		for (int k=0;k<NOP;k++)
		{
			
			Vector p_i=p_i1[k];	
			Vector beta_i=beta_i1[k];
			F[k]=MutualField(i, k, t_h[k]);
			E_h[k]=field->EField(t_h[k],x_h[k])+get<0>(F[k]);
			B_h[k]=field->BField(t_h[k],x_h[k])+get<1>(F[k]);
			//E_h[k] = field->EField(t_h[k],x_h[k]) ;
    			//B_h[k]=field->BField(t_h[k],x_h[k]);
			dp_dt[k]=(cross(beta_i,B_h[k])+E_h[k]/SpeedOfLight)*qm;
			p_h[k]=p_i+dp_dt[k]*t2;			
			Vector p_prime=p_h[k]+E_h[k]/SpeedOfLight * qmt2;			
			double gamma_prime=sqrt(p_prime.abs2nd()+1.0);			
			Vector tau=B_h[k]*qmt2;			
			double u_star=dot(p_prime,tau);			
			double tau2nd=tau.abs2nd();			
			double sigma=gamma_prime*gamma_prime-tau2nd;
			gamma_i1[k]=sqrt(0.5*(sigma+sqrt(sigma*sigma+4.0*(tau2nd+u_star*u_star))));			
			Vector T=tau/gamma_i1[k];
			p_i1[k]=(p_prime+T*dot(p_prime,T)+cross(p_prime,T))/(1.0+T.abs2nd());
			beta_i1[k]=p_i1[k]/gamma_i1[k];
			int stepnumber = i+1;
			b[k].setTrajTime(stepnumber,t_h[k]);
			b[k].setTrajPoint(stepnumber,x_h[k]);
			b[k].setTrajMomentum(stepnumber,p_h[k]);
			b[k].setTrajAcceleration(stepnumber,dp_dt[k]);
			x_h[k]=x_h[k]+beta_i1[k]*SpeedOfLight*tstep;
			t_h[k]=t_h[k]+tstep;
			b[k].setTrajTime(0,b[k].getInitialTime());
			b[k].setTrajMomentum(0,b[k].getInitialMomentum());
		}
	}


	double end=omp_get_wtime();
	cout<<"\033[1;31m Tracking Completed in: \033[0m"<<(end-start)<<"\033[1;31m seconds\033[0m\n"<<endl;

}
void Bunch::WriteSDDS()
{
	SDDS_DATASET data;
	int32_t i;
	int32_t k;
	char buffer[100];
	int Initialize = SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,"trajectory.sdds");
	double charge,mass,nots,nop;
	charge=(double)getCharge();
	mass = (double)getMass();
	nots =(double) getNOTS();
	nop =(double) getNOP();
	double TotalTime=(double)nots*TIMESTEP;
	//nop =(double)10;
	//nots = (double)10;
	double Qtot = nop*charge*ElementaryCharge;
	double Mtot = nop*mass*9.1e-31;
	if(Initialize!=1)
	{
		cout<<"Error Initializing Output\n";
	}
	else
	{
		fprintf(stdout,"output initialized\n");
	}

	if(
			SDDS_DefineSimpleParameter(&data,"NumberOfParticles",NULL, SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleParameter(&data,"TotalTime","seconds", SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleParameter(&data,"BunchTotalCharge","Coloumbs", SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleParameter(&data,"BunchTotalMass","kg", SDDS_DOUBLE)!=1
	  )
	{
		fprintf(stdout,"error in defining parameters\n");
		
	}
	else{
	fprintf(stdout,"parameters defined \n");
	}
	
	if(
			SDDS_DefineSimpleColumn(&data,"Number",NULL, SDDS_DOUBLE)!=1 ||
 			SDDS_DefineSimpleColumn(&data,"t","sec", SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleColumn(&data,"x","m", SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleColumn(&data,"y","m", SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleColumn(&data,"z","m", SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleColumn(&data,"px","gammabetax", SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleColumn(&data,"py","gammabetay",SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleColumn(&data,"pz","gammabetaz", SDDS_DOUBLE) !=1 ||
			SDDS_DefineSimpleColumn(&data,"Ax","msec2", SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleColumn(&data,"Ay","msec2",SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleColumn(&data,"Az","msec2", SDDS_DOUBLE)!=1 
	  )
	{
		fprintf(stdout,"error in defining columns\n");
	}
	else{
	fprintf(stdout,"Columns defined \n");	
	}

	if (SDDS_WriteLayout(&data)!=1)
	{
		fprintf(stdout,"error in writing layout\n");
	}
	else{
	fprintf(stdout,"layout written \n");
	}

	for (i = 0;i<(int32_t)nop;i++)
	{
		if (SDDS_StartPage(&data,(int32_t)nots)!=1)
			{
				fprintf(stdout,"error in starting page\n");
			}
		if (
			SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			"NumberOfParticles",nop,
			"TotalTime",TotalTime,
			"BunchTotalCharge",Qtot,
			"BunchTotalMass",Mtot,
			 NULL)!=1)
			{
				fprintf(stdout,"error in defining parameter\n");
			}
			
		for( k =0;k<(int32_t)nots;k++)
		{
			
			
			if(SDDS_SetRowValues(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,k,
						"Number",(double)(k),
						"t",(double)(b[i].TrajTime(k)),
						"x",(double)((b[i].TrajPoint(k)).x),
						"y",(double)((b[i].TrajPoint(k)).y),
						"z",(double)((b[i].TrajPoint(k)).z),
						"px",(double)((b[i].TrajMomentum(k)).x),
						"py",(double)((b[i].TrajMomentum(k)).y),
						"pz",(double)((b[i].TrajMomentum(k)).z),
						"Ax",(double)((b[i].TrajAcceleration(k)).x),
						"Ay",(double)((b[i].TrajAcceleration(k)).y),
						"Az",(double)((b[i].TrajAcceleration(k)).z),
						NULL)!=1)
			{
				fprintf(stdout,"error in writing columns\n");
			}
	
		}
		if(SDDS_WritePage(&data)!=1)
			{
				fprintf(stdout,"error in writing page\n");
			}

	}
	
	if (SDDS_Terminate(&data)!=1)
	{
		fprintf(stdout,"error terminating data\n");
	}	

}

//write sdds file with time as a parameter. each page has x,y,z,px,py,pz at different time steps
void Bunch::WriteSDDS_Time()
{
	SDDS_DATASET data;
	int32_t i;
	int32_t k;
	char buffer[100];
	int Initialize = SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,"time-trajectory.sdds");
	double charge,mass,nots,nop;
	charge=(double)(getCharge());
	mass = (double)(getMass());
	nots =(double) (getNOTS());
	nop =(double) (getNOP());
	//nop =(double)10;
	//nots = (double)10;
	double Qtot = nop*charge*ElementaryCharge;
	double Mtot = nop*mass*9.1e-31;
	double TotalTime=(double)NT*TIMESTEP;
	if(Initialize!=1)
	{
		cout<<"Error Initializing Output\n";
	}
	else
	{
		fprintf(stdout,"output initialized\n");
	}

	if(
			SDDS_DefineSimpleParameter(&data,"NumberOfParticles",NULL, SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleParameter(&data,"TotalTime","seconds", SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleParameter(&data,"BunchTotalCharge","Coloumbs", SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleParameter(&data,"BunchTotalMass","kg", SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleParameter(&data,"Time","seconds", SDDS_DOUBLE)!=1
	  )
	{
		fprintf(stdout,"error in defining parameters\n");
		
	}
	else{
	fprintf(stdout,"parameters defined \n");
	}
	
	if(
			SDDS_DefineSimpleColumn(&data,"Number",NULL, SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleColumn(&data,"x","m", SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleColumn(&data,"y","m", SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleColumn(&data,"z","m", SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleColumn(&data,"px","gammabetax", SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleColumn(&data,"py","gammabetay",SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleColumn(&data,"pz","gammabetaz", SDDS_DOUBLE) !=1 ||
			SDDS_DefineSimpleColumn(&data,"Ax","msec2", SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleColumn(&data,"Ay","msec2",SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleColumn(&data,"Az","msec2", SDDS_DOUBLE)!=1 
	  )
	{
		fprintf(stdout,"error in defining columns\n");
	}
	else{
	fprintf(stdout,"Columns defined \n");	
	}

	if (SDDS_WriteLayout(&data)!=1)
	{
		fprintf(stdout,"error in writing layout\n");
	}
	else{
	fprintf(stdout,"layout written \n");
	}

	for (i = 0;i<(int32_t)nots;i++)
	{
		if (SDDS_StartPage(&data,(int32_t)nop)!=1)
			{
				fprintf(stdout,"error in starting page\n");
			}
		if (
			SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			"NumberOfParticles",nop,
			"TotalTime",TotalTime,
			"BunchTotalCharge",Qtot,
			"BunchTotalMass",Mtot,
			"Time",(double)(b[0].TrajTime(i)),
			 NULL)!=1)
			{
				fprintf(stdout,"error in defining parameter\n");
			}
			
		for( k =0;k<(int32_t)nop;k++)
		{
			
			
			if(SDDS_SetRowValues(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,k,
						"Number",(double)(k),
						"x",(double)((b[k].TrajPoint(i)).x),
						"y",(double)((b[k].TrajPoint(i)).y),
						"z",(double)((b[k].TrajPoint(i)).z),
						"px",(double)((b[k].TrajMomentum(i)).x),
						"py",(double)((b[k].TrajMomentum(i)).y),
						"pz",(double)((b[k].TrajMomentum(i)).z),
						"Ax",(double)((b[k].TrajAcceleration(i)).x),
						"Ay",(double)((b[k].TrajAcceleration(i)).y),
						"Az",(double)((b[k].TrajAcceleration(i)).z),
						NULL)!=1)
			{
				fprintf(stdout,"error in writing columns\n");
			}
	
		}
		if(SDDS_WritePage(&data)!=1)
			{
				fprintf(stdout,"error in writing page\n");
			}

	}
	
	if (SDDS_Terminate(&data)!=1)
	{
		fprintf(stdout,"error terminating data\n");
	}	

}
void Bunch::InitializeTrajectory(int NOTS)
{
	for (int i=0; i<NOP; i++)
	{
		b[i].setNP(NOTS);
	}

}

void Bunch::LoadBeamProfile(const char *filename, const ChargedParticle *part)
{

	//Generate a stream to the Beam Profile file
	ifstream BunchProfile;
	file=filename;
	//Check if file has correct number of rows and coumns
	//if ok, then stream the values.
	//else print error message and exit properly
	if (FileCheck(file,NOP)==1)
	{		
		ifstream BeamProfile(filename);
		printf("\033[7;31m Loading Beam Profile....\n\033[0m\n");
		for(int i=0;i<NOP;i++)
		{
			double t;			
			Vector position;
			Vector momentum;
			BeamProfile>>t;
			BeamProfile>>position.x;
			BeamProfile>>position.y;
			BeamProfile>>position.z;
			BeamProfile>>momentum.x;
			BeamProfile>>momentum.y;
			BeamProfile>>momentum.z;
			b[i].setInitialPosition(position);
			b[i].setInitialMomentum(momentum);
			b[i].setInitialTime(t);
			b[i].setCharge(Charge);
			b[i].setMass(Mass);
			b[i].setParticleID(i);
			
		}
		printf("\033[7;31m Beam Profile Loaded....\n\033[0m\n");
		BeamProfile.close();
	}
	else
	{
		printf("\033[7;31m Error Loading Beam Profile....\n\033[0m\n");
		printf("\033[7;31m Exiting Program\n\033[0m\n");
		exit (EXIT_FAILURE);
	}
}


int Bunch::FileCheck(const char *filename, int NP)
{
	NOP=NP;
	ifstream InFile(filename);
	ifstream InFile2(filename);
	int rows,columns,tabs;
	if (InFile.is_open() && InFile2.is_open())
	{
	  rows=count(istreambuf_iterator<char>(InFile),istreambuf_iterator<char>(),'\n');
	  tabs=count(istreambuf_iterator<char>(InFile2),istreambuf_iterator<char>(),'\t');
	  columns=tabs/rows+1;
	  if(rows!=NOP || columns!=7)
		{
			printf("\033[1;31m Incorrect Phase Space Distribution\n Number of Particles for which Info Provided: %d \n Number of phase space parameters provided: %d \n \033[0m\n", rows, columns);
			return 0;
		}
	  else {return 1;}
	  
	}
	else {
		return 0;	
		printf("Error opening the file: %p", filename);};
	InFile.close();
	InFile2.close();
}




  

