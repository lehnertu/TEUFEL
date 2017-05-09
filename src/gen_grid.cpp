
#include "gen_grid.h"
#include <math.h>
#include "global.h"
#include <iostream>
#include "particle.h"
using namespace std;


GenGrid::GenGrid()
{

}

void GenGrid::GenPlanarGrid(double x,double y,double z,int NumGridPoints)
{
    Numpoints=NumGridPoints;
    double zz=z;
    TruePoints=NumGridPoints*NumGridPoints;
    double dx=x/((double)Numpoints);
    double dy=y/((double)Numpoints);
    GridLength = x;
    GridWidth = y;
    MeshArea=dx*dy;
    int counter=0;
    r=new Vector[TruePoints];
    double xx=-x/2.0;
    for (int i=0;i<Numpoints;i++)
    {
        xx=xx+dx;
        double yy=-y/2.0;
//#pragma omp parallel for private(yy,xx,zz)
        for (int j=0;j<Numpoints;j++)
        {
            yy=yy+dy;
            r[counter]=Vector(xx,yy,zz);
            counter=counter+1;

        }
    }

}

void GenGrid::GenXZGrid(Vector position,double x,double y,double z,int NumGridPoints)
{
    Numpoints=NumGridPoints;
    TruePoints=NumGridPoints*NumGridPoints;
    double dx=x/((double)Numpoints);
    double dz=z/((double)Numpoints);
    MeshArea=dx*dz;
    int counter=0;
    r=new Vector[TruePoints];
    double xx= position.x - x/2.0;
    for (int i=0;i<Numpoints;i++)
    {
        xx=xx+dx;
        double zz= position.z -z/2.0;
//#pragma omp parallel for private(yy,xx,zz)
        for (int j=0;j<Numpoints;j++)
        {
            zz=zz+dz;
            r[counter]=Vector(xx,position.y,zz);
            counter=counter+1;

        }
    }

}


Vector GenGrid::GetGridPoints(int i)
{
    return r[i];
}

int GenGrid::GetGridSize()
{
    return TruePoints;
}

double GenGrid::GetMeshArea()
{
    return MeshArea;
}

Vector GenGrid::GetNormalVector()
{
    Vector N=Vector(0.0,0.0,1);
    return N;
}


int GenGrid::SDDSRadiationAtGrid(Bunch *bunch, double time_begin, double time_end, int NumberOfPoints, const char* filename )
{
	double dt = (time_end-time_begin)/(double)NumberOfPoints;
	SDDS_DATASET data;
	char buffer[100];
	int Initialize = SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,filename);
	if(Initialize!=1)
	{
			cout<<"Error Initializing Output\n";
			return 1;
	}
	else
	{
			cout<<"output initialized\n";
	}

	if(
			SDDS_DefineSimpleParameter(&data,"NumberOfTiles",NULL, SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleParameter(&data,"TileArea","m2", SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleParameter(&data,"time","s", SDDS_DOUBLE)!=1
	  )
	{
			cout<<"error in defining parameters\n";
			return 2;;
	}
	else
	{
			cout<<"parameters defined \n";
	}
	
	if(
			SDDS_DefineColumn(&data,"x\0","x\0","m","Xposition\0",NULL, SDDS_DOUBLE,0) ==-1 ||
 			SDDS_DefineColumn(&data,"y\0","y\0","m\0","Yposition\0",NULL, SDDS_DOUBLE,0)   ==-1 || 
			SDDS_DefineColumn(&data,"z\0","z\0","m\0","Zposition\0",NULL, SDDS_DOUBLE,0)   ==-1 || 
			SDDS_DefineColumn(&data,"Ex\0","Ex\0","V/m\0","ElectricFieldInX\0",NULL, SDDS_DOUBLE,0) == -1 ||
			SDDS_DefineColumn(&data,"Ey\0","Ey\0","V/m\0","ElectricFieldInY\0",NULL, SDDS_DOUBLE,0) == -1 ||
			SDDS_DefineColumn(&data,"Ez\0","Ez\0","V/m\0","ElectricFieldInZ\0",NULL, SDDS_DOUBLE,0) == -1 || 
			SDDS_DefineColumn(&data,"Bx\0","Bx\0","Tesla\0","MagneticFieldInX\0",NULL, SDDS_DOUBLE,0)== -1 || 
			SDDS_DefineColumn(&data,"By\0","By\0","Tesla\0","MagneticFieldInY\0",NULL,SDDS_DOUBLE,0) == -1 ||
			SDDS_DefineColumn(&data,"Bz\0","Bz\0","Tesla\0","MagneticFieldInZ\0",NULL,SDDS_DOUBLE,0)==-1   ||
			SDDS_DefineColumn(&data,"PoyntingVector\0","|S|\0","dp/da\0","MagnitudeOfPoyntingVector\0",NULL,SDDS_DOUBLE,0) == -1
	  )
	{
			cout<<"error in defining columns\n";
			return 3;
	}
	else
	{
			cout<<"Columns defined \n";	
	}

	if (SDDS_WriteLayout(&data)!=1)
	{
			cout<<"error in writing layout\n";
			return 4;
	}
	else
	{
			cout<<"layout written \n";
	}


	if (SDDS_StartPage(&data,(int32_t)TruePoints)!=1)
	{
		cout<<"error in starting page\n";
		return 5;
	}
	{
		cout<<"pages initialized\n";
	}

	
	for (int32_t i = 0;i<(int32_t)NumberOfPoints;i++)
	{

		double TimeObs = time_begin+i*dt;

		if (
			SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			"NumberOfTiles",(double)TruePoints,
			"TileArea",MeshArea,
			"time",time_begin+i*dt,
			 NULL)!=1
	   	   )
			{
				cout<<"error in defining parameter\n";
				return 6;
			}

		pair<Vector,Vector>*Field = new pair<Vector,Vector>[TruePoints]; 
#pragma omp parallel for 
		for (int tim=0;tim<TruePoints;tim++)
		{
			Vector Robs = r[tim];
			Field[tim] = bunch -> RadiationField(Robs,TimeObs);
		}
			

		for (int tim=0;tim<TruePoints;tim++)
		{
			Vector Robs = r[tim];		
			Vector E = (Field[tim].first);
			Vector B = (Field[tim].second);
			Vector Poynting = cross(E,B/MuNull);
			if(SDDS_SetRowValues(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,tim,
				"x",(double)(Robs.x),
				"y",(double)(Robs.y),
				"z", (double)(Robs.z),
				"Ex",(double)(E.x),
				"Ey",(double)(E.y),
				"Ez",(double)(E.z),
				"Bx",(double)(B.x),
				"By",(double)(B.y),
				"Bz",(double)(B.z),
				"PoyntingVector",(double)(Poynting.norm()),
				NULL)!=1
				)

				{
					cout<<"error in writing columns\n";
					return 7;
				}
		}

		if(SDDS_WritePage(&data)!=1)	
		{
			cout<<"error in writing page\n";
			return 8;
		}

	}

	
	
	if (SDDS_Terminate(&data)!=1)
	{
		cout<<"error terminating data\n";
		return 9;
	}	

	return 0;





}

int GenGrid::SDDSFilteredRadiationAtGrid(double RequiredFrequency, const char* filename, Bunch *bunch, double time_begin, double time_end, int NumberOfPoints )
{
	double dt = (time_end-time_begin)/(double)NumberOfPoints;
	SDDS_DATASET data;
	char buffer[100];
	int Initialize = SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,filename);
	if(Initialize!=1)
	{
			cout<<"Error Initializing Output\n";
			return 1;
	}
	else
	{
			cout<<"output initialized\n";
	}

	if(
			SDDS_DefineSimpleParameter(&data,"NumberOfTiles",NULL, SDDS_DOUBLE)!=1 || 
			SDDS_DefineSimpleParameter(&data,"TileArea","m2", SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleParameter(&data,"Length","m", SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleParameter(&data,"Breadth","m", SDDS_DOUBLE)!=1 ||
			SDDS_DefineSimpleParameter(&data,"RequiredFrequency","Hz", SDDS_DOUBLE)!=1
	  )
	{
			cout<<"error in defining parameters\n";
			return 2;;
	}
	else
	{
			cout<<"parameters defined \n";
	}
	
	if(
			SDDS_DefineColumn(&data,"x\0","x\0","m","Xposition\0",NULL, SDDS_DOUBLE,0) ==-1 ||
 			SDDS_DefineColumn(&data,"y\0","y\0","m\0","Yposition\0",NULL, SDDS_DOUBLE,0)   ==-1 || 
			SDDS_DefineColumn(&data,"z\0","z\0","m\0","Zposition\0",NULL, SDDS_DOUBLE,0)   ==-1 || 
			SDDS_DefineColumn(&data,"ERealNorm\0","|ERe|\0","|V/m|\0","ElectricFieldMagnitude\0",NULL, SDDS_DOUBLE,0) == -1 ||
			SDDS_DefineColumn(&data,"EImagNorm\0","|EImg|\0","|V/m|\0","ElectricFieldMagnitude\0",NULL, SDDS_DOUBLE,0) == -1 ||
			SDDS_DefineColumn(&data,"EAbs\0","|EAbs|\0","|V/m|\0","ElectricFieldMagnitude\0",NULL, SDDS_DOUBLE,0) == -1 ||
			SDDS_DefineColumn(&data,"Phase\0","Phase\0","/theta\0","PhaseOfElectricField\0",NULL,SDDS_DOUBLE,0)==-1
	  )
	{
			cout<<"error in defining columns\n";
			return 3;
	}
	else
	{
			cout<<"Columns defined \n";	
	}

	if (SDDS_WriteLayout(&data)!=1)
	{
			cout<<"error in writing layout\n";
			return 4;
	}
	else
	{
			cout<<"layout written \n";
	}


	if (SDDS_StartPage(&data,(int32_t)TruePoints)!=1)
	{
		cout<<"error in starting page\n";
		return 5;
	}
	{
		cout<<"pages initialized\n";
	}

	if (
			SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			"NumberOfTiles",(double)TruePoints,
			"TileArea",MeshArea,
			"Length",GridLength,
			"Breadth",GridWidth,
			"RequiredFrequency",RequiredFrequency,
			 NULL)!=1
	   	   )
			{
				cout<<"error in defining parameter\n";
				return 6;
			}
	

	for (int32_t i = 0;i<(int32_t)TruePoints;i++)
	{

		
		Vector Robs = r[i];
		Vector EReal =Vector(0,0,0);
		Vector EImag = Vector(0,0,0);
		double ERx=0;
		double ERy=0;
		double ERz=0;
		double EImx=0;
		double EImy=0;
		double EImz=0;
#pragma omp parallel for reduction(+:ERx,ERy,ERz,EImx,EImy,EImz)
		for (int tim=0;tim<NumberOfPoints;tim++)
		{
			double Ex,Ey,Ez;
			double TimeObs = time_begin+tim*dt;
			double OmegaT = RequiredFrequency*2.0*Pi*TimeObs;
			pair<Vector,Vector>Field; 
			Field = bunch -> RadiationField(Robs,TimeObs);
			ERx += (Field.first).x*cos(OmegaT);
			ERy += (Field.first).y*cos(OmegaT);
			ERz += (Field.first).z*cos(OmegaT);
			EImx += (Field.first).x*sin(OmegaT);
			EImy += (Field.first).y*sin(OmegaT);
			EImz += (Field.first).z*sin(OmegaT);
			
		
		}
		EReal = Vector(ERx,ERy,ERz)/NumberOfPoints;
		EImag = Vector(EImx,EImy,EImz)/NumberOfPoints;
			
		if(SDDS_SetRowValues(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,
			"x",(double)(Robs.x),
			"y",(double)(Robs.y),
			"z",(double)(Robs.z),
			"ERealNorm",(double)(EReal.norm()),
			"EImagNorm",(double)(EImag.norm()),
			"EAbs",(double)sqrt((EReal.abs2nd()+EImag.abs2nd())),
			"Phase",(double)(atan(EImag.norm()/EReal.norm())),
			NULL)!=1
			)
			{
				cout<<"error in writing columns\n";
				return 7;
			}
		

		

	}

	
	if(SDDS_WritePage(&data)!=1)	
		{
			cout<<"error in writing page\n";
			return 8;
		}
	if (SDDS_Terminate(&data)!=1)
	{
		cout<<"error terminating data\n";
		return 9;
	}	

	return 0;





}

