
#include "gen_grid.h"
#include <math.h>
#include "global.h"
#include <iostream>
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
