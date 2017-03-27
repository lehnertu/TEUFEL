#include <math.h>
#include "SDDStypes.h"
#include "SDDS.h"
#include <stdio.h>

using namespace std;
int WriteSDDSOutput();
int ReadSDDSInput()
{
	return 1;
};


int WriteSDDSOutput()
{
	SDDS_DATASET data;
	int32_t i;
	char buffer[100];
	if (SDDS_InitializeOutput(&data,SDDS_BINARY,1,NULL,NULL,"demo.sdds")!=1)
	{
		return (1);
	}
	fprintf(stdout,"output initialized\n");
	if(SDDS_DefineSimpleParameter(&data,"Undulator Length","meters", SDDS_DOUBLE)!=1 || SDDS_DefineSimpleParameter(&data,"Undulator Period Length","meters", SDDS_FLOAT)!=1 || SDDS_DefineSimpleParameter(&data,"Peak Field Variation","T", SDDS_LONG)!=1 ||
SDDS_DefineSimpleParameter(&data,"Number of Undulator Periods",NULL, SDDS_SHORT)!=1 ||
SDDS_DefineSimpleParameter(&data,"Test Case of Undulator","meters", SDDS_STRING)!=1)
	{
		return 1;
	}
	fprintf(stdout,"parameters defined \n");

	if(SDDS_DefineSimpleColumn(&data,"Number",NULL, SDDS_LONG)!=1 || SDDS_DefineSimpleColumn(&data,"s","s", SDDS_FLOAT)!=1 || SDDS_DefineSimpleColumn(&data,"element",NULL, SDDS_STRING)!=1 ||
SDDS_DefineSimpleColumn(&data,"x","m", SDDS_DOUBLE)!=1 ||SDDS_DefineSimpleColumn(&data,"y","m", SDDS_DOUBLE)!=1 ||SDDS_DefineSimpleColumn(&data,"z","m", SDDS_DOUBLE)!=1 ||SDDS_DefineSimpleColumn(&data,"px","gammabetax", SDDS_DOUBLE)!=1 ||SDDS_DefineSimpleColumn(&data,"py","gammabetay", SDDS_DOUBLE)!=1 ||SDDS_DefineSimpleColumn(&data,"pz","gammabetaz", SDDS_DOUBLE)!=1)
	{
		return 1;
	}
	fprintf(stdout,"Columns defined \n");

	if (SDDS_WriteLayout(&data)!=1)
	{
		return (1);
	}
	fprintf(stdout,"layout written \n");

	if (SDDS_StartPage(&data,5)!=1)
	{
		return (1);
	}

	if (SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Undulator Length","meters",1.5,"Undulator Period Length",0.45,"Peak Field Variation",0.1239999994,"Number of Undulator Periods",34,"Test Case of Undulator","Radiation From Electron passing through undulator",NULL)!=1)
	{
		return (1);
	}

	for(i=0;i<5;i++)
	{
		sprintf(buffer,"element %d",i);
		if(SDDS_SetRowValues(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,

						"s",(double)(1.0*i+0.1),
						"element",buffer,
						"x",(double)sin(i),
						"y",(double)cos(i),
						"z",(double)sin(i)*cos(i),
						"px",(double)(i+2.0),
						"py",(double)(i+3.0),
						"pz",(double)(i+4.0),
						NULL)!=1)
							{
								return (1);
							}



	}

	fprintf(stdout,"data filled in \n");

	if(SDDS_WritePage(&data)!=1)
	{
		return 1;
	}
	
	fprintf(stdout,"page written out \n");

	if(SDDS_StartPage(&data,6)!=1)
	{
		return 1;
	}

	if (SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Undulator Length","meters",1.5,"Undulator Period Length",0.45,"Peak Field Variation",0.1239999994,"Number of Undulator Periods",34,"Test Case of Undulator","Radiation From Electron passing through undulator",NULL)!=1)
	{
		return (1);
	}

	for(i=0;i<5;i++)
	{
		sprintf(buffer,"element %d",i);
		if(SDDS_SetRowValues(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,

						"s",(double)(1.0*i+0.1),
						"element",buffer,
						"x",(double)sin(i),
						"y",(double)cos(i),
						"z",(double)sin(i)*cos(i),
						"px",(double)(i+2.0),
						"py",(double)(i+3.0),
						"pz",(double)(i+4.0),
						NULL)!=1){
								return (1);
							}



	}

	fprintf(stdout,"data filled in \n");

	if(SDDS_WritePage(&data)!=1)
	{
		return 1;
	}
	
	fprintf(stdout,"page written out \n");
	if (SDDS_Terminate(&data)!=1)
	{
		return 1;
	}

	return 0;
}

int main(int argc, char **argv)
{
	if(WriteSDDSOutput() !=0)
	{
		SDDS_PrintErrors(stderr,SDDS_VERBOSE_PrintErrors);
		return 1;
	}
	if(ReadSDDSInput()!=0)
	{
		SDDS_PrintErrors(stderr,SDDS_VERBOSE_PrintErrors);
		return 1;
	}
	return 0;
}

