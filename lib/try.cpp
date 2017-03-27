#include <math.h>
#include "SDDS.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;


int WriteSDDSOutput();
int ReadSDDSInput();



int WriteSDDSOutput()
{
	SDDS_DATASET data;
	int32_t i;
	char buffer[100];
	if (SDDS_InitializeOutput(&data,SDDS_ASCII,1,NULL,NULL,"try.sdds")!=1)
	{
		return (1);
	}
	fprintf(stdout,"output initialized\n");
	if(SDDS_DefineSimpleParameter(&data,"UndulatorLength","meters", SDDS_DOUBLE)!=1 || SDDS_DefineSimpleParameter(&data,"UndulatorPeriodLength","meters", SDDS_FLOAT)!=1 || SDDS_DefineSimpleParameter(&data,"PeakFieldVariation","T", SDDS_LONG)!=1 ||
SDDS_DefineSimpleParameter(&data,"NumberofUndulatorPeriods",NULL, SDDS_SHORT)!=1 ||
SDDS_DefineSimpleParameter(&data,"TestCaseofUndulator","meters", SDDS_STRING)!=1)
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
	
	fprintf(stdout,"Creating Page\n");
	if (SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "UndulatorLength",1.5,"UndulatorPeriodLength",0.45,"PeakFieldVariation",0.123,"NumberofUndulatorPeriods",34.0,"TestCaseofUndulator","RadiationFromElectronpassingthroughundulator",NULL)!=1)
	{
		return (1);
	}
	
	fprintf(stdout,"Parameters set \n");
	for(i=0;i<5;i++)
	{
		sprintf(buffer,"element %d",i);
		if(SDDS_SetRowValues(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,
						"Number",(long)(i),
						"s",(double)(1.0*i+0.1),
						"element",buffer,
						"x",(double)(i*3.0/2.0),
						"y",(double)(2*i),
						"z",(double)(100),
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

	if(SDDS_StartPage(&data,5)!=1)
	{
		return 1;
	}

	if (SDDS_SetParameters(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "UndulatorLength",1.5,"UndulatorPeriodLength",0.45,"PeakFieldVariation",0.1239999994,"NumberofUndulatorPeriods",34.0,"TestCaseofUndulator","RadiationFromElectronpassingthroughundulator",NULL)!=1)
	{
		return (1);
	}

	for(i=0;i<5;i++)
	{
		sprintf(buffer,"element %d",i);
		if(SDDS_SetRowValues(&data,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,i,
						"Number",(long)(i),
						"s",(double)(1.0*i+0.1),
						"element",buffer,
						"x",(double)(i*3.0/2.0),
						"y",(double)(2*i),
						"z",(double)(100),
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
	if (SDDS_Terminate(&data)!=1)
	{
		return 1;
	}

	return 0;
}


int ReadSDDSInput()
{
	SDDS_DATASET data;
	int32_t i;
	char buffer[100];
	int32_t numberOfParameters;
	int32_t numberOfColumns;
	char **parameterNames;
	char **columnNames;
	int32_t parameterDataType;
	int32_t columnDataType;
	void *parameterData;
	void *columnData;
	long page,rows,row;
	if (SDDS_InitializeInput(&data, "try.sdds") != 1)
	{
		return (1);
	}
	fprintf(stdout,"input initialized\n");
	parameterNames = SDDS_GetParameterNames(&data,&numberOfParameters);
	columnNames = SDDS_GetColumnNames(&data,&numberOfColumns);
	page = SDDS_ReadPage(&data);
	while(page>=1)
	{
		fprintf(stdout,"page %1d\n",page);
		for(i=0;i<numberOfParameters;i++)
		{
			parameterDataType=SDDS_GetParameterType(&data,i);
			parameterData = SDDS_GetParameterByIndex(&data,i,NULL);
			fprintf(stdout,"parameter %s=\n",parameterNames[i]);
			switch(parameterDataType)
			{
				case SDDS_DOUBLE:
					fprintf(stdout,"%21.15e\n",*((double*)parameterData));
					break;
				case SDDS_FLOAT:
					fprintf(stdout,"%15.8e\n",*((float*)parameterData));
					break;
				case SDDS_ULONG:
					fprintf(stdout,"%"PRIu32"\n",*((uint32_t*)parameterData));
					break;
				case SDDS_LONG:
					fprintf(stdout,"%"PRId32"\n",*((int32_t*)parameterData));
					break;
				case SDDS_USHORT:
					fprintf(stdout,"%hu\n",*((unsigned short*)parameterData));
					break;
				case SDDS_SHORT:
					fprintf(stdout,"%hd\n",*((short*)parameterData));
					break;
				case SDDS_CHARACTER:
					fprintf(stdout,"%c\n",*((char*)parameterData));
					break;
				case SDDS_STRING:
					fprintf(stdout,"%s\n",*((char**)parameterData));
					break;
				default:
					fprintf(stdout,"\n");	


			}
			if (parameterDataType == SDDS_STRING)
			{
				if(*((char**)parameterData)) free(*((char**)parameterData));
			}
			if(parameterData) free(parameterData);
		}

		rows = SDDS_RowCount(&data);
		for(i=0;i<numberOfColumns;i++)
		{
			columnDataType=SDDS_GetColumnType(&data,i);
			columnData = SDDS_GetColumn(&data,columnNames[i]);
			fprintf(stdout,"column %s=\n",columnNames[i]);
			switch(columnDataType)
			{
				case SDDS_DOUBLE:
				    for(row=0;row<rows;row++){
					fprintf(stdout,"%21.15e\n",((double*)columnData)[row]);
					}
					break;
				case SDDS_FLOAT:
				    for(row=0;row<rows;row++){
					fprintf(stdout,"%15.8e\n",((float*)columnData)[row]);
					}
					break;
				case SDDS_ULONG:
				     for(row=0;row<rows;row++){
					fprintf(stdout,"%"PRIu32"\n",((uint32_t*)columnData)[row]);
					}
					break;
				case SDDS_LONG:
				     for(row=0;row<rows;row++){
					fprintf(stdout,"%"PRId32"\n",((int32_t*)columnData)[row]);
					}
					break;
				case SDDS_USHORT:
				     for(row=0;row<rows;row++){
					fprintf(stdout,"%hu\n",((unsigned short*)columnData)[row]);
					}
					break;
				case SDDS_SHORT:
				     for(row=0;row<rows;row++){
					fprintf(stdout,"%hd\n",((short*)columnData)[row]);
					}
					break;
				case SDDS_CHARACTER:
				     for(row=0;row<rows;row++){
					fprintf(stdout,"%c\n",((char*)columnData)[row]);
					}
					break;
				case SDDS_STRING:
				     for(row=0;row<rows;row++){
					fprintf(stdout,"%s\n",((char**)columnData)[row]);
					}
					break;
				default:
					fprintf(stdout,"\n");	
			}
			if (columnDataType == SDDS_STRING)
			{
			    for(row=0;row<rows;row++){
				if(((char**)columnData)[row]) free(((char**)columnData)[row]);
				}
			}
			if(columnData) free(columnData);
		}

		page = SDDS_ReadPage(&data);
	}

	for (i=0; i<numberOfParameters;i++){
		free(parameterNames[i]);
		}
	if(parameterNames) free (parameterNames);

	for (i=0; i<numberOfColumns;i++){
		free(columnNames[i]);
		}
	if(columnNames) free (columnNames);
	if(SDDS_Terminate(&data)!=1){
		return(1);
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

