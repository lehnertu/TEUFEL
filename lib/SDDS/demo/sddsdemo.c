#include "mdb.h"
#include "SDDS.h"

int WriteSDDSOutput();
int ReadSDDSInput();

int main(int argc, char **argv)
{
  if (WriteSDDSOutput() != 0) {
  /* SDDS_PrintErrors
      arguments:
       1) FILE stream such as stderr or stdout
       2) mode
        SDDS_VERBOSE_PrintErrors - print all errors
        SDDS_EXIT_PrintErrors - exit after printing errors
  */
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  if (ReadSDDSInput() != 0) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    return(1);
  }
  return(0);
}

int WriteSDDSOutput() {
  SDDS_DATASET SDDS_dataset;
  int32_t i;
  char buffer[100];

  /* SDDS_InitializeOutput
      arguments:
       1) *SDDS_DATASET
       2) data mode = SDDS_ASCII or SDDS_BINARY
       3) lines per row of ascii data
       4) description string
       5) contents string
       6) output file name
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_InitializeOutput(&SDDS_dataset, SDDS_BINARY, 1,
                            NULL, NULL, "demo.sdds" ) != 1) {
    return(1);
  }
  fprintf(stdout, "output initialized\n");

  /* SDDS_DefineSimpleParameter
      arguments:
       1) *SDDS_DATASET
       2) parameter name
       3) parameter units
       4) SDDS datatype
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_DefineSimpleParameter(&SDDS_dataset, "adouble", "meters",  SDDS_DOUBLE)!=1 ||
      SDDS_DefineSimpleParameter(&SDDS_dataset, "afloat",  "C",       SDDS_FLOAT)!=1 ||
      SDDS_DefineSimpleParameter(&SDDS_dataset, "along",   "seconds", SDDS_LONG)!=1 ||
      SDDS_DefineSimpleParameter(&SDDS_dataset, "ashort",  NULL,      SDDS_SHORT)!=1 ||
      SDDS_DefineSimpleParameter(&SDDS_dataset, "astring", NULL,      SDDS_STRING)!=1) {
    return(1);
  }
  fprintf(stdout, "parameters defined\n");

  /* SDDS_DefineSimpleColumn
      arguments:
       1) *SDDS_DATASET
       2) column name
       3) column units
       4) SDDS datatype
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_DefineSimpleColumn(&SDDS_dataset, "number",  NULL, SDDS_LONG)!=1 ||
      SDDS_DefineSimpleColumn(&SDDS_dataset, "s",       "s",  SDDS_DOUBLE)!=1 ||
      SDDS_DefineSimpleColumn(&SDDS_dataset, "element", NULL, SDDS_STRING)!=1 ||
      SDDS_DefineSimpleColumn(&SDDS_dataset, "x",       "x",  SDDS_DOUBLE)!=1 ||
      SDDS_DefineSimpleColumn(&SDDS_dataset, "xp",      "x'", SDDS_DOUBLE)!=1 ||
      SDDS_DefineSimpleColumn(&SDDS_dataset, "y",       "y",  SDDS_DOUBLE)!=1 ||
      SDDS_DefineSimpleColumn(&SDDS_dataset, "yp",      "y'", SDDS_DOUBLE)!=1) {
    return(1);
  }
  fprintf(stdout, "columns defined\n");

  /* SDDS_WriteLayout
      arguments:
       1) *SDDS_DATASET
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_WriteLayout(&SDDS_dataset)!=1) {
    return(1);
  }
  fprintf(stdout, "layout written\n");

  /* SDDS_StartPage
      arguments:
       1) *SDDS_DATASET
       2) expected number of rows
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_StartPage(&SDDS_dataset, 5)!=1) {
    return(1);
  }

  /* SDDS_SetParameters
      arguments:
       1) *SDDS_DATASET
       2) mode = SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE
        3) char *name1, value1, char *name2, value2, ..., NULL
       2) mode = SDDS_SET_BY_NAME | SDDS_PASS_BY_REFERENCE
        3) char *name1, void *data1, char *name2, void *data2, ..., NULL
       2) mode = SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE
        3) long index1, value1, long index2, value2, ..., -1
       2) mode = SDDS_SET_BY_INDEX | SDDS_PASS_BY_REFERENCE
        3) long index1, void *data1, long index2, void *data2, ..., -1
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_SetParameters(&SDDS_dataset, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                         "adouble", 250.0, 
                         "afloat",  250.1, 
                         "along",   1234567891, 
                         "ashort",  12345, 
                         "astring", "This is a string", 
                         NULL)!=1) {
    return(1);
  }

  /* SDDS_SetRowValues
      arguments:
       1) *SDDS_DATASET
       2) mode = SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE
        3) row
        4) char *name1, value1, char *name2, value2, ..., NULL
       2) mode = SDDS_SET_BY_NAME | SDDS_PASS_BY_REFERENCE
        3) row
        4) char *name1, void *data1, char *name2, void *data2, ..., NULL
       2) mode = SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE
        3) row
        4) long index1, value1, long index2, value2, ..., -1
       2) mode = SDDS_SET_BY_INDEX | SDDS_PASS_BY_REFERENCE
        3) row
        4) long index1, void *data1, long index2, void *data2, ..., -1
      return:
       1 on success
       0 on failure 
  */
  for (i=0; i<5; i++) {
    sprintf(buffer, "element %d", i);
    if (SDDS_SetRowValues(&SDDS_dataset, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, i, 
                           "number", i,
                           "x",      (double)(1.0*i+0.1), 
                           "xp",     (double)(1.0*i+0.2), 
                           "y",      (double)(1.0*i+0.3), 
                           "yp",     (double)(1.0*i+0.4), 
                           "s",      (double)(1.0*i), 
                           "element", buffer, 
                           NULL)!=1) {
      return(1);
    }
  }
  fprintf(stdout, "data filled in\n");

  /* SDDS_WritePage
      arguments:
       1) *SDDS_DATASET
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_WritePage(&SDDS_dataset)!=1) {
    return(1);
  }
  fprintf(stdout,"page written out\n");

  if (SDDS_StartPage(&SDDS_dataset, 6)!=1) {
    return(1);
  }
  if (SDDS_SetParameters(&SDDS_dataset, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                         "adouble", 451.0, 
                         "afloat",  451.1, 
                         "along",   234567891, 
                         "ashort",  2345, 
                         "astring", "this is a string", 
                         NULL)!=1) {
    return(1);
  }
  for (i=0; i<6; i++) {
    sprintf(buffer, "element %d", (i+50));
    if (SDDS_SetRowValues(&SDDS_dataset, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, i, 
                           "number", (i+10),
                           "x",      (double)(1.0*(i+50)+0.1), 
                           "xp",     (double)(1.0*(i+50)+0.2), 
                           "y",      (double)(1.0*(i+50)+0.3), 
                           "yp",     (double)(1.0*(i+50)+0.4), 
                           "s",      (double)(1.0*(i+50)), 
                           "element", buffer, 
                           NULL)!=1) {
      return(1);
    }
  }
  fprintf(stdout, "data filled in\n");
  if (SDDS_WritePage(&SDDS_dataset)!=1) {
    return(1);
  }
  fprintf(stdout, "page written out\n");

  /* SDDS_Terminate
      arguments:
       1) *SDDS_DATASET
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_Terminate(&SDDS_dataset)!=1) {
    return(1);
  }
  return(0);
}



int ReadSDDSInput() {
  SDDS_DATASET SDDS_dataset;
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
  long page, rows, row;

  /* SDDS_InitializeInput
      arguments:
       1) *SDDS_DATASET
       2) input file name
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_InitializeInput(&SDDS_dataset, "demo.sdds" ) != 1) {
    return(1);
  }
  fprintf(stdout, "input initialized\n");

  /* SDDS_GetParameterNames
      arguments:
       1) *SDDS_DATASET
       2) int32_t address location to place the number of parameters
      return:
       char** of parameter names on success
       NULL on failure
  */
  parameterNames = SDDS_GetParameterNames(&SDDS_dataset, &numberOfParameters);

  /* SDDS_GetColumnNames
      arguments:
       1) *SDDS_DATASET
       2) int32_t address location to place the number of columns
      return:
       char** of column names on success
       NULL on failure
  */
  columnNames = SDDS_GetColumnNames(&SDDS_dataset, &numberOfColumns);

  /* SDDS_ReadPage
      arguments:
       1) *SDDS_DATASET
      return:
       page number starting with 1 on success
       -1 if no more pages to read
       0 on failure
  */
  page=SDDS_ReadPage(&SDDS_dataset);
  while (page >= 1) {
    fprintf(stdout, "page %ld\n", page);

    for (i = 0; i < numberOfParameters; i++) {

  /* SDDS_GetParameterType
      arguments:
       1) *SDDS_DATASET
       2) parameter index
      return:
       data type on success
       -1 on failure
  */
      parameterDataType = SDDS_GetParameterType(&SDDS_dataset, i);

  /* SDDS_GetParameterByIndex
      arguments:
       1) *SDDS_DATASET
       2) parameter index
       3) pre-allocated memory storage or NULL
      return:
       pointer to data on success
       NULL on failure
  */
      parameterData = SDDS_GetParameterByIndex(&SDDS_dataset, i, NULL);
      fprintf(stdout, "parameter %s=", parameterNames[i]);
      
      switch (parameterDataType) {
      case SDDS_DOUBLE:
        fprintf(stdout, "%21.15e\n", *((double*)parameterData));
        break;
      case SDDS_FLOAT:
        fprintf(stdout, "%15.8e\n", *((float*)parameterData));
        break;
      case SDDS_ULONG:
        fprintf(stdout, "%" PRIu32 "\n", *((uint32_t*)parameterData));
        break;
      case SDDS_LONG:
        fprintf(stdout, "%" PRId32 "\n", *((int32_t*)parameterData));
        break;
      case SDDS_USHORT:
        fprintf(stdout, "%hu\n", *((unsigned short*)parameterData));
        break;
      case SDDS_SHORT:
        fprintf(stdout, "%hd\n", *((short*)parameterData));
        break;
       case SDDS_CHARACTER:
        fprintf(stdout, "%c\n", *((char*)parameterData));
        break;
       case SDDS_STRING:
        fprintf(stdout, "%s\n", *((char**)parameterData));
        break;
     default:
        fprintf(stdout, "\n");
      }

      if (parameterDataType == SDDS_STRING) {
        if (*((char**)parameterData)) free(*((char**)parameterData));
      }
      if (parameterData) free(parameterData);
    }

  /* SDDS_RowCount
      arguments:
       1) *SDDS_DATASET
      return:
       number of rows in page
  */
    rows = SDDS_RowCount(&SDDS_dataset);
    for (i = 0; i < numberOfColumns; i++) {

  /* SDDS_GetColumnType
      arguments:
       1) *SDDS_DATASET
       2) column index
      return:
       data type on success
       -1 on failure
  */
      columnDataType = SDDS_GetColumnType(&SDDS_dataset, i);

  /* SDDS_GetColumn
      arguments:
       1) *SDDS_DATASET
       2) parameter name
      return:
       pointer to data on success
       NULL on failure
  */
      columnData = SDDS_GetColumn(&SDDS_dataset, columnNames[i]);
      fprintf(stdout, "column %s=\n", columnNames[i]);
      
      switch (columnDataType) {
      case SDDS_DOUBLE:
        for (row = 0; row < rows; row++) {
          fprintf(stdout, "%21.15e\n", ((double*)columnData)[row]);
        }
        break;
      case SDDS_FLOAT:
        for (row = 0; row < rows; row++) {
          fprintf(stdout, "%15.8e\n", ((float*)columnData)[row]);
        }
        break;
      case SDDS_ULONG:
        for (row = 0; row < rows; row++) {
          fprintf(stdout, "%" PRIu32 "\n", ((uint32_t*)columnData)[row]);
        }
        break;
      case SDDS_LONG:
        for (row = 0; row < rows; row++) {
          fprintf(stdout, "%" PRId32 "\n", ((int32_t*)columnData)[row]);
        }
        break;
      case SDDS_USHORT:
        for (row = 0; row < rows; row++) {
          fprintf(stdout, "%hu\n", ((unsigned short*)columnData)[row]);
        }
        break;
      case SDDS_SHORT:
        for (row = 0; row < rows; row++) {
          fprintf(stdout, "%hd\n", ((short*)columnData)[row]);
        }
        break;
      case SDDS_CHARACTER:
        for (row = 0; row < rows; row++) {
          fprintf(stdout, "%c\n", ((char*)columnData)[row]);
        }
        break;
      case SDDS_STRING:
        for (row = 0; row < rows; row++) {
          fprintf(stdout, "%s\n", ((char**)columnData)[row]);
        }
        break;
      default:
        fprintf(stdout, "\n");
      }

      if (columnDataType == SDDS_STRING) {
        for (row = 0; row < rows; row++) {
          if (((char**)columnData)[row]) free(((char**)columnData)[row]);
        }
      }
      if (columnData) free(columnData);
    }
    page=SDDS_ReadPage(&SDDS_dataset);
  }

  for (i = 0 ; i < numberOfParameters; i++) {
    free(parameterNames[i]);
  }
  if (parameterNames) free(parameterNames);

  for (i = 0 ; i < numberOfColumns; i++) {
    free(columnNames[i]);
  }
  if (columnNames) free(columnNames);

  if (SDDS_Terminate(&SDDS_dataset)!=1) {
    return(1);
  }

  return(0);
}
