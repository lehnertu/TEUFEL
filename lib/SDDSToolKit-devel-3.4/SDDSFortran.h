C     FORTRAN
C     The following functions do not return any values and thus do not need to be declared.
C     SDDS_SetTerminateMode_F
C     SDDS_DeferSavingLayout_F
C     SDDS_ClearErrors_F
C     SDDS_SetError_F
C     SDDS_Bomb_F
C     SDDS_Warning_F
C     SDDS_RegisterProgramName_F
C     SDDS_PrintErrors_F

C     The following functions return character arrays. 
C     The size of the array may need to be changed if it is too small
      character*100 SDDS_GetDescriptionText_F
      character*100 SDDS_GetDescriptionContents_F
      character*10 SDDS_GetTypeName_F
      character*50 SDDS_GetColumnNameFromIndex_F
      character*50 SDDS_GetParNameFromIndex_F

C     The following functions all return integers.
      integer SDDS_InitializeAppend_F
      integer SDDS_InitializeOutput_F
      integer SDDS_Terminate_F
      integer SDDS_DefineParameter_F
      integer SDDS_DefineParameter1_F
      integer SDDS_DefineColumn_F
      integer SDDS_IsValidName_F
      integer SDDS_SetNameValidityFlags_F
      integer SDDS_DefineSimpleColumn_F
      integer SDDS_DefineSimpleParameter_F
      integer SDDS_WriteLayout_F
      integer SDDS_EraseData_F
      integer SDDS_ProcessColumnString_F
      integer SDDS_ProcessParameterString_F
      integer SDDS_InitializeCopy_F
      integer SDDS_CopyLayout_F
      integer SDDS_AppendLayout_F
      integer SDDS_CopyPage_F
      integer SDDS_CopyParameters_F
      integer SDDS_CopyColumns_F
      integer SDDS_CopyRow_F
      integer SDDS_CopyRowDirect_F
      integer SDDS_CopyAdditionalRows_F
      integer SDDS_SaveLayout_F
      integer SDDS_RestoreLayout_F
      integer SDDS_StartPage_F
      integer SDDS_ClearPage_F
      integer SDDS_LengthenTable_F
      integer SDDS_WritePage_F
      integer SDDS_InitializeInput_F
      integer SDDS_InitHeaderlessInput_F
      integer SDDS_ReadPage_F
      integer SDDS_ReadPageSparse_F
      integer SDDS_RowCount_F
      integer SDDS_SetColumnFlags_F
      integer SDDS_SetRowFlags_F
      integer SDDS_GetRowFlag_F
      integer SDDS_DeleteColumn_F
      integer SDDS_DeleteParameter_F
      integer SDDS_DeleteUnsetColumns_F
      integer SDDS_DeleteUnsetRows_F
      integer SDDS_ColumnCount_F
      integer SDDS_ParameterCount_F
      integer SDDS_GetColumn_F
      integer SDDS_GetInternalColumn_F
      integer SDDS_GetColumnInDoubles_F
      integer SDDS_GetColumnInLong_F
      integer SDDS_GetNumericColumn_F
      integer SDDS_GetValue_F
      integer SDDS_GetValueByIndex_F 
      integer SDDS_GetValueByAbsIndex_F
      integer SDDS_GetParameter_F
      integer SDDS_GetParameterAsDouble_F
      integer SDDS_GetFixedValueParameter_F
      integer SDDS_NumberOfErrors_F
      integer SDDS_TransferColumnDefinition_F
      integer SDDS_DefineColumnLikePar_F
      integer SDDS_TransferParDefinition_F
      integer SDDS_DefineParLikeColumn_F
      integer SDDS_TransferAllParDefs_F
      integer SDDS_GetColumnIndex_F
      integer SDDS_GetParameterIndex_F
      integer SDDS_GetColumnType_F
      integer SDDS_GetNamedColumnType_F
      integer SDDS_GetParameterType_F
      integer SDDS_GetNamedParameterType_F
      integer SDDS_GetTypeSize_F
      integer SDDS_IdentifyType_F
      integer SDDS_CheckColumn_F
      integer SDDS_CheckParameter_F
      integer SDDS_HasWhitespace_F
      integer SDDS_SprintTypedValue_F
      integer SDDS_StringIsBlank_F
      integer SDDS_ApplyFactorToParameter_F
      integer SDDS_ApplyFactorToColumn_F
      integer SDDS_DeleteParFixedValues_F
      integer SDDS_SetDataMode_F
      integer SDDS_CheckDataset_F
      integer SDDS_SetAutoCheckMode_F
      integer SDDS_SetParameterByName_F
      integer SDDS_Set2ParametersByName_F
      integer SDDS_Set3ParametersByName_F
      integer SDDS_Set4ParametersByName_F
      integer SDDS_Set5ParametersByName_F
      integer SDDS_Set6ParametersByName_F
      integer SDDS_SetParameterStringByName_F
      integer SDDS_SetParameterByIndex_F
      integer SDDS_SetColumnByName_F
      integer SDDS_SetColumnByIndex_F
      integer SDDS_GetValuesFrom_C_To_F
      integer SDDS_GetVoidFrom_C_To_F
      integer SDDS_GetParameterData_F
      integer SDDS_GetColumnData_F
      integer SDDS_SetDoubleRowValues_F
      integer SDDS_Set3DoubleRowValues_F
      integer SDDS_Set4DoubleRowValues_F
      integer SDDS_Set5DoubleRowValues_F
C     Here I have declared some variables that are 
C     constants in the SDDS library.
      integer SDDS_VERBOSE_PrintErrors
      integer SDDS_CHECK_NONEXISTENT
      integer SDDS_CHECK_WRONGUNITS
      integer SDDS_DOUBLE
      integer SDDS_FLOAT
      integer SDDS_LONG
      integer SDDS_ULONG
      integer SDDS_SHORT
      integer SDDS_USHORT
      integer SDDS_STRING
      integer SDDS_CHARACTER
      integer SDDS_NUM_TYPES
      integer SDDS_BINARY
      integer SDDS_ASCII
      parameter (SDDS_VERBOSE_PrintErrors = 1)
      parameter (SDDS_CHECK_NONEXISTENT = 1)
      parameter (SDDS_CHECK_WRONGUNITS = 3)
      parameter (SDDS_DOUBLE = 1)
      parameter (SDDS_FLOAT = 2)
      parameter (SDDS_LONG = 3)
      parameter (SDDS_ULONG = 4)
      parameter (SDDS_SHORT = 5)
      parameter (SDDS_USHORT = 6)
      parameter (SDDS_STRING = 7)
      parameter (SDDS_CHARACTER = 8)
      parameter (SDDS_NUM_TYPES = 8)
      parameter (SDDS_BINARY = 1)
      parameter (SDDS_ASCII = 2)
