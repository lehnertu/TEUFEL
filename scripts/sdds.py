import sddsdata, sys, time

class SDDS:
     """This class implements SDDS datasets."""

     def __init__(self, index):
          #define common SDDS definitions
          self.SDDS_VERBOSE_PrintErrors = 1
          self.SDDS_EXIT_PrintErrors = 2
          self.SDDS_CHECK_OKAY = 0
          self.SDDS_CHECK_NONEXISTENT = 1
          self.SDDS_CHECK_WRONGTYPE = 2
          self.SDDS_CHECK_WRONGUNITS = 3
          self.SDDS_DOUBLE = 1
          self.SDDS_REAL64 = 1
          self.SDDS_FLOAT = 2
          self.SDDS_REAL32 = 2
          self.SDDS_LONG = 3
          self.SDDS_INT32 = 3
          self.SDDS_ULONG = 4
          self.SDDS_UINT32 = 4
          self.SDDS_SHORT = 5
          self.SDDS_INT16 = 5
          self.SDDS_USHORT = 6
          self.SDDS_UINT16 = 6
          self.SDDS_STRING = 7
          self.SDDS_CHARACTER = 8
          self.SDDS_NUM_TYPES = 8
          self.SDDS_BINARY = 1
          self.SDDS_ASCII = 2
          #only indexes of 0 through 19 are allowed
          if index >= 0 and index < 20:
               self.index = index
          else:
               self.index = 0
          #initialize data storage variables
          self.description = ["", ""]
          self.parameterName = []
          self.columnName = []
          self.parameterDefinition = []
          self.columnDefinition = []
          self.parameterData = []
          self.columnData = []
          self.mode = self.SDDS_ASCII

     def load(self, input):
          """Load an SDDS file into an SDDS class."""

          try:
               #open SDDS file
               if sddsdata.InitializeInput(self.index, input) != 1:
                    time.sleep(2)
                    if sddsdata.InitializeInput(self.index, input) != 1:
                         sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

               #get data storage mode (SDDS_ASCII or SDDS_BINARY)
               self.mode = sddsdata.GetMode(self.index)

               #get description text and contents
               self.description = sddsdata.GetDescription(self.index)

               #get parameter names
               self.parameterName = sddsdata.GetParameterNames(self.index)
               numberOfParameters = len(self.parameterName)

               #get column names
               self.columnName = sddsdata.GetColumnNames(self.index)
               numberOfColumns = len(self.columnName)

               #get parameter definitions
               self.parameterDefinition = list(range(numberOfParameters))
               for i in range(numberOfParameters):
                    self.parameterDefinition[i] = sddsdata.GetParameterDefinition(self.index, self.parameterName[i])

               #get column definitions
               self.columnDefinition = list(range(numberOfColumns))
               for i in range(numberOfColumns):
                    self.columnDefinition[i] = sddsdata.GetColumnDefinition(self.index, self.columnName[i])

               #initialize parameter and column data
               self.parameterData = list(range(numberOfParameters))
               self.columnData = list(range(numberOfColumns))
               for i in range(numberOfParameters):
                    self.parameterData[i] = []
               for i in range(numberOfColumns):
                    self.columnData[i] = []

               #read in SDDS data
               page = sddsdata.ReadPage(self.index)
               if page != 1:
                    sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
               while page > 0:
                    for i in range(numberOfParameters):
                         self.parameterData[i].append(sddsdata.GetParameter(self.index,i))
                    for i in range(numberOfColumns):
                         self.columnData[i].append(sddsdata.GetColumn(self.index,i))
                    page = sddsdata.ReadPage(self.index)

               #close SDDS file
               if sddsdata.Terminate(self.index) != 1:
                    sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

          except:
               sddsdata.PrintErrors(self.SDDS_VERBOSE_PrintErrors)
               raise

     def listParameters(self):
        names = self.parameterName
        defs = self.parameterDefinition
        val = self.parameterData
        numberOfParameters = len(self.parameterName)
        for i in range(numberOfParameters):
            print "%20s" % names[i], "%20s" % str(val[i][0]), "%10s" % ("[%s]" % defs[i][1]), defs[i][2]

     def listColumns(self):
        names = self.columnName
        defs = self.columnDefinition
        numberOfColumns = len(self.columnName)
        for i in range(numberOfColumns):
            print "%20s" % names[i], "%10s" % ("[%s]" % defs[i][1]), defs[i][2]

     def getColumnData(self, name):
        """get the data of a named column in an SDDS file."""
        # see if we can find that column
        index = self.columnName.index(name)
        # print "index = ", index
        data = self.columnData[index]
        return data[0]

     def save(self, output):
          """Save an SDDS class to an SDDS file."""

          try:
               #check for invalid SDDS data
               numberOfParameters = len(self.parameterName)
               numberOfColumns = len(self.columnName)
               pages = 0
               if numberOfParameters != len(self.parameterData):
                    raise Exception("unmatched parameterName and parameterData")
               if numberOfColumns != len(self.columnData):
                    raise Exception("unmatched columnName and columnData")
               if numberOfParameters != len(self.parameterDefinition):
                    raise Exception("unmatched parameterName and parameterDefinition")
               if numberOfColumns != len(self.columnDefinition):
                    raise Exception("unmatched columnName and columnDefinition")
               if numberOfParameters > 0:
                    pages = len(self.parameterData[0])
               elif numberOfColumns > 0:
                    pages = len(self.columnData[0])
               for i in range(numberOfParameters):
                    if pages != len(self.parameterData[i]):
                         raise Exception("unequal number of pages in parameter data")
               for i in range(numberOfColumns):
                    if pages != len(self.columnData[i]):
                         raise Exception("unequal number of pages in column data")
               for page in range(pages):
                    rows = 0
                    if numberOfColumns > 0:
                         rows = len(self.columnData[0][page])
                    for i in range(numberOfColumns):
                         if rows != len(self.columnData[i][page]):
                              raise Exception("unequal number of rows in column data")

               #open SDDS output file
               if sddsdata.InitializeOutput(self.index, self.mode, 1, self.description[0], self.description[1], output) != 1:
                    sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

               #define parameters and columns
               for i in range(numberOfParameters):
                    if sddsdata.DefineParameter(self.index, self.parameterName[i],
                                            self.parameterDefinition[i][0],
                                            self.parameterDefinition[i][1],
                                            self.parameterDefinition[i][2],
                                            self.parameterDefinition[i][3],
                                            self.parameterDefinition[i][4],
                                            self.parameterDefinition[i][5]) == -1:
                         sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
               for i in range(numberOfColumns):
                    if sddsdata.DefineColumn(self.index, self.columnName[i],
                                            self.columnDefinition[i][0],
                                            self.columnDefinition[i][1],
                                            self.columnDefinition[i][2],
                                            self.columnDefinition[i][3],
                                            self.columnDefinition[i][4],
                                            self.columnDefinition[i][5]) == -1:
                         sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

               #write SDDS header
               if sddsdata.WriteLayout(self.index) != 1:
                    sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

               #write SDDS data
               for page in range(pages):
                    rows = 0
                    if numberOfColumns > 0:
                         rows = len(self.columnData[0][page])
                    if sddsdata.StartPage(self.index, rows) != 1:
                         sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
                    for i in range(numberOfParameters):
                         if sddsdata.SetParameter(self.index, i, self.parameterData[i][page]) != 1:
                              sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
                    for i in range(numberOfColumns):
                         if sddsdata.SetColumn(self.index, i, self.columnData[i][page]) != 1:
                              sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
                    if sddsdata.WritePage(self.index) != 1:
                         sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)

               #close SDDS output file
               if sddsdata.Terminate(self.index) != 1:
                    sddsdata.PrintErrors(self.SDDS_EXIT_PrintErrors)
          except:
               sddsdata.PrintErrors(self.SDDS_VERBOSE_PrintErrors)
               raise


     def setDescription(self, text, contents):
          self.description = [text, contents]

     def defineParameter(self, name, symbol, units, description, formatString, type, fixedValue):
          self.parameterName.append(name)
          self.parameterDefinition.append([symbol, units, description, formatString, type, fixedValue])
          self.parameterData.append([])

     def defineSimpleParameter(self, name, type):
          self.parameterName.append(name)
          self.parameterDefinition.append(["", "", "", "", type, ""])
          self.parameterData.append([])

     def defineColumn(self, name, symbol, units, description, formatString, type, fieldLength):
          self.columnName.append(name)
          self.columnDefinition.append([symbol, units, description, formatString, type, fieldLength])
          self.columnData.append([])

     def defineSimpleColumn(self, name, type):
          self.columnName.append(name)
          self.columnDefinition.append(["", "", "", "", type, 0])
          self.columnData.append([])

     def setParameterValueList(self, name, valueList):
          numberOfParameters = len(self.parameterName)
          for i in range(numberOfParameters):
               if self.parameterName[i] == name:
                    self.parameterData[i] = valueList
                    return
          msg = "invalid parameter name " + name
          raise Exception(msg)

     def setParameterValue(self, name, value, page):
          page = page - 1
          numberOfParameters = len(self.parameterName)
          for i in range(numberOfParameters):
               if self.parameterName[i] == name:
                    if len(self.parameterData[i]) == page:
                         self.parameterData[i][page:] = [value]
                    elif len(self.parameterData[i]) < page or page < 0:
                         msg = "invalid page " + str(page+1)
                         raise Exception(msg)
                    else:
                         self.parameterData[i][page] = [value]
                    return
          msg = "invalid parameter name " + name
          raise Exception(msg)

     def setColumnValueLists(self, name, valueList):
          numberOfColumns = len(self.columnName)
          for i in range(numberOfColumns):
               if self.columnName[i] == name:
                    self.columnData[i] = valueList
                    return
          msg = "invalid column name " + name
          raise Exception(msg)

     def setColumnValueList(self, name, valueList, page):
          page = page - 1
          numberOfColumns = len(self.columnName)
          for i in range(numberOfColumns):
               if self.columnName[i] == name:
                    if len(self.columnData[i]) == page:
                         self.columnData[i][page:] = [valueList]
                    elif len(self.columnData[i]) < page or page < 0:
                         msg = "invalid page " + str(page+1)
                         raise Exception(msg)
                    else:
                         self.columnData[i][page] = [valueList]
                    return
          msg = "invalid column name " + name
          raise Exception(msg)

     def setColumnValue(self, name, value, page, row):
          page = page - 1
          row = row - 1
          numberOfColumns = len(self.columnName)
          for i in range(numberOfColumns):
               if self.columnName[i] == name:
                    if len(self.columnData[i]) == page:
                         if row == 0:
                              self.columnData[i][page:] = [[value]]
                         else:
                              msg = "invalid row " + str(row+1)
                              raise Exception(msg)
                    elif len(self.columnData[i]) < page or page < 0:
                         msg = "invalid page " + str(page+1)
                         raise Exception(msg)
                    else:
                         if len(self.columnData[i][page]) == row:
                              self.columnData[i][page][row:] = [value]
                         elif len(self.columnData[i][page]) < row or row < 0:
                              msg = "invalid row " + str(row+1)
                              raise Exception(msg)
                         else:
                              self.columnData[i][page][row] = [value]
                    return
          msg = "invalid column name " + name
          raise Exception(msg)

def demo(output):
     """Save an demo SDDS file using the SDDS class."""

     x = SDDS(0)
     x.description[0] = "text"
     x.description[1] = "contents"
     x.parameterName = ["ShortP", "LongP", "FloatP", "DoubleP", "StringP", "CharacterP"]
     x.parameterData = [[1, 6], [2, 7], [3.3, 8.8], [4.4, 9.8], ["five", "ten"], ["a", "b"]]
     x.parameterDefinition = [["","","","",x.SDDS_SHORT,""],
                                 ["","","","",x.SDDS_LONG,""],
                                 ["","","","",x.SDDS_FLOAT,""],
                                 ["","","","",x.SDDS_DOUBLE,""],
                                 ["","","","",x.SDDS_STRING,""],
                                 ["","","","",x.SDDS_CHARACTER, ""]]
     x.columnName = ["ShortC", "LongC", "FloatC", "DoubleC", "StringC", "CharacterC"]
     x.columnData = [[[1, 2, 3], [-1, -2, -3.6, -4.4]],
                        [[1, 2, 3], [-1, -2, -3.6, -4.4]],
                        [[1, 2, 3], [-1, -2, -3.6, -4.4]],
                        [[1, 2, 3], [-1, -2, -3.6, -4.4]],
                        [["row 1", "row 2", "row 3"], ["row 1", "row 2", "row 3", "row 4"]],
                        [["x", "y", "z"], ["i", "j", "k", "l"]]]
     x.columnDefinition = [["","","","",x.SDDS_SHORT,0],
                              ["","","","",x.SDDS_LONG,0],
                              ["","","","",x.SDDS_FLOAT,0],
                              ["","","","",x.SDDS_DOUBLE,0],
                              ["","","","",x.SDDS_STRING,0],
                              ["","","","",x.SDDS_CHARACTER,0]]
     x.save(output)
     del x

def demo2(output):
     """Save an demo SDDS file using the SDDS class."""

     x = SDDS(0)
     x.setDescription("text", "contents")
     names = ["Short", "Long", "Float", "Double", "String", "Character"]
     types = [x.SDDS_SHORT, x.SDDS_LONG, x.SDDS_FLOAT, x.SDDS_DOUBLE, x.SDDS_STRING, x.SDDS_CHARACTER]
     for i in range(6):
          x.defineSimpleParameter(names[i] + "P", types[i])
          x.defineSimpleColumn(names[i] + "C", types[i])
     parameterData = [[1, 6], [2, 7], [3.3, 8.8], [4.4, 9.8], ["five", "ten"], ["a", "b"]]
     for i in range(6):
          x.setParameterValueList(names[i] + "P", parameterData[i])
     columnData = [[[1, 2, 3], [-1, -2, -3.6, -4.4]],
                   [[1, 2, 3], [-1, -2, -3.6, -4.4]],
                   [[1, 2, 3], [-1, -2, -3.6, -4.4]],
                   [[1, 2, 3], [-1, -2, -3.6, -4.4]],
                   [["row 1", "row 2", "row 3"], ["row 1", "row 2", "row 3", "row 4"]],
                   [["x", "y", "z"], ["i", "j", "k", "l"]]]
     for i in range(6):
          x.setColumnValueLists(names[i] + "C", columnData[i])
     x.save(output)
