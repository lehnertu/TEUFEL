/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/*
 $Log: not supported by cvs2svn $
 Revision 1.18  2004/03/16 23:26:35  borland
 Added SCANITEMLIST_IGNORE_VALUELESS macro.

 Revision 1.17  2003/07/22 20:01:18  soliday
 Added support for Kylix.

 Revision 1.16  2002/08/14 15:40:17  soliday
 Added Open License

 Revision 1.15  2002/03/22 22:53:18  soliday
 Replaced free_scanargs with free_scanargs2

 Revision 1.14  2002/03/21 23:10:47  soliday
 Added free_scanargs2

 Revision 1.13  2002/03/07 01:18:53  soliday
 Added parse_string

 Revision 1.12  2002/01/28 16:49:49  soliday
 Added free_scanargs.

 Revision 1.11  2000/07/19 16:10:10  soliday
 Added ability to call from C++

 Revision 1.10  2000/04/11 16:19:32  soliday
 Modified prototypes to work with new mdbcommon library.

 Revision 1.9  2000/01/18 19:59:55  soliday
 Added support for ZLIB.

 Revision 1.8  1999/09/14 18:06:21  soliday
 Added export commands for WIN32 dll files.

 Revision 1.7  1999/07/22 16:22:06  soliday
 Added contains_keyword_phrase

 Revision 1.6  1996/05/29 21:44:48  borland
 Added mode flags for scanItemLists().

 * Revision 1.5  1996/02/14  01:02:08  borland
 * Added prototype for scanItemList().
 *
 * Revision 1.4  1996/01/21  00:15:54  borland
 * Added bit flag definitions and prototypes for new versions of unpacking
 * routines.
 *
 * Revision 1.3  1996/01/19  00:18:11  borland
 * SDDS.h: Added popenUsed to SDDS_LAYOUT structure.
 * scan.h: Added prototypes for unpacking functions.
 *
 * Revision 1.2  1995/09/05  21:15:39  saunders
 * First test release of the SDDS1.5 package.
 *
*/

/* define structure for use with scanargs(), scanlist() */
#if !defined(SCAN_INCLUDED)
#define SCAN_INCLUDED 1

#undef epicsShareFuncMDBLIB
#undef epicsShareFuncMDBCOMMON
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_MDBLIB)
#define epicsShareFuncMDBLIB  __declspec(dllexport)
#else
#define epicsShareFuncMDBLIB
#endif
#if defined(EXPORT_MDBCOMMON)
#define epicsShareFuncMDBCOMMON  __declspec(dllexport)
#else
#define epicsShareFuncMDBCOMMON
#endif
#else
#define epicsShareFuncMDBLIB
#define epicsShareFuncMDBCOMMON
#endif

#if defined(zLib)
#include "zlib.h"
#endif

#ifdef __cplusplus 
extern "C" {
#endif

typedef struct {
    long arg_type;	/* type of argument */
    long n_items;	/* number of items in list */
    char **list;	/* the list */
    } SCANNED_ARG;

/* possible values for arg_type */
#define OPTION 		1
#define A_LIST 		2

epicsShareFuncMDBCOMMON extern int scanargs(SCANNED_ARG **scanned, int argc, char **argv);
  /*
epicsShareFuncMDBCOMMON extern void free_scanargs(SCANNED_ARG *scanned, int argc);
  */
epicsShareFuncMDBCOMMON extern void free_scanargs(SCANNED_ARG **scanned, int argc);
epicsShareFuncMDBCOMMON extern int scanargsg(SCANNED_ARG **scanned, int argc, char **argv);
extern int parseList(char ***list, char *string);
extern int parse_string(char ***list, char *string);
epicsShareFuncMDBCOMMON long processPipeOption(char **item, long items, unsigned long *flags);
#define USE_STDIN      0x0001UL
#define USE_STDOUT     0x0002UL
#define DEFAULT_STDIN  0x0004UL
#define DEFAULT_STDOUT 0x0008UL

epicsShareFuncMDBCOMMON void processFilenames(char *programName, char **input, char **output, unsigned long pipeFlags, long noWarnings, long *tmpOutputUsed);


#include <stdio.h>
#define UNPACK_REQUIRE_SDDS 0x00000001UL
#define UNPACK_USE_PIPE     0x00000002UL
long PackSuffixType(char *filename, char **unpackedName,  unsigned long mode);
epicsShareFuncMDBLIB FILE *UnpackFopen(char *filename, unsigned long mode, short *popenUsed, char **tmpFileUsed);
#if defined(zLib)
epicsShareFuncMDBLIB gzFile UnpackGZipOpen(char *filename);
#endif

#include "SDDStypes.h"

epicsShareFuncMDBLIB extern long scanItemList(unsigned long *flags, char **item, long *items, unsigned long mode, ...);
epicsShareFuncMDBLIB extern long scanItemListLong(unsigned long long *flags, char **item, long *items, unsigned long mode, ...);
/* usage: scanItemList(&flags, item, &items, mode,
               <keyword>, <SDDS-type>, <pointer>, <number-required>, <set-flag>, etc.
               NULL)
 */
#define SCANITEMLIST_UNKNOWN_VALUE_OK    0x00000001UL
#define SCANITEMLIST_UNKNOWN_KEYVALUE_OK 0x00000002UL
#define SCANITEMLIST_REMOVE_USED_ITEMS   0x00000004UL
#define SCANITEMLIST_IGNORE_VALUELESS    0x00000008UL

epicsShareFuncMDBLIB extern long contains_keyword_phrase(char *string);
/* Obsolete: */
epicsShareFuncMDBLIB extern long scan_item_list(unsigned long *flags, char **item, long *items, ...);

/* usage: scan_item_list(&flags, item, &items, 
               <keyword>, <SDDS-type>, <pointer>, <number-required>, <set-flag>, etc.
               NULL)
 */

#ifdef __cplusplus
}
#endif

#endif

