/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: match_string.h
 * purpose: flag definitions for use with match_string
 *
 * Michael Borland, 1988
 $Log: not supported by cvs2svn $
 Revision 1.8  2004/12/02 23:02:00  soliday
 The section of code inside the _MATCH_STRING_ ifdef in both the match_string.h
 and mdb.h is now the same. So it will not matter which one is included first.

 Revision 1.7  2003/12/02 20:29:00  soliday
 Exported all the procedures.

 Revision 1.6  2003/07/22 20:01:17  soliday
 Added support for Kylix.

 Revision 1.5  2002/08/14 15:40:15  soliday
 Added Open License

 Revision 1.4  2002/07/24 21:19:28  borland
 Added prototypes for strncmp_case_insensitive() and strcmp_case_insensitive().

 Revision 1.3  1999/09/14 18:04:38  soliday
 Added export commands for WIN32 dll files.

 Revision 1.2  1995/09/05 21:15:10  saunders
 First test release of the SDDS1.5 package.

 */

#undef epicsShareFuncMDBLIB
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_MDBLIB)
#define epicsShareFuncMDBLIB  __declspec(dllexport)
#else
#define epicsShareFuncMDBLIB
#endif
#else
#define epicsShareFuncMDBLIB
#endif

#if !defined(_MATCH_STRING_)
epicsShareFuncMDBLIB extern long match_string(char *string, char **option_list, long n_options,
                        long match_mode_flags);
epicsShareFuncMDBLIB extern int strncmp_case_insensitive(char *s1, char *s2, long n);
epicsShareFuncMDBLIB extern int strcmp_case_insensitive(char *s1, char *s2);
#if defined(_WIN32)
#if defined(__BORLANDC__)
#define strcasecmp(s, t) stricmp(s, t)
#define strncasecmp(s, t, n) strnicmp(s, t, n)
#else
#define strcasecmp(s, t) _stricmp(s, t)
#define strncasecmp(s, t, n) _strnicmp(s, t, n)
#endif
#endif
#define _MATCH_STRING_ 1

#define DCL_STYLE_MATCH 0
#define UNIQUE_MATCH DCL_STYLE_MATCH
#define CASE_SENSITIVE 1
#define MATCH_WHOLE_STRING 2
#define RETURN_FIRST_MATCH 8
#define WILDCARD_MATCH 16
#define EXACT_MATCH (CASE_SENSITIVE|MATCH_WHOLE_STRING|RETURN_FIRST_MATCH)

#endif


