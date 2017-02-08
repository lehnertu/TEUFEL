/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file   : mdbsunos4.h
 * purpose: definitions for general use, for mdblib, and for mdbmth.
 *          These definitions are specific to the SUNOS4.
 *
 * Michael Borland, 1988
 $Log: not supported by cvs2svn $
 Revision 1.2  1995/09/05 21:15:20  saunders
 First test release of the SDDS1.5 package.

 */
#ifndef _MDBSUNOS4_
#define _MDBSUNOS4_ 1

/* escape sequences for VT100/200 series terminals */
#define CLEAR_SCREEN "\033[2J"
#define HOME_SCREEN  "\033[H"
#define HOME_CLEAR "\033[2J\033[H"

/* commonly needed terminal control characters */
#define CARRIAGE_RETURN '\015'
#define LINE_FEED '\012'
#define FORM_FEED '\014'
#define BACK_SPACE '\010'
#define ESCAPE '\033'

/************************ mathematical stuff *******************************/
#ifndef _MDBSUNOS4_MTH_
#define _MDBSUNOS4_MTH_ 1

#define MAX_LONG LONG_MAX

#if !defined(HUGE)
#define HUGE HUGE_VAL
#endif
#define FHUGE ((float)HUGE)

#endif  /* _MDBSUNOS4_MTH_ */

#endif  /* _MDBSUNOS4_ */

/*   -- Run-time statistics: */
extern long cpu_time(void);
extern long page_faults(void);
extern long dio_count(void);
extern long bio_count(void);
extern long init_stats(void);
extern char *elapsed_time();
extern long memory_count(void);
extern void report_stats(FILE *fp, char *label);
extern double delapsed_time();

extern char *mtime(void);

