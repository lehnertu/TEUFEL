/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: matlib.h
 * purpose: definitions for C matrix routines 
 *
 * Michael Borland, 1989
 $Log: not supported by cvs2svn $
 Revision 1.10  2008/01/05 03:35:32  borland
 Adde m_subtract().

 Revision 1.9  2003/07/22 20:01:18  soliday
 Added support for Kylix.

 Revision 1.8  2002/08/14 15:40:15  soliday
 Added Open License

 Revision 1.7  2002/08/07 18:52:39  borland
 Added prototypes for fmat_* and macros for fm_*  routines, which do single-precision
 matrix computations.

 Revision 1.6  1999/09/14 18:04:48  soliday
 Added export commands for WIN32 dll files.

 Revision 1.5  1998/06/11 22:21:31  borland
 Fixed typo with p_merror macro.

 Revision 1.4  1998/04/21 21:26:05  borland
 New version that provides compatibility with Meschach library in the
 same executable (but not in the same .c file).

 * Revision 1.2  1995/09/05  21:15:12  saunders
 * First test release of the SDDS1.5 package.
 *
 */
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MATLIB 

#ifdef MESCHACH
#error "The Meschach header is already included by this file, which just included matlib.h.  MATLIB and Meschach use common function names!"
#endif

#undef epicsShareFuncMATLIB
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_MATLIB)
#define epicsShareFuncMATLIB  __declspec(dllexport)
#else
#define epicsShareFuncMATLIB
#endif
#else
#define epicsShareFuncMATLIB
#endif

typedef struct {
	double **a;
	int n, m;
	} MATRIX;

typedef struct {
	float **a;
	int n, m;
	} FMATRIX;

epicsShareFuncMATLIB extern void mat_error(char *message);
#define m_error(message) mat_error(message)
 
epicsShareFuncMATLIB extern int p_materror(char *message);
#define p_merror(message) p_materror(message)
 
epicsShareFuncMATLIB extern int mat_add(MATRIX *C, MATRIX *A, MATRIX *B);
#define m_add(C,A,B) mat_add(C,A,B)
 
epicsShareFuncMATLIB extern int mat_subtract(MATRIX *C, MATRIX *A, MATRIX *B);
#define m_subtract(C,A,B) mat_subtract(C,A,B)
 
epicsShareFuncMATLIB extern void mat_alloc(MATRIX **A, int n, int m);
#define m_alloc(A,n,m) mat_alloc(A,n,m)
 
epicsShareFuncMATLIB extern void mat_alloc1(MATRIX **A, int n, int m);
#define m_alloc1(A,n,m) mat_alloc1(A,n,m)
 
epicsShareFuncMATLIB extern int mat_copy(MATRIX *A, MATRIX *B);
#define m_copy(A,B) mat_copy(A,B)
 
epicsShareFuncMATLIB extern double mat_det(MATRIX *D);
#define m_det(D) mat_det(D)
 
epicsShareFuncMATLIB extern void mat_free(MATRIX **A);
#define m_free(A) mat_free(A)
 
epicsShareFuncMATLIB extern int mat_invert(MATRIX *A, MATRIX *B);
#define m_invert(A,B) mat_invert(A,B)
 
epicsShareFuncMATLIB extern int mat_mult(MATRIX *C, MATRIX *A, MATRIX *B);
#define m_mult(C,A,B) mat_mult(C,A,B)
 
epicsShareFuncMATLIB extern int mat_scmul(MATRIX *B, MATRIX *A, double a);
#define m_scmul(B,A,a) mat_scmul(B,A,a)
 
epicsShareFuncMATLIB extern void mat_show(MATRIX *A, char *format, char *label, FILE *fp);
#define m_show(A,format,label,fp) mat_show(A,format,label,fp)
 
epicsShareFuncMATLIB extern void mat_get(MATRIX *A);
#define m_get(A) mat_get(A)
 
epicsShareFuncMATLIB extern void mat_rand(MATRIX *A, double lo, double hi);
#define m_rand(A,lo,hi) mat_rand(A,lo,hi)
 
epicsShareFuncMATLIB extern int mat_trans(MATRIX *B, MATRIX *A);
#define m_trans(B,A) mat_trans(B,A)
 
epicsShareFuncMATLIB extern void mat_zero(MATRIX *A);
#define m_zero(A) mat_zero(A)
 
epicsShareFuncMATLIB extern void mat_identity(MATRIX *A);
#define m_identity(A) mat_identity(A)
 
epicsShareFuncMATLIB extern int mat_check(MATRIX *A);
#define m_check(A) mat_check(A)

/* float versions  */

epicsShareFuncMATLIB extern int fmat_add(FMATRIX *C, FMATRIX *A, FMATRIX *B);
#define fm_add(C,A,B) fmat_add(C,A,B)
 
epicsShareFuncMATLIB extern void fmat_alloc(FMATRIX **A, int n, int m);
#define fm_alloc(A,n,m) fmat_alloc(A,n,m)
 
epicsShareFuncMATLIB extern void fmat_alloc1(FMATRIX **A, int n, int m);
#define fm_alloc1(A,n,m) fmat_alloc1(A,n,m)
 
epicsShareFuncMATLIB extern int fmat_copy(FMATRIX *A, FMATRIX *B);
#define fm_copy(A,B) fmat_copy(A,B)
 
epicsShareFuncMATLIB extern float fmat_det(FMATRIX *D);
#define fm_det(D) fmat_det(D)
 
epicsShareFuncMATLIB extern void fmat_free(FMATRIX **A);
#define fm_free(A) fmat_free(A)
 
epicsShareFuncMATLIB extern int fmat_invert(FMATRIX *A, FMATRIX *B);
#define fm_invert(A,B) fmat_invert(A,B)
 
epicsShareFuncMATLIB extern int fmat_mult(FMATRIX *C, FMATRIX *A, FMATRIX *B);
#define fm_mult(C,A,B) fmat_mult(C,A,B)
 
epicsShareFuncMATLIB extern int fmat_scmul(FMATRIX *B, FMATRIX *A, float a);
#define fm_scmul(B,A,a) fmat_scmul(B,A,a)
 
epicsShareFuncMATLIB extern void fmat_show(FMATRIX *A, char *format, char *label, FILE *fp);
#define fm_show(A,format,label,fp) fmat_show(A,format,label,fp)
 
epicsShareFuncMATLIB extern void fmat_get(FMATRIX *A);
#define fm_get(A) fmat_get(A)
 
epicsShareFuncMATLIB extern void fmat_rand(FMATRIX *A, float lo, float hi);
#define fm_rand(A,lo,hi) fmat_rand(A,lo,hi)
 
epicsShareFuncMATLIB extern int fmat_trans(FMATRIX *B, FMATRIX *A);
#define fm_trans(B,A) fmat_trans(B,A)
 
epicsShareFuncMATLIB extern void fmat_zero(FMATRIX *A);
#define fm_zero(A) fmat_zero(A)
 
epicsShareFuncMATLIB extern void fmat_identity(FMATRIX *A);
#define fm_identity(A) fmat_identity(A)
 
epicsShareFuncMATLIB extern int fmat_check(FMATRIX *A);
#define fm_check(A) fmat_check(A)

#define MATLIB 1
#endif

#ifdef __cplusplus
}
#endif
