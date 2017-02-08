/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: rpn.h
 * Michael Borland, 1988
 $Log: not supported by cvs2svn $
 Revision 1.33  2010/02/24 23:59:38  borland
 Added InvFq function to RPN.

 Revision 1.32  2010/02/05 16:35:41  soliday
 Fixed an issue with creating sharied libraries.

 Revision 1.31  2010/02/04 23:42:34  soliday
 Updated so that it can be used by c++

 Revision 1.30  2009/03/10 14:11:28  shang
 added strlen function

 Revision 1.29  2008/11/11 21:05:18  soliday
 Updated to fix an issue on vxWorks.

 Revision 1.28  2006/08/31 15:06:56  soliday
 Updated to work with SDDS2

 Revision 1.27  2005/11/04 22:47:00  soliday
 Updated code to be compiled by a 64 bit processor.

 Revision 1.26  2005/01/13 16:52:23  shang
 added string related functions.

 Revision 1.25  2005/01/12 22:03:51  shang
 added is_string argument to is_memory()

 Revision 1.24  2005/01/10 20:59:13  shang
 added store_in_str_mem for storing string values in memory and modified
 some memory functions to make them work for string type data.

 Revision 1.23  2004/02/09 15:26:57  soliday
 Updated to export definitions of push_string and pop_string

 Revision 1.22  2003/10/30 21:21:59  soliday
 Changed the function name "store" to "store_in_mem" to avoid a system function name.

 Revision 1.21  2003/10/30 19:06:53  soliday
 Added epicsShareExtern definition

 Revision 1.20  2003/10/13 21:11:14  soliday
 Fixed problem with missing extern statement.

 Revision 1.19  2003/10/03 14:55:29  soliday
 Moved some statements from rpn_internal.h to rpn.h

 Revision 1.18  2003/07/22 20:01:18  soliday
 Added support for Kylix.

 Revision 1.17  2003/06/19 16:41:10  shang
 added push_long and pop_long functions

 Revision 1.16  2003/06/15 17:21:50  borland
 Added prototypes for internal functions for modified bessel functins.

 Revision 1.15  2003/03/17 22:36:38  borland
 Added definition of rpn_mudf() function.

 Revision 1.14  2002/08/14 15:40:17  soliday
 Added Open License

 Revision 1.13  2001/04/09 20:38:27  soliday
 Added if2pf.

 Revision 1.12  2000/06/15 15:53:32  soliday
 Renamed help to rpn_help so that it does not conflict with the IOC.

 Revision 1.11  2000/05/09 14:27:45  borland
 Added rpn_srnd() and rpn_grndl() prototypes.

 Revision 1.10  2000/04/06 22:25:46  soliday
 Fixed prototypes by adding void definitions.

 Revision 1.9  2000/03/27 20:26:14  borland
 Added prototypes for rpn_execs() and rpn_execn().

 Revision 1.8  1999/09/16 22:03:56  soliday
 push_num is now exported to the WIN32 dll

 Revision 1.7  1999/09/14 18:06:12  soliday
 Added export commands for WIN32 dll files.

 Revision 1.6  1999/06/01 14:35:45  soliday
 Removed warnings when compiled under Linux

 Revision 1.5  1998/08/21 19:46:33  borland
 New version per R. Soliday for his optimized version of rpn.

 Revision 1.4  1996/10/22 18:48:18  borland
 Added prototypes for poisson statistics significance level routines.

 * Revision 1.3  1996/02/12  17:18:25  borland
 * Added prototype for rpn_quick_store().
 *
 * Revision 1.2  1995/09/05  21:15:37  saunders
 * First test release of the SDDS1.5 package.
 *
 */

#include "SDDStypes.h"
#if defined(_WIN32) && !defined(_MINGW)
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
#define PRId32 "ld"
#define SCNd32 "ld"
#define PRIu32 "lu"
#define SCNu32 "lu"
#define INT32_MAX (2147483647)
#else
#if defined(vxWorks)
#define PRId32 "ld"
#define SCNd32 "ld"
#define PRIu32 "lu"
#define SCNu32 "lu"
#define INT32_MAX (2147483647)
#else
#include <inttypes.h>
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif


#undef epicsShareFuncRPNLIB
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_RPNLIB)
#define epicsShareFuncRPNLIB  __declspec(dllexport)
#define epicsShareExtern extern __declspec(dllexport)
#else
#define epicsShareFuncRPNLIB
#if !defined(EXPORT_SDDS)
#define epicsShareExtern extern __declspec(dllimport)
#endif
#endif
#else
#undef epicsShareExtern
#define epicsShareFuncRPNLIB
#define epicsShareExtern extern
#endif

/* function call for programs that use rpn: */
epicsShareFuncRPNLIB double rpn(char *expression);

/* prototypes for code in file array.c */
void rpn_alloc(void);
void rref(void);
void sref(void);
epicsShareFuncRPNLIB long rpn_createarray(long size);
epicsShareFuncRPNLIB double *rpn_getarraypointer(long memory_number, int32_t *length);
epicsShareFuncRPNLIB long rpn_resizearray(long arraynum, long size);
void udf_createarray(short type, short index, double data, char *rpn, long i_udf);
void udf_cond_createarray(long colon, long i);
void udf_modarray(short type, short index, double data, long i);
void udf_id_createarray(long start_index_value, long end_index_value);
void udf_create_unknown_array(char *ptr, long index);



/* prototypes for code in file conditional.c */
void conditional(void);
long dissect_conditional(char **branch, long is_true);
void conditional_udf(long udf_current_step);

/* prototypes for code in file execute.c */
epicsShareFuncRPNLIB long execute_code(void);
void set_ptrs(char **text, char **buffer, char **token);
long is_func(char *string);
void quit(void);
void rpn_help(void);
long stack_test(long stackptr, long numneeded, char *stackname, char *caller);
void stop(void);
void ttrace(void);
void rep_stats(void);
void rpn_sleep(void);
long cycle_through_udf(void);

/* long func_compare(struct FUNCTION *f1, struct FUNCTION *f2); */
epicsShareFuncRPNLIB int func_compare(const void *f1v, const void *f2v);

/* prototypes for code in file get_token_rpn.c */
char *get_token_rpn(char *s, char *buf, long lbuf, long *spos);

/* prototypes for code in file logical.c */
void greater(void);
void less(void);
void greater_equal(void);
void less_equal(void);
void equal(void);
void not_equal(void);
void log_and(void);
void log_or(void);
void log_not(void);
void poplog(void);
void lton(void);
void ntol(void);

/* prototypes for code in file math.c */
void rpn_strlen(void);
void rpn_streq(void);
void rpn_strgt(void);
void rpn_strlt(void);
void rpn_strmatch(void);
void rpn_add(void);
void rpn_sum(void);
void rpn_subtract(void);
void rpn_multiply(void);
void rpn_divide(void);
void rpn_sqrt(void);
void rpn_inverseFq(void);
void rpn_square(void);
void rpn_power(void);
void rpn_sin(void);
void rpn_cos(void);
void rpn_atan(void);
void rpn_asin(void);
void rpn_acos(void);
void rpn_ex(void);
void rpn_ln(void);
void rpn_erf(void);
void rpn_erfc(void);
void rpn_JN(void);
void rpn_YN(void);
void rpn_IN(void);
void rpn_KN(void);
void rpn_FresS(void);
void rpn_FresC(void);
void rpn_G1y(void);
void rpn_int(void);
void rpn_cei1(void);
void rpn_cei2(void);
void rpn_lngam(void);
void rpn_betai(void), rpn_gammaP(void), rpn_gammaQ(void);
void rpn_rnd(void);
void rpn_grnd(void);
void rpn_grndlim(void);
void rpn_srnd(void);
void rpn_atan2(void);
void rpn_push_nan(void);
void rpn_isinf(void);
void rpn_isnan(void);
void rpn_poissonSL(void);
void rpn_simpson(void);
void rpn_isort_stack(void);
void rpn_dsort_stack(void);

/* prototypes for code in file memory.c */
epicsShareFuncRPNLIB long rpn_create_mem(char *name, short is_string);
epicsShareFuncRPNLIB long rpn_store(double value, char *str_value, long memory_number);
epicsShareFuncRPNLIB long rpn_quick_store(double value, char *str_value, long memory_number);
epicsShareFuncRPNLIB double rpn_recall(long memory_number);
epicsShareFuncRPNLIB char *rpn_str_recall(long memory_number);
void store_in_mem(void);
void store_in_str_mem(void);
epicsShareFuncRPNLIB long is_memory(double *val, char **str_value, short *is_string, char *string);
void revmem(void);

/* prototypes for code in file pcode.c */
void gen_pcode(char *s, long i);

/* prototypes for code in file pop_push.c */
double pop_num(void);
long pop_long(void);
epicsShareFuncRPNLIB long push_num(double num);
epicsShareFuncRPNLIB long push_long(long num);
epicsShareFuncRPNLIB long pop_log(int32_t *logical);
long push_log(long logical);
long pop_file(void);
long push_file(char *filename);
epicsShareFuncRPNLIB char *pop_string(void);
epicsShareFuncRPNLIB void push_string(char *s);
void pop_code(void);
void push_code(char *code, long mode);

/* prototypes for code in file prompt.c */
epicsShareFuncRPNLIB long prompt(char *prompt_s, long do_prompt);

/* prototypes for code in file rpn_io.c */
void open_cominp(void);
void open_io(void);
void close_io(void);
void rpn_gets(void);
void scan(void);
void format(void);
void get_format(void);
epicsShareFuncRPNLIB char *choose_format(long flag, double x);
void view(void);
void view_top(void);
void tsci(void);
void viewlog(void);
void fprf(void);
void view_str(void);
void rpn_puts(void);
void sprf(void);

/* prototypes for code in file stack.c */
void swap(void);
void duplicate(void);
void nduplicate(void);
void stack_lev(void);
void pop(void);
epicsShareFuncRPNLIB void rpn_clear(void);
void pops(void);
void rup(void);
void rdn(void);
void dup_str(void);
void exe_str(void);

/* prototypes for code in file udf.c */
long find_udf(char *udf_name);
long find_udf_mod(char *udf_name);
short get_udf(long number);
void get_udf_indexes(long number);
void make_udf(void);
void rpn_mudf(void);
epicsShareFuncRPNLIB void create_udf(char *name, char *function);
epicsShareFuncRPNLIB void link_udfs(void);
void insert_udf(char *instr, char *udf_string);
void revudf(void);

/* prototypes for code in file rpn_csh.c */
void rpn_csh(void);
void rpn_csh_str(void);
void rpn_execs(void);
void rpn_execn(void);

/* prototypes for code in file rpn_draw.c */
void rpn_draw(void);

/* prototypes for code in file rpn_error.c */
void rpn_set_error(void);
epicsShareFuncRPNLIB long rpn_check_error(void);
epicsShareFuncRPNLIB void rpn_clear_error(void);

/* prototypes for code in file infixtopostfix.c */
#define IFPF_BUF_SIZE    1024
epicsShareFuncRPNLIB int if2pf(char *pfix, char *ifix, size_t size_of_pfix);

#define STACKSIZE 5000
epicsShareExtern long dstackptr;
epicsShareExtern long sstackptr;

#ifdef __cplusplus
}
#endif
