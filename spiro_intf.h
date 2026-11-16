/* spiro_intf.h
libspiro - A sharable library of Spiro formula and functions.
*/
#ifndef _SPIRO_INTF_H
#define _SPIRO_INTF_H

#include "spirobatch.h"

/* Generic helper definitions for shared library support */
/* as explained in - https://gcc.gnu.org/wiki/Visibility */
#if defined (_WIN32) || defined (__CYGWIN__) || (_WIN64) || defined (__CYGWIN64__)
  #ifndef __clang__
    #define LS_DLL_IMPORT __declspec(dllimport)
    #define LS_DLL_EXPORT __declspec(dllexport)
    #define LS_DLL_LOCAL
  #else
    /* CLANG already appears to add this. Turn this off. */
    #define LS_DLL_IMPORT
    #define LS_DLL_EXPORT
    #define LS_DLL_LOCAL
  #endif
#else
  #if __GNUC__ >= 4
    #define LS_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define LS_DLL_EXPORT __attribute__ ((visibility ("default")))
    /* 'internal' to avoid library conflicts. 'hidden' for debugging */
    /* #define UN_DLL_LOCAL  __attribute__ ((visibility ("hidden"))) */
    #define LS_DLL_LOCAL  __attribute__ ((visibility ("internal")))
  #else
    #define LS_DLL_IMPORT
    #define LS_DLL_EXPORT
    #define LS_DLL_LOCAL
  #endif
#endif

int spiroreverse(spiro_cp *src, int n); /* internal use */

typedef int (*spiro_tol_emit)(void *ctx, const spiro_batch_segment *segment);

int spiro_to_bpath_tol0(const spiro_cp *src, const spiro_seg *s,
                        double *dm, int ncq, int n, double tolerance,
                        spiro_tol_emit emit, void *ctx);

#endif
