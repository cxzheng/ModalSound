/******************************************************************************
 *  The definitions of some utility macros
 *
 *  Changxi Zheng (cxzheng@cs.cornell.edu)
 *
 ******************************************************************************/

/*
 * The returned value should be zero, if an operation is succeed.
 * otherwise, non-zero is returned.
 */
#ifndef MACRO_DEF_H
#   define MACRO_DEF_H

#ifndef ERROR_RETURN
#   define ERROR_RETURN        -1
#endif

#ifndef SUCC_RETURN
#   define SUCC_RETURN          0
#endif

#ifndef _FAILED
#   define _FAILED(p)      ((p) < 0)
#endif

#ifndef _SUCCEEDED
#   define _SUCCEEDED(p)  ((p) >= 0)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*
#ifdef HAVE_CONFIG_H
#   include "conts-config.h"
#endif
*/
#ifndef __ASSERT_VOID_CAST
#   define __ASSERT_VOID_CAST   (void)
#endif

#define SHOULD_NEVER_HAPPEN(x)                                              \
        {                                                                   \
            fprintf(stderr, "ERROR: SHOULD NEVER ARRIVE HERE at %s:L.%d\n", \
                    __FILE__, __LINE__);                                    \
            fflush(stderr);                                                 \
            exit(x);                                                        \
        }

#define ABORT_HERE(x)                                                       \
        {                                                                   \
            fprintf(stdout, "Abort execution at %s:L.%d\n",                 \
                    __FILE__, __LINE__);                                    \
            fflush(stdout);                                                 \
            exit(x);                                                        \
        }

#if 0 // FIXME: This conflicts with Eigen
#if defined(DEBUG) | defined(_DEBUG)
#   ifndef V
#       define V(x)                                                     \
        {                                                               \
            int hr = x;                                                 \
            if( _FAILED(hr) ) {                                         \
                fprintf(stderr, "DEBUG Failed: line %d, file \"%s\"\n", \
                        __LINE__, __FILE__);                            \
                fflush(stderr);                                         \
            }                                                           \
        }
#   endif
#   ifndef V_RETURN
#       define V_RETURN(x)                                              \
        {                                                               \
            int hr = x;                                                 \
            if( _FAILED(hr) ) {                                         \
                fprintf(stderr, "DEBUG Failed: line %d, file \"%s\"\n", \
                        __LINE__, __FILE__);                            \
                fflush(stderr);                                         \
                return hr;                                              \
            }                                                           \
        }
#   endif
#   ifndef MSG_ASSERT
#     if defined(__APPLE__) || defined(MACOSX)
#       define MSG_ASSERT(x, ...)                                       \
        {                                                               \
            if ( !(x) ) {                                               \
                fprintf(stderr, __VA_ARGS__);                           \
                fprintf(stderr, "Assertion Failed: %s at line %d, file \"%s\"\n", \
                        __STRING(x), __LINE__, __FILE__);               \
                abort();                                                \
            }                                                           \
        }
#     else      // define for mac
#       define MSG_ASSERT(x, ...)                                       \
        {                                                               \
            if ( !(x) ) {                                               \
                fprintf(stderr, __VA_ARGS__);                           \
                __assert_fail(__STRING(x), __FILE__, __LINE__, __ASSERT_FUNCTION); \
            }                                                           \
        }
#     endif
#   endif
#else
#   ifndef V
#       define V(x)         { x; }
#   endif
#   ifndef V_RETURN
#       define V_RETURN(x)  { int hr = x; if( _FAILED(hr) ) { return hr; } }
#   endif
#   ifndef MSG_ASSERT
#       define MSG_ASSERT(x, ...)   (__ASSERT_VOID_CAST (0))
#   endif
#endif
#endif 

/*
 * Safely delete/release a reference
 */
#ifdef __cplusplus

#ifndef SAFE_DELETE
#   define SAFE_DELETE(p)       { if(p) { delete (p);     (p)=NULL; } }
#endif    
#ifndef SAFE_DELETE_ARRAY
#   define SAFE_DELETE_ARRAY(p) { if(p) { delete[] (p);   (p)=NULL; } }
#endif    
#ifndef SAFE_RELEASE
#   define SAFE_RELEASE(p)      { if(p) { (p)->Release(); (p)=NULL; } }
#endif

#endif

// Using similar macro from Linux kernel code
#if defined(__GNUC__)
#if defined(likely) | defined(unlikely)
#   error likely or unlikely is already defined
#else
#   define likely(x)    __builtin_expect(!!(x), 1)
#   define unlikely(x)  __builtin_expect(!!(x), 0)
#endif
#elif defined(_MSC_VER)
#ifndef likely
#   define likely(x)        (x)
#endif
#ifndef unlikely
#   define unlikely(x)      (x)
#endif
#endif

#endif
