/*
 * =====================================================================================
 *
 *       Filename:  IOEndian.h
 *
 *        Version:  1.0
 *        Created:  11/02/2012 23:39:26
 *       Revision:  none
 *       Compiler:  icpc/gcc
 *
 *         Author:  Changxi Zheng (cz), cxz@cs.columbia.edu
 *                  Columbia University
 *
 * =====================================================================================
 */
#ifndef IO_ENDIAN_INC
#   define IO_ENDIAN_INC

#undef IO_BIG_ENDIAN
#undef IO_LITTLE_ENDIAN

#if defined(__linux)

#   include <endian.h>

#   if __BYTE_ORDER == __LITTLE_ENDIAN
#       define IO_LITTLE_ENDIAN
#   else
#       define IO_BIG_ENDIAN
#   endif

#   if !defined(le32toh) || !defined(be32toh) || !defined(le64toh) || !defined(be64toh)
#       include <byteswap.h>
#   endif 

#   ifndef le32toh
#      if __BYTE_ORDER == __LITTLE_ENDIAN
#          define le32toh(x) (x)
#      else
#          define le32toh(x) __bswap_32 (x)
#      endif
#   endif

#   ifndef be32toh
#      if __BYTE_ORDER == __LITTLE_ENDIAN 
#          define be32toh(x) __bswap_32 (x)
#      else
#          define be32toh(x) (x)
#      endif       
#   endif

#   ifndef le64toh
#      if __BYTE_ORDER == __LITTLE_ENDIAN 
#          define le64toh(x) (x)
#      else
#          define le64toh(x) __bswap_64 (x)
#      endif       
#   endif

#   ifndef be64toh
#      if __BYTE_ORDER == __LITTLE_ENDIAN 
#          define be64toh(x) __bswap_64 (x)
#      else
#          define be64toh(x) (x)
#      endif       
#   endif

#   ifndef  htole16
#      if __BYTE_ORDER == __LITTLE_ENDIAN 
#          define htole16(x) (x)
#      else
#          define htole16(x) __bswap_16 (x)
#      endif       
#   endif

#elif defined(__APPLE__) || defined(MACOSX)

#   include <libkern/OSByteOrder.h>

#   if defined(__LITTLE_ENDIAN__)
#       define IO_LITTLE_ENDIAN
#   elif defined(__BIG_ENDIAN__)
#       define IO_BIG_ENDIAN
#   else
#       error --- Cannot determine OS byte order ---
#   endif

#   ifndef le32toh
#       define le32toh(x)   OSSwapLittleToHostInt32(x)
#   endif

#   ifndef be32toh
#       define be32toh(x)   OSSwapBigToHostInt32(x)
#   endif

#   ifndef le64toh
#       define le64toh(x)   OSSwapLittleToHostInt64(x)
#   endif

#   ifndef be64toh
#       define be64toh(x)   OSSwapBigToHostInt64(x)
#   endif

#   ifndef htole16
#       define htole16(x)   OSSwapHostToLittleInt16(x)
#   endif
#endif

#endif
