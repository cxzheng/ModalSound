/******************************************************************************
 *  File: nano_timer.h
 *
 *  This file is part of isostuffer
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
/*
 * Changxi Zheng (cxzheng@cs.cornell.edu)
 */

#ifndef NANO_TIMER_H
#   define NANO_TIMER_H

#undef __BEGIN_DECLS
#undef __END_DECLS

#ifdef	__cplusplus
#   define __BEGIN_DECLS    extern "C" {
#   define __END_DECLS      }
#else
#   define __BEGIN_DECLS
#   define __END_DECLS
#endif

#undef __p
#if defined (__STDC__) || defined (_AIX) \
    || (defined (__mips) && defined (_SYSTYPE_SVR4)) \
    || defined(WIN32) || defined(__cplusplus)
#   define __p(protos) protos
#else
#   define __p(protos) ()
#endif

__BEGIN_DECLS

double    GetNanoTimed  __p((void));
long long GetNanoTimei  __p((void));
double    GetMilliTimed __p((void));
long long GetMilliTimei __p((void));

__END_DECLS

#endif
