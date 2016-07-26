/******************************************************************************
 *  File: nano_timer.c
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
#include "nano_timer.h"

#include <time.h>
#include <sys/time.h>
#include <assert.h>

double GetNanoTimed(void)
{
    struct timespec ts;
#ifndef NDEBUG
    int h = clock_gettime(CLOCK_MONOTONIC, &ts);
    assert(!h);
#else
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif
    //int h = clock_gettime(CLOCK_REALTIME, &ts);
    return (double)ts.tv_sec+(double)ts.tv_nsec*1e-9;
}

long long GetNanoTimei(void)
{
    struct timespec ts;
#ifndef NDEBUG
    int h = clock_gettime(CLOCK_MONOTONIC, &ts);
    assert(!h);
#else
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif
    //int h = clock_gettime(CLOCK_REALTIME, &ts);
    assert(!h);
    return (long long)ts.tv_sec*1000000000+ts.tv_nsec;
}

double GetMilliTimed(void)
{
    struct timeval tv;
#ifndef NDEBUG
    int h = gettimeofday(&tv, NULL);
    assert(!h);
#else
    gettimeofday(&tv, NULL);
#endif
    return (double)tv.tv_sec+(double)tv.tv_usec*1e-6;
}

long long GetMilliTimei(void)
{
    struct timeval tv;
#ifndef NDEBUG
    int h = gettimeofday(&tv, NULL);
    assert(!h);
#else
    gettimeofday(&tv, NULL);
#endif
    return (long long)tv.tv_sec*1000000+tv.tv_usec;
}
