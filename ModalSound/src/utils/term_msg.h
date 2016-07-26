/*
 * =====================================================================================
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
 * 
 * -------------------------------------------------------------------------------------
 *
 *       Filename:  term_msg.h
 *
 *    Description:  print message in terminal
 *
 *        Version:  1.0
 *        Created:  06/18/12 18:59:52
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef TERM_MSG_H
#   define TERM_MSG_H

#define TERM_COLOR_BEGIN     "\033["
#define TERM_COLOR_END       "\033[0m"

#define TERM_COLOR_BOLD      1
#define TERM_COLOR_BLINK     5

#define TERM_COLOR_RED       31
#define TERM_COLOR_GREEN     32
#define TERM_COLOR_YELLOW    33
#define TERM_COLOR_BLUE      34
#define TERM_COLOR_MAGENTA   35
#define TERM_COLOR_WHITE     37

#define TERM_COLOR_BLACK_BG  40
#define TERM_COLOR_RED_BG    41
#define TERM_COLOR_BLUE_BG   44

#include <stdio.h>

#define PRINT_WARNING(...)                                  \
        {                                                   \
            fprintf(stderr, "%s%d;%dmWARNING:%s %s%dm",     \
                    TERM_COLOR_BEGIN, TERM_COLOR_BOLD,      \
                    TERM_COLOR_YELLOW, TERM_COLOR_END,      \
                    TERM_COLOR_BEGIN, TERM_COLOR_WHITE);    \
            fprintf(stderr, __VA_ARGS__);                   \
            fprintf(stderr, TERM_COLOR_END);                \
        }

#define PRINT_ERROR(...)                                    \
        {                                                   \
            fprintf(stderr, "%s%d;%dmERROR:%s %s%dm",       \
                    TERM_COLOR_BEGIN, TERM_COLOR_BOLD,      \
                    TERM_COLOR_RED, TERM_COLOR_END,         \
                    TERM_COLOR_BEGIN, TERM_COLOR_YELLOW);   \
            fprintf(stderr, __VA_ARGS__);                   \
            fprintf(stderr, TERM_COLOR_END);                \
        }

#define PRINT_MSG(...)                                     \
        {                                                  \
            printf("%s%d;%dmMESSAGE:%s %s%dm",             \
                   TERM_COLOR_BEGIN, TERM_COLOR_BOLD,      \
                   TERM_COLOR_GREEN, TERM_COLOR_END,       \
                   TERM_COLOR_BEGIN, TERM_COLOR_WHITE);    \
            printf(__VA_ARGS__);                           \
            printf(TERM_COLOR_END);                        \
        }

#endif
