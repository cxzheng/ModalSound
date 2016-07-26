/******************************************************************************
 *  File: strings.cpp
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
#include "strings.hpp"
#include <cctype>
#include <cstring>

namespace carbine
{
    
/*!
 * Helper method to remove the space characters from
 * the head and the tail of the string.
 */
extern std::string trim(const std::string text) 
{   
    size_t i;
    for (i = 0; i < text.length() && isspace(text[i]); ++i);
    if ( i == text.length() ) return "";    

    std::string ret = text.substr(i);
    for (i = ret.length(); isspace(ret[i-1]); --i);
    ret = ret.substr(0, i);

    return ret;
}

/*!
 * check if the given string starts with the given prefix
 */
extern bool start_with(const char* pre, const char* str)
{
    for(;*str != '\0' && *pre != '\0' && *pre == *str;++ pre, ++ str);
    return *pre == '\0';
}

extern bool end_with(const char* suf, const char* str)
{
    int i;
    int lsuf = strlen(suf);
    int lstr = strlen(str);

    if ( lsuf > lstr ) return false;

    suf += (lsuf-1);
    str += (lstr-1);
    for(i = 0;i < lsuf && *suf == *str;-- suf, -- str, ++ i);
    return i == lsuf;
}

/*
template <class T>
std::vector< std::basic_string<T> > 
tokenize(const std::basic_string<T> &s, const std::basic_string<T> &delim) 
{
    std::vector<std::basic_string<T> > ret(0);
    for(int b,e=0;;ret.push_back(s.substr(b,e-b)))
        if ( (b=s.find_first_not_of(delim,e))==(e=s.find_first_of(delim,b)) )
            return ret;
}
*/

} // end namespace

