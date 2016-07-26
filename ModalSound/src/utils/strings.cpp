#include "strings.hpp"
#include <cctype>
#include <cstring>
#include <algorithm>

namespace sploosh
{
    
/*!
 * Helper method to remove the space characters from
 * the head and the tail of the string.
 */
extern std::string trim_copy(const std::string& text) 
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

} // end namespace

