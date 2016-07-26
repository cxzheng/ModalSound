#ifndef CARBINE_STRING_UTILS_H
#   define CARBINE_STRING_UTILS_H

#include <vector>
#include <string>
#include <sstream>

//! string to int function
#define S2D(s,d)  {std::istringstream(s)>>d;}

namespace sploosh
{

//! Trim the string
std::string trim_copy(const std::string& text);

//! Return whether the string a starts with string b
bool start_with(const char* pre, const char* str);

bool end_with(const char* suf, const char* str);

template<typename T> 
inline int Int(const T &t)
{
    int r;
    std::stringstream s;
    s << t;
    s >> r;
    return r;
}

template<typename T> 
inline double Double(const T &t)
{
    double r;
    std::stringstream s;
    s << t;
    s >> r;
    return r;
}

template<typename T>
inline float Float(const T& t)
{
    float r;
    std::stringstream s;
    s << t;
    s >> r;
    return r;
}

//parsing routine
#if 1
template <class T>
static std::vector< std::basic_string<T> > 
tokenize(const std::basic_string<T> &s, const std::basic_string<T> &delim = " \t\n")
{
    std::vector<std::basic_string<T> > ret(0);
    for(int b,e=0;;ret.push_back(s.substr(b,e-b)))
    {
        b=s.find_first_not_of(delim,e);
        e=s.find_first_of(delim,b);
        if ( b==e ) return ret;
    }
    return ret;
}
#endif

template <class T>
static std::vector< std::basic_string<T> >
split(const std::basic_string<T>& s, char delim)
{
    std::vector< std::basic_string<T> > ret(0);
    for(int b = 0, e = 0;b < s.size();b = e+1)
        if ( (e = s.find_first_of(delim, b))==std::string::npos )
        {
            ret.push_back(s.substr(b));
            return ret;
        }
        else
        {
            ret.push_back(s.substr(b, e-b));
        }
    return ret;
}

}

#endif
