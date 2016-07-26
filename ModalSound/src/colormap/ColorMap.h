#ifndef COLOR_MAP_INC
#   define COLOR_MAP_INC

#include "sc/Tuple3.hpp"
#include "utils/macros.h"
#include "JET.h"
#include "AUTUMN.h"
#include "BONE.h"
#include "UM4.h"
#include <assert.h>

//! Map a scalar (int) value to a tuple which represent a color.
/*!
 * Similar with the colormap in matlab
 */
class ColorMap
{
    public:
        // =============== Constructor =============
        ColorMap(int num):m_num(num), m_min(0), m_max(1), 
                m_step(1./num), mp_colors(NULL)
        { }

        virtual ~ColorMap()
        {
            SAFE_DELETE_ARRAY(mp_colors);
        }

        // =========================================
        Tuple3f get_color(int s) const
        {
            assert(s >= 0 && s < m_num);
            return mp_colors[s];
        }

        void get_color(int s, Tuple3f& color) const
        {
            assert(s >= 0 && s < m_num);
            color.set(mp_colors[s].r, mp_colors[s].g, mp_colors[s].b);
        }

        void get_color(int s, Tuple3f* color) const
        {
            assert(s >= 0 && s < m_num);
            color->set(mp_colors[s].r, mp_colors[s].g, mp_colors[s].b);
        }

        Tuple3f get_interpolated_color(double v) const
        {
            if ( v >= m_max ) return mp_colors[m_num - 1];
            if ( v <= m_min ) return mp_colors[0];

            double re  = (v - m_min)/m_step;
            int    idx = (int)re;

            if ( idx == m_num - 1 ) // the last one
                return mp_colors[idx];

            re -= idx;
            return mp_colors[idx]*(1.-re) + mp_colors[idx+1]*re;        // linear interpolation
        }

        int  color_num() const { return m_num; }
        void set_interpolation_range(double minv, double maxv)
        {
            assert(minv < maxv);
            m_min = minv;
            m_max = maxv;
            m_step = (maxv - minv) / m_num;
        }

    protected:
        int         m_num;
        double      m_min, m_max;
        double      m_step;
        Tuple3f*   mp_colors;
};

class LinearColorMap : public ColorMap
{
    public:
        LinearColorMap(const Tuple3f& minc, const Tuple3f& maxc, int num):ColorMap(num)
        {
            mp_colors = new Tuple3f[num];
            Tuple3f dc = (maxc - minc) / (num - 1);
            for(int i = 0;i < num;++ i)
                mp_colors[i] = minc + dc * i;
        }
};

class JetColorMap : public ColorMap
{
    public:
        JetColorMap():ColorMap(256)
        {
            mp_colors = (Tuple3f*)colormap_jet;
        }

        ~JetColorMap()
        {
            mp_colors = NULL;
        }
};

class AutumnColorMap : public ColorMap
{
    public:
        AutumnColorMap():ColorMap(256)
        {
            mp_colors = (Tuple3f*)colormap_autumn;
        }
        ~AutumnColorMap()
        {
            mp_colors = NULL;
        }
};

class BoneColorMap : public ColorMap
{
    public:
        BoneColorMap():ColorMap(256)
        {
            mp_colors = (Tuple3f*)colormap_bone;
        }
        ~BoneColorMap()
        {
            mp_colors = NULL;
        }
};

class Um4ColorMap : public ColorMap
{
    public:
        Um4ColorMap():ColorMap(1024)
        {
            mp_colors = (Tuple3f*)colormap_um4;
        }
        ~Um4ColorMap()
        {
            mp_colors = NULL;
        }
};

#endif
