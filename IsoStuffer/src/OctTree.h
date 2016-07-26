/******************************************************************************
 *  File: OctTree.h
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
#ifndef ISOSTUFFING_OCTTREE_H
#   define ISOSTUFFING_OCTTREE_H

#include "linearalgebra/Tuple3.hpp"
#include "utils/arrays.hpp"

namespace IsoStuffing
{

/*
 * The bitmask for children node: (in x,y,z order)
 * (0,0,0)      0
 * (1,0,0)      1
 * (1,0,1)      2
 * (0,0,1)      3
 * (0,1,0)      4
 * (1,1,0)      5
 * (1,1,1)      6
 * (0,1,1)      7
 *
 * face center mask
 * -x           0
 * +x           1
 * -y           2
 * +y           3
 * -z           4
 * +z           5
 *
 * edge middle mask
 *                         E11
 *             4--------------------------5
 *            /.                         /|
 *           / .                        / |
 *        E5/  .                     E6/  |
 *         /   .E0                    /   |E1
 *        /    .                     /    |    y^
 *       7--------------------------6     |     |
 *       |     .    E10             |     |     |
 *       |     .                    |     |     |
 *       |     .                    |     |     |
 *       |    0.....................|.....1     0---------> x
 *       |E3   .          E8        |E2   /    /
 *       |    .                     |    /    /
 *       |   .                      |   /    /
 *       |  .E4                     |  /E7  z
 *       | .                        | /
 *       |.                         |/
 *       3--------------------------2
 *                  E9
 * 
 */
struct Octant
{
    bool    isLeaf;
    int     childrenMask;       // label what children this octant has
    int     faceCenterMask;
    int     edgeMidMask;

    Octant(bool leaf):isLeaf(leaf), childrenMask(0),
            faceCenterMask(0), edgeMidMask(0) { }

    Octant(bool leaf, int cbm):isLeaf(leaf), childrenMask(cbm),
            faceCenterMask(0), edgeMidMask(0) { }
};

template <class DistFunc> class IsoStuffer;

class OctTree
{
    template <class DistFunc> friend class IsoStuffer;

    public:
        typedef boost::multi_array<Octant*, 3>      TOctLevel;

        OctTree(const Tuple3i& highestRes, int nlevel);
        ~OctTree();

        /*
         * create a leaf octant at the finest level
         */
        void create_leaf(int ix, int iy, int iz);

        void create_octree();

        void weak_balance();

        void update_bitmasks();

        /* resolution at given level */
        const Tuple3i& resolution(int l) const
        {   return levelRes_[l]; }

        const TOctLevel& level(int l) const
        {   return levels_[l]; }

        const std::vector<Tuple3i>& resolutions() const
        {   return levelRes_; }

        /*
         * Check if the children mask is consistent with the
         * actual children. USED ONLY FOR TESTING
         */
        void check_children_mask() const;

    private:
        bool check_for_weak_balance(int lid, int ix, int iy, int iz) const;

    private:
        const int               nLevels_;           // # of octree levels

        /*
         * level-0 is the lowest level (finest)
         * level-n is the higer level  (coarser)
         */
        std::vector<TOctLevel>  levels_;
        std::vector<Tuple3i>    levelRes_;        // highest resolution of the oct tree
};

}

#endif
