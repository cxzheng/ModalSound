/******************************************************************************
 *  File: OctTree.cpp
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
#include "OctTree.h"
#include "tables.h"
#include "utils/arrays.hpp"
#include "utils/macros.h"

namespace IsoStuffing
{

OctTree::OctTree(const Tuple3i& highestRes, int nlevel):nLevels_(nlevel)
{
    levels_.resize(nlevel);
    levelRes_.resize(nlevel);

    levelRes_[0] = highestRes;
    levels_[0].resize(boost::extents[highestRes.z][highestRes.y][highestRes.x]);
    zero_multi_array(levels_[0]);

    for(int lid = 1;lid < nlevel;++ lid)
    {
        if ( levelRes_[lid-1].x % 2 ||
             levelRes_[lid-1].y % 2 ||
             levelRes_[lid-1].z % 2 ||
             !levelRes_[lid-1].x ||
             !levelRes_[lid-1].y ||
             !levelRes_[lid-1].z )
        {
            fprintf(stderr, "wrong configuration for resolution and levels\n");
            SHOULD_NEVER_HAPPEN(1);
        }

        levelRes_[lid] = levelRes_[lid-1] / 2;
        levels_[lid].resize(boost::extents[levelRes_[lid].z][levelRes_[lid].y][levelRes_[lid].x]);
        zero_multi_array(levels_[lid]);
    }
}

OctTree::~OctTree()
{
    for(int lid = nLevels_-1;lid >= 0;-- lid)
    {
        for(int iz = 0;iz < levelRes_[lid].z;++ iz)
        for(int iy = 0;iy < levelRes_[lid].y;++ iy)
        for(int ix = 0;ix < levelRes_[lid].x;++ ix)
            if ( levels_[lid][iz][iy][ix] ) delete levels_[lid][iz][iy][ix];
    }
}

void OctTree::create_leaf(int ix, int iy, int iz)
{
    assert(ix >= 0 && ix < levelRes_[0].x);
    assert(iy >= 0 && iy < levelRes_[0].y);
    assert(iz >= 0 && iz < levelRes_[0].z);

    if ( !levels_[0][iz][iy][ix] ) 
        levels_[0][iz][iy][ix] = new Octant(true);
}

void OctTree::create_octree()
{
    printf("INFO: OctTree::creating octree ... ");
    for(int lid = 1;lid < nLevels_;++ lid)
    {
        for(int iz = 0;iz < levelRes_[lid].z;++ iz)
        for(int iy = 0;iy < levelRes_[lid].y;++ iy)
        for(int ix = 0;ix < levelRes_[lid].x;++ ix)
        {
            assert(!levels_[lid][iz][iy][ix]);

            // grid of children
            const int xx = ix << 1;
            const int yy = iy << 1;
            const int zz = iz << 1;

            int bitmask = 0;
            for(int i = 0;i < 8;++ i)
                if ( levels_[lid-1][zz+D_CORNERS[i][2]][yy+D_CORNERS[i][1]][xx+D_CORNERS[i][0]] )
                    bitmask |= 1<<i;

            if ( bitmask )
                levels_[lid][iz][iy][ix] = new Octant(false, bitmask);
        }
    }

    printf("OK\n");
}

void OctTree::weak_balance()
{
    printf("INFO: OctTree::weak balancing ... ");
    bool changed = true;

    while ( changed )
    {
        changed = false;

        for(int lid = nLevels_ - 1;lid > 1;-- lid)    // from top to bottom
        {
            for(int iz = 0;iz < levelRes_[lid].z;++ iz)
            for(int iy = 0;iy < levelRes_[lid].y;++ iy)
            for(int ix = 0;ix < levelRes_[lid].x;++ ix)
            {
                if ( !levels_[lid][iz][iy][ix] ) continue;

                // check if we should create any of the child
                for(int i = 0;i < 8;++ i)
                {
                    // already has that child
                    if ( levels_[lid][iz][iy][ix]->childrenMask & (1 << i) ) continue;

                    const int xx = (ix << 1) + D_CORNERS[i][0];
                    const int yy = (iy << 1) + D_CORNERS[i][1];
                    const int zz = (iz << 1) + D_CORNERS[i][2];

                    if ( check_for_weak_balance(lid-1, xx, yy, zz) )
                    {
                        assert(!levels_[lid-1][zz][yy][xx]);
                        levels_[lid-1][zz][yy][xx] = new Octant(false);
                        levels_[lid][iz][iy][ix]->childrenMask |= (1 << i);

                        changed = true;
                        break;
                    }
                } // end for each child
            }
        }
    } // end while 

    printf("OK\n");
}

/*
 * Check to see if there is cell at level lid is adjoining a cell sharing the
 * face/edge/vertex
 */
bool OctTree::check_for_weak_balance(int lid, int ix, int iy, int iz) const
{
    // check for 6 faces
    for(int idir = 0;idir < 3;++ idir)
    for(int j = 0, d = -1;j < 2;++ j, d *= -1)
    {
        Tuple3i pos(ix, iy, iz);
        pos[idir] += d;

        if ( pos[idir] >= 0 && pos[idir] < levelRes_[lid][idir] &&
             levels_[lid][pos.z][pos.y][pos.x] )
        {
            const Octant* poct = levels_[lid][pos.z][pos.y][pos.x];
            assert(poct);

            const int opsface = (idir << 1) + ((1 - d) >> 1);    // -1 --> 1  1 --> 0 
            if ( poct->childrenMask & FACE_MASK[opsface] ) return true;
        }
    }

    // 12 edge sharing diagonal cell
    for(int i = 0;i < 12;++ i)
    {
        const Tuple3i pos(ix + EDGE_DIAG_DIR[i][0], 
                          iy + EDGE_DIAG_DIR[i][1], 
                          iz + EDGE_DIAG_DIR[i][2]);
        if ( pos.x >= 0 && pos.x < levelRes_[lid].x &&
             pos.y >= 0 && pos.y < levelRes_[lid].y &&
             pos.z >= 0 && pos.z < levelRes_[lid].z &&
             levels_[lid][pos.z][pos.y][pos.x] &&
             (levels_[lid][pos.z][pos.y][pos.x]->childrenMask & 
                IDX_EDGES_MASK[EDGE_DIAG_EDGE[i]]) )
            return true;
    }

    return false;
}

void OctTree::check_children_mask() const
{
    for(int lid = nLevels_-1;lid > 0;-- lid)
    {
        for(int iz = 0;iz < levelRes_[lid].z;++ iz)
        for(int iy = 0;iy < levelRes_[lid].y;++ iy)
        for(int ix = 0;ix < levelRes_[lid].x;++ ix)
        {
            if ( !levels_[lid][iz][iy][ix] ) continue;

            const int xx = ix << 1;
            const int yy = iy << 1;
            const int zz = iz << 1;
            for(int i = 0;i < 8;++ i)
                if ( levels_[lid][iz][iy][ix]->childrenMask & (1 << i) )
                {
                    if ( !levels_[lid-1][zz+D_CORNERS[i][2]][yy+D_CORNERS[i][1]][xx+D_CORNERS[i][0]] )
                        SHOULD_NEVER_HAPPEN(1);
                }
                else
                {
                    if ( levels_[lid-1][zz+D_CORNERS[i][2]][yy+D_CORNERS[i][1]][xx+D_CORNERS[i][0]] )
                        SHOULD_NEVER_HAPPEN(1);
                }
        }
    }
}

void OctTree::update_bitmasks()
{
    for(int lid = nLevels_-1;lid > 0;-- lid)
    {
        for(int iz = 0;iz < levelRes_[lid].z;++ iz)
        for(int iy = 0;iy < levelRes_[lid].y;++ iy)
        for(int ix = 0;ix < levelRes_[lid].x;++ ix)
        {
            if ( !levels_[lid][iz][iy][ix] ) continue;

            // update the mask based on the children layout
            for(int i = 0;i < 8;++ i)
                if ( levels_[lid][iz][iy][ix]->childrenMask & (1 << i) )
                {   // has i-th child

                    levels_[lid][iz][iy][ix]->faceCenterMask |= CHILD_CENTERED_MASK[i][0];
                    levels_[lid][iz][iy][ix]->edgeMidMask    |= CHILD_CENTERED_MASK[i][1];

                    // in 3 directions sharing the face
                    for(int idir = 0;idir < 3;++ idir)
                    {

                        Tuple3i nowpos(ix, iy, iz);
                        nowpos[idir] += D_COVTX[i][idir];

                        if ( nowpos[idir] >= 0 && nowpos[idir] < levelRes_[lid][idir] &&
                             levels_[lid][nowpos.z][nowpos.y][nowpos.x] )
                        {
                            levels_[lid][nowpos.z][nowpos.y][nowpos.x]->faceCenterMask |= 
                                    CHILD_ADJOIN_CENTERED_MASK[i][idir][0];
                            levels_[lid][nowpos.z][nowpos.y][nowpos.x]->edgeMidMask |=
                                    CHILD_ADJOIN_CENTERED_MASK[i][idir][1];
                        }
                    }

                    // in 3 direction sharing the edge
                    for(int idir = 0;idir < 3;++ idir)
                    {
                        const int teid = CORNER_EDGES[i][idir];
                        const Tuple3i nowpos(ix + EDGE_DIAG_DIR[teid][0],
                                             iy + EDGE_DIAG_DIR[teid][1],
                                             iz + EDGE_DIAG_DIR[teid][2]);
                        const int vj1 = (idir+1) % 3, vj2 = (idir+2) % 3;

                        if ( nowpos[vj1] >= 0 && nowpos[vj1] < levelRes_[lid][vj1] &&
                             nowpos[vj2] >= 0 && nowpos[vj2] < levelRes_[lid][vj2] &&
                             levels_[lid][nowpos.z][nowpos.y][nowpos.x] )
                            levels_[lid][nowpos.z][nowpos.y][nowpos.x]->edgeMidMask |= 
                                    (1 << EDGE_DIAG_EDGE[teid]);
                    }
                }
        }
    } // end level
}

}

