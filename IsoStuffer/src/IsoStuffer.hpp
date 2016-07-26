/******************************************************************************
 *  File: IsoStuffer.hpp
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
#ifndef ISOSTUFFING_STUFFER_HPP
#   define ISOSTUFFING_STUFFER_HPP

#include "config.h"
#include <stdlib.h>
#include <map>
#include <set>
#include <queue>
#include "utils/macros.h"
#include "utils/arrays.hpp"

#ifdef USE_HASH_MAP
#   include <unordered_map>
#   include <unordered_set>
#else
#   include <tr1/unordered_map>
#   include <tr1/unordered_set>
#endif

#include "tables.h"
#include "OctTree.h"
#include "linearalgebra/Tuple4.hpp"
#include "geometry/Point3.hpp"
#include "geometry/TetMesh.hpp"

namespace IsoStuffing
{

template <class DistFunc>
class IsoStuffer
{
    public:
        IsoStuffer(const Tuple3i& highestRes, int nlevel, double sz, 
                const Point3d& minCorner, double al, double as, DistFunc& df):
                gridSize_(sz), gs2_(sz*sz), minPoint_(minCorner), 
                highestRes_(highestRes), octTree_(highestRes, nlevel),
                isoValues_(boost::extents[highestRes.z+1][highestRes.y+1][highestRes.x+1]),
                centerIsoValues_(boost::extents[highestRes.z][highestRes.y][highestRes.x]),
                gridPts_(boost::extents[highestRes.z+1][highestRes.y+1][highestRes.x+1]),
                centerPts_(boost::extents[highestRes.z][highestRes.y][highestRes.x]),
                alphaLong_(al), alphaShort_(as), distfunc_(df)
        { 
            zero_multi_array(isoValues_);
            zero_multi_array(centerIsoValues_);

            memset(gridPts_.data(),   0xFF, sizeof(int)*gridPts_.num_elements());
            memset(centerPts_.data(), 0xFF, sizeof(int)*centerPts_.num_elements());
        }
                
        void create_mesh(TetMesh<double>* pmesh);

        /*
         * step 1: create the octants near the surface
         */
        void create_boundary_octants();


        /*
         * step 2: create an oct tree from bottom to up
         *         call octTree_.create_octree()
         * NOTE: this method is only for testing
         */
        void create_octree()
        {   octTree_.create_octree(); }

        /*
         *
         * step 3: make sure the oct tree satisfies the weak balance condition
         *         call octTree_.weak_balance()
         * NOTE: this method is only for testing
         */
        void weak_balance()
        {   octTree_.weak_balance(); }

        /*
         * step 2.8: go over each cell on the top level of the OctTree, if that
         *           cell is entirely inside of mesh, create an octant if that
         *           cell is empty
         */
        void fill_highest_level();

        /*
         *
         * step 3.8: update the faceCenterMask and edgeMidMask in OctTree
         *         call octTree_.update_bitmasks()
         * NOTE: this method is only for testing
         */
        void update_bitmasks()
        {   octTree_.update_bitmasks(); }


        /*
         * step 4: create the tets in background
         */
        void create_background_tets();

        /*
         * step 5: for each breaking edge, computing the cutting point
         */
        void compute_cutting_pts();

        /*
         * step 6: moving the violated vertex
         */
        void wrap_violated_pts();
        /*
         * step 7: breaking the objects using pattern match
         */
        void extract_final_tets();

        void extract_final_mesh(TetMesh<double>* pmesh) const;

        void background_tet_mesh(TetMesh<double>* pmesh);

        const OctTree& oct_tree() const
        {   return octTree_; }

        double smallest_grid_size() const
        {   return gridSize_; }

        const Point3d& min_corner() const
        {   return minPoint_; }

    private:
        void update_cutting_points();

        void wrap_background_tets();

        void final_triangulate();

        /*
         * get the iso value at the finest cell labeled with (x,y,z)
         */
        inline double get_iso_value(int x, int y, int z);
        inline double get_center_iso_value(const Tuple3i& pos);
        inline double get_center_iso_value(int x, int y, int z);

        /*
         * NOTE: whenever a vertex is added into vtx_, the iso value associated 
         * with that vertex is also added into iso_;
         */
        int grid_vtx_idx(int x, int y, int z);
        int center_vtx_idx(int x, int y, int z);

        void face_center_vtx(int lid, int ix, int iy, int iz, int faceid,
                             std::pair<int, double>& out);
        void cell_center_vtx(int lid, int ix, int iy, int iz, 
                             std::pair<int, double>& out); 
        void edge_mid_vtx(int lid, int ix, int iy, int iz, int edgeid, 
                          std::pair<int, double>& out);

        /* ------------ create background tets ------------- */
        void create_half_pyramids(int idir, Tuple3i& pos, int cid, 
                                  double ciso, const int* vids, 
                                  const double* viso, const int* fids);

        void create_bisected_bbc_tets(int lid, int ix, int iy, int iz, int faceid, 
                                      const std::pair<int, double>& cidiso, 
                                      const std::pair<int, double>& fidiso,
                                      const int* vids, const double* fiso);

        /* create tets cross two octants */
        void create_cross_tets_lowest(int cid1, double iso1, int cid2, double iso2,
                                      const int* vids, const double* fiso);
        void create_cross_tets_higher(int lid, int ix, int iy, int iz, int faceid, 
                                      const std::pair<int, double>& cidiso, 
                                      const std::pair<int, double>& fidiso,
                                      const int* vids, const double* fiso);

        /*
         * do BFS search the finest cell, create the lowest level octants
         */
        void bfs_boundary_octants(int ix, int iy, int iz);

        bool stencil_match(int tetid, int vt, int vl, int vm, int vr);

        inline void bisect_tets(int a, int b, int c, int d, int e);

        inline void create_tet_type1(int a, int b, int c, int d, int e, int wz);
        inline void create_tet_type2(int v1, int v3, int c0, int c1, 
                                     int c2, int c3, bool wz);

        inline bool been_added(int id1, int id2) 
        {
            return id1 < id2 ? crossSet_.count(id1) && crossSet_[id1].count(id2) :
                               crossSet_.count(id2) && crossSet_[id2].count(id1);
        }

        inline void add_cross_pair(int id1, int id2)
        {
            if ( id1 < id2 ) 
                crossSet_[id1].insert(id2);
            else
                crossSet_[id2].insert(id1);
        }

        inline void add_breaking_tet(int vid0, int vid1, int vid2, int vid3)
        {
            Tuple4i tidx(vid0, vid1, vid2, vid3);
            breakingTet_.push_back(tidx);

            //// store the long edge bitmask
            //   short edge should be 0.75*gridSize^2, long edge should be gridSize^2
            int bit = 0;
            for(int i = 0;i < 4;++ i)
            for(int j = i+1;j < 4;++ j)
            {
                const double ls2 = (vtx_[tidx[i]] - vtx_[tidx[j]]).lengthSqr();
                assert(ls2 > 0.74*gs2_ && ls2 < 1.01*gs2_);
                if (  ls2 > 0.8*gs2_ ) bit |= (1 << TET_PVTX_EDGE[i][j]);
            }
            longEdgeMask_.push_back(bit);
        }

        inline void cache_iso_and_vtx(int fvid, int ix, int iy, int iz,
                                      double* fiso, int* vids)
        {
            for(int t = 0;t < 4;++ t)
            {
                int tx = ix + D_CORNERS[IDX_FACES[fvid][t]][0];
                int ty = iy + D_CORNERS[IDX_FACES[fvid][t]][1];
                int tz = iz + D_CORNERS[IDX_FACES[fvid][t]][2];

                fiso[t] = get_iso_value(tx, ty, tz);
                vids[t] = grid_vtx_idx(tx, ty, tz);
            }
        }

        inline void cache_iso_and_vtx(int lid, int fvid, int ix, int iy, int iz,
                                      double* fiso, int* vids)
        {
            const int SZ = 1 << lid;
            ix <<= lid;
            iy <<= lid;
            iz <<= lid;

            for(int t = 0;t < 4;++ t)
            {
                int tx = ix + D_CORNERS[IDX_FACES[fvid][t]][0] * SZ;
                int ty = iy + D_CORNERS[IDX_FACES[fvid][t]][1] * SZ;
                int tz = iz + D_CORNERS[IDX_FACES[fvid][t]][2] * SZ;

                fiso[t] = get_iso_value(tx, ty, tz);
                vids[t] = grid_vtx_idx(tx, ty, tz);
            }
        }

        /*
         * computing the cutting point between vid1 and vid2
         * NOTE: make sure vid1 < vid2, iso[vid1]*iso[vid2] < 0
         */
        inline void add_cutting_point(int vid1, int vid2, bool islongedge) 
        {
            assert(iso_[vid1]*iso_[vid2] < 0);
            if ( vid1 > vid2 ) std::swap(vid1, vid2);

            if ( cuttingPts_.count(vid1) && cuttingPts_[vid1].count(vid2) ) return;
            double t = fabs(iso_[vid1] / (iso_[vid1] - iso_[vid2]));
            Point3d pt = vtx_[vid1];
            pt.scaleAdd(t, vtx_[vid2] - vtx_[vid1]); // vtx[id1] + t*(vtx[id2] - vtx[id1])

            vtx_.push_back(pt);
            cuttingPts_[vid1][vid2] = (int)vtx_.size()-1;
            //m_cuttingPts[vid2][vid2] = std::make_pair((int)vtx_.size()-1, 1.-t);

            // check if this edge is violated
            const double thresh = islongedge ? alphaLong_ : alphaShort_;
            if ( t < thresh )
                violatedEdges_[vid1][vid2] = t;
            else if ( t > 1. - thresh )
                violatedEdges_[vid2][vid1] = 1. - t;
        }

        inline int get_cutting_point_idx(int vid1, int vid2)
        {
            if ( vid1 > vid2 ) std::swap(vid1, vid2);
            assert(cuttingPts_.count(vid1) && cuttingPts_[vid1].count(vid2));
            return cuttingPts_[vid1][vid2];
        }

    private:
        typedef boost::multi_array<double, 3>       TLatticed;
        typedef boost::multi_array<bool, 3>         TLatticeb;
        typedef boost::multi_array<int, 3>          TLatticei;

        const double            gridSize_;
        const double            gs2_;           // square of gridSize_
        const Point3d           minPoint_;
        const Tuple3i           highestRes_;    // the same as OctTree.levelRes_[0]

        OctTree                 octTree_;

        /* 
         * the cached iso values at each node of the (finest) octants
         * isoValue_ are all initialized into zero
         *
         * if a iso value is exactly zero, it means the iso value at the 
         * corresponding points hasn't been cached
         */
        TLatticed               isoValues_;
        TLatticed               centerIsoValues_;

        // the cached iso values at each the center of finest octant
        TLatticei               gridPts_;       // index for the grid points
        TLatticei               centerPts_;     // index for the center points of each grid

        std::vector<Point3d>    vtx_;
        std::vector<double>     iso_;
        std::vector<Tuple4i>    tet_;

        std::vector<Tuple4i>    breakingTet_;
        std::vector<int>        longEdgeMask_;

        /*
         * violatedEdges_[i] maintains the violated edges from i to the adjoining pts,
         * point i is violated.
         */
        std::map< int, std::map<int, double> >      violatedEdges_;

#ifdef USE_HASH_MAP
        std::unordered_map< int, std::unordered_set<int> >        crossSet_;
        std::unordered_map< int, std::unordered_map<int, int> >   cuttingPts_;
#else
        std::tr1::unordered_map< int, std::tr1::unordered_set<int> >        crossSet_;
        std::tr1::unordered_map< int, std::tr1::unordered_map<int, int> >   cuttingPts_;
#endif
        // two parameters alpha_{long} and alpha_{short} (see details in the paper)
        const double            alphaLong_;
        const double            alphaShort_;

        DistFunc&               distfunc_;
};

///////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION
/*
 *            v3
 *           /| \
 *          / |  \
 *         /  |   \  
 *      e3/f0 |    \e5
 *       /    |     \
 *      /     |e4    \  
 *     /  e2  |       \
 *    v0......|........v2
 *     \      |       /
 *       \    | f1  /  
 *      e0 \  |   /e1
 *           \| / 
 *            v1
 */


template <class DistFunc>
void IsoStuffer<DistFunc>::create_boundary_octants()
{
    printf("INFO: creating boundary octants ... ");
    // sweep the grid until we find a cell that is crossing the boundary
    for(int iz = 0;iz < highestRes_.z;++ iz)
    for(int iy = 0;iy < highestRes_.y;++ iy)
    for(int ix = 0;ix < highestRes_.x;++ ix)
    {
        int v = 0;
        for(int i = 0;i < 8;++ i)
        {
            double isoval = get_iso_value(
                    ix+D_CORNERS[i][0], 
                    iy+D_CORNERS[i][1], 
                    iz+D_CORNERS[i][2]);
            int    sign   = M_SIGN(isoval);
            assert(sign);   // sign should not be zero

            if ( v*sign >= 0 ) 
                v += sign;
            else 
            {
                // bfs all the octants cross boundary
                bfs_boundary_octants(ix, iy, iz);
                printf("OK\n");
                return;
            }
        }
    }

    fprintf(stderr, "ERROR: no cell crossing boundary is found!");
    SHOULD_NEVER_HAPPEN(1);
}

/*
 * if a leaf octant o has a corner vertex v whose sign
 * is opposite the sign of o’s center point, or if either sign is
 * zero, then we must create the three leaf octants incident
 * on v that share a square face with o
 */
template <class DistFunc>
void IsoStuffer<DistFunc>::bfs_boundary_octants(int ix, int iy, int iz)
{
    using namespace std;

    queue<Tuple3i> que;
    TLatticeb visited(boost::extents[highestRes_.z][highestRes_.y][highestRes_.x]);
    zero_multi_array(visited);
    
    que.push(Tuple3i(ix, iy, iz));
    visited[iz][iy][ix] = true;
    octTree_.create_leaf(ix, iy, iz);

    double ivs[8];
    bool   cross[12];   // indicate if there is a cutting point on each of 12 edges 

    while ( !que.empty() )
    {
        Tuple3i nowpos = que.front();
        que.pop();

        double isovalue = get_center_iso_value(nowpos);

        for(int i = 0;i < 8;++ i)
        {
            ivs[i]= get_iso_value(nowpos.x+D_CORNERS[i][0], 
                                  nowpos.y+D_CORNERS[i][1], 
                                  nowpos.z+D_CORNERS[i][2]);

            if ( isovalue * ivs[i] <= 0. ) 
            {   
                // NOTE: here nowpos.x+D_COVTX[i][0], nowpos.y+D_.., nowpos.z+D_.. should NOT be out of range
                // create three octants 
                octTree_.create_leaf(nowpos.x+D_COVTX[i][0], nowpos.y, nowpos.z);
                octTree_.create_leaf(nowpos.x, nowpos.y+D_COVTX[i][1], nowpos.z);
                octTree_.create_leaf(nowpos.x, nowpos.y, nowpos.z+D_COVTX[i][2]);
            }
        }

        for(int i = 0;i < 12;++ i)
            cross[i] = ivs[IDX_EDGES[i][0]] * ivs[IDX_EDGES[i][1]] < 0;

        for(int i = 0;i < 3;++ i)
        for(int j = 0, dir = -1;j < 2;++ j, dir *= -1)
        {
            Tuple3i newp = nowpos;
            newp[i] += dir;

            if ( newp[i] >= 0 && newp[i] < highestRes_[i] && !visited[newp.z][newp.y][newp.x] )
            {
                const int fid = (i<<1) + j;
                int k;
                for(k = 0;k < 4;++ k) if ( cross[IDX_EDGE_FACES[fid][k]] ) break;

                if ( k < 4 )
                {
                    visited[newp.z][newp.y][newp.x] = true;
                    octTree_.create_leaf(newp.x, newp.y, newp.z);
                    que.push(newp);
                }
            }
        }
    } // end while 
}

static inline int bitmask(double iso1, double iso2, double iso3, double iso4)
{
    return (iso1 > 0) << 3 | (iso2 > 0) << 2 | (iso3 > 0) << 1 | (iso4 > 0);
}

template<class DistFunc>
void IsoStuffer<DistFunc>::create_background_tets()
{
    double fiso[4]; 
    int    vids[4];

    //// go over the leaf nodes
    for(int iz = 0;iz < highestRes_.z;++ iz)
    for(int iy = 0;iy < highestRes_.y;++ iy)
    for(int ix = 0;ix < highestRes_.x;++ ix)
    {
        if ( !octTree_.levels_[0][iz][iy][ix] ) continue;

        const double iso1 = get_center_iso_value(ix, iy, iz);
        const int    cid1 = center_vtx_idx(ix, iy, iz); 
        // iterate on each face
        for(int idir = 0;idir < 3;++ idir)
        for(int j = 0, dir = -1;j < 2;++ j, dir *= -1)
        {
            const int faceid = (idir << 1) + j;
            Tuple3i nowpos(ix, iy, iz);
            nowpos[idir] += dir;

            if ( nowpos[idir] >= 0 && nowpos[idir] < highestRes_[idir] )
            {
                if ( octTree_.levels_[0][nowpos.z][nowpos.y][nowpos.x] )
                {   // face is shared with an adjoining octant with the same size
                    const int cid2 = center_vtx_idx(nowpos.x, nowpos.y, nowpos.z);
                    if ( been_added(cid1, cid2) ) continue;

                    // cache iso value and vertex ids
                    cache_iso_and_vtx(faceid, ix, iy, iz, fiso, vids);
                    const double iso2 = get_center_iso_value(nowpos);

                    create_cross_tets_lowest(cid1, iso1, cid2, iso2, vids, fiso);
                    add_cross_pair(cid1, cid2);
                }
                else // face is shared with larger octant
                {
                    // cache iso value and vertex ids
                    cache_iso_and_vtx(faceid, ix, iy, iz, fiso, vids);
#if defined(DEBUG) | defined(_DEBUG)
                    // sanity check when debugging
                    int bbit = (iso1<0) << 4 | (fiso[3]<0) << 3 | (fiso[2]<0) << 2 |
                               (fiso[1]<0) << 1 | (fiso[0]<0);
                    assert(bbit==0 || bbit == 31);
#endif
                    create_half_pyramids(idir, nowpos, cid1, iso1, 
                                         vids, fiso, IDX_FACES[faceid]);
                }
            }
        } // end for 6 faces
    }

    ////// Create tets for higher level OctTree nodes
    std::pair<int, double> cidiso, fidiso;
    for(int lid = 1;lid < octTree_.nLevels_;++ lid)
    {
        for(int iz = 0;iz < octTree_.levelRes_[lid].z;++ iz)
        for(int iy = 0;iy < octTree_.levelRes_[lid].y;++ iy)
        for(int ix = 0;ix < octTree_.levelRes_[lid].x;++ ix)
        {
            if ( !octTree_.levels_[lid][iz][iy][ix] ) continue;

            cell_center_vtx(lid, ix, iy, iz, cidiso);
            // iterate on each face
            for(int idir = 0;idir < 3;++ idir)
            for(int j = 0,dir = -1;j < 2;++ j, dir *= -1)
            {
                int faceid = (idir << 1) + j;   // faceid 0...5
                if ( octTree_.levels_[lid][iz][iy][ix]->faceCenterMask & (1<<faceid) ) 
                {   // there is a vertex at the center of this face
                    face_center_vtx(lid, ix, iy, iz, faceid, fidiso);
                    cache_iso_and_vtx(lid, faceid, ix, iy, iz, fiso, vids);
                    create_bisected_bbc_tets(lid, ix, iy, iz, faceid, 
                            cidiso, fidiso, vids, fiso);
                }
                else // no vertex at the center of faces
                {
                    Tuple3i nowpos(ix, iy, iz); // pos of neighbor
                    nowpos[idir] += dir;

                    if ( nowpos[idir] >= 0 && nowpos[idir] < octTree_.levelRes_[lid][idir] )
                    {
                        if ( octTree_.levels_[lid][nowpos.z][nowpos.y][nowpos.x] )
                        {   // face is shared with another octant with the same size
                            std::pair<int, double> cidiso2;
                            cell_center_vtx(lid, nowpos.x, nowpos.y, nowpos.z, cidiso2);
                            if ( been_added(cidiso.first, cidiso2.first) ) continue;

                            cache_iso_and_vtx(lid, faceid, ix, iy, iz, fiso, vids);
                            create_cross_tets_higher(lid, ix, iy, iz, faceid,
                                                     cidiso, cidiso2, vids, fiso);

                            add_cross_pair(cidiso.first, cidiso2.first);
                        }
                        else // face is shared with larger octant
                        {
                            cache_iso_and_vtx(lid, faceid, ix, iy, iz, fiso, vids);
                            create_half_pyramids(idir, nowpos, cidiso.first, cidiso.second,
                                                 vids, fiso, IDX_FACES[faceid]);
                        }
                    }
                }
            } // end for 6 faces
        }
    }
}

/*
 * NOTE: make sure
 * when this method is called, there is no center point at the specified face,
 * and the two center points in the two adjoining cell (cidiso1, cidiso2) have
 * never been connected.
 */
template <class DistFunc>
void IsoStuffer<DistFunc>::create_cross_tets_higher(
        int lid, int ix, int iy, int iz, int faceid,
        const std::pair<int, double>& cidiso1,
        const std::pair<int, double>& cidiso2,
        const int* vids, const double* fiso)
{
    std::pair<int, double> meidiso;
    const Octant* oct = octTree_.levels_[lid][iz][iy][ix];

    for(int t = 0;t < 4;++ t)
        if ( oct->edgeMidMask & (1 << IDX_EDGE_FACES[faceid][t]) )
        {   // there is a midpoint of this edge
            edge_mid_vtx(lid, ix, iy, iz, IDX_EDGE_FACES[faceid][t], meidiso);
            // create two tet
            int bit = bitmask(fiso[t], meidiso.second, cidiso1.second, cidiso2.second);
            if ( bit < 15 )
            {
                if ( !bit )
                    tet_.push_back(Tuple4i(vids[t], meidiso.first, 
                                cidiso1.first, cidiso2.first));
                else
                    ; //SHOULD_NEVER_HAPPEN(1);
            }

            bit = bitmask(meidiso.second, fiso[(t+1)%4], cidiso1.second, cidiso2.second);
            if ( bit < 15 )
            {
                if ( !bit )
                    tet_.push_back(Tuple4i(meidiso.first, vids[(t+1)%4], 
                                cidiso1.first, cidiso2.first));
                else
                    ; //SHOULD_NEVER_HAPPEN(1);
            }
        }
        else
        {
            int bit = bitmask(fiso[t], fiso[(t+1)%4], cidiso1.second, cidiso2.second);
            if ( bit == 15 ) continue;

            if ( !bit ) 
                tet_.push_back(Tuple4i(vids[t], vids[(t+1)%4], 
                            cidiso1.first, cidiso2.first));
            else
                ; //SHOULD_NEVER_HAPPEN(1);
        }
}

template <class DistFunc>
void IsoStuffer<DistFunc>::create_cross_tets_lowest(
        int cid1, double iso1, int cid2, double iso2,
        const int* vids, const double* fiso)
{
    for(int t = 0;t < 4;++ t)
    {
        int bit = bitmask(fiso[t], fiso[(t+1)%4], iso1, iso2);
        if ( bit == 15 ) continue;

        if ( bit > 0 ) 
            //NOTE that the red points are always the vtx[2] and vtx[3]
            add_breaking_tet(vids[t], vids[(t+1)%4], cid1, cid2);
        else
            tet_.push_back(Tuple4i(vids[t], vids[(t+1)%4], cid1, cid2));
    }
}

template <class DistFunc>
void IsoStuffer<DistFunc>::create_bisected_bbc_tets(
        int lid, int ix, int iy, int iz, int faceid, 
        const std::pair<int, double>& cidiso, 
        const std::pair<int, double>& fidiso,
        const int* vids, const double* fiso)
{
    std::pair<int, double> meidiso;
    const Octant* oct = octTree_.levels_[lid][iz][iy][ix]; 

    for(int t = 0;t < 4;++ t)
        if ( oct->edgeMidMask & (1 << IDX_EDGE_FACES[faceid][t]) )
        {   // there is a midpoint of this edge 
            edge_mid_vtx(lid, ix, iy, iz, IDX_EDGE_FACES[faceid][t], meidiso);

            if ( !(oct->childrenMask & (1 << IDX_FACES[faceid][t])) ) // there is no child with vertex t
            {
                // create the first tet (vids[t], fcid, meid, cid)
                int bit = bitmask(fiso[t], fidiso.second, meidiso.second, cidiso.second);
                if ( bit < 15 )
                {
                    if ( !bit )
                        tet_.push_back(Tuple4i(vids[t], fidiso.first, meidiso.first, cidiso.first));
                    else
                        ; //SHOULD_NEVER_HAPPEN(1);
                }
            }

            // create the second tet (vids[t+1], meid, fcid, cid)
            if ( !(oct->childrenMask & (1 << IDX_FACES[faceid][(t+1)%4])) )
            {
                int bit = bitmask(fiso[(t+1)%4], meidiso.second, fidiso.second, cidiso.second);
                if ( bit == 15 ) continue;
                if ( !bit )
                    tet_.push_back(Tuple4i(vids[(t+1)%4], meidiso.first, fidiso.first, cidiso.first));
                else
                    ; //SHOULD_NEVER_HAPPEN(1);
            }
        }
        else    // no vertex at midpoint
        {   // tet vids[t] fcid vids[t+1] cid
            const int bit = bitmask(fiso[t], fidiso.second, fiso[(t+1)%4], cidiso.second);
            if ( bit == 15 ) continue;
            if ( !bit )
                tet_.push_back(Tuple4i(vids[t], fidiso.first, vids[(t+1)%4], cidiso.first));
            else
                ; //SHOULD_NEVER_HAPPEN(1);
        }
}

/*
 * |   \ | /   |
 * |(10)\|/(11)|
 * +-----+-----+
 * |    /|\    |
 * |   / | \   |
 * |  /  |  \  |
 * | /   |   \ |
 * |/(00)|(10)\| 
 * +-----------+
 */
template <class DistFunc>
void IsoStuffer<DistFunc>::create_half_pyramids(
        int idir, Tuple3i& pos, int cid, double ciso,
        const int* vids, const double* viso, const int* fids)
{
    // determine the diagonal bisecting
    const int ids[2] = { (idir+1)%3, (idir+2)%3 };
    const int diagpt = (D_CORNERS[fids[0]][ids[0]] + D_CORNERS[fids[0]][ids[1]]) % 2 != 
                       ((pos[ids[0]] % 2) + (pos[ids[1]] % 2)) % 2;

    // create the first half pyramid (diagpt, diagpt+1, diagpt+2)
    int bit = bitmask(viso[diagpt], viso[diagpt+2], viso[diagpt+1], ciso);
    if ( bit < 15 ) // not all the vertices are outside of the mesh
    {
        if ( bit > 0 )
            add_breaking_tet(vids[diagpt], vids[diagpt+2], vids[diagpt+1], cid);
        else
            tet_.push_back(Tuple4i(vids[diagpt], vids[diagpt+2], vids[diagpt+1], cid));
    }

    // create the second half pyramid (diagpt, diagpt+2, (diagpt+3)%4)
    const int pt3 = (diagpt + 3) % 4;
    bit = bitmask(viso[diagpt], viso[pt3], viso[diagpt+2], ciso);
    if ( bit < 15 )
    {
        if ( bit > 0 )
            add_breaking_tet(vids[diagpt], vids[pt3], vids[diagpt+2], cid);
        else
            tet_.push_back(Tuple4i(vids[diagpt], vids[pt3], vids[diagpt+2], cid));
    }
}

template <class DistFunc>
inline double IsoStuffer<DistFunc>::get_iso_value(int x, int y, int z)
{
    if ( isoValues_[z][y][x] != 0 ) return isoValues_[z][y][x];

    Point3d pts(x, y, z);
    pts *= gridSize_;
    pts += minPoint_;
#if defined(DEBUG) | defined(_DEBUG)
    isoValues_[z][y][x] = distfunc_(pts);
    assert(isoValues_[z][y][x] > 1E-12 || isoValues_[z][y][x] < -1E-12);
    return isoValues_[z][y][x];
#else
    //return (isoValues_[z][y][x] = distfunc_(pts)+0.2);       ///// HERE For TIM
    return (isoValues_[z][y][x] = distfunc_(pts));
#endif
}

template <class DistFunc>
inline double IsoStuffer<DistFunc>::get_center_iso_value(const Tuple3i& pos)
{
    if ( centerIsoValues_[pos.z][pos.y][pos.x] != 0 ) return centerIsoValues_[pos.z][pos.y][pos.x];

    Point3d pts(pos);
    pts += 0.5;
    pts *= gridSize_;
    pts += minPoint_;
#if defined(DEBUG) | defined(_DEBUG)
    centerIsoValues_[pos.z][pos.y][pos.x] = distfunc_(pts);
    assert(centerIsoValues_[pos.z][pos.y][pos.x] > 1E-12 || centerIsoValues_[pos.z][pos.y][pos.x] < -1E-12);
    return centerIsoValues_[pos.z][pos.y][pos.x];
#else
    //return (centerIsoValues_[pos.z][pos.y][pos.x] = distfunc_(pts)+0.2);       ///// HERE For TIM
    return (centerIsoValues_[pos.z][pos.y][pos.x] = distfunc_(pts));
#endif
}

template <class DistFunc>
inline double IsoStuffer<DistFunc>::get_center_iso_value(int x, int y, int z)
{
    if ( centerIsoValues_[z][y][x] != 0 ) return centerIsoValues_[z][y][x];

    Point3d pts(x, y, z);
    pts += 0.5;
    pts *= gridSize_;
    pts += minPoint_;
#if defined(DEBUG) | defined(_DEBUG)
    centerIsoValues_[z][y][x] = distfunc_(pts);
    assert(centerIsoValues_[z][y][x] > 1E-12 || centerIsoValues_[z][y][x] < -1E-12);
    return centerIsoValues_[z][y][x];
#else
    //return (centerIsoValues_[z][y][x] = distfunc_(pts)+0.2);       ///// HERE For TIM
    return (centerIsoValues_[z][y][x] = distfunc_(pts));
#endif
}

template <class DistFunc>
int IsoStuffer<DistFunc>::grid_vtx_idx(int x, int y, int z)
{
    if ( gridPts_[z][y][x] >= 0 ) return gridPts_[z][y][x];

    Point3d pts(x, y, z);
    pts *= gridSize_;
    pts += minPoint_;

    vtx_.push_back(pts);
    iso_.push_back(get_iso_value(x, y, z));
    return ( gridPts_[z][y][x] = (int)vtx_.size() - 1 );
}

template <class DistFunc>
int IsoStuffer<DistFunc>::center_vtx_idx(int x, int y, int z)
{
    if ( centerPts_[z][y][x] >= 0 ) return centerPts_[z][y][x];

    Point3d pts(x, y, z);
    pts += 0.5;
    pts *= gridSize_;
    pts += minPoint_;

    vtx_.push_back(pts);
    iso_.push_back(get_center_iso_value(x, y, z));
    return ( centerPts_[z][y][x] = (int)vtx_.size() - 1 );
}

template <class DistFunc>
void IsoStuffer<DistFunc>::cell_center_vtx(
        int lid, int ix, int iy, int iz, std::pair<int, double>& out)
{
    assert(lid > 0);
    const int SZ = 1 << (lid - 1);

    ix = (ix << lid) + SZ;
    iy = (iy << lid) + SZ;
    iz = (iz << lid) + SZ;
    
    out.first = grid_vtx_idx(ix, iy, iz);
    out.second = get_iso_value(ix, iy, iz);
}

template <class DistFunc>
void IsoStuffer<DistFunc>::face_center_vtx(
        int lid, int ix, int iy, int iz, int faceid,
        std::pair<int, double>& out)
{
    assert(lid > 0);

    const int SZ  = 1 << lid;
    const int HSZ = 1 << (lid - 1);
    const int idir = faceid / 2;
    const int offset = faceid % 2;

    Tuple3i tpos(ix << lid, iy << lid, iz << lid);
    tpos[idir] += offset*SZ;
    tpos[(idir+1)%3] += HSZ;
    tpos[(idir+2)%3] += HSZ;

    out.first = grid_vtx_idx(tpos.x, tpos.y, tpos.z);
    out.second = get_iso_value(tpos.x, tpos.y, tpos.z);
}

template <class DistFunc>
void IsoStuffer<DistFunc>::edge_mid_vtx(
        int lid, int ix, int iy, int iz, int edgeid, 
        std::pair<int, double>& out)
{
    assert(lid > 0);
    const int SZ = 1 << lid;
    ix <<= lid;
    iy <<= lid;
    iz <<= lid;

    Tuple3i tgd0(ix + SZ*D_CORNERS[IDX_EDGES[edgeid][0]][0],
                 iy + SZ*D_CORNERS[IDX_EDGES[edgeid][0]][1],
                 iz + SZ*D_CORNERS[IDX_EDGES[edgeid][0]][2]);

    Tuple3i tgd1(ix + SZ*D_CORNERS[IDX_EDGES[edgeid][1]][0],
                 iy + SZ*D_CORNERS[IDX_EDGES[edgeid][1]][1],
                 iz + SZ*D_CORNERS[IDX_EDGES[edgeid][1]][2]);

    Tuple3i gd = (tgd0 + tgd1)/2;
    out.first = grid_vtx_idx(gd.x, gd.y, gd.z);
    out.second = get_iso_value(gd.x, gd.y, gd.z);
}

template <class DistFunc>
void IsoStuffer<DistFunc>::fill_highest_level()
{
    printf("INFO: filling highest level of oct-tree ... ");
    const int L  = octTree_.nLevels_ - 1;
    const int SZ = 1 << L;

    for(int iz = 0;iz < octTree_.levelRes_[L].z;++ iz)
    for(int iy = 0;iy < octTree_.levelRes_[L].y;++ iy)
    for(int ix = 0;ix < octTree_.levelRes_[L].x;++ ix)
    {
        if ( octTree_.levels_[L][iz][iy][ix] ) continue;

        const int xx = ix << L;
        const int yy = iy << L;
        const int zz = iz << L;
        for(int i = 0;i < 8;++ i)
            if ( get_iso_value(xx + SZ*D_CORNERS[i][0],
                               yy + SZ*D_CORNERS[i][1],
                               zz + SZ*D_CORNERS[i][2]) > 0 )
                goto NEXT_CELL;
        octTree_.levels_[L][iz][iy][ix] = new Octant(false);
NEXT_CELL:
        ;
    }
    printf("OK\n");
}

template <class DistFunc>
void IsoStuffer<DistFunc>::background_tet_mesh(TetMesh<double>* pmesh)
{
    for(size_t i = 0;i < vtx_.size();++ i)
        pmesh->add_vertex(vtx_[i]);
    for(size_t i = 0;i < tet_.size();++ i)
    {
        //if ( vtx_[tet_[i][0]].z < -0. && vtx_[tet_[i][1]].z < -0. &&
        //     vtx_[tet_[i][2]].z < -0. && vtx_[tet_[i][3]].z < -0. )
        pmesh->add_tet(tet_[i][0], tet_[i][1], tet_[i][2], tet_[i][3]);
    }
    for(size_t i = 0;i < breakingTet_.size();++ i)
    {
        //if ( vtx_[breakingTet_[i][0]].z < -0. && vtx_[breakingTet_[i][1]].z < -0. &&
        //     vtx_[breakingTet_[i][2]].z < -0. && vtx_[breakingTet_[i][3]].z < -0. )
        pmesh->add_tet(breakingTet_[i][0], breakingTet_[i][1],
                       breakingTet_[i][2], breakingTet_[i][3]);
    }

    printf("INFO: background tet mesh: %d tets\n", (int)pmesh->num_tets());
    pmesh->init();
}

template <class DistFunc>
void IsoStuffer<DistFunc>::compute_cutting_pts()
{
    printf("INFO: computing cutting points ...");
    for(size_t tid = 0;tid < breakingTet_.size();++ tid)
    {
        for(int i = 0;i < 4;++ i)
        for(int j = i+1;j < 4;++ j)
            if ( iso_[breakingTet_[tid][i]] * iso_[breakingTet_[tid][j]] <= 0. )
                add_cutting_point(breakingTet_[tid][i], breakingTet_[tid][j], 
                        longEdgeMask_[tid] & (1 << TET_PVTX_EDGE[i][j]));
    }
    printf(" OK\n");
}

template <class DistFunc>
void IsoStuffer<DistFunc>::wrap_violated_pts()
{
    typedef std::map< int, std::map<int, double> >::const_iterator  VEIT;
    typedef std::map<int, double>::const_iterator                   SVIT;

    VEIT end = violatedEdges_.end();
    for(VEIT it = violatedEdges_.begin();it != end;++ it)
        if ( iso_[it->first] > 1E-12 )  // positive violated point
        {
            int otherend = -1;
            double tmin = 1E+10;
            SVIT end2 = it->second.end();
            for(SVIT it2 = it->second.begin();it2 != end2;++ it2)
            {
                assert(iso_[it2->first] < -1E-12);
                if ( !violatedEdges_.count(it2->first) && tmin > it2->second )
                {
                    otherend = it2->first;
                    tmin = it2->second;
                }
            }

            if ( otherend >= 0 )
            {
                int v1 = it->first, v2 = otherend;
                if ( v1 > v2 ) std::swap(v1, v2);
                vtx_[it->first] = vtx_[cuttingPts_[v1][v2]];
                iso_[it->first] = 0;
            }
        }

    for(VEIT it = violatedEdges_.begin();it != end;++ it)
        if ( iso_[it->first] < -1E-12 )  // positive violated point
        {
            int otherend = -1;
            double tmin = 1E+10;
            SVIT end2 = it->second.end();
            for(SVIT it2 = it->second.begin();it2 != end2;++ it2)
                if ( iso_[it2->first] > 1E-12 && tmin > it2->second )
                {
                    otherend = it2->first;
                    tmin = it2->second;
                }

            if ( otherend >= 0 )
            {
                int v1 = it->first, v2 = otherend;
                if ( v1 > v2 ) std::swap(v1, v2);
                vtx_[it->first] = vtx_[cuttingPts_[v1][v2]];
                iso_[it->first] = 0;
            }
        }
}

/*
 * positive->2
 * negative->1
 * zero    ->0
 */
static inline int encode_iso(double v)
{
    return v > 1E-12 ? 2 : (v < -1E-12 ? 1 : 0);
}

template <class DistFunc>
bool IsoStuffer<DistFunc>::stencil_match(int tetid, 
        int vt, int vl, int vm, int vr)
{
    const Tuple4i& tidx = breakingTet_[tetid];
    const int encode = encode_iso(iso_[tidx[vt]])*1000 +
                       encode_iso(iso_[tidx[vl]])*100  +
                       encode_iso(iso_[tidx[vm]])*10   +
                       encode_iso(iso_[tidx[vr]]);
    int vida, vidb, vidc;

    switch ( encode )
    {
        case 1111:      // all negative
        case 1110:      // 3 negative + 1 zero
        case 1101:
        case 1011:
        case  111:
        case 1100:      // 2 negative + 2 zero
        case   11:
        case 1010:
        case  101:
        case    1:      // 1 negative + 3 zero
        case   10:
        case  100:
        case 1000:
            tet_.push_back(tidx);   // add the entire tet
            return true;

        case 2220:      // 3 positive + zero
        case 2202:
        case    0:
        case 2022:
        case  222:
        case 2200:
        case 2020:
        case  202:
        case   22:
        case 2000:
        case  200:
        case   20:
        case    2:
            return true;

        case 2100:
            tet_.push_back(Tuple4i(tidx[vl], tidx[vm], tidx[vr], 
                        get_cutting_point_idx(tidx[vl], tidx[vt])));
            return true;

        case 2120:
            tet_.push_back(Tuple4i(tidx[vl], 
                        get_cutting_point_idx(tidx[vl], tidx[vm]),
                        tidx[vr], 
                        get_cutting_point_idx(tidx[vl], tidx[vt])));
            return true;

        case 1222:
            tet_.push_back(Tuple4i(
                        get_cutting_point_idx(tidx[vt], tidx[vl]),
                        get_cutting_point_idx(tidx[vt], tidx[vm]),
                        get_cutting_point_idx(tidx[vt], tidx[vr]),
                        tidx[vt]));
            return true;

        case 2122:
            tet_.push_back(Tuple4i(tidx[vl],
                        get_cutting_point_idx(tidx[vl], tidx[vm]),
                        get_cutting_point_idx(tidx[vl], tidx[vr]),
                        get_cutting_point_idx(tidx[vl], tidx[vt])));
            return true;

        case 2110:
            if ( longEdgeMask_[tetid] & (1 << TET_PVTX_EDGE[vt][vl]) )
            {
                vida = get_cutting_point_idx(tidx[vt], tidx[vl]);
                vidb = get_cutting_point_idx(tidx[vt], tidx[vm]);
                tet_.push_back(Tuple4i(tidx[vl], tidx[vm], tidx[vr], vida));
                tet_.push_back(Tuple4i(vida, tidx[vm], tidx[vr], vidb));
            }
            else if ( longEdgeMask_[tetid] & (1 << TET_PVTX_EDGE[vt][vm]) )
            {
                vida = get_cutting_point_idx(tidx[vt], tidx[vl]);
                vidb = get_cutting_point_idx(tidx[vt], tidx[vm]);
                tet_.push_back(Tuple4i(vida, vidb, tidx[vl], tidx[vr]));
                tet_.push_back(Tuple4i(vidb, tidx[vm], tidx[vl], tidx[vr]));
            }
            else
            {
                // the black edge has to be vl-vm
                vida = get_cutting_point_idx(tidx[vt], tidx[vl]);
                vidb = get_cutting_point_idx(tidx[vt], tidx[vm]);

                assert(((1<<vl)|(1<<vm)) == ((1<<2)|(1<<3)) ||
                       ((1<<vl)|(1<<vm)) == ((1<<1)|1));
                create_tet_type1(tidx[vl], tidx[vm], vidb, vida, tidx[vr],
                        ((1<<vl)|(1<<vm)) != ((1<<2)|(1<<3)));
            }
            return true;
        case 2111:
            if ( longEdgeMask_[tetid] & (1 << TET_PVTX_EDGE[vt][vm]) )
            {
                vidc = get_cutting_point_idx(tidx[vt], tidx[vm]);
                tet_.push_back(Tuple4i(tidx[vl], tidx[vm], tidx[vr], vidc));
                vida = get_cutting_point_idx(tidx[vt], tidx[vl]);
                vidb = get_cutting_point_idx(tidx[vt], tidx[vr]);

                assert(((1<<vr)|(1<<vl)) == ((1<<2)|(1<<3)) ||
                       ((1<<vr)|(1<<vl)) == ((1<<1)|1));
                create_tet_type1(tidx[vr], tidx[vl], vida, vidb, vidc,
                        ((1<<vr)|(1<<vl)) != ((1<<2)|(1<<3)));
                return true;
            }
            return false;

        case 2121:
            if ( longEdgeMask_[tetid] & (1 << TET_PVTX_EDGE[vt][vl]) )
            {
                assert(longEdgeMask_[tetid] & (1 << TET_PVTX_EDGE[vm][vr]));

                vida = get_cutting_point_idx(tidx[vt], tidx[vl]);
                vidb = get_cutting_point_idx(tidx[vm], tidx[vr]);
                tet_.push_back(Tuple4i(tidx[vl], 
                        get_cutting_point_idx(tidx[vl], tidx[vm]),
                        vidb, vida));
                tet_.push_back(Tuple4i(tidx[vl], vidb, tidx[vr], vida));
                tet_.push_back(Tuple4i(vida, vidb, tidx[vr],
                        get_cutting_point_idx(tidx[vt], tidx[vr])));
            }
            else if ( longEdgeMask_[tetid] & (1 << TET_PVTX_EDGE[vt][vr]) )
            {
                assert(longEdgeMask_[tetid] & (1 << TET_PVTX_EDGE[vl][vm]));

                vida = get_cutting_point_idx(tidx[vt], tidx[vr]);
                vidb = get_cutting_point_idx(tidx[vm], tidx[vl]);

                tet_.push_back(Tuple4i(vidb, 
                            get_cutting_point_idx(tidx[vm], tidx[vr]),
                            tidx[vr], vida));
                tet_.push_back(Tuple4i(tidx[vl], vidb, tidx[vr], vida));
                tet_.push_back(Tuple4i(tidx[vl], vidb, vida, 
                            get_cutting_point_idx(tidx[vt], tidx[vl])));
            }
            else
            {
                assert(longEdgeMask_[tetid] & (1 << TET_PVTX_EDGE[vt][vm]));
                assert(longEdgeMask_[tetid] & (1 << TET_PVTX_EDGE[vl][vr]));
                create_tet_type2(tidx[vl], tidx[vr],
                            get_cutting_point_idx(tidx[vt], tidx[vl]),
                            get_cutting_point_idx(tidx[vl], tidx[vm]),
                            get_cutting_point_idx(tidx[vm], tidx[vr]),
                            get_cutting_point_idx(tidx[vr], tidx[vt]),
                            ((1<<vl)|(1<<vr)) != ((1<<2)|(1<<3)));
            }
            return true;

        default:
            return false;
    }
}

static inline int point_greater_number(const Point3d& a, const Point3d& b)
{
    //return (a.x > b.x+1E-10) + (a.y > b.y+1E-10) + (a.z > b.z+1E-10);
    return (a.x > b.x) + (a.y > b.y) + (a.z > b.z);
}

/*
 *        d.
 *       /  \       
 *      /     c     
 *     /      |       
 *    a.......|........e
 *     \      |       /
 *       \    |     /  
 *         \  |   /  
 *           \| / 
 *            b
 *
 * The diagonal is assumed to be AC
 */
template <class DistFunc>
inline void IsoStuffer<DistFunc>::bisect_tets(int a, int b, int c, int d, int e)
{
    tet_.push_back(Tuple4i(a, c, b, e));
    tet_.push_back(Tuple4i(a, d, c, e));
}

/*
 *
 *        d.
 *       /  \       
 *      /     c     
 *     /      |       
 *    a.......|........e
 *     \      |       /
 *       \    |     /  
 *         \  |   /  
 *           \| / 
 *            b
 *
 * <ab> is black edge.
 *
 * \param wz if the black edge a-b is in Z not Z+(1/2,1/2,1/2)
 */
template <class DistFunc>
inline void IsoStuffer<DistFunc>::create_tet_type1(
        int a, int b, int c, int d, int e, int wz)
{
    int agc = point_greater_number(vtx_[a], vtx_[c]) % 2;    // 1: odd, 0: even
    assert(agc + (point_greater_number(vtx_[b], vtx_[d]) % 2) == 1);

    if ( (wz + agc) % 2 )   // in Z(1) and agc is even(0) OR in Z+1/2(0) and agc(1) is odd
        // choose bd as the diagonal
        bisect_tets(b, c, d, a, e);
    else
        // choose ac as the diagonal
        bisect_tets(a, b, c, d, e);
}

template <class DistFunc>
inline void IsoStuffer<DistFunc>::create_tet_type2(
        int v1, int v3, int c0, int c1, int c2, int c3, 
        bool wz)
{
    int agc = point_greater_number(vtx_[v3], vtx_[c0]) % 2;
    assert(agc + (point_greater_number(vtx_[v1], vtx_[c3]) % 2) == 1);

    if ( (wz + agc) % 2 )   // in Z and agc is even OR in Z+1/2 and agc is odd
    {   // choose bd as the diagonal
        tet_.push_back(Tuple4i(v1, c1, v3, c3));
        tet_.push_back(Tuple4i(c1, c2, v3, c3));
        tet_.push_back(Tuple4i(v1, c1, c3, c0));
    }
    else
    {
        tet_.push_back(Tuple4i(v1, c1, c2, c0));
        tet_.push_back(Tuple4i(v1, c2, v3, c0));
        tet_.push_back(Tuple4i(c0, c2, v3, c3));
    }
}

template <class DistFunc>
void IsoStuffer<DistFunc>::extract_final_tets()
{
    int ddd[3];
    for(int i = 0;i < (int)breakingTet_.size();++ i)
    {
        ddd[0] = 0; ddd[1] = 1; ddd[2] = 2;
        for(int j = 0;j < 3;++ j)
            if ( stencil_match(i, 3, ddd[j], ddd[(j+1)%3], ddd[(j+2)%3]) )
                goto NEXT_EXTRACT_TET;

        ddd[0] = 1; ddd[1] = 3; ddd[2] = 2;
        for(int j = 0;j < 3;++ j)
            if ( stencil_match(i, 0, ddd[j], ddd[(j+1)%3], ddd[(j+2)%3]) )
                goto NEXT_EXTRACT_TET;

        ddd[0] = 3; ddd[1] = 0; ddd[2] = 2;
        for(int j = 0;j < 3;++ j)
            if ( stencil_match(i, 1, ddd[j], ddd[(j+1)%3], ddd[(j+2)%3]) )
                goto NEXT_EXTRACT_TET;

        ddd[0] = 3; ddd[1] = 1; ddd[2] = 0;
        for(int j = 0;j < 3;++ j)
            if ( stencil_match(i, 2, ddd[j], ddd[(j+1)%3], ddd[(j+2)%3]) )
                goto NEXT_EXTRACT_TET;

        SHOULD_NEVER_HAPPEN(2);
NEXT_EXTRACT_TET:
        ;
    }
}

template <class DistFunc>
void IsoStuffer<DistFunc>::extract_final_mesh(TetMesh<double>* pmesh) const
{
    std::vector<int> idmap(vtx_.size());
    pmesh->clear();
    memset(&idmap[0], 0xFF, sizeof(int)*idmap.size());

    for(size_t i = 0;i < tet_.size();++ i)
    {
        for(int j = 0;j < 4;++ j)
            if ( idmap[tet_[i][j]] < 0 )
                idmap[tet_[i][j]] = pmesh->add_vertex_rt(vtx_[tet_[i][j]]);
        pmesh->add_tet(idmap[tet_[i][0]], idmap[tet_[i][1]],
                       idmap[tet_[i][2]], idmap[tet_[i][3]]);
    }
    printf("INFO: background tet mesh: %d tets\n", (int)pmesh->num_tets());
    pmesh->init();
}

template <class DistFunc>
void IsoStuffer<DistFunc>::create_mesh(TetMesh<double>* pmesh)
{
    create_boundary_octants();      // step 1
    octTree_.create_octree();       // step 2
    fill_highest_level();           // step 2.8
    octTree_.weak_balance();        // step 3
    octTree_.update_bitmasks();     // step 3.8
    create_background_tets();       // step 4
    compute_cutting_pts();          // step 5
    wrap_violated_pts();            // step 6
    extract_final_tets();           // step 7
    extract_final_mesh(pmesh);
}

}

#endif
