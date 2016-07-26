/******************************************************************************
 *  File: KDTree.hpp
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
#ifndef KD_TREE_HPP
#   define KD_TREE_HPP

#include <stdio.h>
#include <vector>
#include <assert.h>

#ifdef USE_NAMESPACE
namespace carbine
{
#endif


template <typename _Val>
struct _BracketAccessor
{
    typedef typename _Val::value_type   result_type;

    result_type operator()(const _Val& V, const size_t N) const
    {
        return V[N];
    }
};

/*!
 * KD-tree Node
 */
template <int _Dim, typename T>
struct TreeNode
{
    int parent, child[2];   // parent and two children
    T   bdBox[2][_Dim];
    int loIdx, hiIdx;       // index of data in the range

    int size() const 
    {
        return hiIdx - loIdx;
    }

    void set(const T* bdMin, const T* bdMax, int p, int c1, int c2, int lidx, int hidx)
    {
        for(int i = 0;i < _Dim;++ i)
        {
            bdBox[0][i] = bdMin[i];
            bdBox[1][i] = bdMax[i];
        }
        parent = p;
        child[0] = c1;
        child[1] = c2;
        loIdx = lidx;
        hiIdx = hidx;
    }
};

#define LARGE_VAL   1E+30

/*!
 * KD tree for date type T, and dimension DIM.
 */
template <int _Dim, typename _Val,
          typename _DistSqr,
          typename _Acc = _BracketAccessor<_Val>,
          typename _ValAlloc = std::allocator<_Val>,
          typename _NodeAlloc = std::allocator<TreeNode<_Dim, typename _Acc::result_type> > >
class KDTree
{
    private:
        const _Val*                         mptr_data;  
        size_t                              m_N;        // number of points
        TreeNode<_Dim, typename _Acc::result_type>*  m_nodes;    // node on the binary tree
        size_t                              m_numNodes; // number of nodes allocated
        std::vector<int>                    m_ptIdx;
        bool                                m_copied;

        // heap used for finding n-nearest
        std::vector<typename _Acc::result_type>     m_auxHeap;  

        _ValAlloc       _val_alloc;
        _NodeAlloc      _node_alloc;
        _Acc            _acc;
        _DistSqr        _dist_sqr;

    private:
        /*!
         * Find a tree node, where the given point could be placed in. We also 
         * want the tree node maintains at least n points. Also we want the tree
         * node is as low as possible on the kd-tree. 
         *
         * this method is used in the first stage in find_nearest method to narrow down
         * the search space.
         */
        int locate_parent(const _Val& pt, const size_t n) const
        {
            if ( m_numNodes == 0 || m_nodes[0].size() < (int)n )
            {
                fprintf(stderr, "Too many points required\n");
                return -1;
            }

            int retId = 0;
            int dim = 0;
            while ( true )
            {
                typename _Acc::result_type v = _acc(pt, dim);
                if ( m_nodes[retId].child[0] && 
                     v < m_nodes[m_nodes[retId].child[0]].bdBox[1][dim] &&
                     m_nodes[m_nodes[retId].child[0]].size() >= (int)n )
                    retId = m_nodes[retId].child[0];
                else if ( m_nodes[retId].child[1] &&
                          v >= m_nodes[m_nodes[retId].child[1]].bdBox[0][dim] &&
                          m_nodes[m_nodes[retId].child[1]].size() >= (int)n )
                    retId = m_nodes[retId].child[1];
                else
                    break;
                dim = (dim + 1) % _Dim;
            }
            return retId;
        }

        int locate_parent(const _Val& pt) const
        {
            int retId = 0, dim = 0;
            while ( m_nodes[retId].child[0] )
            {
                typename _Acc::result_type v = _acc(pt, dim);   // the the dim-th dimension of data of pt
                if ( v < m_nodes[m_nodes[retId].child[0]].bdBox[1][dim] )
                    retId = m_nodes[retId].child[0];
                else
                    retId = m_nodes[retId].child[1];
                dim = (dim + 1) % _Dim;
            }
            return retId;
        }

        void sift_down(std::vector<int>& pid, const size_t n)
        {
            typename _Acc::result_type a = m_auxHeap[0];
            int ia = pid[0];
            int j = 1, jold = 0;
            while (j < (int)n)
            {
                if ( j < (int)n-1 && m_auxHeap[j] < m_auxHeap[j+1] ) ++ j;
                if ( a >= m_auxHeap[j] ) break;

                m_auxHeap[jold] = m_auxHeap[j];
                pid[jold] = pid[j];

                jold = j;
                j = 2*j + 1;
            }
            m_auxHeap[jold] = a;
            pid[jold] = ia;
        }

        /*!
         * Select the Kth element
         */
        void select_kth(const size_t k, int * pid, const size_t npt, int dim)
        {
            if ( npt < k ) return;

            int l = 0, ir = npt-1;
            for(;;) 
            {
                if ( ir <= l+1 )
                {
                    if ( ir == l+1 && 
                         _acc(mptr_data[pid[ir]], dim) < _acc(mptr_data[pid[l]], dim) )
                        std::swap(pid[l], pid[ir]);
                    return;
                }

                int mid = (l + ir) >> 1;
                std::swap(pid[mid], pid[l+1]);

                if ( _acc(mptr_data[pid[l]], dim) > _acc(mptr_data[pid[ir]], dim) )
                    std::swap(pid[l], pid[ir]);
                if ( _acc(mptr_data[pid[l+1]], dim) > _acc(mptr_data[pid[ir]], dim) )
                    std::swap(pid[l+1], pid[ir]);
                if ( _acc(mptr_data[pid[l]], dim) > _acc(mptr_data[pid[l+1]], dim) )
                    std::swap(pid[l], pid[l+1]);

                int i = l+1, j = ir, ia = pid[l+1];
                typename _Acc::result_type a = _acc(mptr_data[ia], dim);
                for(;;)
                {
                    do ++ i; while (_acc(mptr_data[pid[i]], dim) < a);
                    do -- j; while (_acc(mptr_data[pid[j]], dim) > a);
                    if ( j <= i ) break;
                    std::swap(pid[i], pid[j]);
                }
                pid[l+1] = pid[j];
                pid[j] = ia;
                if ( j >= (int)k ) ir = j - 1;
                if ( j <= (int)k ) l = i;
            }
        }

        typename _Acc::result_type dist_sqr_to_node(int nodeId, const _Val& pt) const
        {
            typename _Acc::result_type ret = 0;
            for(int i = 0;i < _Dim;++ i)
            {
                typename _Acc::result_type v = _acc(pt, i);
                if ( v < m_nodes[nodeId].bdBox[0][i] ) ret += (m_nodes[nodeId].bdBox[0][i] - v)*(m_nodes[nodeId].bdBox[0][i] - v);
                if ( v > m_nodes[nodeId].bdBox[1][i] ) ret += (m_nodes[nodeId].bdBox[1][i] - v)*(m_nodes[nodeId].bdBox[1][i] - v);
            }
            return ret;
        }

    public:
        KDTree():m_N(0), m_numNodes(0), m_copied(false) {}
        ~KDTree()
        {
            if ( m_copied && m_N > 0 ) _val_alloc.deallocate(const_cast<_Val*>(mptr_data), m_N);
            if ( m_numNodes > 0 ) _node_alloc.deallocate(m_nodes, m_numNodes);
        }

        /*!
         * Construct KD-tree from an array of points
         */
        KDTree(const _Val* data, const size_t N, bool copy = false):
            m_N(0), m_numNodes(0)
        {
            reinitialize(data, N, copy);
        }

        const _Val* data() const
        {  return mptr_data; }

        /*!
         * find the nearest point to the given point pt, return its index
         */
        int find_nearest(const _Val& pt) const
        {
            int nodeId = locate_parent(pt);  // which node is pt in

            int nrstId;
            typename _Acc::result_type nrstDist = LARGE_VAL;
            for(int i = m_nodes[nodeId].loIdx;i < m_nodes[nodeId].hiIdx;++ i)
            {
                typename _Acc::result_type d = _dist_sqr(pt, mptr_data[m_ptIdx[i]]);
                if ( d < nrstDist )
                {
                    nrstId   = m_ptIdx[i];
                    nrstDist = d;
                }
            }   // end for

            int task[128];
            int ntask = 1;
            task[1] = 0;
            while ( ntask )
            {
                int k = task[ntask --];
                if ( k == nodeId ) continue;
                if ( dist_sqr_to_node(k, pt) < nrstDist )
                {
                    if ( m_nodes[k].child[0] )
                    {
                        task[++ ntask] = m_nodes[k].child[0];
                        task[++ ntask] = m_nodes[k].child[1];
                    }
                    else
                    {
                        for(int i = m_nodes[k].loIdx;i < m_nodes[k].hiIdx;++ i)
                        {
                            typename _Acc::result_type d = _dist_sqr(pt, mptr_data[m_ptIdx[i]]);
                            if ( d < nrstDist )
                            {
                                nrstId   = m_ptIdx[i];
                                nrstDist = d;
                            }
                        }
                    } // end else 
                }
            } // end while
            return nrstId;
        }

        /*!
         * find n nearest points to the given point pt
         */
        void find_nearest(const _Val& pt, const size_t n, std::vector<int>& ret) const
        {
            ret.clear();
            int nodeId = locate_parent(pt, n);
            if ( nodeId < 0 ) return;

            m_auxHeap.resize(n);
            ret.resize(n);
            std::fill(m_auxHeap.begin(), m_auxHeap.end(), LARGE_VAL);

            // examine the points and save the n closest
            for(int i = m_nodes[nodeId].loIdx, j = 0;i < m_nodes[nodeId].hiIdx;++ i, ++ j)
            {
                typename _Acc::result_type d = _dist_sqr(pt, mptr_data[m_ptIdx[i]]);
                if ( d < m_auxHeap[0] )
                {
                    m_auxHeap[0] = d;
                    ret[0] = m_ptIdx[i];
                    // maintain the heap structure
                    sift_down(ret, n);
                }
            }

            // traverse the tree
            int task[128];
            int ntask = 1;
            task[1] = 0;
            while ( ntask )
            {
                int k = task[ntask --];
                if ( k == nodeId ) continue;
                if ( dist_sqr_to_node(k, pt) < m_auxHeap[0] )
                {
                    if ( m_nodes[k].child[0] )
                    {
                        task[++ ntask] = m_nodes[k].child[0];
                        task[++ ntask] = m_nodes[k].child[1];
                    }
                    else
                    {
                        for(int i = m_nodes[k].loIdx;i < m_nodes[k].hiIdx;++ i)
                        {
                            typename _Acc::result_type d = _dist_sqr(pt, mptr_data[m_ptIdx[i]]);
                            if ( d < m_auxHeap[0] )
                            {
                                m_auxHeap[0] = d;
                                ret[0] = m_ptIdx[i];
                                // maintain the heap structure
                                sift_down(ret, n);
                            }
                        }
                    }
                }
            } // end while 
        }

        /*!
         * Find all the data point whose distance computed by _dist is <= r. The indices of those
         * points are put into vector out
         */
        int find_in_ball(const _Val& pt, typename _Acc::result_type r, std::vector<int>& out) const
        {
            out.clear();
            if ( r < (typename _Acc::result_type)0 ) return 0;

            int nodeId = 0, dim = 0;
            while ( m_nodes[nodeId].child[0] )
            {
                typename _Acc::result_type v = _acc(pt, dim);
                if ( v + r < m_nodes[m_nodes[nodeId].child[0]].bdBox[1][dim] )
                    nodeId = m_nodes[nodeId].child[0];
                else if ( v - r > m_nodes[m_nodes[nodeId].child[1]].bdBox[0][dim] )
                    nodeId = m_nodes[nodeId].child[1];
                else
                    break;
                dim = (dim + 1) % _Dim;
            }

            const typename _Acc::result_type r2 = r*r;
            int task[128];
            int ntask = 1;
            task[1] = nodeId;
            while ( ntask ) 
            {
                int k = task[ntask --];
                if ( dist_sqr_to_node(k, pt) > r2 ) continue;
                if ( m_nodes[k].child[0] )
                {
                    task[++ ntask] = m_nodes[k].child[0];
                    task[++ ntask] = m_nodes[k].child[1];
                }
                else
                {
                    for(int i = m_nodes[k].loIdx;i < m_nodes[k].hiIdx;++ i)
                        if ( _dist_sqr(pt, mptr_data[m_ptIdx[i]]) <= r2 )
                            out.push_back(m_ptIdx[i]);
                } // end if
            } // end while
            return (int)out.size();
        }

        /*!
         * NOTE: if copy == true, make sure there is a valid copy operator for
         *       type _Val
         */
        void reinitialize(const _Val* data, const size_t N, bool copy = false)
        {
            assert(N > 0);

            m_copied = copy;
            m_auxHeap.clear();
            //// keep a reference to the original data points
            if ( copy ) // copy the data locally
            {
                // free previously allocated memory
                if ( m_N > 0 ) _val_alloc.deallocate(const_cast<_Val*>(mptr_data), m_N);

                _Val* ptr = _val_alloc.allocate(N);
                for(size_t i = 0;i < N;++ i)
                    _val_alloc.construct(&ptr[i], data[i]);     // invoke the copy constructor

                mptr_data = ptr;
                m_N = N;
            }
            else
            {
                mptr_data = data;
                m_N = N;
            }

            //// initialize the index of points
            m_ptIdx.resize(N);
            for(int i = 0;i < (int)N;++ i) m_ptIdx[i] = i;

            if ( m_numNodes > 0 ) _node_alloc.deallocate(m_nodes, m_numNodes);
            //// calculate the number of boxes and allocate memory for them
            int m = 1;
            for(int ntmp = (int)N;ntmp;ntmp >>= 1) m <<= 1;
            m_numNodes = 2*N - (m >> 1);
            if ( m < (int)m_numNodes ) m_numNodes = m;
            m_nodes = _node_alloc.allocate(m_numNodes);

            //// initialize the root node
            for(int i = 0;i < _Dim;++ i)
            {
                m_nodes[0].bdBox[0][i] = -LARGE_VAL;
                m_nodes[0].bdBox[1][i] =  LARGE_VAL;
            }

            m_nodes[0].parent = m_nodes[0].child[0] = m_nodes[0].child[1] = 0; 
            m_nodes[0].loIdx = 0; 
            m_nodes[0].hiIdx = (int)N;
            if ( N <= 4 ) return;

            int jbox = 0;
            int taskmom[50], taskdim[50];
            int nowtask = 1;
            taskmom[1] = 0;     // which node
            taskdim[1] = 0;     // dimension for that node
            while ( nowtask ) 
            {
                int tmom = taskmom[nowtask];
                int tdim = taskdim[nowtask --];

                int npt = m_nodes[tmom].size();
                int kk = npt / 2;
                
                //// select the middle point on the current dimension
                select_kth(kk, &m_ptIdx[m_nodes[tmom].loIdx], npt, tdim);

                int kkidx = m_ptIdx[m_nodes[tmom].loIdx + kk];
                typename _Acc::result_type v = _acc(mptr_data[kkidx], tdim);

                m_nodes[++ jbox].set(m_nodes[tmom].bdBox[0], m_nodes[tmom].bdBox[1],
                        tmom, 0, 0, m_nodes[tmom].loIdx, m_nodes[tmom].loIdx+kk);
                m_nodes[jbox].bdBox[1][tdim] = v;

                m_nodes[++ jbox].set(m_nodes[tmom].bdBox[0], m_nodes[tmom].bdBox[1],
                        tmom, 0, 0, m_nodes[tmom].loIdx+kk, m_nodes[tmom].hiIdx);
                m_nodes[jbox].bdBox[0][tdim] = v;

                m_nodes[tmom].child[0] = jbox - 1;
                m_nodes[tmom].child[1] = jbox;

                if ( kk > 4 ) 
                {
                    taskmom[++ nowtask] = jbox - 1;
                    taskdim[nowtask] = (tdim + 1) % _Dim;
                }
                if ( npt - kk > 4 )
                {
                    taskmom[++ nowtask] = jbox;
                    taskdim[nowtask] = (tdim + 1) % _Dim;
                }
            }
        }
};

#ifdef USE_NAMESPACE
}
#endif
#endif
