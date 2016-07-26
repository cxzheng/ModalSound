#ifndef KD_TREE_HPP
#   define KD_TREE_HPP

#include <stdio.h>
#include <assert.h>
#include <vector>
#include <limits>

#ifdef USE_NAMESPACE
namespace sploosh
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
        const _Val*                                 data_;  
        size_t                                      N_;        // number of points
        TreeNode<_Dim, typename _Acc::result_type>* nodes_;    // node on the binary tree
        size_t                                      numNodes_; // number of nodes allocated
        std::vector<int>                            ptIdx_;
        bool                                        copied_;

        // heap used for finding n-nearest
        std::vector<typename _Acc::result_type>     auxHeap_;  

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
            if ( numNodes_ == 0 || nodes_[0].size() < (int)n )
            {
                fprintf(stderr, "Too many points required\n");
                return -1;
            }

            int retId = 0;
            int dim = 0;
            while ( true )
            {
                typename _Acc::result_type v = _acc(pt, dim);
                if ( nodes_[retId].child[0] && 
                     v < nodes_[nodes_[retId].child[0]].bdBox[1][dim] &&
                     nodes_[nodes_[retId].child[0]].size() >= (int)n )
                    retId = nodes_[retId].child[0];
                else if ( nodes_[retId].child[1] &&
                          v >= nodes_[nodes_[retId].child[1]].bdBox[0][dim] &&
                          nodes_[nodes_[retId].child[1]].size() >= (int)n )
                    retId = nodes_[retId].child[1];
                else
                    break;
                dim = (dim + 1) % _Dim;
            }
            return retId;
        }

        int locate_parent(const _Val& pt) const
        {
            int retId = 0, dim = 0;
            while ( nodes_[retId].child[0] )
            {
                typename _Acc::result_type v = _acc(pt, dim);   // the dim-th dimension of data of pt
                if ( v < nodes_[nodes_[retId].child[0]].bdBox[1][dim] )
                    retId = nodes_[retId].child[0];
                else
                    retId = nodes_[retId].child[1];
                dim = (dim + 1) % _Dim;
            }
            return retId;
        }

        void sift_down(std::vector<int>& pid, const size_t n)
        {
            typename _Acc::result_type a = auxHeap_[0];
            int ia = pid[0];
            int j = 1, jold = 0;
            while (j < (int)n)
            {
                if ( j < (int)n-1 && auxHeap_[j] < auxHeap_[j+1] ) ++ j;
                if ( a >= auxHeap_[j] ) break;

                auxHeap_[jold] = auxHeap_[j];
                pid[jold] = pid[j];

                jold = j;
                j = 2*j + 1;
            }
            auxHeap_[jold] = a;
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
                         _acc(data_[pid[ir]], dim) < _acc(data_[pid[l]], dim) )
                        std::swap(pid[l], pid[ir]);
                    return;
                }

                int mid = (l + ir) >> 1;
                std::swap(pid[mid], pid[l+1]);

                if ( _acc(data_[pid[l]], dim) > _acc(data_[pid[ir]], dim) )
                    std::swap(pid[l], pid[ir]);
                if ( _acc(data_[pid[l+1]], dim) > _acc(data_[pid[ir]], dim) )
                    std::swap(pid[l+1], pid[ir]);
                if ( _acc(data_[pid[l]], dim) > _acc(data_[pid[l+1]], dim) )
                    std::swap(pid[l], pid[l+1]);

                int i = l+1, j = ir, ia = pid[l+1];
                typename _Acc::result_type a = _acc(data_[ia], dim);
                for(;;)
                {
                    do ++ i; while (_acc(data_[pid[i]], dim) < a);
                    do -- j; while (_acc(data_[pid[j]], dim) > a);
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
                if ( v < nodes_[nodeId].bdBox[0][i] ) ret += (nodes_[nodeId].bdBox[0][i] - v)*(nodes_[nodeId].bdBox[0][i] - v);
                if ( v > nodes_[nodeId].bdBox[1][i] ) ret += (nodes_[nodeId].bdBox[1][i] - v)*(nodes_[nodeId].bdBox[1][i] - v);
            }
            return ret;
        }

    public:
        KDTree():N_(0), numNodes_(0), copied_(false) {}
        ~KDTree()
        {
            if ( copied_ && N_ > 0 ) _val_alloc.deallocate(const_cast<_Val*>(data_), N_);
            if ( numNodes_ > 0 ) _node_alloc.deallocate(nodes_, numNodes_);
        }

        /*!
         * Construct KD-tree from an array of points
         */
        KDTree(const _Val* data, const size_t N, bool copy = false):
            N_(0), numNodes_(0)
        {
            reinitialize(data, N, copy);
        }

        const _Val* data() const
        {  return data_; }

        /*!
         * find the nearest point to the given point pt, return its index
         */
        int find_nearest(const _Val& pt) const
        {
            int nodeId = locate_parent(pt);  // which node is pt in

            int nrstId;
            typename _Acc::result_type nrstDist = std::numeric_limits<_Acc::result_type>::infinity();
            for(int i = nodes_[nodeId].loIdx;i < nodes_[nodeId].hiIdx;++ i)
            {
                typename _Acc::result_type d = _dist_sqr(pt, data_[ptIdx_[i]]);
                if ( d < nrstDist )
                {
                    nrstId   = ptIdx_[i];
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
                    if ( nodes_[k].child[0] )
                    {
                        task[++ ntask] = nodes_[k].child[0];
                        task[++ ntask] = nodes_[k].child[1];
                    }
                    else
                    {
                        for(int i = nodes_[k].loIdx;i < nodes_[k].hiIdx;++ i)
                        {
                            typename _Acc::result_type d = _dist_sqr(pt, data_[ptIdx_[i]]);
                            if ( d < nrstDist )
                            {
                                nrstId   = ptIdx_[i];
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

            auxHeap_.resize(n);
            ret.resize(n);
            std::fill(auxHeap_.begin(), auxHeap_.end(), std::numeric_limits<_Acc::result_type>::infinity());

            // examine the points and save the n closest
            for(int i = nodes_[nodeId].loIdx, j = 0;i < nodes_[nodeId].hiIdx;++ i, ++ j)
            {
                typename _Acc::result_type d = _dist_sqr(pt, data_[ptIdx_[i]]);
                if ( d < auxHeap_[0] )
                {
                    auxHeap_[0] = d;
                    ret[0] = ptIdx_[i];
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
                if ( dist_sqr_to_node(k, pt) < auxHeap_[0] )
                {
                    if ( nodes_[k].child[0] )
                    {
                        task[++ ntask] = nodes_[k].child[0];
                        task[++ ntask] = nodes_[k].child[1];
                    }
                    else
                    {
                        for(int i = nodes_[k].loIdx;i < nodes_[k].hiIdx;++ i)
                        {
                            typename _Acc::result_type d = _dist_sqr(pt, data_[ptIdx_[i]]);
                            if ( d < auxHeap_[0] )
                            {
                                auxHeap_[0] = d;
                                ret[0] = ptIdx_[i];
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
            while ( nodes_[nodeId].child[0] )
            {
                typename _Acc::result_type v = _acc(pt, dim);
                if ( v + r < nodes_[nodes_[nodeId].child[0]].bdBox[1][dim] )
                    nodeId = nodes_[nodeId].child[0];
                else if ( v - r > nodes_[nodes_[nodeId].child[1]].bdBox[0][dim] )
                    nodeId = nodes_[nodeId].child[1];
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
                if ( nodes_[k].child[0] )
                {
                    task[++ ntask] = nodes_[k].child[0];
                    task[++ ntask] = nodes_[k].child[1];
                }
                else
                {
                    for(int i = nodes_[k].loIdx;i < nodes_[k].hiIdx;++ i)
                        if ( _dist_sqr(pt, data_[ptIdx_[i]]) <= r2 )
                            out.push_back(ptIdx_[i]);
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

            copied_ = copy;
            auxHeap_.clear();
            //// keep a reference to the original data points
            if ( copy ) // copy the data locally
            {
                // free previously allocated memory
                if ( N_ > 0 ) _val_alloc.deallocate(const_cast<_Val*>(data_), N_);

                _Val* ptr = _val_alloc.allocate(N);
                for(size_t i = 0;i < N;++ i)
                    _val_alloc.construct(&ptr[i], data[i]);     // invoke the copy constructor

                data_ = ptr;
                N_    = N;
            }
            else
            {
                data_ = data;
                N_    = N;
            }

            //// initialize the index of points
            ptIdx_.resize(N);
            for(int i = 0;i < (int)N;++ i) ptIdx_[i] = i;

            if ( numNodes_ > 0 ) _node_alloc.deallocate(nodes_, numNodes_);
            //// calculate the number of boxes and allocate memory for them
            int m = 1;
            for(int ntmp = (int)N;ntmp;ntmp >>= 1) m <<= 1;
            numNodes_ = 2*N - (m >> 1);
            if ( m < (int)numNodes_ ) numNodes_ = m;
            nodes_ = _node_alloc.allocate(numNodes_);

            //// initialize the root node
            for(int i = 0;i < _Dim;++ i)
            {
                nodes_[0].bdBox[0][i] = -std::numeric_limits<_Acc::result_type>::infinity();
                nodes_[0].bdBox[1][i] =  std::numeric_limits<_Acc::result_type>::infinity();
            }

            nodes_[0].parent = nodes_[0].child[0] = nodes_[0].child[1] = 0; 
            nodes_[0].loIdx = 0; 
            nodes_[0].hiIdx = (int)N;
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

                int npt = nodes_[tmom].size();
                int kk = npt / 2;
                
                //// select the middle point on the current dimension
                select_kth(kk, &ptIdx_[nodes_[tmom].loIdx], npt, tdim);

                int kkidx = ptIdx_[nodes_[tmom].loIdx + kk];
                typename _Acc::result_type v = _acc(data_[kkidx], tdim);

                nodes_[++ jbox].set(nodes_[tmom].bdBox[0], nodes_[tmom].bdBox[1],
                        tmom, 0, 0, nodes_[tmom].loIdx, nodes_[tmom].loIdx+kk);
                nodes_[jbox].bdBox[1][tdim] = v;

                nodes_[++ jbox].set(nodes_[tmom].bdBox[0], nodes_[tmom].bdBox[1],
                        tmom, 0, 0, nodes_[tmom].loIdx+kk, nodes_[tmom].hiIdx);
                nodes_[jbox].bdBox[0][tdim] = v;

                nodes_[tmom].child[0] = jbox - 1;
                nodes_[tmom].child[1] = jbox;

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
