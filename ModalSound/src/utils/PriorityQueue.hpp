/*
 * PriorityQueue.hpp
 * author: Changxi Zheng (cxzheng@cs.cornell.edu)
 */
#ifndef PRIORITY_QUEUE
#   define PRIORITY_QUEUE

#include <functional>
#include <memory>
#include <assert.h>
#include "macros.h"

#ifdef USE_OPENMP_NA
#   include <omp.h>
#endif

/*!
 * Data should have a field Data.qIdx
 * Compare returns true if the first argument is less than the second one
 */
template<
        typename Data, 
        typename Compare = std::less<Data>,
        typename Alloc = std::allocator<Data*>
        >
class PriorityQueue
{
    public:
        // =========== Constructors ==========
        PriorityQueue():dataQue_(NULL), used_(0) 
        { 
#ifdef USE_OPENMP_NA
            omp_init_lock(&lock_);
#endif
        }
        PriorityQueue(size_t s):size_(s), used_(0)
        {
            dataQue_ = alloc_.allocate(s);
#ifdef USE_OPENMP_NA
            omp_init_lock(&lock_);
#endif
        }
        ~PriorityQueue()
        {
            if ( size_ > 0 ) alloc_.deallocate(dataQue_, size_);
#ifdef USE_OPENMP_NA
            omp_destroy_lock(&lock_);
#endif
        }

        void  resize(size_t size);
        void  push(Data* ptr);
        //! Move up or down 
        void  update_node(Data* val);
        //! Move down along the tree
        void  move_up_node(Data* ptr);
        //! Move up along the tree
        void  move_down_node(Data* ptr);

        Data* pop();
        Data* peek()
        { return used_ ? dataQue_[0] : NULL; }

        bool  empty() const
        {
            return used_ == 0; 
        }
        void clear() 
        {
            used_ = 0;
        }
        size_t size() const { return used_; }

    private:
        /* 
         * These private methods are not thread-safe.
         *
         * The returned boolean indicates whether or not the queue position of 
         * the data has been changed due to the update/move operation.
         */
        bool private_move_up_node(Data* ptr);
        bool private_move_down_node(Data* ptr);

        Data**          dataQue_;
        size_t          size_;
        size_t          used_;

#ifdef USE_OPENMP_NA
        omp_lock_t      lock_;
#endif
        Compare         comp_;
        Alloc           alloc_;
};

// ==================== implementation ===================
template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::resize(size_t size)
{
#ifdef USE_OPENMP_NA
    omp_set_lock(&lock_);
#endif
    if ( size_ > 0 ) alloc_.deallocate(dataQue_, size_);
    size_ = size;
    used_ = 0;
    dataQue_ = alloc_.allocate(size_);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&lock_);
#endif
}

template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::push(Data* ptr)
{
#ifdef USE_OPENMP_NA
    omp_set_lock(&lock_);
#endif
    assert(used_ < size_ && ptr != NULL);
    dataQue_[used_] = ptr;
    ptr->qIdx = used_ ++;
    private_move_up_node(ptr);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&lock_);
#endif
}

template<typename Data, typename Compare, typename Alloc>
bool PriorityQueue<Data, Compare, Alloc>::private_move_up_node(Data* ptr)
{
    bool ret = false;
    int parent = (ptr->qIdx - 1) / 2;
    while ( parent >= 0 && comp_(*ptr, *dataQue_[parent]) )
    {
        dataQue_[parent]->qIdx = ptr->qIdx;
        std::swap(dataQue_[parent], dataQue_[ptr->qIdx]);
        ptr->qIdx = parent;
        parent = (ptr->qIdx - 1) / 2;
        ret = true;
    }
    return ret;
}

template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::update_node(Data* ptr)
{
#ifdef USE_OPENMP_NA
    omp_set_lock(&lock_);
#endif
    if ( !private_move_down_node(ptr) )
        private_move_up_node(ptr);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&lock_);
#endif
}

template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::move_up_node(Data* ptr)
{
#ifdef USE_OPENMP_NA
    omp_set_lock(&lock_);
#endif
    private_move_up_node(ptr);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&lock_);
#endif
}

template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::move_down_node(Data* ptr)
{
#ifdef USE_OPENMP_NA
    omp_set_lock(&lock_);
#endif
    private_move_down_node(ptr);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&lock_);
#endif
}

template<typename Data, typename Compare, typename Alloc>
bool PriorityQueue<Data, Compare, Alloc>::private_move_down_node(Data* ptr)
{
    bool ret = false;
    int child = ptr->qIdx * 2 + 1;
    child += (child+1 < used_ && 
            comp_(*dataQue_[child+1], *dataQue_[child]));
    while ( child < used_ && comp_(*dataQue_[child], *ptr) )
    {
        dataQue_[child]->qIdx = ptr->qIdx;
        std::swap(dataQue_[child], dataQue_[ptr->qIdx]);
        ptr->qIdx = child;
        child = child * 2 + 1;
        child += (child+1 < used_ && comp_(*dataQue_[child+1], *dataQue_[child]));
        ret = true;
    }
    return ret;
}

template<typename Data, typename Compare, typename Alloc>
Data* PriorityQueue<Data, Compare, Alloc>::pop()
{
    if ( !used_ ) return NULL;

#ifdef USE_OPENMP_NA
    omp_set_lock(&lock_);
#endif
    Data* ret = dataQue_[0];
    dataQue_[0] = dataQue_[-- used_];
    dataQue_[0]->qIdx = 0;

    Data* cur = dataQue_[0];
    private_move_down_node(cur);
#ifdef USE_OPENMP_NA
    omp_unset_lock(&lock_);
#endif
    return ret;
}

#endif
