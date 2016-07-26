#ifndef DEFORMABLE_FEM_HPP
#   define DEFORMABLE_FEM_HPP

#include <vector>
#include "sc/PardisoMatrix.hpp"
#include "geometry/FixedVtxTetMesh.hpp"

namespace DeformableFEM
{

/*!
 * Compute the mass matrix for a tetrahedron based FEM model
 */
template <typename T>
static void mass_mat(const FixedVtxTetMesh<T>* pmesh, PardisoMatrix<T>& M)
{
    using namespace std;

    const int numFixedVtx = pmesh->num_fixed_vertices();
    const int numFreeVtx  = pmesh->num_vertices() - numFixedVtx;
    M.resize(numFreeVtx*3, true, true);  // symmetric pos-def. matrix

    const vector< Tet<REAL> >& tets = pmesh->tets();
    const vector<TetMesh<REAL>::TetIdx>& idx = pmesh->tet_indices();

    // construct the sparsity structure for mass matrix
    for(size_t tetid = 0;tetid < tets.size();++ tetid)  // each tet
    for(size_t nid = 0;nid < 4;++ nid)                  // each node
    {
        if ( pmesh->is_fixed_vertex(idx[tetid][nid]) ) continue;

        int rowId = (idx[tetid][nid] - numFixedVtx) * 3;
        for(size_t n2id = 0;n2id < 4;++ n2id)           // each node in the same tet
        {
            if ( pmesh->is_fixed_vertex(idx[tetid][n2id]) ) continue;
            int colId = (idx[tetid][n2id] - numFixedVtx) * 3;

            for(size_t dir = 0;dir < 3;++ dir)              // x,y,z component
                M.set_nonzero(rowId+dir, colId+dir);
        }
    }
    M.generate_pattern();

    // construct the mass matrix
    M.zeros();
    for(size_t tetid = 0;tetid < tets.size();++ tetid)  // each tet
    {
        const REAL vol = tets[tetid].volume();
        for(size_t nid = 0;nid < 4;++ nid)              // each node (i)
        {
            if ( pmesh->is_fixed_vertex(idx[tetid][nid]) ) continue;

            int rowId = (idx[tetid][nid] - numFixedVtx) * 3;
            for(size_t n2id = 0;n2id < 4;++ n2id)       // each node in the same tet (j)
            {
                if ( pmesh->is_fixed_vertex(idx[tetid][n2id]) || 
                     idx[tetid][nid] > idx[tetid][n2id] ) continue;

                int colId = (idx[tetid][n2id] - numFixedVtx) * 3;
                for(size_t dir = 0;dir < 3;++ dir)      // (a)
                    M.add(rowId+dir, colId+dir, vol*0.05*(1 + (nid == n2id)));
            }
        }
    }  // end for tetid
}

/*!
 * Compute the stiffness matrix for a tetrahedron based FEM model 
 * with the given material
 */
template <typename T, class TMaterial>
static void stiffness_mat(const FixedVtxTetMesh<T>* pmesh, 
        const TMaterial* pmaterial, PardisoMatrix<T>& K)
{
    using namespace std;

    const int numFixedVtx = pmesh->num_fixed_vertices();
    const int numFreeVtx  = pmesh->num_vertices() - numFixedVtx;
    K.resize(numFreeVtx*3, true);  // symmetric matrix
    
    const vector< Tet<REAL> >& tets = pmesh->tets();
    const vector<TetMesh<REAL>::TetIdx>& idx = pmesh->tet_indices();

    // construct the sparsity structure for stiffness matrix
    for(size_t tetid = 0;tetid < tets.size();++ tetid)  // each tet
    for(size_t nid = 0;nid < 4;++ nid)                  // each node
    {
        if ( pmesh->is_fixed_vertex(idx[tetid][nid]) ) continue;

        int rowId = (idx[tetid][nid] - numFixedVtx) * 3;
        for(size_t dir = 0;dir < 3;++ dir)              // x,y,z component
        for(size_t n2id = 0;n2id < 4;++ n2id)           // each node in the same tet
        {
            if ( pmesh->is_fixed_vertex(idx[tetid][n2id]) ) continue;
            int colId = (idx[tetid][n2id] - numFixedVtx) * 3;

            for(size_t dir2 = 0;dir2 < 3;++ dir2)
                K.set_nonzero(rowId+dir, colId+dir2);
        }
    } // end for
    K.generate_pattern();

    // construct the stiffness matrix
    K.zeros();
    Matrix<REAL> stiff(12, true);
    for(size_t tetid = 0;tetid < tets.size();++ tetid)  // each tet
    {
        pmaterial->stiffness_matrix(tets[tetid], stiff);
        for(size_t nid = 0;nid < 4;++ nid)              // each node
        {
            if ( pmesh->is_fixed_vertex(idx[tetid][nid]) ) continue;

            int rowId = (idx[tetid][nid] - numFixedVtx) * 3;
            for(size_t n2id = 0;n2id < 4;++ n2id)       // each node in the same tet
            {
                if ( pmesh->is_fixed_vertex(idx[tetid][n2id]) || 
                     idx[tetid][nid] > idx[tetid][n2id] ) continue;

                int colId = (idx[tetid][n2id] - numFixedVtx) * 3;
                for(size_t dir = 0;dir < 3;++ dir)
                for(size_t dir2 = 0;dir2 < 3;++ dir2)
                {
                    int srow = nid*3 + dir;
                    int scol = n2id*3 + dir2;

                    /* 
                     * -stiff[srow][scol] because the equation is actually
                     * Mu'' = F + Ku
                     */
                    K.add(rowId+dir, colId+dir2, 
                            srow <= scol ? -stiff[srow][scol] : -stiff[scol][srow]);
                }
            }
        }
    }  // end for tetid
}

}

#endif
