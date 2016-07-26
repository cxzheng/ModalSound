#include "linear.h"
#include "linearalgebra/mat4inv.h"

void LinearMaterial::rest_stiffness_matrix(
        const Tet<REAL>& tet, Matrix<REAL>& out) const
{
    REAL vol = tet.volume();

    REAL mb[16];    // mb^T
    REAL beta[16];

    for(int i = 0;i < 4;++ i)
    {
        mb[i*4]   = tet[i].x;
        mb[i*4+1] = tet[i].y;
        mb[i*4+2] = tet[i].z;
        mb[i*4+3] = 1;
    }
    mat4d_invert(mb, beta);     // beta^T

    for(int i = 0;i < 4;++ i)   // each tet
    for(int j = 0;j < 4;++ j)
    {
        if ( i > j ) continue;

        for(int a = 0;a < 3;++ a)
        for(int b = 0;b < 3;++ b)
        {
            const int nrow = i*3 + a;
            const int ncol = j*3 + b;
            if ( nrow > ncol ) continue;

            const int ia = a*4 + i;
            const int jb = b*4 + j;
            const int ib = b*4 + i;
            const int ja = a*4 + j;
            out[nrow][ncol] = (lambda_*beta[ia]*beta[jb] + mu_*beta[ib]*beta[ja]);
            if ( a == b )
            {
                REAL tsum = 0;
                for(int k = 0;k < 3;++ k)
                {
                    const int ik = k*4 + i;
                    const int jk = k*4 + j;
                    tsum += beta[ik]*beta[jk];
                }
                out[nrow][ncol] += mu_*tsum;
            }

            out[nrow][ncol] *= (-vol * 0.5);
        }
    }
}
