/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#if GMX_NBNXN_SIMD_BITWIDTH == 128
#define GMX_MM128_HERE
#else
#if GMX_NBNXN_SIMD_BITWIDTH == 256
#define GMX_MM256_HERE
#else
#error "unsupported GMX_NBNXN_SIMD_BITWIDTH"
#endif
#endif
#include "gmx_simd_macros.h"

#if GMX_SIMD_WIDTH_HERE >= 2*NBNXN_CPU_CLUSTER_I_SIZE
#define STRIDE_S  (GMX_SIMD_WIDTH_HERE/2)
#else
#define STRIDE_S  NBNXN_CPU_CLUSTER_I_SIZE
#endif

static gmx_inline gmx_mm_pr gmx_load_hpr_hilo_pr(const real *a)
{
    gmx_mm_hpr a_SSE;

    a_SSE = _mm_load_ps(a);

    return gmx_2hpr_to_pr(a_SSE, a_SSE);
}

static gmx_inline gmx_mm_pr gmx_set_2real_shift_pr(const real *a, real shift)
{
    gmx_mm_hpr a0, a1;

    a0 = _mm_set1_ps(a[0] + shift);
    a1 = _mm_set1_ps(a[1] + shift);

    return gmx_2hpr_to_pr(a1, a0);
}

/* Copies PBC shifted i-cell packed atom coordinates to working array */
static gmx_inline void
icell_set_x_simd_2xnn(int ci,
                      real shx, real shy, real shz,
                      int na_c,
                      int stride, const real *x,
                      nbnxn_list_work_t *work)
{
    int                     ia;
    nbnxn_x_ci_simd_2xnn_t *x_ci;

    x_ci = work->x_ci_simd_2xnn;

    ia = X_IND_CI_SIMD_2XNN(ci);

    x_ci->ix_SSE0 = gmx_set_2real_shift_pr(x + ia + 0*STRIDE_S + 0, shx);
    x_ci->iy_SSE0 = gmx_set_2real_shift_pr(x + ia + 1*STRIDE_S + 0, shy);
    x_ci->iz_SSE0 = gmx_set_2real_shift_pr(x + ia + 2*STRIDE_S + 0, shz);
    x_ci->ix_SSE2 = gmx_set_2real_shift_pr(x + ia + 0*STRIDE_S + 2, shx);
    x_ci->iy_SSE2 = gmx_set_2real_shift_pr(x + ia + 1*STRIDE_S + 2, shy);
    x_ci->iz_SSE2 = gmx_set_2real_shift_pr(x + ia + 2*STRIDE_S + 2, shz);
}

/* SIMD code for making a pair list of cell ci vs cell cjf-cjl
 * for coordinates in packed format.
 * Checks bouding box distances and possibly atom pair distances.
 * This is an accelerated version of make_cluster_list_simple.
 */
static gmx_inline void
make_cluster_list_simd_2xnn(const nbnxn_grid_t *gridj,
                            nbnxn_pairlist_t *nbl,
                            int ci, int cjf, int cjl,
                            gmx_bool remove_sub_diag,
                            const real *x_j,
                            real rl2, float rbb2,
                            int *ndistc)
{
    const nbnxn_x_ci_simd_2xnn_t *work;
    const float                  *bb_ci;

    gmx_mm_pr                     jx_SSE, jy_SSE, jz_SSE;

    gmx_mm_pr                     dx_SSE0, dy_SSE0, dz_SSE0;
    gmx_mm_pr                     dx_SSE2, dy_SSE2, dz_SSE2;

    gmx_mm_pr                     rsq_SSE0;
    gmx_mm_pr                     rsq_SSE2;

    gmx_mm_pr                     wco_SSE0;
    gmx_mm_pr                     wco_SSE2;
    gmx_mm_pr                     wco_any_SSE;

    gmx_mm_pr                     rc2_SSE;

    gmx_bool                      InRange;
    float                         d2;
    int                           xind_f, xind_l, cj;

    cjf = CI_TO_CJ_SIMD_2XNN(cjf);
    cjl = CI_TO_CJ_SIMD_2XNN(cjl+1) - 1;

    work = nbl->work->x_ci_simd_2xnn;

    bb_ci = nbl->work->bb_ci;

    rc2_SSE   = gmx_set1_pr(rl2);

    InRange = FALSE;
    while (!InRange && cjf <= cjl)
    {
        d2       = subc_bb_dist2_sse(4, 0, bb_ci, cjf, gridj->bbj);
        *ndistc += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */
        if (d2 < rbb2)
        {
            InRange = TRUE;
        }
        else if (d2 < rl2)
        {
            xind_f  = X_IND_CJ_SIMD_2XNN(CI_TO_CJ_SIMD_2XNN(gridj->cell0) + cjf);

            jx_SSE  = gmx_load_hpr_hilo_pr(x_j+xind_f+0*STRIDE_S);
            jy_SSE  = gmx_load_hpr_hilo_pr(x_j+xind_f+1*STRIDE_S);
            jz_SSE  = gmx_load_hpr_hilo_pr(x_j+xind_f+2*STRIDE_S);

            /* Calculate distance */
            dx_SSE0            = gmx_sub_pr(work->ix_SSE0, jx_SSE);
            dy_SSE0            = gmx_sub_pr(work->iy_SSE0, jy_SSE);
            dz_SSE0            = gmx_sub_pr(work->iz_SSE0, jz_SSE);
            dx_SSE2            = gmx_sub_pr(work->ix_SSE2, jx_SSE);
            dy_SSE2            = gmx_sub_pr(work->iy_SSE2, jy_SSE);
            dz_SSE2            = gmx_sub_pr(work->iz_SSE2, jz_SSE);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_calc_rsq_pr(dx_SSE0, dy_SSE0, dz_SSE0);
            rsq_SSE2           = gmx_calc_rsq_pr(dx_SSE2, dy_SSE2, dz_SSE2);

            wco_SSE0           = gmx_cmplt_pr(rsq_SSE0, rc2_SSE);
            wco_SSE2           = gmx_cmplt_pr(rsq_SSE2, rc2_SSE);

            wco_any_SSE        = gmx_or_pr(wco_SSE0, wco_SSE2);

            InRange            = gmx_movemask_pr(wco_any_SSE);

            *ndistc += 2*GMX_SIMD_WIDTH_HERE;
        }
        if (!InRange)
        {
            cjf++;
        }
    }
    if (!InRange)
    {
        return;
    }

    InRange = FALSE;
    while (!InRange && cjl > cjf)
    {
        d2       = subc_bb_dist2_sse(4, 0, bb_ci, cjl, gridj->bbj);
        *ndistc += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */
        if (d2 < rbb2)
        {
            InRange = TRUE;
        }
        else if (d2 < rl2)
        {
            xind_l  = X_IND_CJ_SIMD_2XNN(CI_TO_CJ_SIMD_2XNN(gridj->cell0) + cjl);

            jx_SSE  = gmx_load_hpr_hilo_pr(x_j+xind_l+0*STRIDE_S);
            jy_SSE  = gmx_load_hpr_hilo_pr(x_j+xind_l+1*STRIDE_S);
            jz_SSE  = gmx_load_hpr_hilo_pr(x_j+xind_l+2*STRIDE_S);

            /* Calculate distance */
            dx_SSE0            = gmx_sub_pr(work->ix_SSE0, jx_SSE);
            dy_SSE0            = gmx_sub_pr(work->iy_SSE0, jy_SSE);
            dz_SSE0            = gmx_sub_pr(work->iz_SSE0, jz_SSE);
            dx_SSE2            = gmx_sub_pr(work->ix_SSE2, jx_SSE);
            dy_SSE2            = gmx_sub_pr(work->iy_SSE2, jy_SSE);
            dz_SSE2            = gmx_sub_pr(work->iz_SSE2, jz_SSE);

            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_calc_rsq_pr(dx_SSE0, dy_SSE0, dz_SSE0);
            rsq_SSE2           = gmx_calc_rsq_pr(dx_SSE2, dy_SSE2, dz_SSE2);

            wco_SSE0           = gmx_cmplt_pr(rsq_SSE0, rc2_SSE);
            wco_SSE2           = gmx_cmplt_pr(rsq_SSE2, rc2_SSE);

            wco_any_SSE        = gmx_or_pr(wco_SSE0, wco_SSE2);

            InRange            = gmx_movemask_pr(wco_any_SSE);

            *ndistc += 2*GMX_SIMD_WIDTH_HERE;
        }
        if (!InRange)
        {
            cjl--;
        }
    }

    if (cjf <= cjl)
    {
        for (cj = cjf; cj <= cjl; cj++)
        {
            /* Store cj and the interaction mask */
            nbl->cj[nbl->ncj].cj   = CI_TO_CJ_SIMD_2XNN(gridj->cell0) + cj;
            nbl->cj[nbl->ncj].excl = get_imask_x86_simd_2xnn(remove_sub_diag, ci, cj);
            nbl->ncj++;
        }
        /* Increase the closing index in i super-cell list */
        nbl->ci[nbl->nci].cj_ind_end = nbl->ncj;
    }
}

#undef STRIDE_S
#undef GMX_MM128_HERE
#undef GMX_MM256_HERE
