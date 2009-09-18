/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>

#include "typedefs.h"
#include "smalloc.h"
#include "genborn.h"
#include "vec.h"
#include "grompp.h"
#include "pdbio.h"
#include "names.h"
#include "physics.h"
#include "partdec.h"
#include "domdec.h"
#include "network.h"
#include "gmx_fatal.h"
#include "mtop_util.h"
#include "genborn.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif


/* Only compile this file if SSE intrinsics are available */
#if ( (defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_SSE2)) && !defined(GMX_DOUBLE) )

#include <gmx_sse2_single.h>
#include <xmmintrin.h>
#include <emmintrin.h>


int 
calc_gb_rad_still_sse(t_commrec *cr, t_forcerec *fr,int natoms, gmx_localtop_t *top,
					  const t_atomtypes *atype, float *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md)
{
	int i,k,n,ii,is3,ii3,nj0,nj1,offset;
    int n0,n1;
	int jnrA,jnrB,jnrC,jnrD,j3A,j3B,j3C,j3D;
	int jnrE,jnrF,jnrG,jnrH,j3E,j3F,j3G,j3H;
	int shift;
    int *mdtype;
	real shX,shY,shZ;
    int *jjnr;
    real *shiftvec;

	float gpi_ai,gpi2;
	float factor;
	float *gb_radius;
    float *vsolv;
    float *work;
    float *dadx;
    
	__m128 ix,iy,iz;
	__m128 jx,jy,jz;
	__m128 dx,dy,dz;
	__m128 tx,ty,tz;
	__m128 jxB,jyB,jzB;
	__m128 dxB,dyB,dzB;
	__m128 txB,tyB,tzB;
	__m128 rsq,rinv,rinv2,rinv4,rinv6;
	__m128 rsqB,rinvB,rinv2B,rinv4B,rinv6B;
	__m128 ratio,gpi,rai,raj,vai,vaj,rvdw;
	__m128 ratioB,rajB,vajB,rvdwB;
	__m128 ccf,dccf,theta,cosq,term,sinq,res,prod,prod_ai,tmp;
	__m128 ccfB,dccfB,thetaB,cosqB,termB,sinqB,resB,prodB;
	__m128 mask,icf4,icf6,mask_cmp;
	__m128 icf4B,icf6B,mask_cmpB;
	
    __m128   mask1 = _mm_castsi128_ps( _mm_set_epi32(0, 0, 0, 0xffffffff) );
	__m128   mask2 = _mm_castsi128_ps( _mm_set_epi32(0, 0, 0xffffffff, 0xffffffff) );
	__m128   mask3 = _mm_castsi128_ps( _mm_set_epi32(0, 0xffffffff, 0xffffffff, 0xffffffff) );
    
	const __m128 half   = {0.5f , 0.5f , 0.5f , 0.5f };
	const __m128 three  = {3.0f , 3.0f , 3.0f , 3.0f };
	const __m128 one    = {1.0f,  1.0f , 1.0f , 1.0f };
	const __m128 two    = {2.0f , 2.0f , 2.0f,  2.0f };
	const __m128 zero   = {0.0f , 0.0f , 0.0f , 0.0f };
	const __m128 four   = {4.0f , 4.0f , 4.0f , 4.0f };
	
	const __m128 still_p5inv  = {STILL_P5INV, STILL_P5INV, STILL_P5INV, STILL_P5INV};
	const __m128 still_pip5   = {STILL_PIP5,  STILL_PIP5,  STILL_PIP5,  STILL_PIP5};
	const __m128 still_p4     = {STILL_P4,    STILL_P4,    STILL_P4,    STILL_P4};
		
	factor  = 0.5 * ONE_4PI_EPS0;
		
    gb_radius = born->gb_radius;
    vsolv     = born->vsolv;
    work      = born->gpol_still_work;
	jjnr      = nl->jjnr;
    shiftvec  = fr->shift_vec[0];
    dadx      = fr->dadx;
    
	jnrA = jnrB = jnrC = jnrD = 0;
    jx = _mm_setzero_ps();
    jy = _mm_setzero_ps();
    jz = _mm_setzero_ps();
    
	n = 0;
    
    n0 = md->start;
    n1 = md->start+md->homenr+natoms/2+1;
		
	for(i=0;i<natoms;i++)
	{
		work[i]=0;
	}

	for(i=0;i<nl->nri;i++)
	{
        ii     = nl->iinr[i];
		ii3	   = ii*3;
        is3    = 3*nl->shift[i];     
        shX    = shiftvec[is3];  
        shY    = shiftvec[is3+1];
        shZ    = shiftvec[is3+2];
        nj0    = nl->jindex[i];      
        nj1    = nl->jindex[i+1];    
        
        ix     = _mm_set1_ps(shX+x[ii3+0]);
		iy     = _mm_set1_ps(shY+x[ii3+1]);
		iz     = _mm_set1_ps(shZ+x[ii3+2]);
		
		offset = (nj1-nj0)%4;
		
		/* Polarization energy for atom ai */
		gpi    = _mm_setzero_ps();
		
        rai     = _mm_load1_ps(gb_radius+ii);
        prod_ai = _mm_set1_ps(STILL_P4*vsolv[ii]);
        
		for(k=nj0;k<nj1-4-offset;k+=8)
		{
			jnrA        = jjnr[k];   
			jnrB        = jjnr[k+1];
			jnrC        = jjnr[k+2];
			jnrD        = jjnr[k+3];
			jnrE        = jjnr[k+4];   
			jnrF        = jjnr[k+5];
			jnrG        = jjnr[k+6];
			jnrH        = jjnr[k+7];
            
            j3A         = 3*jnrA;  
			j3B         = 3*jnrB;
			j3C         = 3*jnrC;
			j3D         = 3*jnrD;
            j3E         = 3*jnrE;  
			j3F         = 3*jnrF;
			j3G         = 3*jnrG;
			j3H         = 3*jnrH;
            
            GMX_MM_LOAD_4JCOORD_1ATOM_PS(x+j3A,x+j3B,x+j3C,x+j3D,jx,jy,jz);
            GMX_MM_LOAD_4JCOORD_1ATOM_PS(x+j3E,x+j3F,x+j3G,x+j3H,jxB,jyB,jzB);
            
            GMX_MM_LOAD_4VALUES_PS(gb_radius+jnrA,gb_radius+jnrB,gb_radius+jnrC,gb_radius+jnrD,raj);
            GMX_MM_LOAD_4VALUES_PS(gb_radius+jnrE,gb_radius+jnrF,gb_radius+jnrG,gb_radius+jnrH,rajB);
			GMX_MM_LOAD_4VALUES_PS(vsolv+jnrA,vsolv+jnrB,vsolv+jnrC,vsolv+jnrD,vaj);
			GMX_MM_LOAD_4VALUES_PS(vsolv+jnrE,vsolv+jnrF,vsolv+jnrG,vsolv+jnrH,vajB);

			dx          = _mm_sub_ps(ix,jx);
			dy          = _mm_sub_ps(iy,jy);
			dz          = _mm_sub_ps(iz,jz);
			dxB         = _mm_sub_ps(ix,jxB);
			dyB         = _mm_sub_ps(iy,jyB);
			dzB         = _mm_sub_ps(iz,jzB);
            
            rsq         = gmx_mm_calc_rsq(dx,dy,dz);
            rsqB        = gmx_mm_calc_rsq(dxB,dyB,dzB);
            rinv        = gmx_mm_invsqrt_ps(rsq);
            rinvB       = gmx_mm_invsqrt_ps(rsqB);
            rinv2       = _mm_mul_ps(rinv,rinv);
            rinv2B      = _mm_mul_ps(rinvB,rinvB);
            rinv4       = _mm_mul_ps(rinv2,rinv2);
            rinv4B      = _mm_mul_ps(rinv2B,rinv2B);
            rinv6       = _mm_mul_ps(rinv4,rinv2);
            rinv6B      = _mm_mul_ps(rinv4B,rinv2B);
            
            rvdw        = _mm_add_ps(rai,raj);
            rvdwB       = _mm_add_ps(rai,rajB);
            ratio       = _mm_mul_ps(rsq, gmx_mm_inv_ps( _mm_mul_ps(rvdw,rvdw)));
            ratioB      = _mm_mul_ps(rsqB, gmx_mm_inv_ps( _mm_mul_ps(rvdwB,rvdwB)));

            mask_cmp    = _mm_cmple_ps(ratio,still_p5inv);
            mask_cmpB   = _mm_cmple_ps(ratioB,still_p5inv);
            
            /* gmx_mm_sincos_ps() is quite expensive, so avoid calculating it if we can! */
            if( 0 == _mm_movemask_ps(mask_cmp) )
            {
                /* if ratio>still_p5inv for ALL elements */
                ccf         = one;
                dccf        = _mm_setzero_ps();
            }
            else 
            {
                ratio       = _mm_min_ps(ratio,still_p5inv);
                theta       = _mm_mul_ps(ratio,still_pip5);
                gmx_mm_sincos_ps(theta,&sinq,&cosq);            
                term        = _mm_mul_ps(half,_mm_sub_ps(one,cosq));
                ccf         = _mm_mul_ps(term,term);
                dccf        = _mm_mul_ps(_mm_mul_ps(two,term),
                                         _mm_mul_ps(sinq,theta));
            }
            if( 0 == _mm_movemask_ps(mask_cmpB) )
            {
                /* if ratio>still_p5inv for ALL elements */
                ccfB        = one;
                dccfB       = _mm_setzero_ps();
            }
            else 
            {
                ratioB      = _mm_min_ps(ratioB,still_p5inv);
                thetaB      = _mm_mul_ps(ratioB,still_pip5);
                gmx_mm_sincos_ps(thetaB,&sinqB,&cosqB);            
                termB       = _mm_mul_ps(half,_mm_sub_ps(one,cosqB));
                ccfB        = _mm_mul_ps(termB,termB);
                dccfB       = _mm_mul_ps(_mm_mul_ps(two,termB),
                                         _mm_mul_ps(sinqB,thetaB));
            }
            
            prod        = _mm_mul_ps(still_p4,vaj);
            prodB       = _mm_mul_ps(still_p4,vajB);
            icf4        = _mm_mul_ps(ccf,rinv4);
            icf4B       = _mm_mul_ps(ccfB,rinv4B);
            icf6        = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four,ccf),dccf), rinv6);
            icf6B       = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four,ccfB),dccfB), rinv6B);
            
            GMX_MM_INCREMENT_4VALUES_PS(work+jnrA,work+jnrB,work+jnrC,work+jnrD,_mm_mul_ps(prod_ai,icf4));
            GMX_MM_INCREMENT_4VALUES_PS(work+jnrE,work+jnrF,work+jnrG,work+jnrH,_mm_mul_ps(prod_ai,icf4B));
            
            gpi           = _mm_add_ps(gpi, _mm_add_ps( _mm_mul_ps(prod,icf4) , _mm_mul_ps(prodB,icf4B) ) );
            
            _mm_store_ps(dadx,_mm_mul_ps(prod,icf6));
            dadx+=4;
            _mm_store_ps(dadx,_mm_mul_ps(prod_ai,icf6));
            dadx+=4;
            _mm_store_ps(dadx,_mm_mul_ps(prodB,icf6B));
            dadx+=4;
            _mm_store_ps(dadx,_mm_mul_ps(prod_ai,icf6B));
            dadx+=4;
		} 
 
        for(;k<nj1-offset;k+=4)
		{
			jnrA        = jjnr[k];   
			jnrB        = jjnr[k+1];
			jnrC        = jjnr[k+2];
			jnrD        = jjnr[k+3];
            
            j3A         = 3*jnrA;  
			j3B         = 3*jnrB;
			j3C         = 3*jnrC;
			j3D         = 3*jnrD;
            
            GMX_MM_LOAD_4JCOORD_1ATOM_PS(x+j3A,x+j3B,x+j3C,x+j3D,jx,jy,jz);
            
            GMX_MM_LOAD_4VALUES_PS(gb_radius+jnrA,gb_radius+jnrB,gb_radius+jnrC,gb_radius+jnrD,raj);
			GMX_MM_LOAD_4VALUES_PS(vsolv+jnrA,vsolv+jnrB,vsolv+jnrC,vsolv+jnrD,vaj);
            
			dx          = _mm_sub_ps(ix,jx);
			dy          = _mm_sub_ps(iy,jy);
			dz          = _mm_sub_ps(iz,jz);
            
            rsq         = gmx_mm_calc_rsq(dx,dy,dz);
            rinv        = gmx_mm_invsqrt_ps(rsq);
            rinv2       = _mm_mul_ps(rinv,rinv);
            rinv4       = _mm_mul_ps(rinv2,rinv2);
            rinv6       = _mm_mul_ps(rinv4,rinv2);
            
            rvdw        = _mm_add_ps(rai,raj);
            ratio       = _mm_mul_ps(rsq, gmx_mm_inv_ps( _mm_mul_ps(rvdw,rvdw)));
            
            mask_cmp    = _mm_cmple_ps(ratio,still_p5inv);
            
            /* gmx_mm_sincos_ps() is quite expensive, so avoid calculating it if we can! */
            if(0 == _mm_movemask_ps(mask_cmp))
            {
                /* if ratio>still_p5inv for ALL elements */
                ccf         = one;
                dccf        = _mm_setzero_ps();
            }
            else 
            {
                ratio       = _mm_min_ps(ratio,still_p5inv);
                theta       = _mm_mul_ps(ratio,still_pip5);
                gmx_mm_sincos_ps(theta,&sinq,&cosq);            
                term        = _mm_mul_ps(half,_mm_sub_ps(one,cosq));
                ccf         = _mm_mul_ps(term,term);
                dccf        = _mm_mul_ps(_mm_mul_ps(two,term),
                                         _mm_mul_ps(sinq,theta));
            }
            
            prod        = _mm_mul_ps(still_p4,vaj);
            icf4        = _mm_mul_ps(ccf,rinv4);
            icf6        = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four,ccf),dccf), rinv6);
            
            
            GMX_MM_INCREMENT_4VALUES_PS(work+jnrA,work+jnrB,work+jnrC,work+jnrD,_mm_mul_ps(prod_ai,icf4));
            
            gpi           = _mm_add_ps(gpi, _mm_mul_ps(prod,icf4));
            
            _mm_store_ps(dadx,_mm_mul_ps(prod,icf6));
            dadx+=4;
            _mm_store_ps(dadx,_mm_mul_ps(prod_ai,icf6));
            dadx+=4;
		} 
        
        if(offset!=0)
        {
            if(offset==1)
            {
                jnrA        = jjnr[k];   
                j3A         = 3*jnrA;  
                GMX_MM_LOAD_1JCOORD_1ATOM_PS(x+j3A,jx,jy,jz);
                GMX_MM_LOAD_1VALUE_PS(gb_radius+jnrA,raj);
                GMX_MM_LOAD_1VALUE_PS(vsolv+jnrA,vaj);
                mask        = mask1;
            } 
            else if(offset==2)
            {
                jnrA        = jjnr[k];   
                jnrB        = jjnr[k+1];
                j3A         = 3*jnrA;  
                j3B         = 3*jnrB;
                GMX_MM_LOAD_2JCOORD_1ATOM_PS(x+j3A,x+j3B,jx,jy,jz);
                GMX_MM_LOAD_2VALUES_PS(gb_radius+jnrA,gb_radius+jnrB,raj);
                GMX_MM_LOAD_2VALUES_PS(vsolv+jnrA,vsolv+jnrB,vaj);
                mask        = mask2;
            }
            else
            {
                /* offset must be 3 */   
                jnrA        = jjnr[k];   
                jnrB        = jjnr[k+1];
                jnrC        = jjnr[k+2];
                j3A         = 3*jnrA;  
                j3B         = 3*jnrB;
                j3C         = 3*jnrC;
                GMX_MM_LOAD_3JCOORD_1ATOM_PS(x+j3A,x+j3B,x+j3C,jx,jy,jz);
                GMX_MM_LOAD_3VALUES_PS(gb_radius+jnrA,gb_radius+jnrB,gb_radius+jnrC,raj);
                GMX_MM_LOAD_3VALUES_PS(vsolv+jnrA,vsolv+jnrB,vsolv+jnrC,vaj);
                mask        = mask3;
            }

			dx          = _mm_sub_ps(ix,jx);
			dy          = _mm_sub_ps(iy,jy);
			dz          = _mm_sub_ps(iz,jz);
            
            rsq         = gmx_mm_calc_rsq(dx,dy,dz);
            rinv        = gmx_mm_invsqrt_ps(rsq);
            rinv2       = _mm_mul_ps(rinv,rinv);
            rinv4       = _mm_mul_ps(rinv2,rinv2);
            rinv6       = _mm_mul_ps(rinv4,rinv2);
            
            rvdw        = _mm_add_ps(rai,raj);
            ratio       = _mm_mul_ps(rsq, gmx_mm_inv_ps( _mm_mul_ps(rvdw,rvdw)));
            
            mask_cmp    = _mm_cmple_ps(ratio,still_p5inv);
            
            if(0 == _mm_movemask_ps(mask_cmp))
            {
                /* if ratio>still_p5inv for ALL elements */
                ccf         = one;
                dccf        = _mm_setzero_ps();
            }
            else 
            {
                ratio       = _mm_min_ps(ratio,still_p5inv);
                theta       = _mm_mul_ps(ratio,still_pip5);
                gmx_mm_sincos_ps(theta,&sinq,&cosq);            
                term        = _mm_mul_ps(half,_mm_sub_ps(one,cosq));
                ccf         = _mm_mul_ps(term,term);
                dccf        = _mm_mul_ps(_mm_mul_ps(two,term),
                                         _mm_mul_ps(sinq,theta));
            }

            prod        = _mm_mul_ps(still_p4,vaj);
            icf4        = _mm_mul_ps(ccf,rinv4);
            icf6        = _mm_mul_ps( _mm_sub_ps( _mm_mul_ps(four,ccf),dccf), rinv6);
            
            gpi           = _mm_add_ps(gpi, _mm_mul_ps(prod,icf4));
            
            _mm_store_ps(dadx,_mm_mul_ps(prod,icf6));
            dadx+=4;
            _mm_store_ps(dadx,_mm_mul_ps(prod_ai,icf6));
            dadx+=4;
            
            tmp = _mm_mul_ps(prod_ai,icf4);
            
            if(offset==1)
            {
                GMX_MM_INCREMENT_1VALUE_PS(work+jnrA,tmp);
            } 
            else if(offset==2)
            {
                GMX_MM_INCREMENT_2VALUES_PS(work+jnrA,work+jnrB,tmp);
            }
            else
            {
                /* offset must be 3 */
                GMX_MM_INCREMENT_3VALUES_PS(work+jnrA,work+jnrB,work+jnrC,tmp);
            }
        }
        gmx_mm_update_1pot_ps(gpi,work+ii);
	}
	
	/* Sum up the polarization energy from other nodes */
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms, work, cr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd, work);
	}
	
	/* Compute the radii */
	for(i=0;i<nl->nri;i++)
	{		
		ii               = nl->iinr[i];
		
		if(born->use[ii] != 0)
		{
			gpi_ai           = born->gpol[ii] + work[ii]; /* add gpi to the initial pol energy gpi_ai*/
			gpi2             = gpi_ai * gpi_ai;
			born->bRad[ii]   = factor*gmx_invsqrt(gpi2);
			fr->invsqrta[ii] = gmx_invsqrt(born->bRad[ii]);
		}
	}
		
	/* Extra (local) communication required for DD */
	if(DOMAINDECOMP(cr))
	{
		dd_atom_spread_real(cr->dd, born->bRad);
		dd_atom_spread_real(cr->dd, fr->invsqrta);
	}
    
	return 0;	
}


int 
calc_gb_rad_hct_obc_sse(t_commrec *cr, t_forcerec * fr, int natoms, gmx_localtop_t *top,
                        const t_atomtypes *atype, float *x, t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md,int gb_algorithm)
{
	int i,ai,k,n,ii,ii3,is3,nj0,nj1,at0,at1,offset;
    int jnrA,jnrB,jnrC,jnrD;
    int j3A,j3B,j3C,j3D;
    int jnrE,jnrF,jnrG,jnrH;
    int j3E,j3F,j3G,j3H;
	float shX,shY,shZ;
	float rr,rr_inv,rr_inv2,sum_tmp,sum,sum2,sum3,gbr;
	float sum_ai2, sum_ai3,tsum,tchain,doffset;
	float *obc_param;
    float *gb_radius;
    float *work;
    int *  jjnr;
    float *dadx;
    float *shiftvec;
    float min_rad,rad;
    
	__m128 ix,iy,iz,jx,jy,jz;
	__m128 dx,dy,dz,t1,t2,t3,t4;
	__m128 rsq,rinv,r;
	__m128 rai,rai_inv,raj, raj_inv,rai_inv2,sk,sk2,lij,dlij,duij;
	__m128 uij,lij2,uij2,lij3,uij3,diff2;
	__m128 lij_inv,sk2_inv,prod,log_term,tmp,tmp_sum;
	__m128 sum_ai, tmp_ai,sk_ai,sk_aj,sk2_ai,sk2_aj,sk2_rinv;
	__m128 dadx1,dadx2;
    __m128 logterm;
	__m128 mask;
	__m128 obc_mask1,obc_mask2,obc_mask3;
    __m128 jxB,jyB,jzB,t1B,t2B,t3B,t4B;
    __m128 dxB,dyB,dzB,rsqB,rinvB,rB;
	__m128 rajB, raj_invB,rai_inv2B,sk2B,lijB,dlijB,duijB;
	__m128 uijB,lij2B,uij2B,lij3B,uij3B,diff2B;
	__m128 lij_invB,sk2_invB,prodB;
	__m128 sk_ajB,sk2_ajB,sk2_rinvB;
	__m128 dadx1B,dadx2B;
    __m128 logtermB;
    __m128 obc_mask1B,obc_mask2B,obc_mask3B;

    __m128   mask1 = _mm_castsi128_ps( _mm_set_epi32(0, 0, 0, 0xffffffff) );
	__m128   mask2 = _mm_castsi128_ps( _mm_set_epi32(0, 0, 0xffffffff, 0xffffffff) );
	__m128   mask3 = _mm_castsi128_ps( _mm_set_epi32(0, 0xffffffff, 0xffffffff, 0xffffffff) );
        
    __m128 oneeighth   = _mm_set1_ps(0.125);
    __m128 onefourth   = _mm_set1_ps(0.25);

	const __m128 neg   = {-1.0f , -1.0f , -1.0f , -1.0f };
	const __m128 zero  = {0.0f , 0.0f , 0.0f , 0.0f };
	const __m128 half  = {0.5f , 0.5f , 0.5f , 0.5f };
	const __m128 one   = {1.0f , 1.0f , 1.0f , 1.0f };
	const __m128 two   = {2.0f , 2.0f , 2.0f , 2.0f };
	const __m128 three = {3.0f , 3.0f , 3.0f , 3.0f };
	
	/* Set the dielectric offset */
	doffset   = born->gb_doffset;
	gb_radius = born->gb_radius;
    obc_param = born->param;
    work      = born->gpol_hct_work;
    jjnr      = nl->jjnr;
    dadx      = fr->dadx;
    shiftvec  = fr->shift_vec[0];
    
    jx        = _mm_setzero_ps();
    jy        = _mm_setzero_ps();
    jz        = _mm_setzero_ps();
    
    jnrA = jnrB = jnrC = jnrD = 0;
    
	for(i=0;i<born->nr;i++)
	{
		work[i] = 0;
	}
	
	for(i=0;i<nl->nri;i++)
	{
        ii     = nl->iinr[i];
		ii3	   = ii*3;
        is3    = 3*nl->shift[i];     
        shX    = shiftvec[is3];  
        shY    = shiftvec[is3+1];
        shZ    = shiftvec[is3+2];
        nj0    = nl->jindex[i];      
        nj1    = nl->jindex[i+1];    
        
        ix     = _mm_set1_ps(shX+x[ii3+0]);
		iy     = _mm_set1_ps(shY+x[ii3+1]);
		iz     = _mm_set1_ps(shZ+x[ii3+2]);
		
		offset = (nj1-nj0)%4;

		rai    = _mm_load1_ps(gb_radius+ii);
		rai_inv= gmx_mm_inv_ps(rai);
				
		sum_ai = _mm_setzero_ps();
		
		sk_ai  = _mm_load1_ps(born->param+ii);
		sk2_ai = _mm_mul_ps(sk_ai,sk_ai);
				
		for(k=nj0;k<nj1-4-offset;k+=8)
		{
			jnrA        = jjnr[k];   
			jnrB        = jjnr[k+1];
			jnrC        = jjnr[k+2];
			jnrD        = jjnr[k+3];
			jnrE        = jjnr[k+4];   
			jnrF        = jjnr[k+5];
			jnrG        = jjnr[k+6];
			jnrH        = jjnr[k+7];
			
            j3A         = 3*jnrA;  
			j3B         = 3*jnrB;
			j3C         = 3*jnrC;
			j3D         = 3*jnrD;
            j3E         = 3*jnrE;  
			j3F         = 3*jnrF;
			j3G         = 3*jnrG;
			j3H         = 3*jnrH;
            
            GMX_MM_LOAD_4JCOORD_1ATOM_PS(x+j3A,x+j3B,x+j3C,x+j3D,jx,jy,jz);
            GMX_MM_LOAD_4JCOORD_1ATOM_PS(x+j3E,x+j3F,x+j3G,x+j3H,jxB,jyB,jzB);
            GMX_MM_LOAD_4VALUES_PS(gb_radius+jnrA,gb_radius+jnrB,gb_radius+jnrC,gb_radius+jnrD,raj);
            GMX_MM_LOAD_4VALUES_PS(gb_radius+jnrE,gb_radius+jnrF,gb_radius+jnrG,gb_radius+jnrH,rajB);
            GMX_MM_LOAD_4VALUES_PS(obc_param+jnrA,obc_param+jnrB,obc_param+jnrC,obc_param+jnrD,sk_aj);
            GMX_MM_LOAD_4VALUES_PS(obc_param+jnrE,obc_param+jnrF,obc_param+jnrG,obc_param+jnrH,sk_ajB);
			
            dx    = _mm_sub_ps(ix, jx);
			dy    = _mm_sub_ps(iy, jy);
			dz    = _mm_sub_ps(iz, jz);
            dxB   = _mm_sub_ps(ix, jxB);
			dyB   = _mm_sub_ps(iy, jyB);
			dzB   = _mm_sub_ps(iz, jzB);
			
            rsq         = gmx_mm_calc_rsq(dx,dy,dz);
            rsqB        = gmx_mm_calc_rsq(dxB,dyB,dzB);
			            
            rinv        = gmx_mm_invsqrt_ps(rsq);
            r           = _mm_mul_ps(rsq,rinv);
            rinvB       = gmx_mm_invsqrt_ps(rsqB);
            rB          = _mm_mul_ps(rsqB,rinvB);

			/* Compute raj_inv aj1-4 */
            raj_inv     = gmx_mm_inv_ps(raj);
            raj_invB    = gmx_mm_inv_ps(rajB);

            /* Evaluate influence of atom aj -> ai */
            t1            = _mm_add_ps(r,sk_aj);
            t2            = _mm_sub_ps(r,sk_aj);
            t3            = _mm_sub_ps(sk_aj,r);
            t1B           = _mm_add_ps(rB,sk_ajB);
            t2B           = _mm_sub_ps(rB,sk_ajB);
            t3B           = _mm_sub_ps(sk_ajB,rB);
            obc_mask1     = _mm_cmplt_ps(rai, t1);
            obc_mask2     = _mm_cmplt_ps(rai, t2);
            obc_mask3     = _mm_cmplt_ps(rai, t3);
            obc_mask1B    = _mm_cmplt_ps(rai, t1B);
            obc_mask2B    = _mm_cmplt_ps(rai, t2B);
            obc_mask3B    = _mm_cmplt_ps(rai, t3B);

            uij           = gmx_mm_inv_ps(t1);
            lij           = _mm_or_ps(   _mm_and_ps(obc_mask2,gmx_mm_inv_ps(t2)),
                                      _mm_andnot_ps(obc_mask2,rai_inv));
            dlij          = _mm_and_ps(one,obc_mask2);
            uij2          = _mm_mul_ps(uij, uij);
            uij3          = _mm_mul_ps(uij2,uij);
            lij2          = _mm_mul_ps(lij, lij);
            lij3          = _mm_mul_ps(lij2,lij);

            uijB          = gmx_mm_inv_ps(t1B);
            lijB          = _mm_or_ps(   _mm_and_ps(obc_mask2B,gmx_mm_inv_ps(t2B)),
                                      _mm_andnot_ps(obc_mask2B,rai_inv));
            dlijB         = _mm_and_ps(one,obc_mask2B);
            uij2B         = _mm_mul_ps(uijB, uijB);
            uij3B         = _mm_mul_ps(uij2B,uijB);
            lij2B         = _mm_mul_ps(lijB, lijB);
            lij3B         = _mm_mul_ps(lij2B,lijB);

            diff2         = _mm_sub_ps(uij2,lij2);
            lij_inv       = gmx_mm_invsqrt_ps(lij2);
            sk2_aj        = _mm_mul_ps(sk_aj,sk_aj);
            sk2_rinv      = _mm_mul_ps(sk2_aj,rinv);
            prod          = _mm_mul_ps(onefourth,sk2_rinv);

            diff2B        = _mm_sub_ps(uij2B,lij2B);
            lij_invB      = gmx_mm_invsqrt_ps(lij2B);
            sk2_ajB       = _mm_mul_ps(sk_ajB,sk_ajB);
            sk2_rinvB     = _mm_mul_ps(sk2_ajB,rinvB);
            prodB         = _mm_mul_ps(onefourth,sk2_rinvB);

            logterm       = gmx_mm_log_ps(_mm_mul_ps(uij,lij_inv));
            logtermB      = gmx_mm_log_ps(_mm_mul_ps(uijB,lij_invB));
            
            t1            = _mm_sub_ps(lij,uij);
            t2            = _mm_mul_ps(diff2,
                                       _mm_sub_ps(_mm_mul_ps(onefourth,r),
                                                  prod));
            t3            = _mm_mul_ps(half,_mm_mul_ps(rinv,logterm));
            t1            = _mm_add_ps(t1,_mm_add_ps(t2,t3));
            t4            = _mm_mul_ps(two,_mm_sub_ps(rai_inv,lij));
            t4            = _mm_and_ps(t4,obc_mask3);
            t1            = _mm_mul_ps(half,_mm_add_ps(t1,t4));
            
            t1B           = _mm_sub_ps(lijB,uijB);
            t2B           = _mm_mul_ps(diff2B,
                                       _mm_sub_ps(_mm_mul_ps(onefourth,rB),
                                                  prodB));
            t3B           = _mm_mul_ps(half,_mm_mul_ps(rinvB,logtermB));
            t1B           = _mm_add_ps(t1B,_mm_add_ps(t2B,t3B));
            t4B           = _mm_mul_ps(two,_mm_sub_ps(rai_inv,lijB));
            t4B           = _mm_and_ps(t4B,obc_mask3B);
            t1B           = _mm_mul_ps(half,_mm_add_ps(t1B,t4B));
            
            sum_ai        = _mm_add_ps(sum_ai, _mm_add_ps( _mm_and_ps(t1,obc_mask1), _mm_and_ps(t1B,obc_mask1B) ));
            
            t1            = _mm_add_ps(_mm_mul_ps(half,lij2),
                                       _mm_mul_ps(prod,lij3));
            t1            = _mm_sub_ps(t1,
                                       _mm_mul_ps(onefourth,
                                                  _mm_add_ps(_mm_mul_ps(lij,rinv),
                                                             _mm_mul_ps(lij3,r))));
            t2            = _mm_mul_ps(onefourth,
                                       _mm_add_ps(_mm_mul_ps(uij,rinv),
                                                  _mm_mul_ps(uij3,r)));
            t2            = _mm_sub_ps(t2,
                                       _mm_add_ps(_mm_mul_ps(half,uij2),
                                                      _mm_mul_ps(prod,uij3)));
            t3            = _mm_mul_ps(_mm_mul_ps(onefourth,logterm),
                                       _mm_mul_ps(rinv,rinv));
            t3            = _mm_sub_ps(t3,
                                       _mm_mul_ps(_mm_mul_ps(diff2,oneeighth),
                                                  _mm_add_ps(one,
                                                             _mm_mul_ps(sk2_rinv,rinv))));
            t1            = _mm_mul_ps(rinv,
                                       _mm_add_ps(_mm_mul_ps(dlij,t1),
                                                  _mm_add_ps(t2,t3)));
            

            
            t1B           = _mm_add_ps(_mm_mul_ps(half,lij2B),
                                       _mm_mul_ps(prodB,lij3B));
            t1B           = _mm_sub_ps(t1B,
                                       _mm_mul_ps(onefourth,
                                                  _mm_add_ps(_mm_mul_ps(lijB,rinvB),
                                                             _mm_mul_ps(lij3B,rB))));
            t2B           = _mm_mul_ps(onefourth,
                                       _mm_add_ps(_mm_mul_ps(uijB,rinvB),
                                                  _mm_mul_ps(uij3B,rB)));
            t2B           = _mm_sub_ps(t2B,
                                       _mm_add_ps(_mm_mul_ps(half,uij2B),
                                                  _mm_mul_ps(prodB,uij3B)));
            t3B           = _mm_mul_ps(_mm_mul_ps(onefourth,logtermB),
                                       _mm_mul_ps(rinvB,rinvB));
            t3B           = _mm_sub_ps(t3B,
                                       _mm_mul_ps(_mm_mul_ps(diff2B,oneeighth),
                                                  _mm_add_ps(one,
                                                             _mm_mul_ps(sk2_rinvB,rinvB))));
            t1B           = _mm_mul_ps(rinvB,
                                       _mm_add_ps(_mm_mul_ps(dlijB,t1B),
                                                  _mm_add_ps(t2B,t3B)));
            
            dadx1         = _mm_and_ps(t1,obc_mask1);
            dadx1B        = _mm_and_ps(t1B,obc_mask1B);


            /* Evaluate influence of atom ai -> aj */
            t1            = _mm_add_ps(r,sk_ai);
            t2            = _mm_sub_ps(r,sk_ai);
            t3            = _mm_sub_ps(sk_ai,r);
            t1B           = _mm_add_ps(rB,sk_ai);
            t2B           = _mm_sub_ps(rB,sk_ai);
            t3B           = _mm_sub_ps(sk_ai,rB);
            obc_mask1     = _mm_cmplt_ps(raj, t1);
            obc_mask2     = _mm_cmplt_ps(raj, t2);
            obc_mask3     = _mm_cmplt_ps(raj, t3);
            obc_mask1B    = _mm_cmplt_ps(rajB, t1B);
            obc_mask2B    = _mm_cmplt_ps(rajB, t2B);
            obc_mask3B    = _mm_cmplt_ps(rajB, t3B);
            
            uij           = gmx_mm_inv_ps(t1);
            lij           = _mm_or_ps(   _mm_and_ps(obc_mask2,gmx_mm_inv_ps(t2)),
                                      _mm_andnot_ps(obc_mask2,raj_inv));
            dlij          = _mm_and_ps(one,obc_mask2);
            uij2          = _mm_mul_ps(uij, uij);
            uij3          = _mm_mul_ps(uij2,uij);
            lij2          = _mm_mul_ps(lij, lij);
            lij3          = _mm_mul_ps(lij2,lij);

            uijB          = gmx_mm_inv_ps(t1B);
            lijB          = _mm_or_ps(   _mm_and_ps(obc_mask2B,gmx_mm_inv_ps(t2B)),
                                      _mm_andnot_ps(obc_mask2B,raj_invB));
            dlijB         = _mm_and_ps(one,obc_mask2B);
            uij2B         = _mm_mul_ps(uijB, uijB);
            uij3B         = _mm_mul_ps(uij2B,uijB);
            lij2B         = _mm_mul_ps(lijB, lijB);
            lij3B         = _mm_mul_ps(lij2B,lijB);

            diff2         = _mm_sub_ps(uij2,lij2);
            lij_inv       = gmx_mm_invsqrt_ps(lij2);
            sk2_rinv      = _mm_mul_ps(sk2_ai,rinv);
            prod          = _mm_mul_ps(onefourth,sk2_rinv);

            diff2B        = _mm_sub_ps(uij2B,lij2B);
            lij_invB      = gmx_mm_invsqrt_ps(lij2B);
            sk2_rinvB     = _mm_mul_ps(sk2_ai,rinvB);
            prodB         = _mm_mul_ps(onefourth,sk2_rinvB);

            logterm       = gmx_mm_log_ps(_mm_mul_ps(uij,lij_inv));
            logtermB      = gmx_mm_log_ps(_mm_mul_ps(uijB,lij_invB));

            t1            = _mm_sub_ps(lij,uij);
            t2            = _mm_mul_ps(diff2,
                                       _mm_sub_ps(_mm_mul_ps(onefourth,r),
                                                  prod));
            t3            = _mm_mul_ps(half,_mm_mul_ps(rinv,logterm));
            t1            = _mm_add_ps(t1,_mm_add_ps(t2,t3));
            t4            = _mm_mul_ps(two,_mm_sub_ps(raj_inv,lij));
            t4            = _mm_and_ps(t4,obc_mask3);
            t1            = _mm_mul_ps(half,_mm_add_ps(t1,t4));

            t1B           = _mm_sub_ps(lijB,uijB);
            t2B           = _mm_mul_ps(diff2B,
                                       _mm_sub_ps(_mm_mul_ps(onefourth,rB),
                                                  prodB));
            t3B           = _mm_mul_ps(half,_mm_mul_ps(rinvB,logtermB));
            t1B           = _mm_add_ps(t1B,_mm_add_ps(t2B,t3B));
            t4B           = _mm_mul_ps(two,_mm_sub_ps(raj_invB,lijB));
            t4B           = _mm_and_ps(t4B,obc_mask3B);
            t1B           = _mm_mul_ps(half,_mm_add_ps(t1B,t4B));
            
            GMX_MM_INCREMENT_4VALUES_PS(work+jnrA,work+jnrB,work+jnrC,work+jnrD,_mm_and_ps(t1,obc_mask1));
            GMX_MM_INCREMENT_4VALUES_PS(work+jnrE,work+jnrF,work+jnrG,work+jnrH,_mm_and_ps(t1B,obc_mask1B));
            
            t1            = _mm_add_ps(_mm_mul_ps(half,lij2),
                                       _mm_mul_ps(prod,lij3));
            t1            = _mm_sub_ps(t1,
                                       _mm_mul_ps(onefourth,
                                                  _mm_add_ps(_mm_mul_ps(lij,rinv),
                                                             _mm_mul_ps(lij3,r))));
            t2            = _mm_mul_ps(onefourth,
                                       _mm_add_ps(_mm_mul_ps(uij,rinv),
                                                  _mm_mul_ps(uij3,r)));
            t2            = _mm_sub_ps(t2,
                                       _mm_add_ps(_mm_mul_ps(half,uij2),
                                                  _mm_mul_ps(prod,uij3)));
            t3            = _mm_mul_ps(_mm_mul_ps(onefourth,logterm),
                                       _mm_mul_ps(rinv,rinv));
            t3            = _mm_sub_ps(t3,
                                       _mm_mul_ps(_mm_mul_ps(diff2,oneeighth),
                                                  _mm_add_ps(one,
                                                             _mm_mul_ps(sk2_rinv,rinv))));
            t1            = _mm_mul_ps(rinv,
                                       _mm_add_ps(_mm_mul_ps(dlij,t1),
                                                  _mm_add_ps(t2,t3)));

            
            t1B           = _mm_add_ps(_mm_mul_ps(half,lij2B),
                                       _mm_mul_ps(prodB,lij3B));
            t1B           = _mm_sub_ps(t1B,
                                       _mm_mul_ps(onefourth,
                                                  _mm_add_ps(_mm_mul_ps(lijB,rinvB),
                                                             _mm_mul_ps(lij3B,rB))));
            t2B           = _mm_mul_ps(onefourth,
                                       _mm_add_ps(_mm_mul_ps(uijB,rinvB),
                                                  _mm_mul_ps(uij3B,rB)));
            t2B           = _mm_sub_ps(t2B,
                                       _mm_add_ps(_mm_mul_ps(half,uij2B),
                                                  _mm_mul_ps(prodB,uij3B)));
            t3B           = _mm_mul_ps(_mm_mul_ps(onefourth,logtermB),
                                       _mm_mul_ps(rinvB,rinvB));
            t3B           = _mm_sub_ps(t3B,
                                       _mm_mul_ps(_mm_mul_ps(diff2B,oneeighth),
                                                  _mm_add_ps(one,
                                                             _mm_mul_ps(sk2_rinvB,rinvB))));
            t1B           = _mm_mul_ps(rinvB,
                                       _mm_add_ps(_mm_mul_ps(dlijB,t1B),
                                                  _mm_add_ps(t2B,t3B)));

            
            dadx2         = _mm_and_ps(t1,obc_mask1);
            dadx2B        = _mm_and_ps(t1B,obc_mask1B);
            
            _mm_store_ps(dadx,dadx1);
            dadx += 4;
            _mm_store_ps(dadx,dadx2);
            dadx += 4;
            _mm_store_ps(dadx,dadx1B);
            dadx += 4;
            _mm_store_ps(dadx,dadx2B);
            dadx += 4;
            
        } /* end normal inner loop */
        
		for(;k<nj1-offset;k+=4)
		{
			jnrA        = jjnr[k];   
			jnrB        = jjnr[k+1];
			jnrC        = jjnr[k+2];
			jnrD        = jjnr[k+3];
			
            j3A         = 3*jnrA;  
			j3B         = 3*jnrB;
			j3C         = 3*jnrC;
			j3D         = 3*jnrD;
            
            GMX_MM_LOAD_4JCOORD_1ATOM_PS(x+j3A,x+j3B,x+j3C,x+j3D,jx,jy,jz);
            GMX_MM_LOAD_4VALUES_PS(gb_radius+jnrA,gb_radius+jnrB,gb_radius+jnrC,gb_radius+jnrD,raj);
            GMX_MM_LOAD_4VALUES_PS(obc_param+jnrA,obc_param+jnrB,obc_param+jnrC,obc_param+jnrD,sk_aj);
			
            dx    = _mm_sub_ps(ix, jx);
			dy    = _mm_sub_ps(iy, jy);
			dz    = _mm_sub_ps(iz, jz);
			
            rsq         = gmx_mm_calc_rsq(dx,dy,dz);
            
            rinv        = gmx_mm_invsqrt_ps(rsq);
            r           = _mm_mul_ps(rsq,rinv);
            
			/* Compute raj_inv aj1-4 */
            raj_inv     = gmx_mm_inv_ps(raj);
            
            /* Evaluate influence of atom aj -> ai */
            t1            = _mm_add_ps(r,sk_aj);
            obc_mask1     = _mm_cmplt_ps(rai, t1);
            
            if(_mm_movemask_ps(obc_mask1))
            {
                /* If any of the elements has rai<dr+sk, this is executed */
                t2            = _mm_sub_ps(r,sk_aj);
                t3            = _mm_sub_ps(sk_aj,r);
                
                obc_mask2     = _mm_cmplt_ps(rai, t2);
                obc_mask3     = _mm_cmplt_ps(rai, t3);
                
                uij           = gmx_mm_inv_ps(t1);
                lij           = _mm_or_ps(   _mm_and_ps(obc_mask2,gmx_mm_inv_ps(t2)),
                                          _mm_andnot_ps(obc_mask2,rai_inv));
                dlij          = _mm_and_ps(one,obc_mask2);
                uij2          = _mm_mul_ps(uij, uij);
                uij3          = _mm_mul_ps(uij2,uij);
                lij2          = _mm_mul_ps(lij, lij);
                lij3          = _mm_mul_ps(lij2,lij);
                diff2         = _mm_sub_ps(uij2,lij2);
                lij_inv       = gmx_mm_invsqrt_ps(lij2);
                sk2_aj        = _mm_mul_ps(sk_aj,sk_aj);
                sk2_rinv      = _mm_mul_ps(sk2_aj,rinv);
                prod          = _mm_mul_ps(onefourth,sk2_rinv);
                logterm       = gmx_mm_log_ps(_mm_mul_ps(uij,lij_inv));
                t1            = _mm_sub_ps(lij,uij);
                t2            = _mm_mul_ps(diff2,
                                           _mm_sub_ps(_mm_mul_ps(onefourth,r),
                                                      prod));
                t3            = _mm_mul_ps(half,_mm_mul_ps(rinv,logterm));
                t1            = _mm_add_ps(t1,_mm_add_ps(t2,t3));
                t4            = _mm_mul_ps(two,_mm_sub_ps(rai_inv,lij));
                t4            = _mm_and_ps(t4,obc_mask3);
                t1            = _mm_mul_ps(half,_mm_add_ps(t1,t4));
                sum_ai        = _mm_add_ps(sum_ai,_mm_and_ps(t1,obc_mask1));
                t1            = _mm_add_ps(_mm_mul_ps(half,lij2),
                                           _mm_mul_ps(prod,lij3));
                t1            = _mm_sub_ps(t1,
                                           _mm_mul_ps(onefourth,
                                                      _mm_add_ps(_mm_mul_ps(lij,rinv),
                                                                 _mm_mul_ps(lij3,r))));
                t2            = _mm_mul_ps(onefourth,
                                           _mm_add_ps(_mm_mul_ps(uij,rinv),
                                                      _mm_mul_ps(uij3,r)));
                t2            = _mm_sub_ps(t2,
                                           _mm_add_ps(_mm_mul_ps(half,uij2),
                                                      _mm_mul_ps(prod,uij3)));
                t3            = _mm_mul_ps(_mm_mul_ps(onefourth,logterm),
                                           _mm_mul_ps(rinv,rinv));
                t3            = _mm_sub_ps(t3,
                                           _mm_mul_ps(_mm_mul_ps(diff2,oneeighth),
                                                      _mm_add_ps(one,
                                                                 _mm_mul_ps(sk2_rinv,rinv))));
                t1            = _mm_mul_ps(rinv,
                                           _mm_add_ps(_mm_mul_ps(dlij,t1),
                                                      _mm_add_ps(t2,t3)));
                
                dadx1         = _mm_and_ps(t1,obc_mask1);
            }
            else 
            {
                dadx1         = _mm_setzero_ps();
            }
            
            /* Evaluate influence of atom ai -> aj */
            t1            = _mm_add_ps(r,sk_ai);
            obc_mask1     = _mm_cmplt_ps(raj, t1);
            
            if(_mm_movemask_ps(obc_mask1))
            {
                t2            = _mm_sub_ps(r,sk_ai);
                t3            = _mm_sub_ps(sk_ai,r);
                obc_mask2     = _mm_cmplt_ps(raj, t2);
                obc_mask3     = _mm_cmplt_ps(raj, t3);
                
                uij           = gmx_mm_inv_ps(t1);
                lij           = _mm_or_ps(   _mm_and_ps(obc_mask2,gmx_mm_inv_ps(t2)),
                                          _mm_andnot_ps(obc_mask2,raj_inv));
                dlij          = _mm_and_ps(one,obc_mask2);
                uij2          = _mm_mul_ps(uij, uij);
                uij3          = _mm_mul_ps(uij2,uij);
                lij2          = _mm_mul_ps(lij, lij);
                lij3          = _mm_mul_ps(lij2,lij);
                diff2         = _mm_sub_ps(uij2,lij2);
                lij_inv       = gmx_mm_invsqrt_ps(lij2);
                sk2_rinv      = _mm_mul_ps(sk2_ai,rinv);
                prod          = _mm_mul_ps(onefourth,sk2_rinv);
                logterm       = gmx_mm_log_ps(_mm_mul_ps(uij,lij_inv));
                t1            = _mm_sub_ps(lij,uij);
                t2            = _mm_mul_ps(diff2,
                                           _mm_sub_ps(_mm_mul_ps(onefourth,r),
                                                      prod));
                t3            = _mm_mul_ps(half,_mm_mul_ps(rinv,logterm));
                t1            = _mm_add_ps(t1,_mm_add_ps(t2,t3));
                t4            = _mm_mul_ps(two,_mm_sub_ps(raj_inv,lij));
                t4            = _mm_and_ps(t4,obc_mask3);
                t1            = _mm_mul_ps(half,_mm_add_ps(t1,t4));
                
                GMX_MM_INCREMENT_4VALUES_PS(work+jnrA,work+jnrB,work+jnrC,work+jnrD,_mm_and_ps(t1,obc_mask1));
                
                t1            = _mm_add_ps(_mm_mul_ps(half,lij2),
                                           _mm_mul_ps(prod,lij3));
                t1            = _mm_sub_ps(t1,
                                           _mm_mul_ps(onefourth,
                                                      _mm_add_ps(_mm_mul_ps(lij,rinv),
                                                                 _mm_mul_ps(lij3,r))));
                t2            = _mm_mul_ps(onefourth,
                                           _mm_add_ps(_mm_mul_ps(uij,rinv),
                                                      _mm_mul_ps(uij3,r)));
                t2            = _mm_sub_ps(t2,
                                           _mm_add_ps(_mm_mul_ps(half,uij2),
                                                      _mm_mul_ps(prod,uij3)));
                t3            = _mm_mul_ps(_mm_mul_ps(onefourth,logterm),
                                           _mm_mul_ps(rinv,rinv));
                t3            = _mm_sub_ps(t3,
                                           _mm_mul_ps(_mm_mul_ps(diff2,oneeighth),
                                                      _mm_add_ps(one,
                                                                 _mm_mul_ps(sk2_rinv,rinv))));
                t1            = _mm_mul_ps(rinv,
                                           _mm_add_ps(_mm_mul_ps(dlij,t1),
                                                      _mm_add_ps(t2,t3)));
                dadx2         = _mm_and_ps(t1,obc_mask1);
            }
            else
            {
                dadx2         = _mm_setzero_ps();
            }
            
            _mm_store_ps(dadx,dadx1);
            dadx += 4;
            _mm_store_ps(dadx,dadx2);
            dadx += 4;            
        } /* end normal inner loop */
        
        if(offset!=0)
        {
            if(offset==1)
            {
                jnrA        = jjnr[k];   
                j3A         = 3*jnrA;  
                GMX_MM_LOAD_1JCOORD_1ATOM_PS(x+j3A,jx,jy,jz);
                GMX_MM_LOAD_1VALUE_PS(gb_radius+jnrA,raj);
                GMX_MM_LOAD_1VALUE_PS(obc_param+jnrA,sk_aj);
                mask        = mask1;
            } 
            else if(offset==2)
            {
                jnrA        = jjnr[k];   
                jnrB        = jjnr[k+1];
                j3A         = 3*jnrA;  
                j3B         = 3*jnrB;
                GMX_MM_LOAD_2JCOORD_1ATOM_PS(x+j3A,x+j3B,jx,jy,jz);
                GMX_MM_LOAD_2VALUES_PS(gb_radius+jnrA,gb_radius+jnrB,raj);
                GMX_MM_LOAD_2VALUES_PS(obc_param+jnrA,obc_param+jnrB,sk_aj);
                mask        = mask2;
            }
            else
            {
                /* offset must be 3 */   
                jnrA        = jjnr[k];   
                jnrB        = jjnr[k+1];
                jnrC        = jjnr[k+2];
                j3A         = 3*jnrA;  
                j3B         = 3*jnrB;
                j3C         = 3*jnrC;
                GMX_MM_LOAD_3JCOORD_1ATOM_PS(x+j3A,x+j3B,x+j3C,jx,jy,jz);
                GMX_MM_LOAD_3VALUES_PS(gb_radius+jnrA,gb_radius+jnrB,gb_radius+jnrC,raj);
                GMX_MM_LOAD_3VALUES_PS(obc_param+jnrA,obc_param+jnrB,obc_param+jnrC,sk_aj);
                mask        = mask3;
            }

            dx    = _mm_sub_ps(ix, jx);
			dy    = _mm_sub_ps(iy, jy);
			dz    = _mm_sub_ps(iz, jz);
			
            rsq         = gmx_mm_calc_rsq(dx,dy,dz);
            
            rinv        = gmx_mm_invsqrt_ps(rsq);
            r           = _mm_mul_ps(rsq,rinv);
            
			/* Compute raj_inv aj1-4 */
            raj_inv     = gmx_mm_inv_ps(raj);
            
            /* Evaluate influence of atom aj -> ai */
            t1            = _mm_add_ps(r,sk_aj);
            obc_mask1     = _mm_cmplt_ps(rai, t1);
            obc_mask1     = _mm_and_ps(obc_mask1,mask);

            if(_mm_movemask_ps(obc_mask1))
            {
                t2            = _mm_sub_ps(r,sk_aj);
                t3            = _mm_sub_ps(sk_aj,r);
                obc_mask2     = _mm_cmplt_ps(rai, t2);
                obc_mask3     = _mm_cmplt_ps(rai, t3);
                
                uij           = gmx_mm_inv_ps(t1);
                lij           = _mm_or_ps(   _mm_and_ps(obc_mask2,gmx_mm_inv_ps(t2)),
                                          _mm_andnot_ps(obc_mask2,rai_inv));
                dlij          = _mm_and_ps(one,obc_mask2);
                uij2          = _mm_mul_ps(uij, uij);
                uij3          = _mm_mul_ps(uij2,uij);
                lij2          = _mm_mul_ps(lij, lij);
                lij3          = _mm_mul_ps(lij2,lij);
                diff2         = _mm_sub_ps(uij2,lij2);
                lij_inv       = gmx_mm_invsqrt_ps(lij2);
                sk2_aj         = _mm_mul_ps(sk_aj,sk_aj);
                sk2_rinv      = _mm_mul_ps(sk2_aj,rinv);
                prod          = _mm_mul_ps(onefourth,sk2_rinv);
                logterm       = gmx_mm_log_ps(_mm_mul_ps(uij,lij_inv));
                t1            = _mm_sub_ps(lij,uij);
                t2            = _mm_mul_ps(diff2,
                                           _mm_sub_ps(_mm_mul_ps(onefourth,r),
                                                      prod));
                t3            = _mm_mul_ps(half,_mm_mul_ps(rinv,logterm));
                t1            = _mm_add_ps(t1,_mm_add_ps(t2,t3));
                t4            = _mm_mul_ps(two,_mm_sub_ps(rai_inv,lij));
                t4            = _mm_and_ps(t4,obc_mask3);
                t1            = _mm_mul_ps(half,_mm_add_ps(t1,t4));
                sum_ai        = _mm_add_ps(sum_ai,_mm_and_ps(t1,obc_mask1));
                t1            = _mm_add_ps(_mm_mul_ps(half,lij2),
                                           _mm_mul_ps(prod,lij3));
                t1            = _mm_sub_ps(t1,
                                           _mm_mul_ps(onefourth,
                                                      _mm_add_ps(_mm_mul_ps(lij,rinv),
                                                                 _mm_mul_ps(lij3,r))));
                t2            = _mm_mul_ps(onefourth,
                                           _mm_add_ps(_mm_mul_ps(uij,rinv),
                                                      _mm_mul_ps(uij3,r)));
                t2            = _mm_sub_ps(t2,
                                           _mm_add_ps(_mm_mul_ps(half,uij2),
                                                      _mm_mul_ps(prod,uij3)));
                t3            = _mm_mul_ps(_mm_mul_ps(onefourth,logterm),
                                           _mm_mul_ps(rinv,rinv));
                t3            = _mm_sub_ps(t3,
                                           _mm_mul_ps(_mm_mul_ps(diff2,oneeighth),
                                                      _mm_add_ps(one,
                                                                 _mm_mul_ps(sk2_rinv,rinv))));
                t1            = _mm_mul_ps(rinv,
                                           _mm_add_ps(_mm_mul_ps(dlij,t1),
                                                      _mm_add_ps(t2,t3)));
                dadx1         = _mm_and_ps(t1,obc_mask1);
            }
            else
            {
                dadx1         = _mm_setzero_ps();
            }

                /* Evaluate influence of atom ai -> aj */
            t1            = _mm_add_ps(r,sk_ai);
            obc_mask1     = _mm_cmplt_ps(raj, t1);
            obc_mask1     = _mm_and_ps(obc_mask1,mask);
            
            if(_mm_movemask_ps(obc_mask1))
            {
                t2            = _mm_sub_ps(r,sk_ai);
                t3            = _mm_sub_ps(sk_ai,r);
                obc_mask2     = _mm_cmplt_ps(raj, t2);
                obc_mask3     = _mm_cmplt_ps(raj, t3);
            
                uij           = gmx_mm_inv_ps(t1);
                lij           = _mm_or_ps(_mm_and_ps(obc_mask2,gmx_mm_inv_ps(t2)),
                                          _mm_andnot_ps(obc_mask2,raj_inv));
                dlij          = _mm_and_ps(one,obc_mask2);
                uij2          = _mm_mul_ps(uij, uij);
                uij3          = _mm_mul_ps(uij2,uij);
                lij2          = _mm_mul_ps(lij, lij);
                lij3          = _mm_mul_ps(lij2,lij);
                diff2         = _mm_sub_ps(uij2,lij2);
                lij_inv       = gmx_mm_invsqrt_ps(lij2);
                sk2_rinv      = _mm_mul_ps(sk2_ai,rinv);
                prod          = _mm_mul_ps(onefourth,sk2_rinv);
                logterm       = gmx_mm_log_ps(_mm_mul_ps(uij,lij_inv));
                t1            = _mm_sub_ps(lij,uij);
                t2            = _mm_mul_ps(diff2,
                                           _mm_sub_ps(_mm_mul_ps(onefourth,r),
                                                      prod));
                t3            = _mm_mul_ps(half,_mm_mul_ps(rinv,logterm));
                t1            = _mm_add_ps(t1,_mm_add_ps(t2,t3));
                t4            = _mm_mul_ps(two,_mm_sub_ps(raj_inv,lij));
                t4            = _mm_and_ps(t4,obc_mask3);
                t1            = _mm_mul_ps(half,_mm_add_ps(t1,t4));
                
                tmp           = _mm_and_ps(t1,obc_mask1);
                
                t1            = _mm_add_ps(_mm_mul_ps(half,lij2),
                                           _mm_mul_ps(prod,lij3));
                t1            = _mm_sub_ps(t1,
                                           _mm_mul_ps(onefourth,
                                                      _mm_add_ps(_mm_mul_ps(lij,rinv),
                                                                 _mm_mul_ps(lij3,r))));
                t2            = _mm_mul_ps(onefourth,
                                           _mm_add_ps(_mm_mul_ps(uij,rinv),
                                                      _mm_mul_ps(uij3,r)));
                t2            = _mm_sub_ps(t2,
                                           _mm_add_ps(_mm_mul_ps(half,uij2),
                                                      _mm_mul_ps(prod,uij3)));
                t3            = _mm_mul_ps(_mm_mul_ps(onefourth,logterm),
                                           _mm_mul_ps(rinv,rinv));
                t3            = _mm_sub_ps(t3,
                                           _mm_mul_ps(_mm_mul_ps(diff2,oneeighth),
                                                      _mm_add_ps(one,
                                                                 _mm_mul_ps(sk2_rinv,rinv))));
                t1            = _mm_mul_ps(rinv,
                                           _mm_add_ps(_mm_mul_ps(dlij,t1),
                                                      _mm_add_ps(t2,t3)));
                dadx2         = _mm_and_ps(t1,obc_mask1);
            }
            else
            {
                dadx2         = _mm_setzero_ps();
                tmp           = _mm_setzero_ps();
            }
            
            _mm_store_ps(dadx,dadx1);
            dadx += 4;
            _mm_store_ps(dadx,dadx2);
            dadx += 4;
            
            if(offset==1)
            {
                GMX_MM_INCREMENT_1VALUE_PS(work+jnrA,tmp);
            } 
            else if(offset==2)
            {
                GMX_MM_INCREMENT_2VALUES_PS(work+jnrA,work+jnrB,tmp);
            }
            else
            {
                /* offset must be 3 */
                GMX_MM_INCREMENT_3VALUES_PS(work+jnrA,work+jnrB,work+jnrC,tmp);
            }
            
        }
        gmx_mm_update_1pot_ps(sum_ai,work+ii);
        
	}
	
	/* Parallel summations */
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms, work, cr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd, work);
	}
	
    if(gb_algorithm==egbHCT)
    {
        /* HCT */
        for(i=0;i<nl->nri;i++)
        {
            ai      = nl->iinr[i];
            
            if(born->use[ai] != 0)
            {
                rr      = top->atomtypes.gb_radius[md->typeA[ai]]-doffset; 
                sum     = 1.0/rr - work[ai];
                min_rad = rr + doffset;
                rad     = 1.0/sum; 
                
                born->bRad[ai]   = rad > min_rad ? rad : min_rad;
                fr->invsqrta[ai] = gmx_invsqrt(born->bRad[ai]);
            }
        }
        
        /* Extra communication required for DD */
        if(DOMAINDECOMP(cr))
        {
            dd_atom_spread_real(cr->dd, born->bRad);
            dd_atom_spread_real(cr->dd, fr->invsqrta);
        }
    }
    else
    {
        /* OBC */
        for(i=0;i<nl->nri;i++)
        {
            ii      = nl->iinr[i];
            
            if(born->use[ii] != 0)
            {
                rr      = top->atomtypes.gb_radius[md->typeA[ii]];
                rr_inv2 = 1.0/rr;
                rr      = rr-doffset; 
                rr_inv  = 1.0/rr;
                sum     = rr * work[ii];
                sum2    = sum  * sum;
                sum3    = sum2 * sum;
                
                tsum    = tanh(born->obc_alpha*sum-born->obc_beta*sum2+born->obc_gamma*sum3);
                born->bRad[ii] = rr_inv - tsum*rr_inv2;
                born->bRad[ii] = 1.0 / born->bRad[ii];
                
                fr->invsqrta[ii]=gmx_invsqrt(born->bRad[ii]);
                
                tchain  = rr * (born->obc_alpha-2*born->obc_beta*sum+3*born->obc_gamma*sum2);
                born->drobc[ii] = (1.0-tsum*tsum)*tchain*rr_inv2;
            }
        }
        /* Extra (local) communication required for DD */
        if(DOMAINDECOMP(cr))
        {
            dd_atom_spread_real(cr->dd, born->bRad);
            dd_atom_spread_real(cr->dd, fr->invsqrta);
            dd_atom_spread_real(cr->dd, born->drobc);
        }
    }

	
	
	return 0;
}



float calc_gb_chainrule_sse(int natoms, t_nblist *nl, float *dadx, float *dvda, 
							float *x, float *f, float *fshift, float *shiftvec,
							int gb_algorithm, gmx_genborn_t *born, t_mdatoms *md)						
{
	int    i,k,n,ii,jnr,ii3,is3,nj0,nj1,offset,n0,n1;
	int	   jnrA,jnrB,jnrC,jnrD;
    int    j3A,j3B,j3C,j3D;
	int	   jnrE,jnrF,jnrG,jnrH;
    int    j3E,j3F,j3G,j3H;
	int *  jjnr;
    
	float   rbi,shX,shY,shZ;
	float   *rb;
    
	__m128 ix,iy,iz;
	__m128 jx,jy,jz;
	__m128 jxB,jyB,jzB;
	__m128 fix,fiy,fiz;
	__m128 dx,dy,dz;
    __m128 tx,ty,tz;
	__m128 dxB,dyB,dzB;
    __m128 txB,tyB,tzB;

	__m128 rbai,rbaj,rbajB, f_gb, f_gb_ai,f_gbB,f_gb_aiB;
	__m128 xmm1,xmm2,xmm3;
	
	const __m128 two = {2.0f , 2.0f , 2.0f , 2.0f };
    
	rb     = born->work; 
			
    jjnr   = nl->jjnr;
    
	/* Loop to get the proper form for the Born radius term, sse style */
	offset=natoms%4;
	
    n0 = md->start;
    n1 = md->start+md->homenr+1+natoms/2;
    
	if(gb_algorithm==egbSTILL) 
	{
		for(i=n0;i<n1;i++)
		{
            k = i % natoms;
			rbi   = born->bRad[k];
			rb[k] = (2 * rbi * rbi * dvda[k])/ONE_4PI_EPS0;
		}
	}
	else if(gb_algorithm==egbHCT) 
	{
		for(i=n0;i<n1;i++)
		{
            k = i % natoms;
			rbi   = born->bRad[k];
			rb[k] = rbi * rbi * dvda[k];
		}
	}
	else if(gb_algorithm==egbOBC) 
	{
		for(i=n0;i<n1;i++)
		{
            k = i % natoms;
			rbi   = born->bRad[k];
			rb[k] = rbi * rbi * born->drobc[k] * dvda[k];
		}
	}
    
    jz = _mm_setzero_ps();
    
    n = j3A = j3B = j3C = j3D = 0;
    
	for(i=0;i<nl->nri;i++)
	{
        ii     = nl->iinr[i];
		ii3	   = ii*3;
        is3    = 3*nl->shift[i];     
        shX    = shiftvec[is3];  
        shY    = shiftvec[is3+1];
        shZ    = shiftvec[is3+2];
        nj0    = nl->jindex[i];      
        nj1    = nl->jindex[i+1];    

        ix     = _mm_set1_ps(shX+x[ii3+0]);
		iy     = _mm_set1_ps(shY+x[ii3+1]);
		iz     = _mm_set1_ps(shZ+x[ii3+2]);
		
		offset = (nj1-nj0)%4;
		
		rbai   = _mm_load1_ps(rb+ii);			
		fix    = _mm_setzero_ps();
		fiy    = _mm_setzero_ps();
		fiz    = _mm_setzero_ps();	
				

        for(k=nj0;k<nj1-offset;k+=4)
		{
			jnrA        = jjnr[k];   
			jnrB        = jjnr[k+1];
			jnrC        = jjnr[k+2];
			jnrD        = jjnr[k+3];
            
            j3A         = 3*jnrA;  
			j3B         = 3*jnrB;
			j3C         = 3*jnrC;
			j3D         = 3*jnrD;
            
            GMX_MM_LOAD_4JCOORD_1ATOM_PS(x+j3A,x+j3B,x+j3C,x+j3D,jx,jy,jz);
            
			dx          = _mm_sub_ps(ix,jx);
			dy          = _mm_sub_ps(iy,jy);
			dz          = _mm_sub_ps(iz,jz);
            
            GMX_MM_LOAD_4VALUES_PS(rb+jnrA,rb+jnrB,rb+jnrC,rb+jnrD,rbaj);
            
			/* load chain rule terms for j1-4 */
			f_gb        = _mm_load_ps(dadx);
			dadx += 4;
			f_gb_ai     = _mm_load_ps(dadx);
			dadx += 4;
			
            /* calculate scalar force */
            f_gb    = _mm_mul_ps(f_gb,rbai); 
            f_gb_ai = _mm_mul_ps(f_gb_ai,rbaj);
            f_gb    = _mm_add_ps(f_gb,f_gb_ai);
            
            tx     = _mm_mul_ps(f_gb,dx);
            ty     = _mm_mul_ps(f_gb,dy);
            tz     = _mm_mul_ps(f_gb,dz);
            
            fix    = _mm_add_ps(fix,tx);
            fiy    = _mm_add_ps(fiy,ty);
            fiz    = _mm_add_ps(fiz,tz);
            
            GMX_MM_DECREMENT_4JCOORD_1ATOM_PS(f+j3A,f+j3B,f+j3C,f+j3D,tx,ty,tz);
		}
        
		/*deal with odd elements */
		if(offset!=0) 
        {
            if(offset==1)
            {
                jnrA        = jjnr[k];   
                j3A         = 3*jnrA; 
                GMX_MM_LOAD_1JCOORD_1ATOM_PS(x+j3A,jx,jy,jz);
                GMX_MM_LOAD_1VALUE_PS(rb+jnrA,rbaj);
            } 
            else if(offset==2)
            {
                jnrA        = jjnr[k];   
                jnrB        = jjnr[k+1];
                j3A         = 3*jnrA;  
                j3B         = 3*jnrB;
                GMX_MM_LOAD_2JCOORD_1ATOM_PS(x+j3A,x+j3B,jx,jy,jz);
                GMX_MM_LOAD_2VALUES_PS(rb+jnrA,rb+jnrB,rbaj);
            }
            else
            {
                /* offset must be 3 */   
                jnrA        = jjnr[k];   
                jnrB        = jjnr[k+1];
                jnrC        = jjnr[k+2];
                j3A         = 3*jnrA;  
                j3B         = 3*jnrB;
                j3C         = 3*jnrC;
                GMX_MM_LOAD_3JCOORD_1ATOM_PS(x+j3A,x+j3B,x+j3C,jx,jy,jz);
                GMX_MM_LOAD_3VALUES_PS(rb+jnrA,rb+jnrB,rb+jnrC,rbaj);
            }
            
            dx          = _mm_sub_ps(ix,jx);
            dy          = _mm_sub_ps(iy,jy);
            dz          = _mm_sub_ps(iz,jz);
            
            /* load chain rule terms for j1-4 */
            f_gb        = _mm_load_ps(dadx);
            dadx += 4;
            f_gb_ai     = _mm_load_ps(dadx);
            dadx += 4;
            
            /* calculate scalar force */
            f_gb    = _mm_mul_ps(f_gb,rbai); 
            f_gb_ai = _mm_mul_ps(f_gb_ai,rbaj);
            f_gb    = _mm_add_ps(f_gb,f_gb_ai);
            
            tx     = _mm_mul_ps(f_gb,dx);
            ty     = _mm_mul_ps(f_gb,dy);
            tz     = _mm_mul_ps(f_gb,dz);
            
            fix    = _mm_add_ps(fix,tx);
            fiy    = _mm_add_ps(fiy,ty);
            fiz    = _mm_add_ps(fiz,tz);
            
            if(offset==1)
            {
                GMX_MM_DECREMENT_1JCOORD_1ATOM_PS(f+j3A,tx,ty,tz);
            } 
            else if(offset==2)
            {
                GMX_MM_DECREMENT_2JCOORD_1ATOM_PS(f+j3A,f+j3B,tx,ty,tz);
            }
            else
            {
                /* offset must be 3 */
                GMX_MM_DECREMENT_3JCOORD_1ATOM_PS(f+j3A,f+j3B,f+j3C,tx,ty,tz);
            }
        } 
        
		/* fix/fiy/fiz now contain four partial force terms, that all should be
         * added to the i particle forces and shift forces. 
         */
 		gmx_mm_update_iforce_1atom_ps(fix,fiy,fiz,f+ii3,fshift+is3);
	}	
    
	return 0;	
}


float gb_bonds_analytic(real *x, real *f, real *charge, real *bRad, real *dvda, 
					   t_idef *idef, real gb_epsilon_solvent, real facel)
{
	
	int i,j,nral,nri;
	int type1,type2,type3,type4,ai1,ai2,ai3,ai4,aj1,aj2,aj3,aj4;
	int ai13,ai23,ai33,ai43,aj13,aj23,aj33,aj43;
	int offset;
	
	float vctot;
	
	__m128 ix,iy,iz,jx,jy,jz,dx,dy,dz;
	__m128 rsq11,rinv,r,t1,t2,t3;
	__m128 isai,isaj,isaprod,inv_isaprod,expterm;
	__m128 ci,cj,qq,fac,dva,dva_i,dva_j,di,dj;
	__m128 vgb,vgb_tot,fgb,fgb2,fijC,fscal;
	__m128 tx,ty,tz,fix1,fiy1,fiz1;
	__m128 xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8;
	__m128 mask;
	
	const __m128 neg   = {-1.0f , -1.0f , -1.0f , -1.0f };
	const __m128 qrtr  = {0.25f , 0.25f , 0.25f , 0.25f };
	const __m128 eigth = {0.125f , 0.125f , 0.125f , 0.125f };
	const __m128 half  = {0.5f , 0.5f , 0.5f , 0.5f };
	const __m128 three = {3.0f , 3.0f , 3.0f , 3.0f };
	
	t_iatom *forceatoms;
	
	/* Keep the compiler happy */
	vgb_tot = _mm_setzero_ps();
	vgb     = _mm_setzero_ps();
	xmm1    = _mm_setzero_ps();
	xmm2    = _mm_setzero_ps();
	xmm3    = _mm_setzero_ps();
	xmm4    = _mm_setzero_ps();
	ai1  = ai2  = ai3  = ai4  = 0;
	aj1  = aj2  = aj3  = aj4  = 0;
	ai13 = ai23 = ai33 = ai43 = 0;
	aj13 = aj23 = aj33 = aj43 = 0;
	
	/* Scale the electrostatics by gb_epsilon_solvent */
	facel = (-1.0) * facel * (1.0 - 1.0/gb_epsilon_solvent);
	fac        = _mm_load1_ps(&facel);
	
	vctot = 0.0;
	
	for(j=F_GB12;j<=F_GB14;j++)
	{
		forceatoms = idef->il[j].iatoms;
		
		/* Number of atoms involved in each GB interaction plus the interaction type */
		nral       = NRAL(j)+1; 
		nri        = idef->il[j].nr/nral;
		offset     = nri%4;
		
		/* In the for loop, i is updated for every element in the forceatom list, so we have to stop
		 * the loop at offset*nral before the maximum idef->il[j].nr number, instead of just offset
		 */
		for(i=0;i<idef->il[j].nr-(offset*nral); )
		{
			/* Load everything separately for now, test with sse load and shuffles later
			 * Also, to avoid reading in the interaction type, we just increment i to pass over
			 * the types in the forceatoms array, this saves some memory accesses.
			 */
			i++;
			ai1            = forceatoms[i++];
			aj1            = forceatoms[i++];
		
			i++;
			ai2            = forceatoms[i++];
			aj2            = forceatoms[i++];
		
			i++;
			ai3            = forceatoms[i++];
			aj3            = forceatoms[i++];
		
			i++;
			ai4            = forceatoms[i++];
			aj4            = forceatoms[i++];
			
			ai13 = ai1*3;
			ai23 = ai2*3;
			ai33 = ai3*3;
			ai43 = ai4*3;
			
			aj13 = aj1*3;
			aj23 = aj2*3;
			aj33 = aj3*3;
			aj43 = aj4*3;
			
			/* Load particle ai1-4 and transpose */
			xmm1 = _mm_loadh_pi(xmm1,(__m64 *) (x+ai13));
			xmm2 = _mm_loadh_pi(xmm2,(__m64 *) (x+ai23));
			xmm3 = _mm_loadh_pi(xmm3,(__m64 *) (x+ai33));
			xmm4 = _mm_loadh_pi(xmm4,(__m64 *) (x+ai43));
			
			xmm5    = _mm_load1_ps(x+ai13+2);  
			xmm6    = _mm_load1_ps(x+ai23+2); 
			xmm7    = _mm_load1_ps(x+ai33+2); 
			xmm8    = _mm_load1_ps(x+ai43+2);
			
			xmm5    = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(0,0,0,0));
			xmm6    = _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(0,0,0,0));
			iz      = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(2,0,2,0));
			
			xmm1    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,2,3,2));
			xmm2    = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(3,2,3,2));
			ix      = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(2,0,2,0));
			iy      = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,1,3,1));
		
			/* Load particle aj1-4 and transpose */
			xmm1 = _mm_loadh_pi(xmm1,(__m64 *) (x+aj13));
			xmm2 = _mm_loadh_pi(xmm2,(__m64 *) (x+aj23));
			xmm3 = _mm_loadh_pi(xmm3,(__m64 *) (x+aj33));
			xmm4 = _mm_loadh_pi(xmm4,(__m64 *) (x+aj43));
			
			xmm5    = _mm_load1_ps(x+aj13+2);  
			xmm6    = _mm_load1_ps(x+aj23+2); 
			xmm7    = _mm_load1_ps(x+aj33+2); 
			xmm8    = _mm_load1_ps(x+aj43+2);
			
			xmm5    = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(0,0,0,0));
			xmm6    = _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(0,0,0,0));
			jz      = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(2,0,2,0));
			
			xmm1    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,2,3,2));
			xmm2    = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(3,2,3,2));
			jx      = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(2,0,2,0));
			jy      = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,1,3,1));
			
			/* Distances */
			dx    = _mm_sub_ps(ix, jx);
			dy    = _mm_sub_ps(iy, jy);
			dz    = _mm_sub_ps(iz, jz);
			
			rsq11 = _mm_add_ps( _mm_add_ps( _mm_mul_ps(dx,dx) , _mm_mul_ps(dy,dy) ) , _mm_mul_ps(dz,dz) );
			
			rinv  = gmx_mm_inv_ps(rsq11);
			r     = _mm_mul_ps(rsq11,rinv);
		
			/* Load Born radii for ai's and aj's */
			xmm1 = _mm_load_ss(bRad+ai1); 
			xmm2 = _mm_load_ss(bRad+ai2); 
			xmm3 = _mm_load_ss(bRad+ai3); 
			xmm4 = _mm_load_ss(bRad+ai4);
			
			xmm1 = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
			xmm3 = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(0,0,0,0)); 
			isai  = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0)); 
			
			xmm1 = _mm_load_ss(bRad+aj1); 
			xmm2 = _mm_load_ss(bRad+aj2); 
			xmm3 = _mm_load_ss(bRad+aj3); 
			xmm4 = _mm_load_ss(bRad+aj4);
			
			xmm1 = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
			xmm3 = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(0,0,0,0)); 
			isaj  = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
			
			isaprod = _mm_mul_ps(isai,isaj); // rb2 in tinker
			inv_isaprod = _mm_mul_ps(isaprod,isaprod);
			inv_isaprod = gmx_mm_inv_ps(inv_isaprod); //1/rb2 in tinker
			
			/* Load charges for ai's and aj's */
			xmm1 = _mm_load_ss(charge+ai1); 
			xmm2 = _mm_load_ss(charge+ai2); 
			xmm3 = _mm_load_ss(charge+ai3); 
			xmm4 = _mm_load_ss(charge+ai4);
			
			xmm1 = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
			xmm3 = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(0,0,0,0)); 
			ci   = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
			
			xmm1 = _mm_load_ss(charge+aj1); 
			xmm2 = _mm_load_ss(charge+aj2); 
			xmm3 = _mm_load_ss(charge+aj3); 
			xmm4 = _mm_load_ss(charge+aj4);
			
			xmm1 = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
			xmm3 = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(0,0,0,0)); 
			cj   = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
			
			qq   = _mm_mul_ps(ci,fac);
			qq   = _mm_mul_ps(cj,qq);
	 
			/* Calculate scalar GB force */
			expterm = _mm_mul_ps(rsq11,inv_isaprod);
			expterm = _mm_mul_ps(expterm,qrtr);
			expterm = _mm_mul_ps(expterm,neg);
			expterm = gmx_mm_exp_ps(expterm);
	
			fgb2    = _mm_mul_ps(isaprod,expterm);
			fgb2    = _mm_add_ps(fgb2,rsq11);
			
			fgb     = gmx_mm_inv_ps(fgb2);
			
			/* Potential energy */
			vgb     = _mm_mul_ps(qq,fgb);
			vgb_tot = _mm_add_ps(vgb_tot,vgb);
			
			fijC    = _mm_mul_ps(qrtr,r);
			fijC   = _mm_mul_ps(fijC,expterm);
			fijC   = _mm_sub_ps(r,fijC);
			fijC   = _mm_mul_ps(fijC,vgb);
			fijC   = _mm_mul_ps(fijC,neg);
			fijC   = _mm_mul_ps(fijC,fgb);
			fijC   = _mm_mul_ps(fijC,fgb); /* fijC = fijC in tab code, de */
			
			/* Chain rule terms */
			dva     = _mm_mul_ps(rsq11,inv_isaprod);
			dva     = _mm_mul_ps(dva,eigth);
			dva     = _mm_add_ps(dva,half);
			dva     = _mm_mul_ps(dva,expterm);
			dva     = _mm_mul_ps(dva,vgb);
			dva     = _mm_mul_ps(dva,neg);
			dva     = _mm_mul_ps(dva,fgb);
			dva     = _mm_mul_ps(dva,fgb); /* dva * isaj = dvatmp*isai*isai */
			
			/* Calculate vectorial force */
			fscal   = _mm_mul_ps(fijC,rinv);
			fscal   = _mm_mul_ps(fscal,neg);
			
			tx      = _mm_mul_ps(fscal,dx);
			ty      = _mm_mul_ps(fscal,dy);
			tz      = _mm_mul_ps(fscal,dz);
			
			/* Load, update and store chain rule terms for ai and aj atoms */
			dva_i   = _mm_mul_ps(dva,isaj);			
			dva_j   = _mm_mul_ps(dva,isai);
						
			xmm1 = _mm_load_ss(dvda+ai1); 
			xmm2 = _mm_load_ss(dvda+aj1);
			xmm1 = _mm_add_ss(xmm1,dva_i);
			xmm2 = _mm_add_ss(xmm2,dva_j);
			_mm_store_ss(dvda+ai1,xmm1);
			_mm_store_ss(dvda+aj1,xmm2);
		
			/* Need to shuffle dva_i and dva_j */
			dva_i    = _mm_shuffle_ps(dva_i,dva_i,_MM_SHUFFLE(0,3,2,1));
			dva_j    = _mm_shuffle_ps(dva_j,dva_j,_MM_SHUFFLE(0,3,2,1));
			
			xmm1 = _mm_load_ss(dvda+ai2); 
			xmm2 = _mm_load_ss(dvda+aj2);
			xmm1 = _mm_add_ss(xmm1,dva_i);
			xmm2 = _mm_add_ss(xmm2,dva_j);
			_mm_store_ss(dvda+ai2,xmm1);
			_mm_store_ss(dvda+aj2,xmm2);
			
			/* Need to shuffle dva_i and dva_j */
			dva_i    = _mm_shuffle_ps(dva_i,dva_i,_MM_SHUFFLE(0,3,2,1));
			dva_j    = _mm_shuffle_ps(dva_j,dva_j,_MM_SHUFFLE(0,3,2,1));
			
			xmm1 = _mm_load_ss(dvda+ai3); 
			xmm2 = _mm_load_ss(dvda+aj3);
			xmm1 = _mm_add_ss(xmm1,dva_i);
			xmm2 = _mm_add_ss(xmm2,dva_j);
			_mm_store_ss(dvda+ai3,xmm1);
			_mm_store_ss(dvda+aj3,xmm2);
		
			/* Need to shuffle dva_i and dva_j */
			dva_i    = _mm_shuffle_ps(dva_i,dva_i,_MM_SHUFFLE(0,3,2,1));
			dva_j    = _mm_shuffle_ps(dva_j,dva_j,_MM_SHUFFLE(0,3,2,1));
			
			xmm1 = _mm_load_ss(dvda+ai4); 
			xmm2 = _mm_load_ss(dvda+aj4);
			xmm1 = _mm_add_ss(xmm1,dva_i);
			xmm2 = _mm_add_ss(xmm2,dva_j);
			_mm_store_ss(dvda+ai4,xmm1);
			_mm_store_ss(dvda+aj4,xmm2);
		
			/* Load, update and store partial forces on ai and ai atoms 
			 * This needs to be done in the order ai1+aj1,ai2+aj2 etc...
			 * since loading all four values at once will mean the force
			 * is not updated properly
			 */
			/* ai1 and aj1 */
			xmm1 = _mm_load_ss(f+ai13);
			xmm2 = _mm_load_ss(f+ai13+1);
			xmm3 = _mm_load_ss(f+ai13+2);
			
			xmm4 = _mm_load_ss(f+aj13);
			xmm5 = _mm_load_ss(f+aj13+1);
			xmm6 = _mm_load_ss(f+aj13+2);
			
			xmm1 = _mm_add_ss(xmm1,tx);
			xmm2 = _mm_add_ss(xmm2,ty);
			xmm3 = _mm_add_ss(xmm3,tz);
			
			xmm4 = _mm_sub_ss(xmm4,tx);
			xmm5 = _mm_sub_ss(xmm5,ty);
			xmm6 = _mm_sub_ss(xmm6,tz);
			
			_mm_store_ss(f+ai13,xmm1);
			_mm_store_ss(f+ai13+1,xmm2);
			_mm_store_ss(f+ai13+2,xmm3);
			
			_mm_store_ss(f+aj13,xmm4);
			_mm_store_ss(f+aj13+1,xmm5);
			_mm_store_ss(f+aj13+2,xmm6);
			
			/* ai2 and aj2 */
			xmm1 = _mm_load_ss(f+ai23);
			xmm2 = _mm_load_ss(f+ai23+1);
			xmm3 = _mm_load_ss(f+ai23+2);
			
			xmm4 = _mm_load_ss(f+aj23);
			xmm5 = _mm_load_ss(f+aj23+1);
			xmm6 = _mm_load_ss(f+aj23+2);
			
			/* Need to shuffle tx/ty/tz */
			tx    = _mm_shuffle_ps(tx,tx,_MM_SHUFFLE(0,3,2,1));
			ty    = _mm_shuffle_ps(ty,ty,_MM_SHUFFLE(0,3,2,1));
			tz    = _mm_shuffle_ps(tz,tz,_MM_SHUFFLE(0,3,2,1));
			
			xmm1 = _mm_add_ss(xmm1,tx);
			xmm2 = _mm_add_ss(xmm2,ty);
			xmm3 = _mm_add_ss(xmm3,tz);
			
			xmm4 = _mm_sub_ss(xmm4,tx);
			xmm5 = _mm_sub_ss(xmm5,ty);
			xmm6 = _mm_sub_ss(xmm6,tz);
			
			_mm_store_ss(f+ai23,xmm1);
			_mm_store_ss(f+ai23+1,xmm2);
			_mm_store_ss(f+ai23+2,xmm3);
			
			_mm_store_ss(f+aj23,xmm4);
			_mm_store_ss(f+aj23+1,xmm5);
			_mm_store_ss(f+aj23+2,xmm6);
			
			/* ai3 and aj3 */
			xmm1 = _mm_load_ss(f+ai33);
			xmm2 = _mm_load_ss(f+ai33+1);
			xmm3 = _mm_load_ss(f+ai33+2);
			
			xmm4 = _mm_load_ss(f+aj33);
			xmm5 = _mm_load_ss(f+aj33+1);
			xmm6 = _mm_load_ss(f+aj33+2);
			
			/* Need to shuffle tx/ty/tz */
			tx    = _mm_shuffle_ps(tx,tx,_MM_SHUFFLE(0,3,2,1));
			ty    = _mm_shuffle_ps(ty,ty,_MM_SHUFFLE(0,3,2,1));
			tz    = _mm_shuffle_ps(tz,tz,_MM_SHUFFLE(0,3,2,1));
			
			xmm1 = _mm_add_ss(xmm1,tx);
			xmm2 = _mm_add_ss(xmm2,ty);
			xmm3 = _mm_add_ss(xmm3,tz);
			
			xmm4 = _mm_sub_ss(xmm4,tx);
			xmm5 = _mm_sub_ss(xmm5,ty);
			xmm6 = _mm_sub_ss(xmm6,tz);
			
			_mm_store_ss(f+ai33,xmm1);
			_mm_store_ss(f+ai33+1,xmm2);
			_mm_store_ss(f+ai33+2,xmm3);
			
			_mm_store_ss(f+aj33,xmm4);
			_mm_store_ss(f+aj33+1,xmm5);
			_mm_store_ss(f+aj33+2,xmm6);
			
			/* ai4 and aj4 */
			xmm1 = _mm_load_ss(f+ai43);
			xmm2 = _mm_load_ss(f+ai43+1);
			xmm3 = _mm_load_ss(f+ai43+2);
			
			xmm4 = _mm_load_ss(f+aj43);
			xmm5 = _mm_load_ss(f+aj43+1);
			xmm6 = _mm_load_ss(f+aj43+2);
			
			/* Need to shuffle tx/ty/tz */
			tx    = _mm_shuffle_ps(tx,tx,_MM_SHUFFLE(0,3,2,1));
			ty    = _mm_shuffle_ps(ty,ty,_MM_SHUFFLE(0,3,2,1));
			tz    = _mm_shuffle_ps(tz,tz,_MM_SHUFFLE(0,3,2,1));
			
			xmm1 = _mm_add_ss(xmm1,tx);
			xmm2 = _mm_add_ss(xmm2,ty);
			xmm3 = _mm_add_ss(xmm3,tz);
			
			xmm4 = _mm_sub_ss(xmm4,tx);
			xmm5 = _mm_sub_ss(xmm5,ty);
			xmm6 = _mm_sub_ss(xmm6,tz);
			
			_mm_store_ss(f+ai43,xmm1);
			_mm_store_ss(f+ai43+1,xmm2);
			_mm_store_ss(f+ai43+2,xmm3);
			
			_mm_store_ss(f+aj43,xmm4);
			_mm_store_ss(f+aj43+1,xmm5);
			_mm_store_ss(f+aj43+2,xmm6);
		}	
		
		if(offset!=0)
		{
			if(offset==1)
			{
				i++;
				ai1            = forceatoms[i++];
				aj1            = forceatoms[i++];
				
				ai13           = ai1*3;
				aj13           = aj1*3;
				
				/* Load ai and aj coordinates */
				xmm1 = _mm_loadl_pi(xmm1,(__m64 *) (x+ai13)); 
				iz   = _mm_load1_ps(x+ai13+2); 
				
				ix   = _mm_shuffle_ps(xmm1, xmm1, _MM_SHUFFLE(0,0,0,0)); 
				iy   = _mm_shuffle_ps(xmm1, xmm1, _MM_SHUFFLE(1,1,1,1)); 
				
				xmm1 = _mm_loadl_pi(xmm1,(__m64 *) (x+aj13)); 
				jz   = _mm_load1_ps(x+aj13+2); 
				
				jx   = _mm_shuffle_ps(xmm1, xmm1, _MM_SHUFFLE(0,0,0,0)); 
				jy   = _mm_shuffle_ps(xmm1, xmm1, _MM_SHUFFLE(1,1,1,1)); 
				
				/* Load Born radii for ai1 and aj1 */
				isai   = _mm_set_ps(0.0f, 0.0f, 0.0f, bRad[ai1]);	
				isaj   = _mm_set_ps(0.0f, 0.0f, 0.0f, bRad[aj1]);
				
				/* Load charges for ai1 and aj1 */
				ci     = _mm_set_ps(0.0f, 0.0f, 0.0f, charge[ai1]);
				cj     = _mm_set_ps(0.0f, 0.0f, 0.0f, charge[aj1]);
				
				/* Load chain rule terms for ai1 and aj1 */
				di     = _mm_set_ps(0.0f,0.0f,0.0f,dvda[ai1]);
				dj     = _mm_set_ps(0.0f,0.0f,0.0f,dvda[aj1]);
				
				mask = gmx_castsi128_ps( _mm_set_epi32(0,0,0,0xffffffff) );
			}
			else if(offset==2)
			{
				i++;
				ai1            = forceatoms[i++];
				aj1            = forceatoms[i++];
				
				i++;
				ai2            = forceatoms[i++];
				aj2            = forceatoms[i++];
				
				ai13 = ai1 * 3; 
				ai23 = ai2 * 3;
				aj13 = aj1 * 3; 
				aj23 = aj2 * 3;
				
				/* Load ai1-2 and aj1-2 coordinates */
				xmm1 = _mm_loadh_pi(xmm1,(__m64 *) (x+ai13)); 
				xmm2 = _mm_loadh_pi(xmm2,(__m64 *) (x+ai23)); 
				
				xmm5   = _mm_load1_ps(x+ai13+2); 
				xmm6   = _mm_load1_ps(x+ai23+2); 
				
				xmm5   = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(0,0,0,0)); 
				iz     = _mm_shuffle_ps(xmm5,xmm5,_MM_SHUFFLE(2,0,2,0)); 
				
				xmm1   = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,2,3,2)); 
				ix     = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0)); 
				iy     = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(3,1,3,1)); 	
				
				xmm1 = _mm_loadh_pi(xmm1,(__m64 *) (x+aj13)); 
				xmm2 = _mm_loadh_pi(xmm2,(__m64 *) (x+aj23)); 
				
				xmm5   = _mm_load1_ps(x+aj13+2); 
				xmm6   = _mm_load1_ps(x+aj23+2); 
				
				xmm5   = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(0,0,0,0)); 
				jz     = _mm_shuffle_ps(xmm5,xmm5,_MM_SHUFFLE(2,0,2,0)); 
				
				xmm1   = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,2,3,2)); 
				jx     = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0)); 
				jy     = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(3,1,3,1)); 	
				
				/* Load Born radii for ai1-2 and aj1-2 */
				isai   = _mm_set_ps(0.0f, 0.0f, bRad[ai2], bRad[ai1]);	
				isaj   = _mm_set_ps(0.0f, 0.0f, bRad[aj2], bRad[aj1]);
				
				/* Load charges for ai1-2 and aj1-2 */
				ci   = _mm_set_ps(0.0f, 0.0f, charge[ai2], charge[ai1]);	
				cj   = _mm_set_ps(0.0f, 0.0f, charge[aj2], charge[aj1]);
				
				/* Load chain rule terms for ai1-2 and aj1-2 */
				di     = _mm_set_ps(0.0f,0.0f,dvda[ai2],dvda[ai1]);
				dj     = _mm_set_ps(0.0f,0.0f,dvda[aj2],dvda[aj1]);
				
				mask   = gmx_castsi128_ps( _mm_set_epi32(0,0,0xffffffff,0xffffffff) );
			}
			else
			{
				i++;
				ai1            = forceatoms[i++];
				aj1            = forceatoms[i++];
				
				i++;
				ai2            = forceatoms[i++];
				aj2            = forceatoms[i++];
				
				i++;
				ai3            = forceatoms[i++];
				aj3            = forceatoms[i++];
				
				ai13 = ai1 * 3; 
				ai23 = ai2 * 3;
				ai33 = ai3 * 3;

				aj13 = aj1 * 3; 
				aj23 = aj2 * 3;
				aj33 = aj3 * 3;
				
				/* Load ai1-3 and aj1-3 coordinates */
				xmm1 = _mm_loadh_pi(xmm1,(__m64 *) (x+ai13)); 
				xmm2 = _mm_loadh_pi(xmm2,(__m64 *) (x+ai23)); 
				xmm3 = _mm_loadh_pi(xmm3,(__m64 *) (x+ai33)); 
				
				xmm5 = _mm_load1_ps(x+ai13+2); 
				xmm6 = _mm_load1_ps(x+ai23+2); 
				xmm7 = _mm_load1_ps(x+ai33+2); 
				
				xmm5 = _mm_shuffle_ps(xmm5,xmm6, _MM_SHUFFLE(0,0,0,0)); 
				iz   = _mm_shuffle_ps(xmm5,xmm7, _MM_SHUFFLE(3,1,3,1)); 						
				
				xmm1 = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(3,2,3,2)); 
				xmm2 = _mm_shuffle_ps(xmm3,xmm3, _MM_SHUFFLE(3,2,3,2)); 
				
				ix   = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(2,0,2,0)); 
				iy   = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(3,1,3,1)); 
				
				xmm1 = _mm_loadh_pi(xmm1,(__m64 *) (x+aj13)); 
				xmm2 = _mm_loadh_pi(xmm2,(__m64 *) (x+aj23)); 
				xmm3 = _mm_loadh_pi(xmm3,(__m64 *) (x+aj33)); 
				
				xmm5 = _mm_load1_ps(x+aj13+2); 
				xmm6 = _mm_load1_ps(x+aj23+2); 
				xmm7 = _mm_load1_ps(x+aj33+2); 
				
				xmm5 = _mm_shuffle_ps(xmm5,xmm6, _MM_SHUFFLE(0,0,0,0)); 
				jz   = _mm_shuffle_ps(xmm5,xmm7, _MM_SHUFFLE(3,1,3,1)); 						
				
				xmm1 = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(3,2,3,2)); 
				xmm2 = _mm_shuffle_ps(xmm3,xmm3, _MM_SHUFFLE(3,2,3,2)); 
				
				jx   = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(2,0,2,0)); 
				jy   = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(3,1,3,1)); 
				
				/* Load Born radii for ai1-3 and aj1-3 */
				isai   = _mm_set_ps(0.0f, bRad[ai3], bRad[ai2], bRad[ai1]);	
				isaj   = _mm_set_ps(0.0f, bRad[aj3], bRad[aj2], bRad[aj1]);
				
				/* Load charges for ai1-3 and aj1-3 */
				ci   = _mm_set_ps(0.0f, charge[ai3], charge[ai2], charge[ai1]);	
				cj   = _mm_set_ps(0.0f, charge[aj3], charge[aj2], charge[aj1]);
				
				/* Load chain rule terms for ai1-3 and aj1-3 */
				di     = _mm_set_ps(0.0f,dvda[ai3],dvda[ai2],dvda[ai1]);
				dj     = _mm_set_ps(0.0f,dvda[aj3],dvda[aj2],dvda[aj1]);
				
				mask = gmx_castsi128_ps( _mm_set_epi32(0,0xffffffff,0xffffffff,0xffffffff) );
			}
			
			ix = _mm_and_ps( mask, ix);
			iy = _mm_and_ps( mask, iy);
			iz = _mm_and_ps( mask, iz);
			
			jx = _mm_and_ps( mask, jx);
			jy = _mm_and_ps( mask, jy);
			jz = _mm_and_ps( mask, jz);
			
			/* Distances */
			dx    = _mm_sub_ps(ix, jx);
			dy    = _mm_sub_ps(iy, jy);
			dz    = _mm_sub_ps(iz, jz);
			
			rsq11 = _mm_add_ps( _mm_add_ps( _mm_mul_ps(dx,dx) , _mm_mul_ps(dy,dy) ) , _mm_mul_ps(dz,dz) );
			
			rinv  = gmx_mm_inv_ps(rsq11);
			r     = _mm_mul_ps(rsq11,rinv);
			
			isaprod = _mm_mul_ps(isai,isaj);
			inv_isaprod = _mm_mul_ps(isaprod,isaprod);
			inv_isaprod = gmx_mm_inv_ps(inv_isaprod); 
			
			qq   = _mm_mul_ps(ci,fac);
			qq   = _mm_mul_ps(cj,qq);
			
			/* Calculate scalar GB force */
			expterm = _mm_mul_ps(rsq11,inv_isaprod);
			expterm = _mm_mul_ps(expterm,qrtr);
			expterm = _mm_mul_ps(expterm,neg);
			expterm = gmx_mm_exp_ps(expterm);
			
			fgb2    = _mm_mul_ps(isaprod,expterm);
			fgb2    = _mm_add_ps(fgb2,rsq11);
			
			fgb     = gmx_mm_inv_ps(fgb2);
			
			/* Potential energy */
			vgb     = _mm_mul_ps(qq,fgb);
			vgb     = _mm_and_ps(mask,vgb);
			vgb_tot = _mm_add_ps(vgb_tot,vgb);
			
			fijC    = _mm_mul_ps(qrtr,r);
			fijC   = _mm_mul_ps(fijC,expterm);
			fijC   = _mm_sub_ps(r,fijC);
			fijC   = _mm_mul_ps(fijC,vgb);
			fijC   = _mm_mul_ps(fijC,neg);
			fijC   = _mm_mul_ps(fijC,fgb);
			fijC   = _mm_mul_ps(fijC,fgb); /* fscal = fijC */
						
			/* Chain rule terms */
			dva     = _mm_mul_ps(rsq11,inv_isaprod);
			dva     = _mm_mul_ps(dva,eigth);
			dva     = _mm_add_ps(dva,half);
			dva     = _mm_mul_ps(dva,expterm);
			dva     = _mm_mul_ps(dva,vgb);
			dva     = _mm_mul_ps(dva,neg);
			dva     = _mm_mul_ps(dva,fgb);
			dva     = _mm_mul_ps(dva,fgb);
			
			/* Calculate vectorial force */
			fscal   = _mm_mul_ps(fijC,rinv);
			fscal   = _mm_mul_ps(fscal,neg);
			
			tx      = _mm_mul_ps(fscal,dx);
			ty      = _mm_mul_ps(fscal,dy);
			tz      = _mm_mul_ps(fscal,dz);
			
			/* Calculate chain rule terms */
			dva_i   = _mm_mul_ps(dva,isaj);
			dva_j   = _mm_mul_ps(dva,isai);
			
			if(offset==1)
			{
				/* Load, update and store chain rule terms for ai1 and aj1 */
				xmm1 = _mm_load_ss(dvda+ai1); 
				xmm2 = _mm_load_ss(dvda+aj1);
				xmm1 = _mm_add_ss(xmm1,dva_i);
				xmm2 = _mm_add_ss(xmm2,dva_j);
				_mm_store_ss(dvda+ai1,xmm1);
				_mm_store_ss(dvda+aj1,xmm2);
				
				/* Load, update and store partial forces on ai1 and aj1 */
				xmm1 = _mm_load_ss(f+ai13);
				xmm2 = _mm_load_ss(f+ai13+1);
				xmm3 = _mm_load_ss(f+ai13+2);
				
				xmm4 = _mm_load_ss(f+aj13);
				xmm5 = _mm_load_ss(f+aj13+1);
				xmm6 = _mm_load_ss(f+aj13+2);
				
				xmm1 = _mm_add_ss(xmm1,tx);
				xmm2 = _mm_add_ss(xmm2,ty);
				xmm3 = _mm_add_ss(xmm3,tz);
				
				xmm4 = _mm_sub_ss(xmm4,tx);
				xmm5 = _mm_sub_ss(xmm5,ty);
				xmm6 = _mm_sub_ss(xmm6,tz);
				
				_mm_store_ss(f+ai13,xmm1);
				_mm_store_ss(f+ai13+1,xmm2);
				_mm_store_ss(f+ai13+2,xmm3);
				
				_mm_store_ss(f+aj13,xmm4);
				_mm_store_ss(f+aj13+1,xmm5);
				_mm_store_ss(f+aj13+2,xmm6);			
			}
			else if(offset==2)
			{
				/* Load, update and store chain rule terms for ai1-2 and aj1-2 atoms */
				dva_i   = _mm_mul_ps(dva,isaj);			
				dva_j   = _mm_mul_ps(dva,isai);
				
				xmm1 = _mm_load_ss(dvda+ai1); 
				xmm2 = _mm_load_ss(dvda+aj1);
				xmm1 = _mm_add_ss(xmm1,dva_i);
				xmm2 = _mm_add_ss(xmm2,dva_j);
				_mm_store_ss(dvda+ai1,xmm1);
				_mm_store_ss(dvda+aj1,xmm2);
				
				/* Need to shuffle dva_i and dva_j */
				dva_i    = _mm_shuffle_ps(dva_i,dva_i,_MM_SHUFFLE(0,3,2,1));
				dva_j    = _mm_shuffle_ps(dva_j,dva_j,_MM_SHUFFLE(0,3,2,1));
				
				xmm1 = _mm_load_ss(dvda+ai2); 
				xmm2 = _mm_load_ss(dvda+aj2);
				xmm1 = _mm_add_ss(xmm1,dva_i);
				xmm2 = _mm_add_ss(xmm2,dva_j);
				_mm_store_ss(dvda+ai2,xmm1);
				_mm_store_ss(dvda+aj2,xmm2);
				
				/* Load, update and store partial forces on ai1-2 and aj1-2 */
				xmm1 = _mm_load_ss(f+ai13);
				xmm2 = _mm_load_ss(f+ai13+1);
				xmm3 = _mm_load_ss(f+ai13+2);
				
				xmm4 = _mm_load_ss(f+aj13);
				xmm5 = _mm_load_ss(f+aj13+1);
				xmm6 = _mm_load_ss(f+aj13+2);
				
				xmm1 = _mm_add_ss(xmm1,tx);
				xmm2 = _mm_add_ss(xmm2,ty);
				xmm3 = _mm_add_ss(xmm3,tz);
				
				xmm4 = _mm_sub_ss(xmm4,tx);
				xmm5 = _mm_sub_ss(xmm5,ty);
				xmm6 = _mm_sub_ss(xmm6,tz);
				
				_mm_store_ss(f+ai13,xmm1);
				_mm_store_ss(f+ai13+1,xmm2);
				_mm_store_ss(f+ai13+2,xmm3);
				
				_mm_store_ss(f+aj13,xmm4);
				_mm_store_ss(f+aj13+1,xmm5);
				_mm_store_ss(f+aj13+2,xmm6);
				
				/* ai2 and aj2 */
				xmm1 = _mm_load_ss(f+ai23);
				xmm2 = _mm_load_ss(f+ai23+1);
				xmm3 = _mm_load_ss(f+ai23+2);
				
				xmm4 = _mm_load_ss(f+aj23);
				xmm5 = _mm_load_ss(f+aj23+1);
				xmm6 = _mm_load_ss(f+aj23+2);
				
				/* Need to shuffle tx/ty/tz */
				tx    = _mm_shuffle_ps(tx,tx,_MM_SHUFFLE(0,3,2,1));
				ty    = _mm_shuffle_ps(ty,ty,_MM_SHUFFLE(0,3,2,1));
				tz    = _mm_shuffle_ps(tz,tz,_MM_SHUFFLE(0,3,2,1));
				
				xmm1 = _mm_add_ss(xmm1,tx);
				xmm2 = _mm_add_ss(xmm2,ty);
				xmm3 = _mm_add_ss(xmm3,tz);
				
				xmm4 = _mm_sub_ss(xmm4,tx);
				xmm5 = _mm_sub_ss(xmm5,ty);
				xmm6 = _mm_sub_ss(xmm6,tz);
				
				_mm_store_ss(f+ai23,xmm1);
				_mm_store_ss(f+ai23+1,xmm2);
				_mm_store_ss(f+ai23+2,xmm3);
				
				_mm_store_ss(f+aj23,xmm4);
				_mm_store_ss(f+aj23+1,xmm5);
				_mm_store_ss(f+aj23+2,xmm6);
			}
			else
			{
				/* Load, update and store chain rule terms for ai1-3 and aj1-3 atoms */
				dva_i   = _mm_mul_ps(dva,isaj);			
				dva_j   = _mm_mul_ps(dva,isai);
				
				xmm1 = _mm_load_ss(dvda+ai1); 
				xmm2 = _mm_load_ss(dvda+aj1);
				xmm1 = _mm_add_ss(xmm1,dva_i);
				xmm2 = _mm_add_ss(xmm2,dva_j);
				_mm_store_ss(dvda+ai1,xmm1);
				_mm_store_ss(dvda+aj1,xmm2);
				
				/* Need to shuffle dva_i and dva_j */
				dva_i    = _mm_shuffle_ps(dva_i,dva_i,_MM_SHUFFLE(0,3,2,1));
				dva_j    = _mm_shuffle_ps(dva_j,dva_j,_MM_SHUFFLE(0,3,2,1));
				
				xmm1 = _mm_load_ss(dvda+ai2); 
				xmm2 = _mm_load_ss(dvda+aj2);
				xmm1 = _mm_add_ss(xmm1,dva_i);
				xmm2 = _mm_add_ss(xmm2,dva_j);
				_mm_store_ss(dvda+ai2,xmm1);
				_mm_store_ss(dvda+aj2,xmm2);
				
				/* Need to shuffle dva_i and dva_j */
				dva_i    = _mm_shuffle_ps(dva_i,dva_i,_MM_SHUFFLE(0,3,2,1));
				dva_j    = _mm_shuffle_ps(dva_j,dva_j,_MM_SHUFFLE(0,3,2,1));
				
				xmm1 = _mm_load_ss(dvda+ai3); 
				xmm2 = _mm_load_ss(dvda+aj3);
				xmm1 = _mm_add_ss(xmm1,dva_i);
				xmm2 = _mm_add_ss(xmm2,dva_j);
				_mm_store_ss(dvda+ai3,xmm1);
				_mm_store_ss(dvda+aj3,xmm2);
				
				/* Load, update and store partial forces on ai1-3 and aj1-3 */
				xmm1 = _mm_load_ss(f+ai13);
				xmm2 = _mm_load_ss(f+ai13+1);
				xmm3 = _mm_load_ss(f+ai13+2);
				
				xmm4 = _mm_load_ss(f+aj13);
				xmm5 = _mm_load_ss(f+aj13+1);
				xmm6 = _mm_load_ss(f+aj13+2);
				
				xmm1 = _mm_add_ss(xmm1,tx);
				xmm2 = _mm_add_ss(xmm2,ty);
				xmm3 = _mm_add_ss(xmm3,tz);
				
				xmm4 = _mm_sub_ss(xmm4,tx);
				xmm5 = _mm_sub_ss(xmm5,ty);
				xmm6 = _mm_sub_ss(xmm6,tz);
				
				_mm_store_ss(f+ai13,xmm1);
				_mm_store_ss(f+ai13+1,xmm2);
				_mm_store_ss(f+ai13+2,xmm3);
				
				_mm_store_ss(f+aj13,xmm4);
				_mm_store_ss(f+aj13+1,xmm5);
				_mm_store_ss(f+aj13+2,xmm6);
				
				/* ai2 and aj2 */
				xmm1 = _mm_load_ss(f+ai23);
				xmm2 = _mm_load_ss(f+ai23+1);
				xmm3 = _mm_load_ss(f+ai23+2);
				
				xmm4 = _mm_load_ss(f+aj23);
				xmm5 = _mm_load_ss(f+aj23+1);
				xmm6 = _mm_load_ss(f+aj23+2);
				
				/* Need to shuffle tx/ty/tz */
				tx    = _mm_shuffle_ps(tx,tx,_MM_SHUFFLE(0,3,2,1));
				ty    = _mm_shuffle_ps(ty,ty,_MM_SHUFFLE(0,3,2,1));
				tz    = _mm_shuffle_ps(tz,tz,_MM_SHUFFLE(0,3,2,1));
				
				xmm1 = _mm_add_ss(xmm1,tx);
				xmm2 = _mm_add_ss(xmm2,ty);
				xmm3 = _mm_add_ss(xmm3,tz);
				
				xmm4 = _mm_sub_ss(xmm4,tx);
				xmm5 = _mm_sub_ss(xmm5,ty);
				xmm6 = _mm_sub_ss(xmm6,tz);
				
				_mm_store_ss(f+ai23,xmm1);
				_mm_store_ss(f+ai23+1,xmm2);
				_mm_store_ss(f+ai23+2,xmm3);
				
				_mm_store_ss(f+aj23,xmm4);
				_mm_store_ss(f+aj23+1,xmm5);
				_mm_store_ss(f+aj23+2,xmm6);
				
				/* ai3 and aj3 */
				xmm1 = _mm_load_ss(f+ai33);
				xmm2 = _mm_load_ss(f+ai33+1);
				xmm3 = _mm_load_ss(f+ai33+2);
				
				xmm4 = _mm_load_ss(f+aj33);
				xmm5 = _mm_load_ss(f+aj33+1);
				xmm6 = _mm_load_ss(f+aj33+2);
				
				/* Need to shuffle tx/ty/tz */
				tx    = _mm_shuffle_ps(tx,tx,_MM_SHUFFLE(0,3,2,1));
				ty    = _mm_shuffle_ps(ty,ty,_MM_SHUFFLE(0,3,2,1));
				tz    = _mm_shuffle_ps(tz,tz,_MM_SHUFFLE(0,3,2,1));
				
				xmm1 = _mm_add_ss(xmm1,tx);
				xmm2 = _mm_add_ss(xmm2,ty);
				xmm3 = _mm_add_ss(xmm3,tz);
				
				xmm4 = _mm_sub_ss(xmm4,tx);
				xmm5 = _mm_sub_ss(xmm5,ty);
				xmm6 = _mm_sub_ss(xmm6,tz);
				
				_mm_store_ss(f+ai33,xmm1);
				_mm_store_ss(f+ai33+1,xmm2);
				_mm_store_ss(f+ai33+2,xmm3);
				
				_mm_store_ss(f+aj33,xmm4);
				_mm_store_ss(f+aj33+1,xmm5);
				_mm_store_ss(f+aj33+2,xmm6);
			}
		}
	}
	
	/* End processing for potential energy */
	vgb  = _mm_movehl_ps(vgb,vgb_tot);
	vgb_tot   = _mm_add_ps(vgb_tot,vgb);
	vgb  = _mm_shuffle_ps(vgb_tot,vgb_tot,_MM_SHUFFLE(1,1,1,1));
	vgb_tot   = _mm_add_ss(vgb_tot,vgb);
	
	_mm_store_ss(&vctot,vgb_tot);
	
	return vctot;
}




#else
/* keep compiler happy */
int genborn_sse_dummy;

#endif /* SSE intrinsics available */