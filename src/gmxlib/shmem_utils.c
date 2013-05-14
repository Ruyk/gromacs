/*                                                                                                                   
 * This file is part of the GROMACS molecular simulation package.                                                    
 *                                                                                                                   
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
/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-              
 *                                                                                                                   
 *                                                                                                                   
 *                This source code is part of                                                                        
 *                                                                                                                   
 *                 G   R   O   M   A   C   S                                                                         
 *                                                                                                                   
 *          GROningen MAchine for Chemical Simulations                                                               
 *                                                                                                                   
 * Written by the Gromacs development team under coordination of                                                     
 * David van der Spoel, Berk Hess, and Erik Lindahl.                                                                 
 *                                                                                                                   
 * This library is free software; you can redistribute it and/or                                                     
 * modify it under the terms of the GNU Lesser General Public License                                                
 * as published by the Free Software Foundation; either version 2                                                    
 * of the License, or (at your option) any later version.                                                            
 *                                                                                                                   
 * To help us fund GROMACS development, we humbly ask that you cite                                                  
 * the research papers on the package. Check out http://www.gromacs.org                                              
 *                                                                                                                                                                                                 
 */    

#include "shmem_utils.h"

#ifdef GMX_SHMEM

void init_shmem_buf(gmx_domdec_shmem_buf_t * shmem)
{
	/* static int pWrk[2>_SHMEM_REDUCE_MIN_WRKDATA_SIZE?2:_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
	static long pSync[_SHMEM_REDUCE_SYNC_SIZE]; */

	if (!shmem)
			gmx_fatal(FARGS, "Cannot initialise an empty shmem structure");
    shmem->int_alloc  = 0;
    shmem->real_alloc = 0;
    shmem->rvec_alloc = 0;
    shmem->int_buf  = NULL;
    shmem->real_buf = NULL;
    shmem->rvec_buf = NULL;
    /* Init max_alloc pWrk and pSync arrays */
    sh_snew(shmem->max_alloc_pWrk1, 2>_SHMEM_REDUCE_MIN_WRKDATA_SIZE?2:_SHMEM_REDUCE_MIN_WRKDATA_SIZE);
    sh_snew(shmem->max_alloc_pWrk2, 2>_SHMEM_REDUCE_MIN_WRKDATA_SIZE?2:_SHMEM_REDUCE_MIN_WRKDATA_SIZE);
    sh_snew(shmem->max_alloc_pSync1, _SHMEM_REDUCE_SYNC_SIZE);
    sh_snew(shmem->max_alloc_pSync2, _SHMEM_REDUCE_SYNC_SIZE);

    /* Init p2p sync */
    sh_snew(shmem->post_events, _num_pes());
    shmem_clear_event(shmem->post_events + _my_pe());
    sh_snew(shmem->done_events, _num_pes());
    shmem_clear_event(shmem->done_events + _my_pe());
    /* Init locks */
    sh_snew(shmem->lock, _num_pes());
    shmem_clear_lock(shmem->lock + _my_pe());
    shmem_barrier_all();
}

void done_shmem_buf(gmx_domdec_shmem_buf_t * shmem)
{
	if (!shmem)
		gmx_fatal(FARGS, "Cannot cleanup an empty shmem structure");

	shmem_barrier_all();
	sh_sfree(shmem->int_buf);
	sh_sfree(shmem->real_buf);
	sh_sfree(shmem->byte_buf);

	shmem->int_alloc = 0;
	shmem->real_alloc = 0;
	shmem->byte_alloc = 0;

	sh_sfree(shmem->post_events);
	sh_sfree(shmem->done_events);
	sh_sfree(shmem->lock);

	sh_sfree(shmem->max_alloc_pWrk1);
	sh_sfree(shmem->max_alloc_pWrk2);
	sh_sfree(shmem->max_alloc_pSync1);
	sh_sfree(shmem->max_alloc_pSync2);
}


void shmem_set_post      (gmx_domdec_shmem_buf_t * shmem, int pe)
{
	SHDEBUG(" Setting post on %d \n", pe);
	shmem_quiet();
	shmem_set_event(shmem->post_events + pe);
}

void shmem_clear_post    (gmx_domdec_shmem_buf_t * shmem, int pe)
{
	shmem_quiet();
	SHDEBUG(" Clearing post on %d \n", pe);
	shmem_clear_event(shmem->post_events + pe);
}

void shmem_wait_post     (gmx_domdec_shmem_buf_t * shmem, int pe)
{
	shmem_quiet();
	SHDEBUG(" Waiting for post on %d \n", pe);
	shmem_wait_event(shmem->post_events + pe);
	SHDEBUG(" Posted on %d \n", pe);
}

void shmem_wait_done     (gmx_domdec_shmem_buf_t * shmem, int pe)
{
	SHDEBUG(" Waiting for done on %d \n", pe);
	shmem_quiet();
	shmem_wait_event(shmem->done_events + pe);
	SHDEBUG(" Received done on %d \n", pe);
}

void shmem_set_done     (gmx_domdec_shmem_buf_t * shmem, int pe)
{
	SHDEBUG(" Setting done on %d \n", pe);
	shmem_quiet();
	shmem_set_event(shmem->done_events + pe);
}

void shmem_clear_done  (gmx_domdec_shmem_buf_t * shmem, int pe)
{
	SHDEBUG(" Clearing done on %d \n", pe);
	shmem_clear_event(shmem->done_events + pe);
}

gmx_bool shmem_is_posted(gmx_domdec_shmem_buf_t * shmem, int pe)
{
	return shmem_test_event(shmem->post_events + pe);
}


/* Lock handling
 */


void shmem_lock(gmx_domdec_shmem_buf_t * shmem, int pe)
{
	// SHDEBUG(" Waiting for lock on %d \n", pe);
    shmem_set_lock(shmem->lock + pe);
    // SHDEBUG(" Acquired lock on %d \n", pe);
}

void shmem_unlock(gmx_domdec_shmem_buf_t * shmem, int pe)
{
	// SHDEBUG(" Leaving lock on %d \n", pe);
	shmem_clear_lock(shmem->lock + pe);
}

gmx_bool shmem_is_locked(gmx_domdec_shmem_buf_t * shmem, int pe)
{
	return shmem_test_lock(shmem->lock + pe);
}




int round_to_next_multiple(int nbytes, int type_size)
{
     int size;
	if (nbytes%type_size) {
		size = nbytes + type_size - (nbytes%type_size);
	}
	else
	{
		size = nbytes;
	}
    return size;
}

int get_max_alloc_shmem(int local_value)
{
   static int global_max = 0;
   static int local_max;
   static int pWrk[2>_SHMEM_REDUCE_MIN_WRKDATA_SIZE?2:_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
   static long pSync[_SHMEM_REDUCE_SYNC_SIZE];
   int i;
   for (i = 0; i < _SHMEM_REDUCE_SYNC_SIZE; i++)
   {
	   pSync[i] = _SHMEM_SYNC_VALUE;
   }
   local_max = local_value;
   shmem_barrier_all();
   shmem_int_max_to_all(&global_max, &local_max, 1, 0, 0, _num_pes(), pWrk, pSync);
   return global_max;
}

int get_max_alloc_shmem_dd(gmx_domdec_shmem_buf_t * shmem, int local_value)
{
	static int global_max = 0;
	static int call = 0;
	static int local_max = 0;
	int * pWrk;
	long * pSync;

	if (call%2)
	{
		pWrk = shmem->max_alloc_pWrk1;
		pSync = shmem->max_alloc_pSync1;
	}
	else
	{
		pWrk = shmem->max_alloc_pWrk2;
		pSync = shmem->max_alloc_pSync2;
	}
	local_max = local_value;
	shmem_int_max_to_all(&global_max, &local_max, 1, 0, 0, _num_pes(), pWrk, pSync);
	call++;
	return global_max;
}


void * sh_renew_buf(gmx_domdec_shmem_buf_t * shmem, void * buf, int * alloc, const int new_size, const int elem_size)
{
	void * p;
	int global_max;
	SHDEBUG(" Before get max alloc \n");
	global_max = get_max_alloc_shmem_dd(shmem, over_alloc_dd(new_size));
	SHDEBUG(" After get max alloc \n");
	if (global_max > (*alloc)) {
		SHDEBUG(" Updating alloc (%d) to new global max (%d) with elem size %d \n", (*alloc), global_max, elem_size);
		// BUGGY: sh_srenew(buf, (*alloc));
		(*alloc) = global_max;
		p = shrealloc(buf, global_max * elem_size);
		if (!p){
			SHDEBUG(" shrealloc returned NULL \n")
        			   p = buf;
		}
		SHDEBUG(" After update to global max (%d) new buf ptr is %p \n", global_max, p);
	} else {
		p = buf;
		SHDEBUG(" Not updating, global max (%d) same buf ptr is %p (alloc: %d) \n", global_max, p, global_max * elem_size);
	}

   	return p;
}


#endif /* GMX_SHMEM */
