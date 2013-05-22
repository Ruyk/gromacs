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

#include <shmem_utils.h>

#ifdef GMX_SHMEM

void init_shmem_buf(gmx_domdec_shmem_buf_t * shmem)
{


	if (!shmem)
			gmx_fatal(FARGS, "Cannot initialise an empty shmem structure");

    shmem->int_buf  = NULL;
    shmem->real_buf = NULL;
    shmem->rvec_buf = NULL;
    shmem->byte_buf = NULL;

#ifdef SHMEM_PRE_INITIALIZE
#define SHMEM_STARTING_ELEMS 500
    /* Initial size for the buffers to avoid repeated reallocation */
    shmem->int_buf  = shmalloc(SHMEM_STARTING_ELEMS * sizeof(int));
    shmem->real_buf = shmalloc(SHMEM_STARTING_ELEMS * sizeof(real));
    shmem->rvec_buf = shmalloc(SHMEM_STARTING_ELEMS * sizeof(rvec));
    shmem->byte_buf = shmalloc(SHMEM_STARTING_ELEMS * sizeof(void *));

    shmem->int_alloc  = SHMEM_STARTING_ELEMS;
    shmem->real_alloc = SHMEM_STARTING_ELEMS;
    shmem->rvec_alloc = SHMEM_STARTING_ELEMS;
    shmem->byte_alloc = SHMEM_STARTING_ELEMS;
#else
    shmem->int_alloc  = 0;
    shmem->real_alloc = 0;
    shmem->rvec_alloc = 0;
    shmem->byte_alloc = 0;

    shmem->int_buf  =NULL;
        shmem->real_buf = NULL;
        shmem->rvec_buf = NULL;
        shmem->byte_buf = NULL;
#endif
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
	sh_sfree(shmem->rvec_buf);
	sh_sfree(shmem->byte_buf);

	shmem->int_alloc = 0;
	shmem->real_alloc = 0;
	shmem->rvec_alloc = 0;
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

int over_alloc_shmem(int n)
{
    return SHMEM_OVER_ALLOC_FAC*n + 100;
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
		if (buf == NULL)
		{
			p = shmalloc(global_max * elem_size);
		}
		else
		{
			p = shrealloc(buf, global_max * elem_size);
		}
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

void shmem_sync_pair(gmx_domdec_shmem_buf_t* shmem, int rank_s, int rank_r)
{
	static int token_a = -1;
	static int token_b = -1;

	shmem_int_p(&token_a, 1, rank_s);
	shmem_quiet();
	shmem_int_wait(&token_a, -1);

	shmem_int_p(&token_b, 1, rank_r);
	shmem_quiet();
	token_a = -1;
	shmem_int_wait(&token_b, -1);

    token_b = -1;
}

/* SHMEM MPI replacements without buffering
 * ==============================================
 */
void shmem_real_sendrecv_nobuf(gmx_domdec_shmem_buf_t* shmem, real* buf_s, int n_s,
		int rank_s, real* buf_r, int n_r, int rank_r)
{
	static int remote_size = -1;
	SHDEBUG(" SendRecv (S: %d,R: %d) using SHMEM (n_s %d, n_r %d) \n", rank_s,
			rank_r, n_s, n_r);
	// shrenew(shmem, shmem->real_buf, &(shmem->real_alloc), n_s);
	// shmem_barrier_all();
	shmem_lock(shmem, rank_s);

	shmem_int_p(&remote_size, n_r, rank_r);
	shmem_sync_pair(shmem, rank_s, rank_r);

	if (min(remote_size,n_s)) {
		// Put buf_is in rank_s
		//               T       S     Len   Pe
		shmem_float_put(buf_r, buf_s, min(remote_size,n_s), rank_s);
	}

	shmem_sync_pair(shmem, rank_s, rank_r);
//	shmem_set_post(shmem, rank_s);
//	shmem_wait_post(shmem, _my_pe());
//	/* if (n_r) {
//		SHDEBUG(" Updating reception buffer \n");
//		memcpy(buf_r, shmem->real_buf, n_r * sizeof(real));
//	} */
//	shmem_set_done(shmem, rank_r);
//	shmem_clear_post(shmem, _my_pe());
//	shmem_wait_done(shmem, _my_pe());
//	shmem_clear_done(shmem, _my_pe());
	shmem_unlock(shmem, rank_s);
}

void shmem_int_sendrecv_nobuf(gmx_domdec_shmem_buf_t* shmem, int* buf_s, int n_s,
		int rank_s, int* buf_r, int n_r, int rank_r)
{
	static int remote_size = -1;
	int data_to_send;
	SHDEBUG(" SendRecv (S: %d,R: %d) using SHMEM (n_s %d, n_r %d) \n", rank_s,
			rank_r, n_s, n_r);

	// shmem_barrier_all();
	shmem_lock(shmem, rank_s);

    shmem_int_p(&remote_size, n_r, rank_r);
    shmem_sync_pair(shmem, rank_s, rank_r);


	if (min(remote_size,n_s)) {
		// Put buf_is in rank_s
		//               T       S     Len   Pe
		shmem_int_put(buf_r, buf_s, min(n_s,remote_size), rank_s);
		shmem_quiet();
		SHDEBUG(" Data posted by %d to %d \n", _my_pe(), rank_s);
	}
	shmem_set_post(shmem, rank_s);
	SHDEBUG(" Event posted to %d \n", rank_s);
	// Receive data
	shmem_wait_post(shmem, _my_pe());
	SHDEBUG(" Waiting for event from %d \n", _my_pe());
	/* if (n_r) {
		SHDEBUG(" Updating reception buffer from %d \n", rank_r);
		memcpy(buf_r, shmem->int_buf, n_r * sizeof(int));
		SHDEBUG(" Copied %p into %p, size %ld, from %d \n", shmem->int_buf,
				buf_r, n_r * sizeof(int), rank_r);
	} */
	shmem_set_done(shmem, rank_r);
	shmem_clear_post(shmem, _my_pe());
	SHDEBUG(" Event of %d cleared \n", _my_pe());
	shmem_wait_done(shmem, _my_pe());
	shmem_clear_done(shmem, _my_pe());
	shmem_unlock(shmem, rank_s);
}

void shmem_rvec_sendrecv_nobuf(gmx_domdec_shmem_buf_t* shmem, rvec* buf_s, int n_s,
		int rank_s, rvec* buf_r, int n_r, int rank_r)
{
	static int remote_size = -1;
	shmem_lock(shmem, rank_s);

	shmem_int_p(&remote_size, n_r, rank_r);
	shmem_sync_pair(shmem, rank_s, rank_r);


	if (min(remote_size,n_s)) {
		int i = 0;
		// Put buf_is in rank_s
		//               T       S     Len   Pe
		shmem_float_put((real *) buf_r, (real *) buf_s, min(remote_size,n_s) * DIM,
				rank_s);
		SHDEBUG(" Putting n_s * srvec %ld data to rank_s %d \n",
				n_s * sizeof(rvec), rank_s);
	}
	shmem_set_post(shmem, rank_s);
	shmem_wait_post(shmem, _my_pe());
	/* if (n_r) {
		SHDEBUG(" Updating reception buffer \n");
		memcpy(buf_r, shmem->rvec_buf, n_r * sizeof(rvec));
	} */
	shmem_set_done(shmem, rank_r);
	shmem_clear_post(shmem, _my_pe());
	shmem_wait_done(shmem, _my_pe());
	shmem_clear_done(shmem, _my_pe());
	shmem_unlock(shmem, rank_s);
}


/* SHMEM MPI replacements WITH buffering
 * ==============================================
 */
void shmem_int_sendrecv(gmx_domdec_shmem_buf_t* shmem, int* buf_s, int n_s,
		int rank_s, int* buf_r, int n_r, int rank_r)
{
	SHDEBUG(" SendRecv (S: %d,R: %d) using SHMEM (n_s %d, n_r %d) \n", rank_s,
			rank_r, n_s, n_r);
	shrenew(shmem, shmem->int_buf, &(shmem->int_alloc), n_s);
	// shmem_barrier_all();
	shmem_lock(shmem, rank_s);
	if (n_s) {
		// Put buf_is in rank_s
		//               T       S     Len   Pe
		shmem_int_put(shmem->int_buf, buf_s, n_s, rank_s);
		shmem_quiet();
		SHDEBUG(" Data posted by %d to %d \n", _my_pe(), rank_s);
	}
	shmem_set_post(shmem, rank_s);
	SHDEBUG(" Event posted to %d \n", rank_s);
	// Receive data
	shmem_wait_post(shmem, _my_pe());
	SHDEBUG(" Waiting for event from %d \n", _my_pe());
	if (n_r) {
		SHDEBUG(" Updating reception buffer from %d \n", rank_r);
		memcpy(buf_r, shmem->int_buf, n_r * sizeof(int));
		SHDEBUG(" Copied %p into %p, size %ld, from %d \n", shmem->int_buf,
				buf_r, n_r * sizeof(int), rank_r);
	}
	shmem_set_done(shmem, rank_r);
	shmem_clear_post(shmem, _my_pe());
	SHDEBUG(" Event of %d cleared \n", _my_pe());
	shmem_wait_done(shmem, _my_pe());
	shmem_clear_done(shmem, _my_pe());
	shmem_unlock(shmem, rank_s);
}


void shmem_real_sendrecv(gmx_domdec_shmem_buf_t* shmem, real* buf_s, int n_s,
		int rank_s, real* buf_r, int n_r, int rank_r)
{
	SHDEBUG(" SendRecv (S: %d,R: %d) using SHMEM (n_s %d, n_r %d) \n", rank_s,
			rank_r, n_s, n_r);
	shrenew(shmem, shmem->real_buf, &(shmem->real_alloc), n_s);
	// shmem_barrier_all();
	shmem_lock(shmem, rank_s);
	if (n_s) {
		// Put buf_is in rank_s
		//               T       S     Len   Pe
		shmem_float_put(shmem->real_buf, buf_s, n_s, rank_s);
	}
	shmem_set_post(shmem, rank_s);
	shmem_wait_post(shmem, _my_pe());
	if (n_r) {
		SHDEBUG(" Updating reception buffer \n");
		memcpy(buf_r, shmem->real_buf, n_r * sizeof(real));
	}
	shmem_set_done(shmem, rank_r);
	shmem_clear_post(shmem, _my_pe());
	shmem_wait_done(shmem, _my_pe());
	shmem_clear_done(shmem, _my_pe());
	shmem_unlock(shmem, rank_s);
}


void shmem_rvec_sendrecv(gmx_domdec_shmem_buf_t* shmem, rvec* buf_s, int n_s,
		int rank_s, rvec* buf_r, int n_r, int rank_r)
{
	shrenew(shmem, shmem->rvec_buf, &(shmem->rvec_alloc), n_s * sizeof(rvec));

	shmem_lock(shmem, rank_s);
	if (n_s) {
		int i = 0;
		// Put buf_is in rank_s
		//               T       S     Len   Pe
		shmem_float_put((real *) shmem->rvec_buf, (real *) buf_s, n_s * DIM,
				rank_s);
		SHDEBUG(" Putting n_s * srvec %ld data to rank_s %d \n",
				n_s * sizeof(rvec), rank_s)
	}
	SHDEBUG(" Set flag %d \n", rank_s);
	shmem_set_post(shmem, rank_s);
	SHDEBUG(" Wait flag %d \n", _my_pe());
	shmem_wait_post(shmem, _my_pe());
	if (n_r) {
		SHDEBUG(" Updating reception buffer \n");
		memcpy(buf_r, shmem->rvec_buf, n_r * sizeof(rvec));
	}
	shmem_set_done(shmem, rank_r);
	shmem_clear_post(shmem, _my_pe());
	shmem_wait_done(shmem, _my_pe());
	shmem_clear_done(shmem, _my_pe());
	shmem_unlock(shmem, rank_s);
}



void shmem_void_sendrecv(gmx_domdec_shmem_buf_t* shmem, void* buf_s, int n_s,
		int rank_s, void* buf_r, int n_r, int rank_r)
{
	SHDEBUG(" SendRecv (S: %d,R: %d) using SHMEM (n_s %d, n_r %d) \n", rank_s,
			rank_r, n_s, n_r);
	shrenew(shmem, shmem->byte_buf, &(shmem->byte_alloc), n_s);
	// shmem_barrier_all();
	shmem_lock(shmem, rank_s);
	if (n_s) {
		// Put buf_is in rank_s
		//               T       S     Len   Pe
		shmem_float_put(shmem->byte_buf, buf_s, n_s, rank_s);
	}
	shmem_set_post(shmem, rank_s);
	shmem_wait_post(shmem, _my_pe());
	if (n_r) {
		SHDEBUG(" Updating reception buffer \n");
		memcpy(buf_r, shmem->byte_buf, n_r * sizeof(real));
	}
	shmem_set_done(shmem, rank_r);
	shmem_clear_post(shmem, _my_pe());
	shmem_wait_done(shmem, _my_pe());
	shmem_clear_done(shmem, _my_pe());
	shmem_unlock(shmem, rank_s);
}


#endif /* GMX_SHMEM */
