/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2008
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "domdec_network.h"

#ifdef GMX_SHMEM
#include <shmem.h>
#include "smalloc.h"
#include "typedefs.h"
#include "gmx_fatal.h"

#ifndef MAX_INT
#warning "Using manual MAX_INT setting"
#define MAX_INT 9999999
#endif

const int MAX_BUFF = 5000;

#define shmem_reset_flag(FLAG) {  FLAG = 0; shmem_barrier_all(); }
#define shmem_set_flag(FLAG, TARGET) { shmem_int_p(FLAG, 1, TARGET); }
#define shmem_wait_flag(FLAG) { shmem_int_wait(FLAG, 0); }
#define SHDEBUG(...) { printf("SHMEM(ID:%d) (domdec_network.c,%d)", _my_pe(), __LINE__); printf(__VA_ARGS__); }

/* get_max_alloc
 * ========================
 *   Computes the maximum value of memory requested across PEs
 */
int get_max_alloc(int local_value) {
   static int global_max = 0;
   static int local_max;
   static int pWrk[2>_SHMEM_REDUCE_MIN_WRKDATA_SIZE?2:_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
   static long pSync[_SHMEM_REDUCE_SYNC_SIZE];
   int i;
   for (i = 0; i < _SHMEM_REDUCE_SYNC_SIZE; i++)
   {
	   pSync[i] = _SHMEM_SYNC_VALUE;
   }
   local_max= local_value;
   shmem_barrier_all();
   shmem_int_max_to_all(&global_max, &local_max, 1, 0, 0, _num_pes(), pWrk, pSync);
   return global_max;
}

#define GMX_SHMEM_DEBUG
#ifndef GMX_SHMEM_DEBUG
#define shrenew(PTR, OLD_SIZE, NEW_SIZE) (PTR) = sh_renew_buf((PTR), (OLD_SIZE), (NEW_SIZE), sizeof(*(PTR)))
#else

#define shrenew(PTR, OLD_SIZE, NEW_SIZE) { SHDEBUG(" Before renew , %p size %d \n", PTR, *(OLD_SIZE)); \
 					   (PTR) = sh_renew_buf((PTR), (OLD_SIZE), (NEW_SIZE), sizeof(*(PTR))); \
    					   SHDEBUG(" After renew , %p size %d \n",  PTR, *(OLD_SIZE));\
					 }
#endif


/* renew the sh tmp buffer if required */
void * sh_renew_buf(void * buf, int * alloc, const int new_size, const int elem_size) {
	void * p;
	int global_max;
	global_max = get_max_alloc(over_alloc_dd(new_size));
   	if (global_max > (*alloc)) {
   		SHDEBUG(" Updating alloc (%d) to new global max (%d) with elem size %d \n", (*alloc), global_max, elem_size);
   		// BUGGY: sh_srenew(buf, (*alloc));
   		(*alloc) = global_max;
       	        p = shrealloc(buf, global_max * elem_size);
       	        if (!p){
        	   SHDEBUG(" shrealloc returned NULL \n", p)
        	   p = buf;
                }
               SHDEBUG(" After update to global max (%d) new buf ptr is %p \n", global_max, p);
   	} else {
   		p = buf;
   	}

   	return p;
}



#endif

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREAD_MPI
#include "tmpi.h"
#endif


#define DDMASTERRANK(dd)   (dd->masterrank)



void dd_sendrecv_int(const gmx_domdec_t *dd,
                     int ddimind, int direction,
                     int *buf_s, int n_s,
                     int *buf_r, int n_r)
{
#ifdef GMX_SHMEM
#warning " Replacing SendRecv int"
    int        rank_s, rank_r;
    gmx_domdec_shmem_buf_t * shmem = dd->shmem;
    static int sh_flag;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    SHDEBUG(" SendRecv (S: %d,R: %d) using SHMEM (n_s %d, n_r %d) \n", rank_s, rank_r, n_s, n_r);
    shrenew(shmem->int_buf, &(shmem->int_alloc), n_s);
    shmem_reset_flag(sh_flag);

    if (n_s) {
    	// Put buf_is in rank_s
    	//               T       S     Len   Pe
    	shmem_int_put(shmem->int_buf, buf_s, n_s, rank_s);
    	shmem_quiet();
    	shmem_set_flag(&sh_flag, rank_s);
    }
    SHDEBUG(" After the shmem_int_put (%d => %d)\n", _my_pe(), rank_s);
    if (n_r) {
    	shmem_wait_flag(&sh_flag);
    	SHDEBUG(" Updating reception buffer \n");
    	memcpy(buf_r, shmem->int_buf, n_r * sizeof(int));
    }

    SHDEBUG("After shfree (end of subroutine)\n");

#elif defined(GMX_MPI)
    int        rank_s, rank_r;
    MPI_Status stat;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s, n_s*sizeof(int), MPI_BYTE, rank_s, 0,
                     buf_r, n_r*sizeof(int), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s, n_s*sizeof(int), MPI_BYTE, rank_s, 0,
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
        MPI_Recv(    buf_r, n_r*sizeof(int), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }

#endif
}

void dd_sendrecv_real(const gmx_domdec_t *dd,
                      int ddimind, int direction,
                      real *buf_s, int n_s,
                      real *buf_r, int n_r)
{
#ifdef GMX_SHMEM
#warning " Replacing SendRecv real"
    int        rank_s, rank_r;
    gmx_domdec_shmem_buf_t * shmem = dd->shmem;
    static int sh_flag;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    SHDEBUG(" SendRecv (S: %d,R: %d) using SHMEM (n_s %d, n_r %d) \n", rank_s, rank_r, n_s, n_r);
    shrenew(shmem->real_buf, &(shmem->real_alloc), n_s);
    shmem_reset_flag(sh_flag);


    if (n_s) {
    	// Put buf_is in rank_s
    	//               T       S     Len   Pe
    	shmem_float_put(shmem->real_buf, buf_s, n_s, rank_s);
    	shmem_quiet();
    	shmem_set_flag(&sh_flag, rank_s);
    }
    SHDEBUG(" After the shmem_real_put (%d => %d)\n", _my_pe(), rank_s);
    if (n_r) {
        shmem_wait_flag(&sh_flag);
    	SHDEBUG(" Updating reception buffer \n");
    	memcpy(buf_r, shmem->real_buf, n_r * sizeof(real));
    }

    SHDEBUG("After shfree (end of subroutine)\n");

#elif defined(GMX_MPI)
    int        rank_s, rank_r;
    MPI_Status stat;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s, n_s*sizeof(real), MPI_BYTE, rank_s, 0,
                     buf_r, n_r*sizeof(real), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s, n_s*sizeof(real), MPI_BYTE, rank_s, 0,
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
        MPI_Recv(    buf_r, n_r*sizeof(real), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }

#endif
}

void dd_sendrecv_rvec(const gmx_domdec_t *dd,
                      int ddimind, int direction,
                      rvec *buf_s, int n_s,
                      rvec *buf_r, int n_r)
{
#ifdef GMX_SHMEM
#warning " Replacing SendRecv rvec"
    int        rank_s, rank_r;
    gmx_domdec_shmem_buf_t * shmem = dd->shmem;
    static int sh_flag;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    SHDEBUG(" SendRecv (S: %d,R: %d) using SHMEM (n_s %d, n_r %d) (size rvec: %ld) \n", rank_s, rank_r, n_s, n_r, sizeof(rvec));
    shrenew(shmem->rvec_buf, &(shmem->rvec_alloc), n_s * sizeof(rvec));
    shmem_reset_flag(sh_flag);

    if (n_s) {
    	// Put buf_is in rank_s
    	//               T       S     Len   Pe
    	shmem_float_put((real *) shmem->rvec_buf, (real *) buf_s, n_s * sizeof(rvec), rank_s);
    	shmem_quiet();
    	shmem_set_flag(&sh_flag, rank_s);
    }
    SHDEBUG(" After the shmem_real_put (RVEC) (%d => %d)\n", _my_pe(), rank_s);
    if (n_r) {
        shmem_wait_flag(&sh_flag);
    	SHDEBUG(" Updating reception buffer \n");
    	memcpy(buf_r, shmem->rvec_buf, n_r * sizeof(rvec));
    }

    SHDEBUG("After shfree (end of subroutine)\n");

#elif defined(GMX_MPI)
    int        rank_s, rank_r;
    MPI_Status stat;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s[0], n_s*sizeof(rvec), MPI_BYTE, rank_s, 0,
                     buf_r[0], n_r*sizeof(rvec), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s[0], n_s*sizeof(rvec), MPI_BYTE, rank_s, 0,
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
        MPI_Recv(    buf_r[0], n_r*sizeof(rvec), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }

#endif
}

void dd_sendrecv2_rvec(const gmx_domdec_t *dd,
                       int ddimind,
                       rvec *buf_s_fw, int n_s_fw,
                       rvec *buf_r_fw, int n_r_fw,
                       rvec *buf_s_bw, int n_s_bw,
                       rvec *buf_r_bw, int n_r_bw)
{
#ifdef GMX_MPI
    int         rank_fw, rank_bw, nreq;
    MPI_Request req[4];
    MPI_Status  stat[4];

    rank_fw = dd->neighbor[ddimind][0];
    rank_bw = dd->neighbor[ddimind][1];

    if (!dd->bSendRecv2)
    {
        /* Try to send and receive in two directions simultaneously.
         * Should be faster, especially on machines
         * with full 3D communication networks.
         * However, it could be that communication libraries are
         * optimized for MPI_Sendrecv and non-blocking MPI calls
         * are slower.
         * SendRecv2 can be turned on with the env.var. GMX_DD_SENDRECV2
         */
        nreq = 0;
        if (n_r_fw)
        {
            MPI_Irecv(buf_r_fw[0], n_r_fw*sizeof(rvec), MPI_BYTE,
                      rank_bw, 0, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_r_bw)
        {
            MPI_Irecv(buf_r_bw[0], n_r_bw*sizeof(rvec), MPI_BYTE,
                      rank_fw, 1, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_s_fw)
        {
            MPI_Isend(buf_s_fw[0], n_s_fw*sizeof(rvec), MPI_BYTE,
                      rank_fw, 0, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_s_bw)
        {
            MPI_Isend(buf_s_bw[0], n_s_bw*sizeof(rvec), MPI_BYTE,
                      rank_bw, 1, dd->mpi_comm_all, &req[nreq++]);
        }
        if (nreq)
        {
            MPI_Waitall(nreq, req, stat);
        }
    }
    else
    {
        /* Communicate in two ordered phases.
         * This is slower, even on a dual-core Opteron cluster
         * with a single full-duplex network connection per machine.
         */
        /* Forward */
        MPI_Sendrecv(buf_s_fw[0], n_s_fw*sizeof(rvec), MPI_BYTE, rank_fw, 0,
                     buf_r_fw[0], n_r_fw*sizeof(rvec), MPI_BYTE, rank_bw, 0,
                     dd->mpi_comm_all, &stat[0]);
        /* Backward */
        MPI_Sendrecv(buf_s_bw[0], n_s_bw*sizeof(rvec), MPI_BYTE, rank_bw, 0,
                     buf_r_bw[0], n_r_bw*sizeof(rvec), MPI_BYTE, rank_fw, 0,
                     dd->mpi_comm_all, &stat[0]);
    }
#endif
}

/* IBM's BlueGene(/L) MPI_Bcast dereferences the data pointer
 * even when 0 == nbytes, so we protect calls to it on BlueGene.
 * Fortunately dd_bcast() and dd_bcastc() are only
 * called during DD setup and partition.
 */

void dd_bcast(gmx_domdec_t *dd, int nbytes, void *data)
{
#ifdef GMX_SHMEM
   static long pSync[_SHMEM_BCAST_SYNC_SIZE];
   void * buf;
   int i, size;
   SHDEBUG(" Bcast of %d bytes (base ptr %p) \n", nbytes, data);
   for (i = 0; i < _SHMEM_BCAST_SYNC_SIZE; i++)
   {
   	pSync[i] = _SHMEM_SYNC_VALUE;
   }
   shmem_barrier_all();
   /* Since the dd_bcast receives a number of bytes, but
    * broadcast expects a number of elements, we need to adjust
    * the size of the temporary buffer so it is a multiple of the
    * size of a pointer to void.
    * Otherwhise pointer get corrupted.
    */
   if (nbytes%sizeof(void *)) {
	size = nbytes + sizeof(void *) - (nbytes%sizeof(void *));
   } 
   else 
   {
	size = nbytes;
   }
   buf = shmalloc(size);
   if (DDMASTERRANK(dd) == _my_pe())
   {
   	memcpy(buf, data, nbytes);
	memset(buf + nbytes, 0, size - nbytes);	
   }
   SHDEBUG("  buf ptr %p , masterrank %d  \n", buf, DDMASTERRANK(dd));
   shmem_broadcast(buf, buf, size, DDMASTERRANK(dd), 0, 0, _num_pes(), pSync); 
   SHDEBUG(" After broadcast  \n");
   if (DDMASTERRANK(dd) != _my_pe()) 
   {
	memcpy(data, buf, nbytes);
   }
   shmem_barrier_all();
   shfree(buf);
   SHDEBUG(" End of routine \n");
#elif defined(GMX_MPI)
#ifdef GMX_BLUEGENE
    if (nbytes > 0)
    {
#endif
    MPI_Bcast(data, nbytes, MPI_BYTE,
              DDMASTERRANK(dd), dd->mpi_comm_all);
#ifdef GMX_BLUEGENE
}
#endif
#endif
}

void dd_bcastc(gmx_domdec_t *dd, int nbytes, void *src, void *dest)
{
    if (DDMASTER(dd))
    {
        memcpy(dest, src, nbytes);
    }
#ifdef GMX_MPI
#ifdef GMX_BLUEGENE
    if (nbytes > 0)
    {
#endif
    MPI_Bcast(dest, nbytes, MPI_BYTE,
              DDMASTERRANK(dd), dd->mpi_comm_all);
#ifdef GMX_BLUEGENE
}
#endif
#endif
}

void dd_scatter(gmx_domdec_t *dd, int nbytes, void *src, void *dest)
{
#ifdef GMX_MPI
    MPI_Scatter(src, nbytes, MPI_BYTE,
                dest, nbytes, MPI_BYTE,
                DDMASTERRANK(dd), dd->mpi_comm_all);
#endif
}

void dd_gather(gmx_domdec_t *dd, int nbytes, void *src, void *dest)
{
#ifdef GMX_MPI
    MPI_Gather(src, nbytes, MPI_BYTE,
               dest, nbytes, MPI_BYTE,
               DDMASTERRANK(dd), dd->mpi_comm_all);
#endif
}

void dd_scatterv(gmx_domdec_t *dd,
                 int *scounts, int *disps, void *sbuf,
                 int rcount, void *rbuf)
{
#ifdef GMX_MPI
    int dum;

    if (rcount == 0)
    {
        /* MPI does not allow NULL pointers */
        rbuf = &dum;
    }
    MPI_Scatterv(sbuf, scounts, disps, MPI_BYTE,
                 rbuf, rcount, MPI_BYTE,
                 DDMASTERRANK(dd), dd->mpi_comm_all);
#endif
}

void dd_gatherv(gmx_domdec_t *dd,
                int scount, void *sbuf,
                int *rcounts, int *disps, void *rbuf)
{
#ifdef GMX_MPI
    int dum;

    if (scount == 0)
    {
        /* MPI does not allow NULL pointers */
        sbuf = &dum;
    }
    MPI_Gatherv(sbuf, scount, MPI_BYTE,
                rbuf, rcounts, disps, MPI_BYTE,
                DDMASTERRANK(dd), dd->mpi_comm_all);
#endif
}
