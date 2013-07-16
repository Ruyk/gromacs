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
#include "shmem_utils.h"
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
#ifdef GMX_SHMEM_INT
#warning " Replacing SendRecv int"
    int        rank_s, rank_r;
    gmx_domdec_shmem_buf_t * shmem = dd->shmem;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

	shmem_int_sendrecv(shmem, buf_s, n_s, rank_s, buf_r, n_r, rank_r);
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

#ifdef GMX_SHMEM
void dd_sendrecv_real_off(const gmx_domdec_t *dd,
                      int ddimind, int direction,
                      real *buf_s, int off_s, int n_s,
                      real *buf_r, int off_r, int n_r)
{
    int        rank_s, rank_r;
    static int off_l = -1;
    static int call = 1;
    gmx_domdec_shmem_buf_t * shmem = dd->shmem;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    shmem_float_sendrecv_off(shmem, buf_s, off_s, n_s, rank_s, buf_r, off_r, n_r, rank_r);
}

void dd_sendrecv_int_nobuf(const gmx_domdec_t *dd,
		int ddimind, int direction,
		int *buf_s, int n_s,
		int *buf_r, int n_r)
{
	int        rank_s, rank_r;
	gmx_domdec_shmem_buf_t * shmem = dd->shmem;
	rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
	rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

	if (n_s)
	{
		shmem_int_put(buf_r, buf_s, n_s, rank_s);
	}
	shmem_quiet();
	shmem_fence();

	shmem_set_post(shmem, rank_s);

	shmem_wait_post(shmem, _my_pe());
	shmem_clear_post(shmem, _my_pe());

}

void dd_sendrecv_real_nobuf(const gmx_domdec_t *dd,
		int ddimind, int direction,
		real *buf_s, int n_s,
		real *buf_r, int n_r)
{
	int        rank_s, rank_r;
	gmx_domdec_shmem_buf_t * shmem = dd->shmem;
	rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
	rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

	shmem_sendrecv_nobuf(shmem, buf_s, n_s * sizeof(real), rank_s,
			buf_r, n_r * sizeof(real), rank_r);

}

void dd_sendrecv_rvec_off(const gmx_domdec_t *dd,
                      int ddimind, int direction,
                      rvec *buf_s, int off_s, int n_s,
                      rvec *buf_r, int off_r, int n_r)
{
    int        rank_s, rank_r;
#ifdef GMX_SHMEM_GETMEM
    static int off_l = -1;
    static int call = 1;
#endif
    gmx_domdec_shmem_buf_t * shmem = dd->shmem;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];
#ifdef GMX_SHMEM_GETMEM
    SHDEBUG(" Get from rank_r %d , ptr %p. Put off %d  to proc %d \n", rank_r, buf_s, off_s,  rank_s);
    shmem_wait_for_previous_call(dd->shmem, &call, rank_s);
    shmem_int_p(&off_l, off_s, rank_s);
    shmem_quiet();
    SHDEBUG(" Waiting for fw to be != -1 \n")
    shmem_int_wait(&off_l, -1);
    shmem_wait_for_previous_call(dd->shmem, &call, rank_r);
	if (n_r)
	{
		shmem_getmem( buf_r + off_r, buf_s + off_l, n_r * sizeof(rvec), rank_r );
	}
	shmem_set_done(shmem, rank_r);

	shmem_wait_done(shmem, _my_pe());
	shmem_clear_done(shmem, _my_pe());

    off_l = -1;
    call++;
#else
    shmem_sendrecv_nobuf_put(shmem, *(buf_s + off_s), n_s * sizeof(rvec), rank_s,
    							    *(buf_r + off_r), n_r * sizeof(rvec), rank_r);
#endif

}

void dd_sendrecv_rvec_swap_off(const gmx_domdec_t *dd,
                      int ddimind, int direction,
                      rvec *buf_s, int off_s, int n_s,
                      rvec *buf_r, int off_r, int n_r)
{
    int        rank_s, rank_r;
#ifdef GMX_SHMEM_GETMEM
    static int off_l = -1;
    static int call = 1;
#endif
    gmx_domdec_shmem_buf_t * shmem = dd->shmem;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    /* shmem_float_sendrecv_nobuf(shmem, *(buf_s + off_s), n_s * sizeof(rvec), rank_s,
    							    *(buf_r + off_r), n_r * sizeof(rvec), rank_r); */
    shmem_float_sendrecv_off(shmem, *buf_s, off_s*DIM, n_s*DIM, rank_s,
    		                        *buf_r, off_r*DIM, n_r*DIM, rank_r);

}

void dd_sendrecv2_rvec_off(const gmx_domdec_t *dd,
                       int ddimind,
                       rvec *buf_s_fw, int off_s_fw, int n_s_fw,
                       rvec *buf_r_fw, int off_r_fw, int n_r_fw,
                       rvec *buf_s_bw, int off_s_bw, int n_s_bw,
                       rvec *buf_r_bw, int off_r_bw, int n_r_bw)
{

	int         rank_fw, rank_bw;
	static int off_fw = -1;
	static int off_bw = -1;
	static int rank_fw_ready = -1;
	static int rank_bw_ready = -1;
	int nrcall;
	gmx_domdec_shmem_buf_t * shmem = dd->shmem;
	static int call = 0;
	int done_fw = 0;
	int done_bw = 0;
	int completed = 0;


	rank_fw = dd->neighbor[ddimind][0];
	rank_bw = dd->neighbor[ddimind][1];

	SHDEBUG(" SendRecv2 (S1: %d,R1: %d) using SHMEM (n_s_fw %d, n_r_bw %d) call %d \n",
			rank_fw, rank_bw, n_s_fw, n_r_fw, call);

    /* Send offset to the getting PE */
    shmem_wait_for_previous_call(dd->shmem, &call, rank_fw);
    shmem_int_p(&off_fw, off_s_fw, rank_fw);
    shmem_quiet();
    shmem_wait_for_previous_call(dd->shmem, &call, rank_bw);
    shmem_int_p(&off_bw, off_s_bw, rank_bw);
    /* No need to send the size */
    shmem_quiet();
#ifdef FAST_SENDRECV2

    // SHDEBUG(" Before WHILE \n");
    while(!completed)
    {
    	usleep(SHMEM_SLEEP_TIME);
    	if ( !done_fw &&  (((volatile int) off_fw) != -1) )
    	{
    		// SHDEBUG(" Data appears to be received from the FW rank \n ");
    		if (n_r_fw)
    		{
    			shmem_getmem( buf_r_fw + off_r_fw, buf_s_fw + off_fw, n_r_fw * sizeof(rvec), rank_bw);
    		}
    		shmem_int_p(&rank_bw_ready, 1, rank_bw);
    		shmem_quiet();
    		// SHDEBUG(" Posted ACK to rank_bw \n")
    		done_fw = 1;
    	}
    	if ( !done_bw &&  (((volatile int) off_bw) != -1) )
    	{
    		// SHDEBUG(" Data appears to be received from the BW rank \n");
    		if (n_r_bw)
    		{
    			shmem_getmem( buf_r_bw + off_r_bw, buf_s_fw + off_bw, n_r_bw * sizeof(rvec), rank_fw);
    		}

    		shmem_int_p(&rank_fw_ready, 1, rank_fw);
    		shmem_quiet();
    		// SHDEBUG(" Posted ACK to rank_fw \n");
    		done_bw = 1;
    	}
    	if ( done_bw && done_fw &&
    			( ((volatile int) rank_fw_ready) == 1 )  && (((volatile int) rank_bw_ready) == 1 ) )
    	{
    		completed = 1;
    	}

    }

    // SHDEBUG(" After WHILE \n");

#else
       while ( (((volatile int) off_fw) == -1) )
       		{
       			usleep(SHMEM_SLEEP_TIME);
       		}
    /* Forward */
	  		if (n_r_fw)
	  			{
	  				shmem_getmem( buf_r_fw + off_r_fw, buf_s_fw + off_fw, n_r_fw * sizeof(rvec), rank_bw);
	  			}
	  		shmem_int_p(&rank_bw_ready, 1, rank_bw);
	  		shmem_quiet();
    /* Backward */

	  	while  ( (((volatile int) off_bw) == -1) )
		{
			usleep(SHMEM_SLEEP_TIME);
		}

	    if (n_r_bw)
   	  			{
   	  				shmem_getmem( buf_r_bw + off_r_bw, buf_s_fw + off_bw, n_r_bw * sizeof(rvec), rank_fw);
   	  			}

	    shmem_int_p(&rank_fw_ready, 1, rank_fw);
	    shmem_quiet();


	    while ( ((volatile int) rank_fw_ready) != 1  || (((volatile int) rank_bw_ready) != 1) )
	    	{
	    	usleep(SHMEM_SLEEP_TIME);
	    	}
#endif

       off_fw = -1;
       off_bw = -1;
	   rank_fw_ready = -1;
	   rank_bw_ready = -1;

	   call++;
	   SHDEBUG(" END OF SendRecv2 (S1: %d,R1: %d) using SHMEM (n_s_fw %d, n_r_bw %d) call %d \n",
	   			rank_fw, rank_bw, n_s_fw, n_r_fw, call);

}
#endif /* GMX_SHMEM */

void dd_sendrecv_rvec(const gmx_domdec_t *dd,
                      int ddimind, int direction,
                      rvec *buf_s, int n_s,
                      rvec *buf_r, int n_r)
{
#ifdef GMX_SHMEM_XXX
    int        rank_s, rank_r;
    gmx_domdec_shmem_buf_t * shmem = dd->shmem;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    SHDEBUG(" SendRecv RVEC (S: %d,R: %d) using SHMEM (n_s %d, n_r %d) (size rvec: %ld) \n", rank_s, rank_r, n_s, n_r, sizeof(rvec));

	shmem_rvec_sendrecv(shmem, buf_s, n_s, rank_s, buf_r, n_r, rank_r);
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
#ifdef GMX_SHMEM
	int         rank_fw, rank_bw, nreq;
    gmx_domdec_shmem_buf_t * shmem = dd->shmem;

    rank_fw = dd->neighbor[ddimind][0];
    rank_bw = dd->neighbor[ddimind][1];

    shmem_rvec_sendrecv(shmem, buf_s_fw, n_s_fw, rank_fw, buf_r_fw, n_r_fw, rank_bw);

    shmem_barrier_all();
    SHDEBUG(" Second SendRecv S2:%d R2:%d (n_s_bw %d, n_r_bw %d) \n",
    		rank_bw, rank_fw, n_s_bw, n_r_bw);

    shmem_rvec_sendrecv(shmem, buf_s_bw, n_s_bw, rank_bw, buf_r_bw, n_r_bw, rank_fw);

   SHDEBUG("Done \n");
#elif defined(GMX_MPI)
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

void dd_sendrecv_real(const gmx_domdec_t *dd,
                      int ddimind, int direction,
                      real *buf_s, int n_s,
                      real *buf_r, int n_r)
{
#ifdef GMX_SHMEM
    int        rank_s, rank_r;
    gmx_domdec_shmem_buf_t * shmem = dd->shmem;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

	shmem_real_sendrecv(shmem, buf_s, n_s, rank_s, buf_r, n_r, rank_r);
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


/* IBM's BlueGene(/L) MPI_Bcast dereferences the data pointer
 * even when 0 == nbytes, so we protect calls to it on BlueGene.
 * Fortunately dd_bcast() and dd_bcastc() are only
 * called during DD setup and partition.
 */

void dd_bcast(gmx_domdec_t *dd, int nbytes, void *data)
{
#ifdef GMX_SHMEM_COLECTIVES
	static long pSync[_SHMEM_BCAST_SYNC_SIZE];
	void * buf;
	gmx_domdec_shmem_buf_t * shmem = dd->shmem;
	int i, size;
	SHDEBUG(" Bcast of %d bytes (base ptr %p) \n", nbytes, data);
	for (i = 0; i < _SHMEM_BCAST_SYNC_SIZE; i++)
	{
		pSync[i] = _SHMEM_SYNC_VALUE;
	}
	shmem_barrier_all();
	/* shmem_broadcast expects a number of elements,
	 *  and assumes that each element of the array
	 * is separated by sizeof(void *).
	 * However, nbytes is not necessary a multiple of void *,
	 * so we have to find the nearest multiple to create the
	 * temporary buffer.
	 */
	size = round_to_next_multiple(nbytes, sizeof(void *));
	shrenew(shmem, shmem->byte_buf, &(shmem->byte_alloc), size);
	buf = shmem->byte_buf;
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
#ifdef GMX_SHMEM_COLECTIVES
	static long pSync[_SHMEM_BCAST_SYNC_SIZE];
	void * buf;
	gmx_domdec_shmem_buf_t * shmem = dd->shmem;
	int i, size;
	SHDEBUG(" Bcast of %d bytes (base ptr %p) \n", nbytes, src);
	for (i = 0; i < _SHMEM_BCAST_SYNC_SIZE; i++)
	{
		pSync[i] = _SHMEM_SYNC_VALUE;
	}
	shmem_barrier_all();
	/* shmem_broadcast expects a number of elements,
	 *  and assumes that each element of the array
	 * is separated by sizeof(void *).
	 * However, nbytes is not necessary a multiple of void *,
	 * so we have to find the nearest multiple to create the
	 * temporary buffer.
	 */
	// size = round_to_next_multiple(nbytes, sizeof(void *));
	shrenew(shmem, shmem->int_buf, &(shmem->int_alloc), nbytes);
	buf = shmem->int_buf;
	if (DDMASTERRANK(dd) == _my_pe())
	{
		memcpy(buf, src, nbytes);
	}
	SHDEBUG("  buf ptr %p , masterrank %d  \n", buf, DDMASTERRANK(dd));
	shmem_broadcast(buf, buf, nbytes/sizeof(int), DDMASTERRANK(dd), 0, 0, _num_pes(), pSync);
	memcpy(dest, buf, nbytes);
	SHDEBUG(" End of routine \n");
#else
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
#endif /* GMX_MPI */
#endif /* GMX_SHMEM */
}

void dd_scatter(gmx_domdec_t *dd, int nbytes, void *src, void *dest)
{
#ifdef GMX_SHMEM_COLECTIVES
	int i;
	gmx_domdec_shmem_buf_t * shmem = dd->shmem;

	shrenew(shmem, shmem->byte_buf, &(shmem->byte_alloc), nbytes);

    SHDEBUG(" Scatter %p (nbytes %d) \n", shmem->byte_buf, nbytes);

	if (_my_pe() == DDMASTERRANK(dd))
	{

		memcpy(shmem->byte_buf, src, nbytes * _num_pes());
	}

	shmem_barrier_all();

	if (_my_pe() != DDMASTERRANK(dd))
	{
		shmem_getmem(dest, shmem->byte_buf + (_my_pe() * nbytes), nbytes, DDMASTERRANK(dd));
	}
	else
	{
		memcpy(dest, src + (_my_pe() * nbytes), nbytes);
	}
    SHDEBUG(" Outside scatter \n");
    shmem_barrier_all();
#elif defined(GMX_MPI)
    MPI_Scatter(src, nbytes, MPI_BYTE,
                dest, nbytes, MPI_BYTE,
                DDMASTERRANK(dd), dd->mpi_comm_all);
#endif
}

void dd_gather(gmx_domdec_t *dd, int nbytes, void *src, void *dest)
{
#ifdef GMX_SHMEM_COLECTIVES
	gmx_domdec_shmem_buf_t * shmem = dd->shmem;
	int size;

    size = _num_pes() * nbytes;
	shrenew(shmem, shmem->byte_buf, &(shmem->byte_alloc), size);

    shmem_putmem(shmem->byte_buf + (_my_pe() * nbytes), src, nbytes, DDMASTERRANK(dd));
	shmem_barrier_all();

	if (_my_pe() == DDMASTERRANK(dd))
	{
		memcpy(dest, shmem->byte_buf, nbytes * _num_pes());
	}
#elif defined(GMX_MPI)
    MPI_Gather(src, nbytes, MPI_BYTE,
               dest, nbytes, MPI_BYTE,
               DDMASTERRANK(dd), dd->mpi_comm_all);
#endif
}


void dd_scatterv(gmx_domdec_t *dd,
                 int *scounts, int *disps, void *sbuf,
                 int rcount, void *rbuf)
{
#ifdef GMX_SHMEM_COLECTIVES
	int i, max_count;
	gmx_domdec_shmem_buf_t * shmem = dd->shmem;
	SHDEBUG(" ScatterV %p (rcount %d) \n", shmem->byte_buf, rcount);

    max_count = 0;
    if (_my_pe() == DDMASTERRANK(dd))
    {
    	/* Extract the maximum size of the buffer to ensure
    	 * we are not exceeding the global maximum
    	 */
    	for (i = 0; i < _num_pes(); i++)
    	{
    		if (scounts[i] > max_count)
    		{
    			max_count = scounts[i];
    		}
    	}
    }

	SHDEBUG(" ScatterV %p (size to allocate %d) \n", shmem->byte_buf, max_count);
    shrenew(shmem, shmem->byte_buf, &(shmem->byte_alloc), max_count);

	if (_my_pe() == DDMASTERRANK(dd))
	{
		for (i = 0; i < _num_pes(); i++)
		{
			if (_my_pe() != i )
			{
				// shmem_reset_flag(shmem, _my_pe());
				shmem_putmem(shmem->byte_buf, sbuf + (disps[i]), scounts[i], i);
				// shmem_set_flag(shmem, _my_pe());
				SHDEBUG(" After putmem of %p (rcount %d) to %d \n", shmem->byte_buf, scounts[i], i);
			}
			else
			{
				memcpy(rbuf, sbuf + disps[i], rcount);
			}
		}
	}

	shmem_barrier_all();

	if (_my_pe() != DDMASTERRANK(dd)  )
	{
		// shmem_wait_flag(shmem, DDMASTERRANK(dd));
		memcpy(rbuf, shmem->byte_buf, rcount);
		SHDEBUG(" Received in %p (rcount %d) from root \n", shmem->byte_buf, rcount);
	}
	SHDEBUG(" Outside scatterV \n");

#elif defined(GMX_MPI)
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
#ifdef GMX_SHMEM_COLECTIVES
		int i, max_count;
		gmx_domdec_shmem_buf_t * shmem = dd->shmem;
		const int npes = _num_pes();
		int local_disps[npes];
		static long pSync[_SHMEM_BCAST_SYNC_SIZE];
		for (i = 0; i < _SHMEM_BCAST_SYNC_SIZE; i++)
		{
			pSync[i] = _SHMEM_SYNC_VALUE;
		}
		shmem_barrier_all();

		SHDEBUG(" GatherV %p (rcount %d) \n", shmem->byte_buf, scount);

	    max_count = 0;
	    if (_my_pe() == DDMASTERRANK(dd))
	    {
	    	max_count = rcounts[npes-1] + disps[npes-1];
	    }

		SHDEBUG(" GatterV %p (size to allocate %d) \n", shmem->byte_buf, max_count);
		/* We have to ensure that the size of the buffer is greater or equal than
		 *  the number of process involved, as we use the shmem->byte_buf to communicate
		 *  the displacement arrays
		 */
	    shrenew(shmem, shmem->byte_buf, &(shmem->byte_alloc), max(max_count, _num_pes() * sizeof(int)));
        if (_my_pe() == DDMASTERRANK(dd))
        {
        	memcpy(shmem->byte_buf, disps, npes * sizeof(int));
        	for (i = 0; i < npes; i++)
               	{
               		SHDEBUG("==> %d = %d (num: %d)\n", i, disps[i], rcounts[i])
               	}
        }
       	SHDEBUG("Broadcast displ   ptr %p , masterrank %d  \n", disps, DDMASTERRANK(dd));

       	shmem_broadcast(shmem->byte_buf, shmem->byte_buf, npes * sizeof(int), DDMASTERRANK(dd), 0, 0, npes, pSync);
       	memcpy(local_disps, shmem->byte_buf, npes * sizeof(int));
       	SHDEBUG(" Broadcast done, now putting %d from %p in %p + %d \n",
       			  scount, sbuf, shmem->byte_buf, local_disps[my_pe()]);
    	for (i = 0; i < npes; i++)
           	{
           		SHDEBUG("==> %d = %d \n", i, local_disps[i])
           	}
    	if (scount)
       	    shmem_putmem(shmem->byte_buf + local_disps[_my_pe()], sbuf, scount, DDMASTERRANK(dd));
       	shmem_barrier_all();
       	SHDEBUG(" Buf of 0 wrote, this PE contributed %d bytes \n", scount);
       	if (_my_pe() == DDMASTERRANK(dd))
       	{
       		memcpy(rbuf, shmem->byte_buf, max_count);
       		SHDEBUG("Root has its copy  \n")
       	}

#elif defined(GMX_MPI)
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
