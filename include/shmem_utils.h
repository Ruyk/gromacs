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

#ifndef SHMEM_UTILS_H_
#define SHMEM_UTILS_H_

#ifdef GMX_SHMEM
#include <shmem.h>
#include <macros.h>

#include "smalloc.h"
#include "typedefs.h"
#include "types/commrec.h"
#include "gmx_fatal.h"



#define shmem_reset_flag(FLAG) {  FLAG = 0; shmem_barrier_all(); }
#define shmem_set_flag(FLAG, TARGET) { shmem_int_p(FLAG, 1, TARGET); }
#define shmem_wait_flag(FLAG) { shmem_int_wait(FLAG, 0); }

/* round_to_next_multiple
 * =========================
 *    Round a number of bytes up to the closest value multiple of the type_size.
 *    Example: nbytes: 7 , type_size: 8   , size = 8;
 */
int round_to_next_multiple(int nbytes, int type_size);


/* get_max_alloc
 * ========================
 *   Computes the maximum value of memory requested across PEs
 */
int get_max_alloc(int local_value);

/* sh_renew_buf
 * =========================
 *   Renew a temporary shmem buffer by collectively gathering the maximum
 *   size of the array in each pe and calling the shrealloc function with that value
 */
void * sh_renew_buf(void * buf, int * alloc, const int new_size, const int elem_size);

/* shrenew macro simplifies the usage of the sh_renew_buf function in the same fashion of
 *  those found in the smalloc.{h,c} files.
 */
#ifndef GMX_SHMEM_DEBUG
#define shrenew(PTR, OLD_SIZE, NEW_SIZE) (PTR) = sh_renew_buf((PTR), (OLD_SIZE), (NEW_SIZE), sizeof(*(PTR)))
#else

#define shrenew(PTR, OLD_SIZE, NEW_SIZE) { SHDEBUG(" Before renew , %p size %d \n", PTR, *(OLD_SIZE)); \
 					   (PTR) = sh_renew_buf((PTR), (OLD_SIZE), (NEW_SIZE), sizeof(*(PTR))); \
    					   SHDEBUG(" After renew , %p size %d \n",  PTR, *(OLD_SIZE));\
					 }

/* shmem_cleanup
 * =========================
 * 	  Free the shmem temporary buffers , fatal if shmem is NULL.
 */
void shmem_cleanup(gmx_domdec_shmem_buf_t * shmem);

#endif



#endif /* GMX_SHMEM */

#endif /* SHMEM_UTILS_H_ */
