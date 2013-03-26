#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.
#
# Manage the SHMEM setup, which requires avaiability of MPI support
if(GMX_SHMEM)
  message(STATUS " Trying to build with SHMEM support.")
  # set(SHMEM_BASE "/home/rreyesc/soft/openshmem/")
  # set(SHMEM_INCLUDE "${SHMEM_BASE}/include/")
  # set(SHMEM_LIB     "${SHMEM_BASE}/lib/")
  # include_directories(${SHMEM_INCLUDE}) 
  
  TRY_COMPILE(SHMEM_FOUND ${CMAKE_BINARY_DIR}
    "${CMAKE_SOURCE_DIR}/cmake/TestSHMEM.c"
    OUTPUT_VARIABLE output )

  IF (NOT SHMEM_FOUND)
    message(ERROR " Could not compile SHMEM support. ")
    message(ERROR " ${output} ${SHMEM_INCLUDE}")
  ENDIF()

  set(GMX_SHMEM 1)
  set(PKG_CFLAGS "${PKG_CFLAGS} -DGMX_SHMEM")
  add_definitions("-DGMX_SHMEM=1")

endif(GMX_SHMEM)
