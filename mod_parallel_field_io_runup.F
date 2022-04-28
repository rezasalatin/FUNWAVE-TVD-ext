!------------------------------------------------------------------------------------
!
!      FILE mod_parallel_field_io_runup.F
!
!      This file is part of the FUNWAVE-TVD program under the Simplified BSD license
!
!-------------------------------------------------------------------------------------
!
!    Copyright (c) 2016, FUNWAVE Development Team
!
!    (See http://www.udel.edu/kirby/programs/funwave/funwave.html
!     for Development Team membership)
!
!    All rights reserved.
!
!    FUNWAVE_TVD is free software: you can redistribute it and/or modify
!    it under the terms of the Simplified BSD License as released by
!    the Berkeley Software Distribution (BSD).
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!
!    1. Redistributions of source code must retain the above copyright notice, this
!       list of conditions and the following disclaimer.
!    2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!    The views and conclusions contained in the software and documentation are those
!    of the authors and should not be interpreted as representing official policies,
!    either expressed or implied, of the FreeBSD Project.
!
!-------------------------------------------------------------------------------------
!
!    PARALLEL_FIELD_IO_Runup is subroutine for print-out of field data
!
!    HISTORY:
!      10/26/2020  Chen, Salatin
!-------------------------------------------------------------------------------------

# if defined (PARALLEL)

MODULE PARALLEL_FIELD_IO_RUNUP
USE GLOBAL

CONTAINS

SUBROUTINE PutFileASCII_RUNUP(FILE,PHI)

  IMPLICIT NONE
  CHARACTER(LEN=80), INTENT(IN) :: FILE
  REAL(SP),DIMENSION(Mloc,Nloc),INTENT(IN) :: PHI
  INTEGER :: irank, lenx, leny, lenxy, ireq, jj
  INTEGER, ALLOCATABLE :: iistas(:), iiends(:), jjstas(:), jjends(:)
  INTEGER :: istatus(mpi_status_size)
  INTEGER, DIMENSION(1,NGlob) :: runup_matrix
  REAL(SP), ALLOCATABLE :: LocalPHI(:,:)
  REAL(SP), ALLOCATABLE :: xx(:,:)
  REAL(SP),DIMENSION( MGlob, NGlob ) :: PHIGLOB

  ALLOCATE(LocalPHI(Mloc-2*Nghost, Nloc-2*Nghost))
  LocalPHI=PHI(Ibeg:Iend,Jbeg:Jend)

  allocate( iistas(nprocs), iiends(nprocs), jjstas(nprocs), jjends(nprocs) )

  call MPI_Gather( iista, 1, MPI_INTEGER, iistas, 1, MPI_INTEGER, &
  0, MPI_COMM_WORLD, ier )
  call MPI_Gather( iiend, 1, MPI_INTEGER, iiends, 1, MPI_INTEGER, &
  0, MPI_COMM_WORLD, ier )
  call MPI_Gather( jjsta, 1, MPI_INTEGER, jjstas, 1, MPI_INTEGER, &
  0, MPI_COMM_WORLD, ier )
  call MPI_Gather( jjend, 1, MPI_INTEGER, jjends, 1, MPI_INTEGER, &
  0, MPI_COMM_WORLD, ier )

  if( myid == 0 ) then
    PHIGLOB( iista:iiend, jjsta:jjend ) = LocalPHI
  endif

  do irank=1, px*py-1

    if( myid == 0 ) then
      lenx = iiends(irank+1) - iistas(irank+1) + 1
      leny = jjends(irank+1) - jjstas(irank+1) + 1
      lenxy = lenx*leny
      allocate( xx(lenx, leny) )

      call mpi_irecv( xx, lenxy, mpi_sp, irank, 1, mpi_comm_world, ireq, ier )
      call mpi_wait( ireq, istatus, ier )

      PHIGLOB( iistas(irank+1):iiends(irank+1), jjstas(irank+1):jjends(irank+1) ) = xx
      deallocate( xx )

    elseif( myid == irank ) then

      lenxy = ( iiend-iista+1 )*( jjend-jjsta+1 )
      call mpi_isend( LocalPHI, lenxy, mpi_sp, 0, 1, mpi_comm_world, ireq, ier )
      call mpi_wait( ireq, istatus, ier )

    endif
  enddo

  deallocate(LocalPHI, iistas, iiends, jjstas, jjends)

  if (myid.eq.0) then
    DO J=1,NGlob
      DO I=MGlob,1,-1
        IF(real(PHIGLOB(I,J)).eq.1.0) goto 2001
      ENDDO
2001  IF((I == MGlob).OR.(I == 1)) THEN
      ! do nothing
    ELSE
      runup_matrix(1,J) = I
    ENDIF
  ENDDO

  open(1, FILE=TRIM(FILE), status = 'unknown',POSITION='APPEND')

# if defined(DEBUG)
  DO jj=1,NGlob
    write(1,'(I8)',advance="no") runup_matrix(1,jj)
  ENDDO
# else
  write(1,'(20E16.5E4)',advance="no") time
  DO jj=1,NGlob
    write(1,'(I8)',advance="no") runup_matrix(1,jj)
  ENDDO
# endif
  CLOSE(1)
endif

END SUBROUTINE PutfileASCII_RUNUP

END MODULE PARALLEL_FIELD_IO_RUNUP

# endif
! end parallel
