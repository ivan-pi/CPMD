! External procedures moved out of wv30_utils.mod.F90 (LFortran workaround, report 26)
#include "cpmd_global.h"
SUBROUTINE wr30pot(nw,ierror,potr)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sync,mp_send,mp_recv
  USE system , ONLY: fpar,spar, parap
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  INTEGER                                    :: nw, ierror
  REAL(real_8)                               :: potr(fpar%kr1,*)

  CHARACTER(*), PARAMETER                    :: procedureN = 'wr30pot'

  INTEGER                                    :: i, ierr, ip, ipp, ipx, ir, &
                                                irr, kk, msgid
  REAL(real_8), ALLOCATABLE                  :: pscr(:)

! ==--------------------------------------------------------------==

  kk=fpar%kr2s*fpar%kr3s
  ALLOCATE(pscr(kk),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  CALL mp_sync(parai%allgrp)
  ! For each plane
  DO ir=1,spar%nr1s
     ! Find the processor which has the given plane.
     DO ipp=1,parai%nproc
        IF (ir.GE.parap%nrxpl(ipp-1,1).AND.ir.LE.parap%nrxpl(ipp-1,2)) THEN
           ipx=ipp
           GOTO 10
        ENDIF
     ENDDO
     CALL stopgm('WR30POT','ERROR IN NRXPL',& 
          __LINE__,__FILE__)
10   CONTINUE
     ip=parap%pgroup(ipx)
     ! IP has the plane.
     IF (parai%me.EQ.ip) THEN
        irr=ir-parap%nrxpl(parai%mepos,1)+1
        CALL dcopy(kk,potr(irr,1),fpar%kr1,pscr(1),1)
        IF (.NOT.paral%parent) THEN
           ! IP sends array.
           msgid=1
           !msglen=kk * 8
           CALL mp_send(pscr,kk,parap%pgroup(1),msgid,parai%allgrp)
        ENDIF
     ELSEIF (paral%parent) THEN
        ! PARENT receives array.
        msgid=1
        !msglen=kk * 8
        CALL mp_recv(pscr,kk,ip,msgid,parai%allgrp)
     ENDIF
     IF (paral%parent) THEN
        IF (paral%io_parent)&
             WRITE(nw,iostat=ierror) (pscr(i),i=1,kk)
     ENDIF
     CALL mp_sync(parai%allgrp)
  ENDDO
  DEALLOCATE(pscr,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE wr30pot
! ==================================================================
