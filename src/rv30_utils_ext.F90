! External procedures moved out of rv30_utils.mod.F90 (LFortran workaround, report 26)
! ==================================================================
SUBROUTINE rd30pot(nr,irecord,nr1s0,nr2s0,nr3s0,potr)
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  ! == Read Potential or density array in RESTART file              ==
  ! == Works if the dimensions of reading array are not             ==
  ! == the same as POTR (but no guarantee)                          ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sync,mp_send, mp_recv
  USE system , ONLY:fpar,spar, parap
  USE parac, ONLY : paral,parai
  USE zeroing_utils,                   ONLY: zeroing
  USE loadpa_utils, ONLY: leadim
  IMPLICIT NONE
  INTEGER                                    :: nr, irecord, nr1s0, nr2s0, &
                                                nr3s0
  REAL(real_8) :: potr(fpar%kr1,fpar%kr2s*fpar%kr3s)

  CHARACTER(*), PARAMETER                    :: procedureN = 'rd30pot'

  INTEGER                                    :: i, ierr, ip, ipp, ipx, ir, &
                                                irr, kk, kkmax, kks0, kr1s0, &
                                                kr2s0, kr3s0, msgid
  REAL(real_8), ALLOCATABLE                  :: pscr(:)

! Variables
! (kr2s*kr3s)
! ==--------------------------------------------------------------==

  kk=fpar%kr2s*fpar%kr3s
  ! Added to avoid a bug in T3E (Thierry Deutsch 28/04/97)
  CALL zeroing(potr)!,kr1*kk)
  ! Dimensions for POT in reading file.
  CALL leadim(nr1s0,nr2s0,nr3s0,kr1s0,kr2s0,kr3s0)
  kks0=kr2s0*kr3s0
  kkmax=MAX(kks0,kk)
  ALLOCATE(pscr(kkmax),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
       __LINE__,__FILE__)
  DO ir=1,irecord
     IF (paral%io_parent) READ(nr) (pscr(i),i=1,kks0)
     IF (kks0.LT.kk) CALL zeroing(pscr(kks0+1:kk))!,kk-kks0)
     IF (irecord.GT.spar%nr1s) THEN
        GOTO 20
     ELSE
        DO ipp=1,parai%nproc
           IF (ir.GE.parap%nrxpl(ipp-1,1).AND.ir.LE.parap%nrxpl(ipp-1,2)) THEN
              ipx=ipp
              GOTO 10
           ENDIF
        ENDDO
        CALL stopgm('RD30POT','ERROR IN NRXPL',& 
             __LINE__,__FILE__)
     ENDIF
10   CONTINUE
     ip=parap%pgroup(ipx)
     CALL mp_sync(parai%allgrp)
     IF (ip.NE.parai%source) THEN
        IF (paral%parent) THEN
           msgid=1
           !msglen=kk * 8
           CALL mp_send(pscr,kk,ip,msgid,parai%allgrp)
        ELSEIF (parai%me.EQ.ip) THEN
           msgid=1
           !msglen=kk * 8
           CALL mp_recv(pscr,kk,parap%pgroup(1),msgid,parai%allgrp)
        ENDIF
        CALL mp_sync(parai%allgrp)
     ENDIF
     IF (ip.EQ.parai%me) THEN
        irr=ir-parap%nrxpl(parai%mepos,1)+1
        CALL dcopy(kk,pscr(1),1,potr(irr,1),fpar%kr1)
     ENDIF
20   CONTINUE
  ENDDO
  DEALLOCATE(pscr,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE rd30pot
