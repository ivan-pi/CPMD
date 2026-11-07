! External procedures moved out of md_driver.mod.F90 (LFortran workaround, report 26)
#include "cpmd_global.h"
! ==================================================================
SUBROUTINE extrapwf(infi,c0,gam,cold,nnow,numcold,nstate,m,scold,only_save)
  USE dotp_utils, ONLY: dotp
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum
  USE system , ONLY:cntl,ncpw,nkpt,parm
  USE parac, ONLY : paral,parai
  USE kpts, ONLY : tkpts
  USE kpnt, ONLY : wk
  USE spin, ONLY : spin_mod
  USE elct, ONLY : crge
  USE geq0mod , ONLY:geq0
  USE pslo, ONLY : pslo_com
  USE legendre_p_utils, ONLY : d_binom
  USE ortho_utils, ONLY: ortho,preortho
  USE ovlap_utils, ONLY : ovlap
  USE rnlsm_utils, ONLY: rnlsm
  USE rotate_utils, ONLY: rotate
  USE rgsvan_utils, ONLY: rgsvan
  USE spsi_utils, ONLY: spsi
  USE rnlsm_utils, ONLY: rnlsm
  USE sfac, ONLY: fnl_packed, il_fnl_packed
  USE rhov_utils, ONLY: rhov
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: nstate, infi, nnow, numcold, m
  COMPLEX(real_8)                            :: cold(nkpt%ngwk,nstate,nkpt%nkpnt,*)
  REAL(real_8)                               :: gam(nstate,*)
  COMPLEX(real_8)                            :: c0(nkpt%ngwk,nstate,*)
  COMPLEX(real_8), optional                  :: scold(nkpt%ngwk,nstate,nkpt%nkpnt,*)
  
  LOGICAL, INTENT(IN), OPTIONAL              :: only_save
  logical                                    :: save
  INTEGER                                    :: i, ik, isub, ma, n1, nm, nnow_save, ig, &
                                                nstates(2,1)
  REAL(real_8)                               :: fa, rsum, scalef, rsumv, fac_r
  COMPLEX(real_8)                            :: fac_c
  REAL(real_8), EXTERNAL                     :: ddot
  CHARACTER(*), PARAMETER                    :: procedureN = 'extrapwf'
  CALL tiset(procedureN,isub)

  IF(pslo_com%tivan)THEN
     IF (PRESENT(only_save)) save=only_save
     IF(save)THEN
        nnow_save=MOD(nnow,m)+1
        CALL dcopy(2*nkpt%ngwk*nstate*nkpt%nkpnt,c0,1,cold(1,1,1,nnow_save),1)
        IF(pslo_com%tivan)THEN
           CALL dcopy(2*nkpt%ngwk*nstate*nkpt%nkpnt,c0,1,scold(1,1,1,nnow_save),1)
           CALL spsi(nstate,scold(:,:,1,nnow_save),fnl_packed,.FALSE.)
        END IF
        numcold=MIN(numcold+1,m)
        CALL tihalt(procedureN,isub)
        RETURN
     ELSE
        nnow=MOD(nnow,m)+1        
     END IF
  ELSE
     nnow=MOD(nnow,m)+1
     CALL dcopy(2*nkpt%ngwk*nstate*nkpt%nkpnt,c0,1,cold(1,1,1,nnow),1)
     numcold=MIN(numcold+1,m)
  END IF

  ma=numcold
  IF (ma.GT.1) THEN
     n1=nnow
     fa=-1.0_real_8
     DO i=1,ma
        IF(i.eq.1)THEN
           fac_r=0._real_8
           fac_c=CMPLX(0._real_8,0._real_8)
        ELSE
           fac_r=1._real_8
           fac_c=CMPLX(1._real_8,0._real_8)
        END IF

        nm=nnow-i+1
        IF (nm.LE.0) nm=nm+m
        ! construct extrapolation polynomial coefficient.
        IF (cntl%taspc) THEN
           fa=(-1._real_8)**(i+1)*REAL(i,kind=real_8)*d_binom(2*ma,ma-i)&
                /d_binom(2*ma-2,ma-1)
        ELSE
           fa=-1._real_8*fa*REAL(ma-i+1,kind=real_8)/REAL(i,kind=real_8)
        ENDIF
        IF (tkpts%tkpnt) THEN
           DO ik=1,nkpt%nkpnt
              CALL ovlap_c(nstate,gam,cold(:,:,ik,nm),cold(:,:,ik,n1))
              CALL mp_sum(gam,2*nstate*nstate,parai%allgrp)
              CALL rotate_c(CMPLX(fa,0._real_8,kind=real_8),cold(1,1,ik,nm),&
                   fac_c,c0(1,1,ik),gam,nstate)
           ENDDO
        ELSE
           IF(pslo_com%tivan)THEN
              CALL ovlap(nstate,gam,scold(:,:,1,nm),cold(:,:,1,n1))
           ELSE
              CALL ovlap(nstate,gam,cold(:,:,1,nm),cold(:,:,1,n1))
           END IF
           CALL mp_sum(gam,nstate*nstate,parai%allgrp)
           CALL rotate(fa,cold(:,:,1,nm),fac_r,c0(:,:,1),gam,nstate,2*nkpt%ngwk,&
                cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
        ENDIF
     ENDDO
     ! ..COUNT NUMBER OF ELECTRONS IN EXTRAPOLATED WFN AND RENORMALIZE
     IF(.NOT.cntl%nonort) THEN
        rsum=0._real_8
        IF (tkpts%tkpnt) THEN
           DO ik=1,nkpt%nkpnt
              DO i=1,nstate
                 IF (crge%f(i,ik).NE.0._real_8) THEN
                    rsum=rsum+&
                         crge%f(i,ik)*wk(ik) *&
                         ddot(nkpt%ngwk*2,c0(1,i, ik),1,c0(1,i,ik),1)
                 ENDIF
              ENDDO
           ENDDO
        ELSE
           DO i=1,nstate
              IF (crge%f(i,1).NE.0._real_8) THEN
                 rsum=rsum+crge%f(i,1)*dotp(ncpw%ngw,c0(:,i,1),c0(:,i,1))
              ENDIF
           ENDDO
           IF(pslo_com%tivan)THEN
              CALL rnlsm(c0(:,:,1),nstate,1,1,.FALSE.,unpack_dfnl_fnl=.FALSE.)
              nstates(1,1)=1
              nstates(2,1)=nstate
              CALL rhov(nstates,rsumv,hfx=.FALSE.,dipole=.FALSE.)
              rsum=rsum+parm%omega*rsumv
           END IF
        ENDIF
        CALL mp_sum(rsum,parai%allgrp)
        scalef=crge%nel/rsum
        CALL dscal(2*nkpt%ngwk*nstate*nkpt%nkpnt,scalef,c0,1)
        IF(pslo_com%tivan) THEN
           CALL dscal(nstate*il_fnl_packed(1),scalef,fnl_packed,1)
           CALL rgsvan(c0(:,:,1),nstate,gam,store_nonort=.FALSE.)
        END IF
     END IF
  ENDIF
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE extrapwf
! ==================================================================
