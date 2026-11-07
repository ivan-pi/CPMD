! External procedures moved out of epr_efg_utils.mod.F90 (LFortran workaround, report 26)
! =======================================================
SUBROUTINE calc_gradient (rhor,rhoe,psi,i,nfi)
  ! =======================================================  
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum
  USE system , ONLY:maxsys,ncpw,fpar
  USE parac, ONLY : paral,parai
  USE ions , ONLY:ions0,ions1
  USE spin , ONLY:clsd
  USE epr_efg_utils, ONLY : calcion,apply_fpibg2,apply_p_d,trace,&
       & epr_efg_outp,epr_efg_outp_ev
  USE eicalc_utils, ONLY : eicalc
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  REAL(real_8)                               :: rhor(fpar%nnr1)
  COMPLEX(real_8)                            :: rhoe(ncpw%nhg), &
                                                psi(fpar%nnr1,clsd%nlsd)
  INTEGER                                    :: i, nfi

  CHARACTER(*), PARAMETER                    :: procedureN = 'calc_gradient'

  COMPLEX(real_8), ALLOCATABLE               :: scr(:,:)
  INTEGER                                    :: dir, dir1, iat, ierr, ig, &
                                                info, isa, isp
  INTEGER, SAVE                              :: ifirst = 0
  REAL(real_8)                               :: eigr_dot, work(3*3)
  REAL(real_8), ALLOCATABLE, SAVE            :: efield(:,:,:), &
                                                eigenvalues(:,:,:), &
                                                gradient_matrix(:,:,:,:), &
                                                pot(:,:)

  ALLOCATE(scr(ncpw%nhg, 3),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)

  IF (ifirst.EQ.0) THEN
     ifirst=1
     ! TODO deallocate - ?
     ALLOCATE(gradient_matrix(3,3,maxsys%nax,maxsys%nsx),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(efield(3,maxsys%nax,maxsys%nsx),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(pot(maxsys%nax,maxsys%nsx),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     ALLOCATE(eigenvalues(3,maxsys%nax,maxsys%nsx),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     CALL zeroing(gradient_matrix)!,3*3*maxsys%nax*maxsys%nsx)
     CALL zeroing(efield)!,3*maxsys%nax*maxsys%nsx)
     CALL zeroing(pot)!,maxsys%nax*maxsys%nsx)
     CALL zeroing(eigenvalues)!,3*maxsys%nax*maxsys%nsx)
  ENDIF
  CALL ffttog(rhor,rhoe,psi,ncpw%nhg,.TRUE.)
  IF (i.EQ.2) THEN
     CALL eicalc(scr(:,1),scr(:,2))
     DO ig=1,ncpw%nhg
        rhoe(ig)=rhoe(ig)+scr(ig,2)
     ENDDO
  ENDIF
  isa = 0
  DO isp=1,ions1%nsp
     DO iat=1,ions0%na(isp)
        isa = isa + 1
        DO ig=1,ncpw%nhg
           scr(ig,1)=rhoe(ig)
        ENDDO
        IF (i.EQ.2) THEN
           CALL calcion(scr(1,2),isa,isp)
           DO ig=1,ncpw%nhg
              scr(ig,1)=scr(ig,1)-scr(ig,2)
           ENDDO
        ENDIF
        CALL apply_fpibg2 (scr(1,1))
        IF (i.EQ.1) THEN
           pot(iat,isp)=eigr_dot(isa,rhoe)
        ELSE
           pot(iat,isp)=(-1._real_8)*eigr_dot(isa,scr(1,1))
        ENDIF
        DO dir=1,3
           CALL apply_p_d (scr(1,1),scr(1,2),dir)
           IF (i.EQ.2) THEN
              efield(dir,iat,isp)=eigr_dot(isa,scr(1,2))
           ENDIF
           DO dir1=dir,3
              CALL apply_p_d (scr(1,2),scr(1,3),dir1)
              IF ((i.EQ.2).AND.(dir1.EQ.dir)) THEN
                 efield(dir,iat,isp)=eigr_dot(isa,scr(1,2))
              ENDIF
              gradient_matrix(dir1,dir,iat,isp) &
                   =eigr_dot(isa,scr(1,3))
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  DEALLOCATE(scr,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)

  CALL mp_sum(gradient_matrix,3*3*maxsys%nax*maxsys%nsx,parai%allgrp)
  CALL mp_sum(efield,3*maxsys%nax*maxsys%nsx,parai%allgrp)
  CALL mp_sum(pot,maxsys%nax*maxsys%nsx,parai%allgrp)
  IF (i.EQ.2) CALL trace(gradient_matrix)
  CALL epr_efg_outp(gradient_matrix,efield,pot,i,nfi)
  DO isp=1,ions1%nsp
     DO iat=1,ions0%na(isp)
        CALL dsyev('N','L',3,gradient_matrix(1,1,iat,isp),3,&
             eigenvalues(1,iat,isp),work,9,info)
        IF (paral%io_parent.AND.(info.NE.0)) THEN
           WRITE (6, *)'ERROR IN DSYEV, code ',info
        ENDIF
     ENDDO
  ENDDO
  CALL epr_efg_outp_ev(eigenvalues,i,nfi)
  RETURN
END SUBROUTINE calc_gradient
