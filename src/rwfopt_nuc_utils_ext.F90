! External procedures moved out of rwfopt_nuc_utils.mod.F90 (LFortran workaround, report 26)
! ==================================================================
SUBROUTINE lseprof_dum(c0,rhoe,psi,tau0,fion,nstate)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE fft_maxfft,                      ONLY: maxfftn
  USE system , ONLY:maxsys,ncpw,fpar
  USE parac, ONLY : paral,parai
  USE spin , ONLY:clsd,lspin2,lspin3
  USE ener , ONLY:ener_com
  USE cnst , ONLY:fpi
  USE rscpot_utils, ONLY : rscpot
  IMPLICIT NONE
  REAL(real_8)                               :: rhoe(fpar%nnr1,clsd%nlsd)
  COMPLEX(real_8)                            :: psi(maxfftn,clsd%nlsd)
  REAL(real_8)                               :: tau0(3,maxsys%nax,maxsys%nsx),&
                                                fion(3,maxsys%nax,maxsys%nsx)
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: c0(ncpw%ngw,nstate)

  CHARACTER(*), PARAMETER                    :: procedureN = 'lseprof_dum'

  COMPLEX(real_8), ALLOCATABLE               :: ca(:), cb(:), sca(:)
  INTEGER                                    :: ia, ierr
  LOGICAL                                    :: tccc
  REAL(real_8)                               :: ang, angd, cosa, deee, ezero, &
                                                sina, tpi

  IF (lspin2%tcas22) THEN
     lspin2%troks=.TRUE.
     lspin2%tcas22=.FALSE.
     tccc=.TRUE.
  ELSE
     tccc=.FALSE.
  ENDIF
  tpi=0.5_real_8*fpi

  ALLOCATE(ca(ncpw%ngw),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  ALLOCATE(cb(ncpw%ngw),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)

  CALL dcopy(2*ncpw%ngw,c0(1,clsd%ialpha),1,ca,1)
  CALL dcopy(2*ncpw%ngw,c0(1,clsd%ibeta),1,cb,1)
  IF (paral%parent) THEN
     IF (paral%io_parent)&
          WRITE(6,'(/,1X,64(1H*))')
     IF (paral%io_parent)&
          WRITE(6,*)
     IF (paral%io_parent)&
          WRITE(6,'(T20,A,T60,A)') '  Rotation Angle ','Energy'
  ENDIF
  lspin3%rotab=lspin3%rotab/tpi*360._real_8
  DO ia=0,9
     ang=ia*5._real_8
     angd=ang/360._real_8*tpi
     cosa=COS(angd)
     sina=SIN(angd)
     CALL dcopy(2*ncpw%ngw,ca,1,c0(1,clsd%ialpha),1)
     CALL dscal(2*ncpw%ngw,cosa,c0(1,clsd%ialpha),1)
     CALL daxpy(2*ncpw%ngw,sina,cb,1,c0(1,clsd%ialpha),1)
     CALL dcopy(2*ncpw%ngw,cb,1,c0(1,clsd%ibeta),1)
     CALL dscal(2*ncpw%ngw,cosa,c0(1,clsd%ibeta),1)
     CALL daxpy(2*ncpw%ngw,-sina,ca,1,c0(1,clsd%ibeta),1)

     CALL rscpot(c0,tau0,fion,rhoe,psi,&
          .FALSE.,.FALSE.,nstate,1)
     IF (ia.EQ.0) ezero=ener_com%etot
     deee=ener_com%etot-ezero
     IF (paral%io_parent)&
          WRITE(6,'(T28,F8.0,T50,F16.8)') lspin3%rotab+ang,deee
  ENDDO
  IF (paral%parent) THEN
     IF (paral%io_parent)&
          WRITE(6,'(1X,64(1H*))')
     IF (paral%io_parent)&
          WRITE(6,'(A,T50,F16.8)') "  LAGRANGE MULTIPLIERS : A-B ",lspin3%hablse
     IF (paral%io_parent)&
          WRITE(6,'(A,T50,F16.8)') "  LAGRANGE MULTIPLIERS : A-O ",lspin3%haolse
     IF (paral%io_parent)&
          WRITE(6,'(A,T50,F16.8)') "  LAGRANGE MULTIPLIERS : B-O ",lspin3%hbolse
     IF (paral%io_parent)&
          WRITE(6,'(1X,64(1H*))')
     IF (paral%io_parent)&
          WRITE(6,*)
  ENDIF

  DEALLOCATE(ca,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  DEALLOCATE(cb,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  DEALLOCATE(sca,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)

  IF (tccc) THEN
     lspin2%tcas22=.TRUE.
     lspin2%troks=.FALSE.
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE lseprof_dum
! ==================================================================
