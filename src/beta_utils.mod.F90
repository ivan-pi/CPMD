#include "cpmd_global.h"

MODULE beta_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlps_com,&
                                             nghtol
  USE system,                          ONLY: maxsys

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: build_beta

CONTAINS
  ! ==================================================================
  SUBROUTINE build_beta(na,eigkr,eiscr,ngw,ig_start,ld_eiscr,td_eiscr,igeq0,geq0,tkpnt,twnl_nghtol,twnl_nghtol_gk)
    ! ==--------------------------------------------------------------==
    ! == build derivative of beta projectors according to na          ==
    ! == call from parallel region to get omp parallelism here        ==
    ! == WARNING: do not forget to set geq0 component outside!        ==
    ! == WARNING: no omp barrier at the end of this routine!          ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: na(2,*),ngw,ig_start,ld_eiscr
    COMPLEX(real_8),INTENT(IN)               :: eigkr(ngw,*)
    COMPLEX(real_8),INTENT(OUT)              :: eiscr(ld_eiscr,*)
    REAL(real_8),INTENT(IN),OPTIONAL,TARGET  :: twnl_nghtol(ngw,maxsys%nhxs,ions1%nsp,*)
    REAL(real_8),INTENT(IN),OPTIONAL,TARGET  :: twnl_nghtol_gk(ngw,3,maxsys%nhxs,ions1%nsp,*)
    LOGICAL,INTENT(IN),OPTIONAL              :: geq0,tkpnt
    INTEGER,INTENT(IN),OPTIONAL              :: igeq0,td_eiscr
    INTEGER                                  :: isa0, is, start_ia, end_ia, ia_sum, iv,&
                                                ia, isa, offset, offset0, ig, k, kit, i
    REAL(real_8)                             :: cir
    COMPLEX(real_8)                          :: ci
    LOGICAL                                  :: derivative, clean_geq0
    COMPLEX(real_8), PARAMETER               :: zzero = (0._real_8,0._real_8)
    REAL(real_8), POINTER __CONTIGUOUS       :: twnl(:)

    ! ==--------------------------------------------------------------==
    IF(PRESENT(twnl_nghtol_gk))THEN
       derivative=.TRUE.
       kit=3
    ELSE
       derivative=.FALSE.
       kit=1
    END IF
    IF(PRESENT(geq0))THEN
       clean_geq0=geq0
    ELSE
       clean_geq0=.FALSE.
    END IF
    !$omp parallel private(offset0,isa0,is,start_ia,end_ia,ia_sum,iv,&
    !$omp& k,ci,cir,twnl,ia,isa,offset)
    offset0=0
    isa0=0
    DO is=1,ions1%nsp
       start_ia=na(1,is)
       end_ia=na(2,is)
       ia_sum=end_ia-start_ia+1
       IF(ia_sum.GT.0)THEN
          DO iv=1,nlps_com%ngh(is)
             DO k=1,kit
                IF(derivative) THEN
                   twnl(1:)=>twnl_nghtol_gk(ig_start:,k,iv,is,1)
                   ci=(0.0_real_8,-1.0_real_8)**(nghtol(iv,is)+1)
                   cir=REAL(ci,kind=real_8)
                ELSE
                   ci=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
                   cir=REAL(ci,kind=real_8)
                   twnl(1:)=>twnl_nghtol(ig_start:,iv,is,1)
                END IF
                !$omp do
                DO ia=start_ia,end_ia
                   isa=isa0+ia
                   offset=offset0+ia-start_ia+1
                   IF (ABS(cir).GT.0.5_real_8) THEN
                      CALL build_betaproj_r(eiscr(:,offset),eigkr(ig_start,isa),twnl,ld_eiscr)
                   ELSE
                      CALL build_betaproj_i(eiscr(:,offset),eigkr(ig_start,isa),twnl,ld_eiscr)
                   END IF
                END DO
                !$omp end do nowait
                offset0=offset0+ia_sum
             END DO
          END DO
       END IF
       isa0=isa0+ions0%na(is)
    END DO

    IF(clean_geq0)THEN
       !$omp barrier
       IF(tkpnt)THEN
          !$omp do
          DO i=1,td_eiscr
             eiscr(igeq0,i)=zzero
          END DO
          !$omp end do nowait
       ELSE
          !$omp do
          DO i=1,td_eiscr
             eiscr(igeq0,i)=eiscr(igeq0,i)*0.5_real_8
          END DO
          !$omp end do nowait
       END IF
    END IF
    !$omp end parallel 
    RETURN
  END SUBROUTINE build_beta
  ! ==================================================================
  PURE SUBROUTINE build_betaproj_r(eiscr,eigr_,twnl_,lda)
    INTEGER,INTENT(IN)                       :: lda
    COMPLEX(real_8), INTENT(OUT)             :: eiscr(lda)
    COMPLEX(real_8),INTENT(IN)               :: eigr_(lda,*)
    REAL(real_8),INTENT(IN)                  :: twnl_(lda,*)
    REAL(real_8)                             :: t1, t2
    INTEGER                                  :: ig

    DO ig=1,lda
       t1=REAL(eigr_(ig,1))*twnl_(ig,1)
       t2=AIMAG(eigr_(ig,1))*twnl_(ig,1)
       eiscr(ig)=CMPLX(t1,t2,kind=real_8)
    END DO

  END SUBROUTINE build_betaproj_r
  ! ==================================================================
  PURE SUBROUTINE build_betaproj_i(eiscr,eigr_,twnl_,lda)
    INTEGER,INTENT(IN)                       :: lda
    COMPLEX(real_8),INTENT(OUT)              :: eiscr(lda)
    COMPLEX(real_8),INTENT(IN)               :: eigr_(lda,*)
    REAL(real_8),INTENT(IN)                  :: twnl_(lda,*)
    REAL(real_8)                             :: t1, t2
    INTEGER                                  :: ig

    DO ig=1,lda
       t1=REAL(eigr_(ig,1))*twnl_(ig,1)
       t2=AIMAG(eigr_(ig,1))*twnl_(ig,1)
       eiscr(ig)=CMPLX(-t2,t1,kind=real_8)
    END DO

  END SUBROUTINE build_betaproj_i
  ! ==================================================================
END MODULE beta_utils
