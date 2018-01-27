      SUBROUTINE ANALYSE(vect, dimvect, startvect, endvect, xymat2, norm2)
!     INCLUDE "nbtots.prm"
      IMPLICIT NONE
      INTEGER   dimvect, nbtots, startvect, endvect
      REAL*8,dimension(dimvect)::vect
      INTEGER (kind=selected_int_kind(16))::add,kvect
      INTEGER (kind=selected_int_kind(16))::iaa2,i
      INTEGER ,dimension(natomax)::ideter
      INTEGER ,dimension(natomax)::ideter2
      REAL*8,allocatable       ::xz(:)
      REAL*8::xmat,xymat
      REAL*8::xmat1,xymat1
      REAL*8::xmat2
      REAL*8,INTENT(INOUT)::xymat2
      REAL*8::xmat3,xymat3
      REAL*8::sym,nonsym,proj_trou
      REAL*8,allocatable      ::xalpha1(:)
!     REAL*8,allocatable      ::vect(:)
      REAL*8,INTENT(INOUT)::norm2
      REAL*8::norm,norm1,norm3,proj_2trou
      REAL*8::t1,t2,XS,XS1,XS2,XS3
      REAL*8::resta_mono,resta_one,resta_bi,delta
      INTEGER ::kko,kok,kkio,j,eigen,nigen,count
      INTEGER ::cntrou,countlvect,ndim,iaa
      INTEGER ::ipt_1,ipt_2,iptemp_1,iptemp_2
      INTEGER ::ipt_3,iptemp_3
      INTEGER ::ibougetrou,jstart
      REAL*8 ,allocatable::eigenvectors(:,:)
      REAL*8 ,allocatable::eigenvalues(:)
      REAL*8 ,allocatable::WORK(:)
      REAL*8 ,allocatable::AP(:)
      REAL*8 ,allocatable::densmat(:,:)
      REAL*8 ,allocatable::densmat2(:,:)
      REAL*8 proj_1,extradiag_dmat2,ionic,nonionic
      REAL*8 proj_2,sum,conduction,prob,prob2
      INTEGER INFO,nrow,ncol,mmo,mom,kk,k,omm,okk
      CHARACTER*1 JOBZ,UPLO
      INTEGER::RESTA=0
  
!     allocate(vect(nbtots))
      allocate(xalpha1(natomax))
      allocate(xz  (natom/2))

!     OPEN (unit=59,file='FIL1',form='formatted',status='old')
!     OPEN (unit=217,file='SBOX217',form='formatted',status='REPLACE')
!     REWIND  59
!     READ (59,*)
  
!     print *,' in analyse', startvect, endvect
!     print *,'nalpha=',nalpha,'nbeta=',nbeta
!     PRINT *,natom,ntrou,nbtots,nt1,nt2,isz
!     ndim=3
      ndim=(natom/2)*((natom/2)-1)/2
      allocate(AP((ndim)*((ndim)+1)/2))
      allocate(WORK(3*(ndim)))
      allocate(eigenvectors((ndim)*(ndim),1))
      allocate(eigenvalues((ndim)))
      allocate(densmat(ndim,ndim))
      allocate(densmat2(ndim,ndim))
!     Touch isz maxdet maxial maxlien maxplac nalpha natom natomax nbeta nbtots nt1 nt2 ntrou
!     PRINT *,(vect(j),j=1,30)
  
	IF(RESTA .eq. 1)THEN
	do i=1,natom/2
	if(mod(natom/2,2).eq.0)then
		xz(i)=(((natom/2)/2)-0.5d0)-(i-1.0d0)
		write(6  ,*)i,xz(i)
	else
		xz(i)=((natom/2)-1.0d0)/2.0d0-(i-1.0d0)
		write(6  ,*)i,xz(i)
	endif
	enddo
        ENDIF

!!    PROVIDE det deth
      nigen=1
      DO eigen=1,nigen
!     READ (59,10) (vect(j),j=1,nbtots)
      IF (ntrou.eq.1) THEN
        norm=0.d0
        norm1=0.d0
        norm2=0.d0
        norm3=0.d0
        count=0
        cntrou=0
        iaa=0
        iaa2=0
        countlvect=0
        proj_2trou=0.d0
        resta_bi=0.d0
        resta_mono=0.d0
        resta_one=0.d0
        xymat = 0.0d0
        xymat1 = 0.0d0
        xymat2 = 0.d0
        xymat3 = 0.d0
        proj_1=0d0
        proj_2=0d0
        densmat=0d0
        densmat2=0d0
        extradiag_dmat2=0d0
        conduction=0d0
        nrow=0
        ncol=0

            jstart=1
            ipt_1=0
            ipt_2=0
            ipt_3=0
            iptemp_1=0
            iptemp_2=0
            iptemp_3=0
            nonionic=0.d0
            ionic=0.d0
            prob=0.d0
            prob2=0.d0
            sym=0.d0
            nonsym=0.d0
            ibougetrou=0

            DO kvect=1, endvect-startvect
              CALL getdet(kvect+startvect,ideter)

!!----------------------------------------
!!  RESTA
!!----------------------------------------
!!! mono
!                 proj_trou=vect(kvect)**2
!                 DO i=1,natom/2
!                   IF (ideter(i).eq.3) THEN
!                     delta=0.0d0
!                   ELSE
!                     delta=1.0d0
!                   END IF
!                   resta_mono=resta_mono+delta*xz(i)*xz(i)*proj_trou
!                   resta_one=resta_one+delta*xz(i)*proj_trou
!                 END DO
!!! bi
!                 DO i=1,natom/2
!                   DO j=1,natom/2
!                     IF (ideter(i).eq.3.or.ideter(j).eq.3.or.i.eq.j)           &
!                      THEN
!                       delta=0.0d0
!                     ELSE
!                       delta=1.0d0
!                     END IF
!                     resta_bi=resta_bi+delta*xz(i)*xz(j)*proj_trou
!                   END DO
!                 END DO
!!----------------------------------------


!!----------------------------------------
!! Prob ionic non-ionic
!!----------------------------------------
!             ipt_1=0
!             ipt_2=0
!             ipt_3=0
!             DO kko=1,3
!               IF(ideter(kko).eq.3)THEN
!                 ipt_1=ipt_1+1
!               ENDIF
!             ENDDO

!             IF(ipt_1.eq.1)THEN
!               DO kko=4,6
!                 IF(ideter(kko).eq.3)THEN
!                   ipt_2=ipt_2+1
!                 ENDIF
!               ENDDO
!             ENDIF

!             IF(ipt_2.eq.1)THEN
!               DO kko=7,9
!                 IF(ideter(kko).eq.3)THEN
!                   ipt_3=ipt_3+1
!                 ENDIF
!               ENDDO
!             ENDIF

!             IF(ipt_3 .eq. 1)THEN
!               nonionic=nonionic+vect(kvect)**2
!             ELSE
!               ionic=ionic+vect(kvect)**2
!             ENDIF
!!----------------------------------------
!! S_box
!!----------------------------------------

              xmat=0.0d0
              xmat1=0.0d0
              xmat2=0.0d0
!             IF (.TRUE.)THEN
!!              IF (ideter(6).eq.3 ) THEN
!!                norm=norm+vect(kvect)**2
!!                DO kko=5,7
!!                  DO kok=kko,7
!!                    IF (kok.eq.kko.and.ideter(kok).ne.3) THEN
!!                      xmat=xmat+(3.d0/4.d0)*(vect(kvect)**2)
!!                    ELSE
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.1) THEN
!!                        xmat=xmat+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.2) THEN
!!                        xmat=xmat+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.2) THEN
!!                        xmat=xmat-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=2
!!                        ideter2(kok)=1
!!                        CALL adr(ideter2, iaa2)
!!                        xmat=xmat+vect(kvect)*vect(iaa2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.1) THEN
!!                        xmat=xmat-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=1
!!                        ideter2(kok)=2
!!                        CALL adr(ideter2, iaa2)
!!                        xmat=xmat+vect(kvect)*vect(iaa2)
!!                      END IF
!!                    END IF
!!                  END DO
!!                END DO
!!  
!!                DO kko=16,18
!!                  DO kok=kko,18
!!                    IF (kok.eq.kko.and.ideter(kok).ne.3) THEN
!!                      xmat=xmat+(3.d0/4.d0)*(vect(kvect)**2)
!!                    ELSE
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.1) THEN
!!                        xmat=xmat+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.2) THEN
!!                        xmat=xmat+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.2) THEN
!!                        xmat=xmat-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=2
!!                        ideter2(kok)=1
!!                        CALL adr(ideter2, iaa2)
!!                        xmat=xmat+vect(kvect)*vect(iaa2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.1) THEN
!!                        xmat=xmat-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=1
!!                        ideter2(kok)=2
!!                        CALL adr(ideter2, iaa2)
!!                        xmat=xmat+vect(kvect)*vect(iaa2)
!!                      END IF
!!                    END IF
!!                  END DO
!!                END DO
!!  
!!                DO kko=5,7
!!                  DO kok=16,18
!!                    IF (kok.eq.kko.and.ideter(kok).ne.3) THEN
!!                      xmat=xmat+(3.d0/4.d0)*(vect(kvect)**2)
!!                    ELSE
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.1) THEN
!!                        xmat=xmat+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.2) THEN
!!                        xmat=xmat+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.2) THEN
!!                        xmat=xmat-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=2
!!                        ideter2(kok)=1
!!                        CALL adr(ideter2, iaa2)
!!                        xmat=xmat+vect(kvect)*vect(iaa2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.1) THEN
!!                        xmat=xmat-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=1
!!                        ideter2(kok)=2
!!                        CALL adr(ideter2, iaa2)
!!                        xmat=xmat+vect(kvect)*vect(iaa2)
!!                      END IF
!!                    END IF
!!                  END DO
!!                END DO
!!              END IF
!!!----------------------------------------
!!              xymat=xymat+xmat
!!
!!              IF (ideter(7).eq.3 ) THEN
!!                norm1=norm1+vect(kvect)**2
!!                DO kko=5,9
!!                  DO kok=kko,9
!!                    IF (kok.eq.kko.and.ideter(kok).ne.3) THEN
!!                      xmat1=xmat1+(3.d0/4.d0)*(vect(kvect)**2)
!!                    ELSE
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.1) THEN
!!                        xmat1=xmat1+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.2) THEN
!!                        xmat1=xmat1+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.2) THEN
!!                        xmat1=xmat1-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=2
!!                        ideter2(kok)=1
!!                        CALL adr(ideter2, iaa2)
!!                        xmat1=xmat1+vect(kvect)*vect(iaa2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.1) THEN
!!                        xmat1=xmat1-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=1
!!                        ideter2(kok)=2
!!                        CALL adr(ideter2, iaa2)
!!                        xmat1=xmat1+vect(kvect)*vect(iaa2)
!!                      END IF
!!                    END IF
!!                  END DO
!!                END DO
!!  
!!                DO kko=14,18
!!                  DO kok=kko,18
!!                    IF (kok.eq.kko.and.ideter(kok).ne.3) THEN
!!                      xmat1=xmat1+(3.d0/4.d0)*(vect(kvect)**2)
!!                    ELSE
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.1) THEN
!!                        xmat1=xmat1+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.2) THEN
!!                        xmat1=xmat1+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.2) THEN
!!                        xmat1=xmat1-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=2
!!                        ideter2(kok)=1
!!                        CALL adr(ideter2, iaa2)
!!                        xmat1=xmat1+vect(kvect)*vect(iaa2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.1) THEN
!!                        xmat1=xmat1-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=1
!!                        ideter2(kok)=2
!!                        CALL adr(ideter2, iaa2)
!!                        xmat1=xmat1+vect(kvect)*vect(iaa2)
!!                      END IF
!!                    END IF
!!                  END DO
!!                END DO
!!
!!                DO kko=5,9
!!                  DO kok=14,18
!!                    IF (kok.eq.kko.and.ideter(kok).ne.3) THEN
!!                      xmat1=xmat1+(3.d0/4.d0)*(vect(kvect)**2)
!!                    ELSE
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.1) THEN
!!                        xmat1=xmat1+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.2) THEN
!!                        xmat1=xmat1+(1.d0/2.d0)*(vect(kvect)**2)
!!                      END IF
!!                      IF (ideter(kko).eq.1.and.ideter(kok).eq.2) THEN
!!                        xmat1=xmat1-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=2
!!                        ideter2(kok)=1
!!                        CALL adr(ideter2, iaa2)
!!                        xmat1=xmat1+vect(kvect)*vect(iaa2)
!!                      END IF
!!                      IF (ideter(kko).eq.2.and.ideter(kok).eq.1) THEN
!!                        xmat1=xmat1-(1.d0/2.d0)*(vect(kvect)**2)
!!                        DO kkio=1,natom
!!                          ideter2(kkio)=ideter(kkio)
!!                        END DO
!!                        ideter2(kko)=1
!!                        ideter2(kok)=2
!!                        CALL adr(ideter2, iaa2)
!!                        xmat1=xmat1+vect(kvect)*vect(iaa2)
!!                      END IF
!!                    END IF
!!                  END DO
!!                END DO
!!              END IF
!!!----------------------------------------
!!              xymat1=xymat1+xmat1

              IF (.TRUE.)THEN
                norm2=norm2+vect(kvect)**2
                DO kko=1,natom/2
                  DO kok=kko,natom/2
                    IF (kok.eq.kko.and.ideter(kok).ne.3) THEN
                      xmat2=xmat2+(3.d0/4.d0)*(vect(kvect)**2)
                    ELSE
                      IF (ideter(kko).eq.1.and.ideter(kok).eq.1) THEN
                        xmat2=xmat2+(1.d0/2.d0)*(vect(kvect)**2)
                      END IF
                      IF (ideter(kko).eq.2.and.ideter(kok).eq.2) THEN
                        xmat2=xmat2+(1.d0/2.d0)*(vect(kvect)**2)
                      END IF
                      IF (ideter(kko).eq.1.and.ideter(kok).eq.2) THEN
                        xmat2=xmat2-(1.d0/2.d0)*(vect(kvect)**2)
                        DO kkio=1,natom
                          ideter2(kkio)=ideter(kkio)
                        END DO
                        ideter2(kko)=2
                        ideter2(kok)=1
                        CALL adr(ideter2, iaa2)
                        xmat2=xmat2+vect(kvect)*vect(iaa2-startvect)
                      END IF
                      IF (ideter(kko).eq.2.and.ideter(kok).eq.1) THEN
                        xmat2=xmat2-(1.d0/2.d0)*(vect(kvect)**2)
                        DO kkio=1,natom
                          ideter2(kkio)=ideter(kkio)
                        END DO
                        ideter2(kko)=1
                        ideter2(kok)=2
                        CALL adr(ideter2, iaa2)
                        xmat2=xmat2+vect(kvect)*vect(iaa2-startvect)
                      END IF
                    END IF
                  END DO
                END DO

                DO kko=(natom/2)+1,natom
                  DO kok=kko,natom
                    IF (kok.eq.kko.and.ideter(kok).ne.3) THEN
                      xmat2=xmat2+(3.d0/4.d0)*(vect(kvect)**2)
                    ELSE
                      IF (ideter(kko).eq.1.and.ideter(kok).eq.1) THEN
                        xmat2=xmat2+(1.d0/2.d0)*(vect(kvect)**2)
                      END IF
                      IF (ideter(kko).eq.2.and.ideter(kok).eq.2) THEN
                        xmat2=xmat2+(1.d0/2.d0)*(vect(kvect)**2)
                      END IF
                      IF (ideter(kko).eq.1.and.ideter(kok).eq.2) THEN
                        xmat2=xmat2-(1.d0/2.d0)*(vect(kvect)**2)
                        DO kkio=1,natom
                          ideter2(kkio)=ideter(kkio)
                        END DO
                        ideter2(kko)=2
                        ideter2(kok)=1
                        CALL adr(ideter2, iaa2)
                        xmat2=xmat2+vect(kvect)*vect(iaa2-startvect)
                      END IF
                      IF (ideter(kko).eq.2.and.ideter(kok).eq.1) THEN
                        xmat2=xmat2-(1.d0/2.d0)*(vect(kvect)**2)
                        DO kkio=1,natom
                          ideter2(kkio)=ideter(kkio)
                        END DO
                        ideter2(kko)=1
                        ideter2(kok)=2
                        CALL adr(ideter2, iaa2)
                        xmat2=xmat2+vect(kvect)*vect(iaa2-startvect)
                      END IF
                    END IF
                  END DO
                END DO
  
                DO kko=1,natom/2
                  DO kok=(natom/2)+1,natom
                    IF (kok.eq.kko.and.ideter(kok).ne.3) THEN
                      xmat2=xmat2+(3.d0/4.d0)*(vect(kvect)**2)
                    ELSE
                      IF (ideter(kko).eq.1.and.ideter(kok).eq.1) THEN
                        xmat2=xmat2+(1.d0/2.d0)*(vect(kvect)**2)
                      END IF
                      IF (ideter(kko).eq.2.and.ideter(kok).eq.2) THEN
                        xmat2=xmat2+(1.d0/2.d0)*(vect(kvect)**2)
                      END IF
                      IF (ideter(kko).eq.1.and.ideter(kok).eq.2) THEN
                        xmat2=xmat2-(1.d0/2.d0)*(vect(kvect)**2)
                        DO kkio=1,natom
                          ideter2(kkio)=ideter(kkio)
                        END DO
                        ideter2(kko)=2
                        ideter2(kok)=1
                        CALL adr(ideter2, iaa2)
                        xmat2=xmat2+vect(kvect)*vect(iaa2-startvect)
                      END IF
                      IF (ideter(kko).eq.2.and.ideter(kok).eq.1) THEN
                        xmat2=xmat2-(1.d0/2.d0)*(vect(kvect)**2)
                        DO kkio=1,natom
                          ideter2(kkio)=ideter(kkio)
                        END DO
                        ideter2(kko)=1
                        ideter2(kok)=2
                        CALL adr(ideter2, iaa2)
                        xmat2=xmat2+vect(kvect)*vect(iaa2-startvect)
                      END IF
                    END IF
                  END DO
                END DO
              END IF
!----------------------------------------
!             print *,"norm = ",norm2,"xmat2 = ",xmat2,"vect =",vect(kvect),"natom=",natom
              xymat2=xymat2+xmat2

!           XS=(1.d0/2.d0)*(-1.d0+dsqrt(1.d0+(4.d0*xymat/norm)))
!           XS1=(1.d0/2.d0)*(-1.d0+dsqrt(1.d0+(4.d0*xymat1/norm1)))
            XS2=(1.d0/2.d0)*(-1.d0+dsqrt(1.d0+(4.d0*xymat2/norm2)))
            END DO
!           print *,eigen,xymat2,XS2

!----------------------------------------
! Resta and probabilities of dets
!----------------------------------------
!            DO kvect=1,nbtots
!              ideter2=ideter
!                   CALL getdet(kvect,ideter)
!
!              IF (jstart.eq.1) THEN
!                jstart=0
!                ideter2=ideter
!                DO kko=1,natom
!                  IF (ideter(kko).eq.3) THEN
!                    ipt_1=kko
!                    DO kok=kko+1,natom
!                      IF (ideter(kok).eq.3) THEN
!                        ipt_2=kok
!                          DO okk=kok+1,natom
!                            IF (ideter(okk).eq.3) THEN
!                              ipt_3=okk
!                              EXIT
!                            END IF
!                          END DO
!                        EXIT
!                      END IF
!                    END DO
!                    EXIT
!                  END IF
!                END DO
!              END IF
!  
!
!              DO kko=1,natom
!                IF (ideter(kko).eq.3) THEN
!                  iptemp_1=kko
!                  DO kok=kko+1,natom
!                    IF (ideter(kok).eq.3) THEN
!                      iptemp_2=kok
!                        DO okk=kok+1,natom
!                          IF (ideter(okk).eq.3) THEN
!                            iptemp_3=okk
!                            EXIT
!                          END IF
!                        END DO
!                      EXIT
!                    END IF
!                  END DO
!                  EXIT
!                END IF
!              END DO
!  
!              IF (iptemp_1.ne.ipt_1.or.iptemp_2.ne.ipt_2                        &
!                    .or. iptemp_3.ne.ipt_3) THEN
!                ibougetrou=1
!                ipt_1=iptemp_1
!                ipt_2=iptemp_2
!                ipt_3=iptemp_3
!              ELSE
!                proj_trou=proj_trou+vect(kvect)**2
!                ibougetrou=0
!              END IF
!
!              IF (iptemp_1.eq.(9-iptemp_3+1)                                    &
!                        .and. iptemp_2 .eq. 5 .and. iptemp_1.ne.4) THEN
!                sym=sym+vect(kvect)**2
!              ELSEIF (ideter(4).eq.3 .and. ideter(5).eq.3                       &
!                        .and. ideter(6).eq.3) THEN
!                nonsym=nonsym+vect(kvect)**2
!              ELSE
!!             ELSEIF(ideter(2).eq.3 .and. ideter(3).eq.3                        &
!!                 .and. ideter(6).eq.3)then
!                nonsym=nonsym+vect(kvect)**2
!              END IF
!
!              IF (ideter(2).eq.3 .and. ideter(4).eq.3                           &
!                        .and. ideter(8).eq.3) THEN
!                prob=prob+vect(kvect)**2
!              END IF
!  
!!             IF (ideter(1).eq.3 .and. ideter(2).eq.3                           &
!!                       .and. ideter(3).eq.3) THEN
!!               prob2=prob2+vect(kvect)**2
!!             END IF
!
!              IF (ibougetrou.eq.1.or.kvect.eq.nbtots) THEN
!!----------------------------------------
!! mono
!                DO i=1,natom/2
!                  IF (ideter2(i).eq.3) THEN
!                    delta=1.0d0
!                  ELSE
!                    delta=0.0d0
!                  END IF
!                  resta_mono=resta_mono+delta*xz(i)*xz(i)*proj_trou
!                  resta_one=resta_one+delta*xz(i)*proj_trou
!                END DO
!! bi
!                DO i=1,natom/2
!                  DO j=1,natom/2
!                    IF (ideter2(i).ne.3.or.ideter2(j).ne.3.or.i.eq.j)           &
!                     THEN
!                      delta=0.0d0
!                    ELSE
!                      delta=1.0d0
!                    END IF
!                    resta_bi=resta_bi+delta*xz(i)*xz(j)*proj_trou
!                  END DO
!                END DO
!!----------------------------------------
!                proj_trou=0.d0
!                proj_trou=proj_trou+vect(kvect)**2
!              END IF
!            END DO

!----------------------------------------
! One particle density matrix
!----------------------------------------
!            DO kko=1,3
!              DO kok=1,3
!
!                DO kvect=1,nbtots

!                  CALL getdet(kvect,ideter)
!                  ideter2=ideter
!                  IF (ideter(kko).ne.3) THEN
!                    IF (ideter(kok).eq.3) THEN
!                      ideter2(kok)=ideter(kko)
!                      ideter2(kko)=3
!                      CALL adr(ideter2, iaa2)
!                      densmat(kko,kok)=densmat(kko,kok)+vect(kvect)*        &
!                       vect(iaa2)
!                    END IF
!                  END IF
!                  IF (kko.eq.kok.and.ideter(kko).ne.3) THEN
!                    densmat(kko,kko)=densmat(kko,kko)+vect(kvect)**2
!                  END IF

!                END DO

!              END DO
!            END DO

!----------------------------------------
! two particle density matrix
!----------------------------------------

!           DO kko=1,(natom/2)-1
!             DO kok=kko+1,natom/2

!               nrow=nrow+1
!               ncol=0
!               DO mmo=1,(natom/2)-1
!                 DO mom=mmo+1,natom/2

!                   ncol=ncol+1

!                   DO kvect=1,nbtots
!                     CALL getdet(kvect,ideter)
!                     ideter2=ideter
!                     IF (ideter(kko).eq.3.and.ideter(kok)                      &
!                      .eq.3.and.ideter(mmo).ne.3.and.ideter(mom).ne.3)         &
!                      THEN
!                      if(ideter(kok).ne.3 .and. ideter(mom).ne.3)then
!                       ideter2(kko)=ideter(mmo)
!                       ideter2(mmo)=3
!                       ideter2(kok)=ideter(mom)
!                       ideter2(mom)=3
!                       CALL adr(ideter2, iaa2)
!                       densmat2(nrow,ncol)=densmat2(nrow,ncol)+                &
!                        vect(kvect)*vect(iaa2)
!                           print *,nrow,ncol,kko,kok,mmo,mom,
!    *                     densmat2(nrow,ncol)
!                      endif
!                     END IF


!                     IF (nrow.eq.ncol.and.ideter(mmo)                          &
!                      .ne.3.and.ideter(mom).ne.3) THEN
!                       densmat2(nrow,ncol)=densmat2(nrow,ncol)+                &
!                        vect(kvect)**2

!                     END IF

!                   END DO


!                 END DO
!               END DO

!             END DO
!           END DO

!----------------------------------------

!----------------------------------------
! conduction
!----------------------------------------

!            count=0
!            DO kko=1,(natom/2)-2
!              DO kok=kko+1,(natom/2)-1
!                DO okk=kok+1,natom/2
!  
!                nrow=nrow+1
!                ncol=0
!                DO mmo=kko,kko+1
!                  DO mom=kok,kok+1
!                    DO omm=okk,okk+1
!  
!                    ncol=ncol+1
!                    DO kvect=1,nbtots
!                      CALL getdet(kvect,ideter)
!                      ideter2=ideter
!                      IF (abs(kko-mmo).eq.1.or.abs(kok-mom).eq.1                        &
!                            .or. abs(okk-omm).eq.1) THEN
!                        IF (mmo.le.natom/2.and.mom.le.natom/2 .and.                     &
!                            omm.le.natom/2) THEN
!                          IF (mmo.ne.mom .and. mom.ne.omm) THEN
!                            IF (ideter(kko).eq.3 .and. ideter(kok).eq.3                 &
!                                .and. ideter(okk).eq.3) THEN
!                              ideter2(okk)=ideter2(omm)
!                              ideter2(omm)=3
!                              ideter2(kok)=ideter2(mom)
!                              ideter2(mom)=3
!                              ideter2(kko)=ideter2(mmo)
!                              ideter2(mmo)=3
!                              CALL adr(ideter2, iaa2)
!!                             count=0
!!                             do i=1,natom/2
!!                                   if(ideter2(i).eq.3)then
!!                                     count+=1
!!                                   endif
!!                             enddo
!!                             print *,kko,kok,okk,mmo,mom,omm,iaa2
!                              conduction=conduction+dabs(vect(kvect)*vect(iaa2))
!                            END IF
!                          END IF
!                        END IF
!                      END IF
!                    END DO
!  
!                    END DO
!                  END DO
!                END DO
!  
!                END DO
!              END DO
!            END DO
  
!----------------------------------------

!            DO j=1,ndim
!               write(217,1022)j,(densmat(j,kko),kko=1,ndim)
!            END DO
!----------------------------------------
! diagonalisation de mat
! affiche vecteur
!           JOBZ='V'
! matrice sup
!           UPLO='U'

! matrice en vecteur ligne ...
!           extradiag_dmat2=0d0
!           k=0
!           DO j=1,ndim
!             DO i=1,j-1
!               if(i.ne.j)then
!                   extradiag_dmat2 = extradiag_dmat2 + dabs(densmat2(i,j))
!               endif
!             END DO
!           END DO

! appel subroutine LAPACK de diagonalisation :: double précision !!
!           INFO=0
!           CALL DSPEV (JOBZ, UPLO, ndim, AP, eigenvalues, eigenvectors,            &
!             ndim, WORK, INFO)

!           IF (INFO.ne.0) THEN
!             PRINT *,'SUBROUTINE MATRIX: Error at dspev',info
!             CALL exit (1)
!           END IF

!           proj_2=0.d0
!           sum=0d0
!           DO j=1,ndim
!             proj_2=proj_2-eigenvalues(j)*log(eigenvalues(j))
!             sum+=eigenvalues(j)
!             write(214,*)eigenvalues(j)
!           END DO

!           XS=(1.d0/2.d0)*(-1.d0+dsqrt(1.d0+(4.d0*xymat/norm)))
!           XS2=(1.d0/2.d0)*(-1.d0+dsqrt(1.d0+(4.d0*xymat2/norm2)))
!           XS3=(1.d0/2.d0)*(-1.d0+dsqrt(1.d0+(4.d0*xymat3/norm3)))
!           WRITE (217,*) eigen,XS,norm
      END IF

      END DO



   10 FORMAT (E25.0)
 1022 FORMAT(3x,I3,6(2x,F12.4))
      END
