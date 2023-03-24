        subroutine elem_diag(xmatd)

        implicit none

        integer :: i, j
        real*8,intent(inout)  :: xmatd
        real*8           :: xrep
        logical :: yw

!       write(6,*)'in elem_diag'
!---

!       if(yw)write(6,*)'ideter',(det(iaa,i),i=1,natom)

       yw=.FALSE.
       xmatd=0.d0
       xrep = 0.00d0
       !xrep = 1.25d0
!      testing code
!      xmatd=1000.d0
!      if(.not.redo)write(6,*)'vijayyves'
       do i=1,nlientot
          if(yalt(i))then
!           xmatd=-(xj1+xeneJ(i)*xbJ)+xmatd
            xmatd= -(xjz(i))+xmatd
	          if(yw)write(6,*)xmatd,'xmatd'
	          if(yw)write(6,*)'xjz',xjz(i)
          endif
!         if(yrep1(i))then
!              xmat=xv1+xmat
!       if(yw)write(6,*)iaa,'diag,v1'
!         endif
       enddo
       do i=1, natom
            if(deter(i).ne.3) xmatd = xmatd + E(i)
       enddo
       do i=1, natom-1
          if(deter(i).eq.3 .and. deter(i+1).eq.3) xmatd = xmatd - xrep*xtt(1)
       enddo
       xmatd = xmatd - E(natom+1)
       if(yw)write(6,*)'xmatd=',xmatd

!-----stockage de l element diag

!      imat=iaa
!      imat3(isto3+1)=imat
!      jmat3(isto3+1)=imat
!      dmat3(isto3+1)=xmat/2

    return
    end



