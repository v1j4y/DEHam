!   SUBROUTINE ylogic(ideter,yalt,ytrou,yrep1)
BEGIN_PROVIDER [logical, yalt,(maxlien)]
&BEGIN_PROVIDER [logical, ytrou,(maxlien)]
&BEGIN_PROVIDER [logical, yrep1,(maxlien)]

    implicit none
    BEGIN_DOC
    ! provides the important logical units
    ! for the presence of holes (ytrou) and
    ! an alternance (yalt) or a rep (yrep)
    END_DOC

        integer :: i,il
        integer :: ik1,ik2
        integer :: ntr,nalt,naltmax
        integer :: itypl(maxlien)
        logical :: yw

!       write(6,*)'in elem_diag'

!---mise a zero
            nalt=0
            ntr=0
            yw=.FALSE.
            do il=1,nlientot
               itypl(il)=0
               yalt(il)=.false.
               ytrou(il)=.false.
               yrep1(il)=.false.
            enddo
!---

        if(yw)write(6,*)'ideter',(deter(il),il=1,natom)

        do il=1,nlientot

            ik1=iliatom1(il)
            ik2=iliatom2(il)

            if(yw)write(6,*)'ik1',ik1,'ik2',ik2
	  if(deter(ik1).eq.deter(ik2).and. &
            deter(ik1).eq.3)then
            yrep1(il)=.true.
	  endif

          if(deter(ik1).ne.deter(ik2))then
	     if(deter(ik1).eq.1.and.deter(ik2).eq.2)then
	       yalt(il)=.true.
	       nalt=nalt+1
	       itypl(il)=1
!       goto 32
                CYCLE
	     endif
	     if(deter(ik1).eq.2.and.deter(ik2).eq.1)then
	       yalt(il)=.true.
	       nalt=nalt+1
	       itypl(il)=2
	       goto 32
                CYCLE
	     endif
 32          continue
	     if(deter(ik1).eq.1.and.deter(ik2).eq.3)then
               ytrou(il)=.true.
	       ntr=ntr+1
	       itypl(il)=3
!       goto 33
                CYCLE
             endif
	     if(deter(ik1).eq.3.and.deter(ik2).eq.1)then
               ytrou(il)=.true.
	       ntr=ntr+1
	       itypl(il)=4
!       goto 33
                CYCLE
             endif
	     if(deter(ik1).eq.2.and.deter(ik2).eq.3)then
               ytrou(il)=.true.
	       ntr=ntr+1
	       itypl(il)=5
!       goto 33
                CYCLE
             endif
	     if(deter(ik1).eq.3.and.deter(ik2).eq.2)then
               ytrou(il)=.true.
	       ntr=ntr+1
	       itypl(il)=6
!       goto 33
                CYCLE
             endif
	  endif
!33    continue
        enddo

END_PROVIDER
