subroutine searchdetfull()
    BEGIN_DOC
    ! this subroutine is at the heart of the idea
    ! it will generate all the determinants in a fixed order
    ! then find the posistion of the determinant given and
    ! return it's position in add.
    END_DOC
!   integer(kind=selected_int_kind(16)),INTENT(INOUT)::foundetadr(maxlien,4)
    integer(kind=selected_int_kind(16))::add
!   integer(kind=selected_int_kind(16)),INTENT(INOUT)::deth
    integer(kind=selected_int_kind(16))::dethsh
    integer(kind=selected_int_kind(16))::addh
    integer(kind=selected_int_kind(16))::a
    integer(kind=selected_int_kind(16))::i
    integer::const,count
    i=1
    a=0
    add=0
    const=0
    count=1
    If(ntrou.ge.1)then

            const=0
            a=0
            addh=0
            i=1
            dethsh = ISHFT(foundaddh(count,1),-natom/2)
            do while (i.le.(2*nt1))
                if(a.eq.dethsh)then
                    addh=i-1
                    foundaddh(count,2)=addh
                    count+=1
                    dethsh = ISHFT(foundaddh(count,1),-natom/2)
                    do while(count .le. detfound .and. a .eq. dethsh)
                        foundaddh(count,2)=addh
                        count+=1
                        dethsh = ISHFT(foundaddh(count,1),-natom/2)
                    enddo
                    if(count.gt.detfound)then
                    EXIT
                    endif
                endif

                i+=1
!               const=1
                a+=1
!               do while(popcnt(a).ne.ntrou .or. const==1)
                do while(popcnt(a).ne.ntrou)
                    a+=1
!                   const=0
                enddo
            enddo

            if(a.eq.foundaddh(count,1))then
            addh=i-1
            foundaddh(count,2)=addh
                    count+=1
                    dethsh = ISHFT(foundaddh(count,1),-natom/2)
                    do while(count .le. detfound .and. a .eq. foundaddh(count,1))
                        foundaddh(count,2)=addh
                        count+=1
                        dethsh = ISHFT(foundaddh(count,1),-natom/2)
                    enddo
            endif

    endif

    !C if det=0 then exit
    a=0
    i=1
    count=1
    const=0
    if(a.eq.foundadd(count,1))then
    add=1
               foundadd(count,2)=add
               count+=1
               do while(foundadd(count,1).eq.a .and. count.le.detfound)
                foundadd(count,2)=add
                count+=1
               enddo
    endif

    do while (i.le.(nt2) .and. count .le.detfound)
        if(a.eq.foundadd(count,1))then
            if(a.eq.1)then
               add=i-1
               foundadd(count,2)=add
               count+=1
               do while(foundadd(count,1).eq.a .and. count.le.detfound)
                foundadd(count,2)=add
                count+=1
               enddo
               if(count.eq.detfound)then
               const=-1
               EXIT
               endif
            else
               add=i-1
               foundadd(count,2)=add
               count+=1
               do while(foundadd(count,1).eq.a .and. count.le.detfound)
                foundadd(count,2)=add
                count+=1
               enddo
               if(count.gt.detfound)then
               const=-1
               EXIT
               endif
            endif
        endif

        i+=1
        a+=1
!C      write(6,16)a,a,i-2
        do while(popcnt(a).ne.nbeta)
            a+=1
        enddo
    enddo
!   if(a.eq.foundadd(count,1) .and. foundadd(count,1).eq.0)then
    if(a.eq.foundadd(count,1) .and. foundadd(count,2) .eq. 0)then
        foundadd(count,2) = nt2
        count+=1
               do while(foundadd(count,1).eq.a .and. count.le.detfound)
                foundadd(count,2)=nt2
                count+=1
               enddo
    endif
!   do i=1,detfound
!       write(6,16)foundadd(i,1) ,foundadd(i,2),foundadd(i,3)
!       write(6,16)foundaddh(i,1),foundaddh(i,2),foundaddh(i,3)
!   enddo

10  FORMAT(B64,I8,F8.2)
15  FORMAT(B64,I8,I8,I8)
11  FORMAT(B64,I3,B64)
12  FORMAT(I5,$)
13  FORMAT(B64,B64)
14  FORMAT(B64,I8)
16  FORMAT(B64,I8,I8)
end
