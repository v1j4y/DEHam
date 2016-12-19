subroutine searchdet(det,add,deth,addh)
    BEGIN_DOC
    ! this subroutine is at the heart of the idea
    ! it will generate all the determinants in a fixed order
    ! then find the posistion of the determinant given and
    ! return it's position in add.
    END_DOC
    integer(kind=selected_int_kind(16)),INTENT(INOUT)::det
    integer(kind=selected_int_kind(16)),INTENT(INOUT)::add
    integer(kind=selected_int_kind(16)),INTENT(INOUT)::deth
    integer(kind=selected_int_kind(16)),INTENT(INOUT)::addh
    integer(kind=selected_int_kind(16))::dethsh
    integer(kind=selected_int_kind(16))::a
    integer(kind=selected_int_kind(16))::i
    integer::const
    i=1
    a=0
    add=0
    const=0

    If(ntrou.ge.1)then

            const=0
            dethsh = ISHFT(deth,-natom/2)
            addh=0
!           i=nt1
            do while (i.le.(nt1))
                if(a.eq.dethsh)then
                    addh=i-1
                    EXIT
                endif

                i+=1
                a+=1
                do while(popcnt(a).ne.ntrou)
                    a+=1
                enddo
            enddo
            if(a.eq.dethsh .and. addh.eq.0)then
            addh=nt1
            endif

    endif

    !C if det=0 then exit
    a=0
    i=0
    count=0
    if(a.eq.det)then
    add=1
    Return
    endif

    do while (i.le.(nt2))
        if(a.eq.det)then
            if(a.eq.1)then
               add=i
               count=-1
               EXIT
            else
               add=i
               count=-1
               EXIT
            endif
        endif

        i+=1
        a+=1
!C      write(6,16)a,a,i-2
        do while(popcnt(a).ne.nbeta)
            a+=1
        enddo
    enddo
    if(a.eq.det .and. count.ne.-1)then
    add=i-1
    endif
    

10  FORMAT(B64,I8,F8.2)
15  FORMAT(B64,I8,I8,I8)
11  FORMAT(B64,I3,B64)
12  FORMAT(I5,$)
13  FORMAT(B64,B64)
14  FORMAT(B64,I8)
16  FORMAT(B64,I8,I8)
end
