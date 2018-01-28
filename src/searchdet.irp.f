subroutine searchdet(deti,add,dethi,addh)
    use iso_c_binding
    BEGIN_DOC
    ! this subroutine is at the heart of the idea
    ! it will generate all the determinants in a fixed order
    ! then find the posistion of the determinant given and
    ! return it's position in add.
    END_DOC
    integer(C_SIZE_T),INTENT(INOUT)::deti
    integer(C_SIZE_T),INTENT(INOUT)::add
    integer(C_SIZE_T),INTENT(INOUT)::dethi
    integer(C_SIZE_T),INTENT(INOUT)::addh
    integer(C_SIZE_T)::dethsh
    integer(C_SIZE_T)::a
    integer(C_SIZE_T)::i,j
    integer::count
    logical::found
    
    i=1
    j=nt1
    found=.FALSE.
    do while(.not.found)
        if(deth((i+j)/2,1).eq.dethi)then
            addh=deth((i+j)/2,2)
            found=.TRUE.
            EXIT
        elseif (abs(i-j).eq.1)then
            if(deth(i,1).eq.dethi)then
            addh=deth(i,2)
            elseif(deth(j,1).eq.dethi)then
            addh=deth(j,2)
            endif

            found=.TRUE.
            EXIT
        else
            if(deth((i+j)/2,1).gt.dethi)then
                j=(i+j)/2
            else
                i=(i+j)/2
            endif
        endif
    enddo

    i=1
    j=nt2
    found=.FALSE.
    do while(.not.found)
        if(det((i+j)/2,1).eq.deti)then
            add=det((i+j)/2,2)
            found=.TRUE.
            EXIT
        elseif (abs(i-j).eq.1)then
            if(det(i,1).eq.deti)then
            add=det(i,2)
            elseif(det(j,1).eq.deti)then
            add=det(j,2)
            endif
            found=.TRUE.
            EXIT
        else
            if(det((i+j)/2,1).gt.deti)then
                j=(i+j)/2
            else
                i=(i+j)/2
            endif
        endif
    enddo

10  FORMAT(B64,I8,F8.2)
15  FORMAT(B64,I8,I8,I8)
11  FORMAT(B64,I3,B64)
12  FORMAT(I5,$)
13  FORMAT(B64,B64)
14  FORMAT(B64,I8)
16  FORMAT(B64,I8,I8)
end
