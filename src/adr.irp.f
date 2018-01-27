subroutine adr(ideter,add)
    implicit none
    BEGIN_DOC
    ! this subroutine provides the address of a detrminant 
    ! given in old format.
    ! It searches in a list of generated determinants and
    ! matches the given determinant.
    END_DOC
    integer,INTENT(INOUT)::ideter(natomax)
    integer(kind=selected_int_kind(16)),INTENT(INOUT)::add
    integer(kind=selected_int_kind(16))::deti,dethi,addh,detnew
    integer::count,i,j

    deti=0
    detnew=0
    dethi=0
    count=0
    call conv(ideter,deti,dethi)
    Do i=0,natom-1
        if(BTEST(dethi,i))then
            count=count+1
        endif
        if(BTEST(deti,i))then
            detnew=IBSET(detnew,i-count)
        endif
    enddo
    deti=detnew
    call searchdet(deti,add,dethi,addh)
        add = add + (nt1-addh)*(nt2)


10  FORMAT(B64,I8,F8.2)
15  FORMAT(B64,I8,I8,I8)
11  FORMAT(B64,I3,B64)
12  FORMAT(I5,$)
13  FORMAT(B64,B64)
14  FORMAT(B64,I14)
16  FORMAT(B64,I14,I14)
end
