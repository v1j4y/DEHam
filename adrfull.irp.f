subroutine adrfull()
    implicit none
    BEGIN_DOC
    ! this subroutine provides the address of a detrminant 
    ! given in old format.
    ! It searches in a list of generated determinants and
    ! matches the given determinant.
    END_DOC
    integer,dimension(natomax)::ideter
    integer(kind=selected_int_kind(16))::add
    integer(kind=selected_int_kind(16))::det,deth,addh,detnew
    integer::count,i,j

    det=0
    detnew=0
    deth=0
    do j=1,detfound
    detnew=0
    count=0
    ideter=foundet(:,j)
    call conv(ideter,det,deth)
    Do i=0,natom-1
        if(BTEST(deth,i))then
            count=count+1
        endif
        if(BTEST(det,i))then
            detnew=IBSET(detnew,i-count)
        endif
    enddo
    det=detnew
        foundadd(j,1)=det
        foundadd(j,3)=j
        foundaddh(j,1)=deth
        foundaddh(j,3)=j
    enddo
    call sort()
    call searchdetfull()
    call desort()
    do i=1,detfound
        add = foundadd(i,2)
        addh = foundaddh(i,2)
        add = add + (nt1-addh)*(nt2)
        foundetadr(i)=add
    enddo


10  FORMAT(B64,I8,F8.2)
15  FORMAT(B64,I8,I8,I8)
11  FORMAT(B64,I3,B64)
12  FORMAT(I5,$)
13  FORMAT(B64,B64)
14  FORMAT(B64,I14)
16  FORMAT(B64,I14,I14)
end
