subroutine adrfull()
    implicit none
    use iso_c_binding
    BEGIN_DOC
    ! this subroutine provides the address of a detrminant 
    ! given in old format.
    ! It searches in a list of generated determinants and
    ! matches the given determinant.
    END_DOC
    integer,dimension(natomax)::ideter
    integer(C_SIZE_T)::add
    integer(C_SIZE_T)::deti,dethi,addh,detnew
    integer::count,i,j

    deti=0
    detnew=0
    dethi=0
    do j=1,detfound
        detnew=0
        count=0
        ideter=foundet(:,j)
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
        foundadd(j,1)=deti
        foundadd(j,3)=j
        foundaddh(j,1)=dethi
        foundaddh(j,3)=j
        call searchdet(deti,add,dethi,addh)
!   enddo
!   call sort()
!   call searchdetfull()
!   call desort()
!   do i=1,detfound
!       add = foundadd(i,2)
!       addh = foundaddh(i,2)
        foundadd(j,2) = add
        foundaddh(j,2)= addh
        add = add + (nt1-addh)*(nt2)
        foundetadr(j)=add
    enddo


10  FORMAT(B64,I8,F8.2)
15  FORMAT(B64,I8,I8,I8)
11  FORMAT(B64,I3,B64)
12  FORMAT(I5,$)
13  FORMAT(B64,B64)
14  FORMAT(B64,I14)
16  FORMAT(B64,I14,I14)
end
