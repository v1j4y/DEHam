subroutine getdet(add,ideter)
    implicit none
    BEGIN_DOC
    ! this routing gives the determinant in
    ! the traditional form given its address
    END_DOC
    integer,INTENT(INOUT)::ideter(natomax)
    integer(kind=selected_int_kind(16)),INTENT(IN)::add
    integer(kind=selected_int_kind(16))::deta,detb
    integer::i,const,ia,ib

    ib = MOD(add,nt2)
    if(MOD(add,nt2).eq.0)then
        ib=nt2
    endif
    ia = (add-ib)/nt2
    ia = nt1 - ia
    ideter=1
    const=0
    detb=0
    deta=0
    i=1
    do while (i.le.(ib))
        const=1
        do while(popcnt(detb).ne.nbeta .or. const==1)
            if(nbeta.eq.0)then
                detb=0
                EXIT
            endif
            detb+=1
            const=0
        enddo
        i+=1
!       write(6,14)detb,detb
    enddo
    i=1
    do while (i.le.(ia))
        const=1
        do while(popcnt(deta).ne.ntrou .or. const==1)
            deta+=1
            const=0
        enddo
        i+=1
!       write(6,14)deta,deta
    enddo
    const=0
    do i=0,(natom/2) - 1
        if(BTEST(deta,i))then
            ideter((natom/2)-i)=3
        endif
    enddo
    do i=0,natom-1
        if(ideter(natom-i).eq.1)then
            if(BTEST(detb,const))then
                ideter(natom-i)=2
            endif
            const=const+1
        endif
    enddo

return
10  FORMAT(B64,I8,F8.2)
15  FORMAT(B64,I8,I8,I8)
11  FORMAT(B64,I3,B64)
12  FORMAT(I5,$)
13  FORMAT(B64,B64)
14  FORMAT(B64,I8)
16  FORMAT(B64,I8,I8)
end subroutine getdet
