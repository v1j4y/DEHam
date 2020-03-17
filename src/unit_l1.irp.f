        subroutine unit_l1(tl1,&
        tl2,                   &
        tktyp,                 &
        tistart,               &
        tnrows,                 &
        txjjxy,                &
        txjjz ,                &
        txtt  ,                &
        tE    ,                &
        tcountcol,             &
        tntrou,                &
        tisz,                  &
        tfix_trou1,            &
        tfix_trou2,            &
        tfam1,                 &
        tmpiid,                &
        tcol,tval)
        use iso_c_binding
        implicit none
        integer,INTENT(INOUT)::tistart, tnrows, tmpiid
        integer,INTENT(INOUT)::tntrou, tisz
        integer,INTENT(INOUT)::tfix_trou1, tfix_trou2
        logical*1,INTENT(INOUT)::tfam1
        integer::i
        real*8,INTENT(INOUT)::tval(maxlien)
        integer(C_SIZE_T),INTENT(INOUT)::tcol(maxlien)
        integer(C_SIZE_T),INTENT(INOUT),dimension(tnrows)::tcountcol
        integer(C_SIZE_T),INTENT(INOUT)::tl1(maxlien),tl2(maxlien),tktyp(maxlien)
        real*8,INTENT(INOUT)::txtt(maxlien),txjjz(maxlien),txjjxy(maxlien), tE(maxlien)
   
		nrows = tnrows
		provide nrows
        do i=1,maxlien
        l1(i)=tl1(i)
        l2(i)=tl2(i)
        ktyp(i)=tktyp(i)
        xtt(i)     = txtt(i)
        xjjxy(i)   = txjjxy(i)
        xjjz (i)   = txjjz (i)
        E (i)   = tE (i)
        enddo
        ntrou = tntrou
        isz = tisz
        FAM1 = tfam1
        fix_trou1 = tfix_trou1
        fix_trou2 = tfix_trou2
        tcol=0
        tval=0d0
        mpiid=tmpiid
        provide l1 l2 ktyp xtt xjjxy xjjz ntrou mpiid
!print *,"l1"
!print *,l1
!print *,"xjjz"
!print *,xjjz
!print *,FAM1
        call unit(tistart, tcountcol,tcol,tval)

        end


