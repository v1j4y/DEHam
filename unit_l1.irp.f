        subroutine unit_l1(tl1,&
        tl2,                   &
        tktyp,                 &
        tistart,               &
        tnrows,                 &
        txjjxy,                &
        txjjz ,                &
        txtt  ,                &
        tcountcol,             &
        tntrou,                &
        tisz,                  &
        tcol,tval)
        implicit none
        integer,INTENT(INOUT)::tistart, tnrows, tntrou, tisz
        integer::i
        real*8,INTENT(INOUT)::tval(maxlien)
        integer(kind=selected_int_kind(16)),INTENT(INOUT)::tcol(maxlien)
        integer(kind=selected_int_kind(16)),INTENT(INOUT),dimension(tnrows)::tcountcol
        integer(kind=selected_int_kind(16)),INTENT(INOUT)::tl1(maxlien),tl2(maxlien),tktyp(maxlien)
        real*8,INTENT(INOUT)::txtt(maxlien),txjjz(maxlien),txjjxy(maxlien)
   
		nrows = tnrows
		provide nrows
        do i=1,maxlien
        l1(i)=tl1(i)
        l2(i)=tl2(i)
        ktyp(i)=tktyp(i)
        xtt(i)     = txtt(i)
        xjjxy(i)   = txjjxy(i)
        xjjz (i)   = txjjz (i)
        enddo
        ntrou = tntrou
        isz = tisz
        tcol=0
        tval=0d0
        provide l1 l2 ktyp xtt xjjxy xjjz ntrou
!print *,"l1"
!print *,l1
!print *,"xjjz"
!print *,xjjz
        call unit(tistart, tcountcol,tcol,tval)

        end


