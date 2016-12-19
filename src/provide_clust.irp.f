BEGIN_PROVIDER[integer,l1,   (maxlien)]
&BEGIN_PROVIDER[integer,l2,  (maxlien)]
&BEGIN_PROVIDER[integer,ktyp,(maxlien)]
&BEGIN_PROVIDER [real*8, xtt  ,(maxlien)]
&BEGIN_PROVIDER[real*8, xjjz ,(maxlien)]
&BEGIN_PROVIDER[real*8, xjjxy,(maxlien)]
&BEGIN_PROVIDER [integer, ntrou]
&BEGIN_PROVIDER [integer, isz]
      implicit none
!     integer::i
!       open(unit=11,file="l1.dat",form="formatted")
!       open(unit=12,file="l2.dat",form="formatted")
!       open(unit=13,file="ktyp.dat",form="formatted")
!       rewind(11)
!       rewind(12)
!       rewind(13)
!       do i=1,42
!       read(11,*)l1(i)
!       print *,l1(i)
!       read(12,*)l2(i)
!       read(13,*)ktyp(i)
!       enddo
!       close(11)
!       close(12)
!       close(13)
END_PROVIDER
