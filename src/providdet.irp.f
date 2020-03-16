use iso_c_binding
BEGIN_PROVIDER[integer(C_SIZE_T),det,(nt2,2)]
&BEGIN_PROVIDER[integer(C_SIZE_T),deth,(nt1,2)]
    BEGIN_DOC
    ! provides det and deth array
    END_DOC
    use iso_c_binding
    implicit none
!   integer(kind=selected_int_kind(16))::dethsh
    integer(C_SIZE_t)::a
    integer(C_SIZE_T)::i,count,iik
    integer::const
    i=1
    a=0
    const=0
    count=0

    If(ntrou.ge.1)then

            const=0
            do while (i.le.(nt1))
                i+=1
                a+=1
                do while(popcnt(a).ne.ntrou)
                    a+=1
                enddo
                count+=1
                if(FAM1) then
                    if(fix_trou1 .eq. fix_trou2) then
                        deth(count,1)=ISHFT(a,natom/2)
                    else
                        deth(count,1)=ISHFT(a,natom - (fix_trou2-fix_trou1))
                    endif
                else 
                    deth(count, 1) = a
                endif
                deth(count,2)=i-1
            enddo
    endif

    !C if det=0 then exit
    a=0
    i=0
    count=0
    if(mpiid==0)then
      print *,'nt2=',nt2,'nbeta=',nbeta
    endif
    do while (i.lt.(nt2))

        i+=1
        a+=1
        do while(popcnt(a).ne.nbeta)
            if(nbeta.eq.0)then
                a=0
                count+=1
                det(count,1)=a
                det(count,2)=i
                EXIT
            endif
            a+=1
        enddo
        
        if(nbeta.ne.0)then
        count+=1
        det(count,1)=a
        det(count,2)=i
        endif
!       write(6,16)a,a,i
    enddo
!   do i=0,nt2
!     print *,i," det=",det(i,1)," add=",det(i,2)
!   enddo
    

10  FORMAT(B64,I8,F8.2)
15  FORMAT(B64,I8,I8,I8)
11  FORMAT(B64,I3,B64)
12  FORMAT(I5,$)
13  FORMAT(B64,B64)
14  FORMAT(B64,I8)
16  FORMAT(B64,I8,I8)
END_PROVIDER
