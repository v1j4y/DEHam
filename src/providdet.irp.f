BEGIN_PROVIDER[integer(kind=selected_int_kind(16)),det,(nt2,2)]
&BEGIN_PROVIDER[integer(kind=selected_int_kind(16)),deth,(nt1,2)]
    BEGIN_DOC
    ! provides det and deth array
    END_DOC
    implicit none
!   integer(kind=selected_int_kind(16))::dethsh
    integer(kind=selected_int_kind(16))::a
    integer(kind=selected_int_kind(16))::i,count
    integer::const
    i=1
    a=0
    const=0
    count=0

    If(ntrou.ge.1)then

            const=0
!           dethsh = ISHFT(deth,-natom/2)
!           i=nt1
            do while (i.le.(nt1))
!               if(a.eq.dethsh)then
!                   addh=i-1
!                   EXIT
!               endif

                i+=1
                a+=1
                do while(popcnt(a).ne.ntrou)
                    a+=1
                enddo
                count+=1
                if(FAM1) then
                    deth(count,1)=ISHFT(a,natom/2)
                else 
                    deth(count, 1) = a
                endif
                deth(count,2)=i-1
!       write(6,16)ISHFT(a,natom/2),ISHFT(a,natom/2),i-1
            enddo
!           if(a.eq.dethsh )then
!           count+=1
!           deth(1,1)=ISHFT(a,natom/2)
!           deth(1,2)=nt1
!           endif

    endif

    !C if det=0 then exit
    a=0
    i=0
    count=0
    print *,'nt2=',nt2,'nbeta=',nbeta
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
    

10  FORMAT(B64,I8,F8.2)
15  FORMAT(B64,I8,I8,I8)
11  FORMAT(B64,I3,B64)
12  FORMAT(I5,$)
13  FORMAT(B64,B64)
14  FORMAT(B64,I8)
16  FORMAT(B64,I8,I8)
END_PROVIDER
