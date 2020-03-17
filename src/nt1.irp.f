use iso_c_binding
BEGIN_PROVIDER [integer(C_SIZE_T), nt1]
    BEGIN_DOC
    ! calculates the number of det the 3's moving
    END_DOC
    implicit none
    integer(C_SIZE_T)::natom2

!   call combin(idet1(1,nt1+1),natom,ntrou,nt1,32,jrangmax)
    natom2=natom
    if(FAM1) then
        if(fix_trou1 .eq. fix_trou2) then
            natom2=natom/2
        else
            natom2 = fix_trou2 - fix_trou1
        endif
    endif
    nt1=   nint(gamma(1.0*(natom2+1))/(gamma(1.0*(natom2-ntrou+1))*gamma(1.0*(ntrou+1))),selected_int_kind(16))
    if(mpiid==0)then
      write(6,*)'nt1',nt1
    endif
END_PROVIDER
