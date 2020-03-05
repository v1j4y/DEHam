use iso_c_binding
BEGIN_PROVIDER [integer(C_SIZE_T), nt2]
    BEGIN_DOC
    ! calculates the number of det the 1's moving
    END_DOC


!   call combin(idet2(1,nt2+1),natrest,ial0,nt2,32,jrangmax)
    nt2=   nint(gamma(1.0*(natom-ntrou+1))/((gamma(1.0*(nalpha+1))*gamma(1.0*(nbeta+1)))),selected_int_kind(16))
    print *,"nt2=",nt2
END_PROVIDER
