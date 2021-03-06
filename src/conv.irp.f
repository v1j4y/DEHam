subroutine conv(ideter,deti,dethi)
    use iso_c_binding
    implicit none
    BEGIN_DOC
    ! this routine converts a detrminant in the old
    ! format into the new one and returns the determinant.
    END_DOC
    integer,INTENT(INOUT)::ideter(natomax)
    integer(C_SIZE_T),INTENT(INOUT)::deti
    integer(C_SIZE_T),INTENT(INOUT)::dethi
    integer::i
    deti=0
    dethi=0
    do i=1,natom
        if(ideter(natom-i+1).eq.2 .and. ideter(natom-i+1).ne.3)then
            deti=IBSET(deti,i-1)
        elseif(ideter(natom-i+1).eq.3)then
            dethi=IBSET(dethi,i-1)
        endif
    enddo
end
