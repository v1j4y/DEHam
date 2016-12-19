subroutine conv(ideter,det,deth)
    implicit none
    BEGIN_DOC
    ! this routine converts a detrminant in the old
    ! format into the new one and returns the determinant.
    END_DOC
    integer,INTENT(INOUT)::ideter(natomax)
    integer(kind=selected_int_kind(16)),INTENT(INOUT)::det
    integer(kind=selected_int_kind(16)),INTENT(INOUT)::deth
    integer::i
    det=0
    deth=0
    do i=1,natom
        if(ideter(natom-i+1).eq.2 .and. ideter(natom-i+1).ne.3)then
            det=IBSET(det,i-1)
        elseif(ideter(natom-i+1).eq.3)then
            deth=IBSET(deth,i-1)
        endif
    enddo
end
