subroutine sort()
    implicit none
    integer::i,j,ord,ordh
    integer(kind=selected_int_kind(16))::add,addh,det,deth,addt
    
    do i=1,detfound-1
        do j=i+1,detfound
            if(foundaddh(i,1).gt.foundaddh(j,1))then
                deth = foundaddh(i,1)
                foundaddh(i,1) = foundaddh(j,1)
                foundaddh(j,1) = deth
                addh = foundaddh(i,2)
                foundaddh(i,2) = foundaddh(j,2)
                foundaddh(j,2) = addh
                ordh = foundaddh(i,3)
                foundaddh(i,3) = foundaddh(j,3)
                foundaddh(j,3) = ordh
            endif
            if(foundadd(i,1).gt.foundadd(j,1))then
                det = foundadd(i,1)
                foundadd(i,1) = foundadd(j,1)
                foundadd(j,1) = det
                add = foundadd(i,2)
                foundadd(i,2) = foundadd(j,2)
                foundadd(j,2) = add
                ord = foundadd(i,3)
                foundadd(i,3) = foundadd(j,3)
                foundadd(j,3) = ord
            endif
        enddo
    enddo
end
