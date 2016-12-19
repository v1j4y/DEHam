BEGIN_PROVIDER [integer, rank]
&BEGIN_PROVIDER [integer, rank_16]
    BEGIN_DOC
    ! calculates the rank of matrix
    END_DOC

    implicit none

    rank=nt1*nt2
    if(MOD(rank,16).eq.0)then
        rank_16=rank
    else
        rank_16=rank+16-MOD(rank,16)
    endif
END_PROVIDER
