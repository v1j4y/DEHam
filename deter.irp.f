BEGIN_PROVIDER [integer, deter, (natomax)]

    implicit none
    BEGIN_DOC
    ! provides ideter and iaa
    END_DOC
    integer ::ideter(natomax)
    deter=0
END_PROVIDER
