use iso_c_binding
BEGIN_PROVIDER[integer, foundet,(natomax,maxlien)]
&BEGIN_PROVIDER[integer(C_SIZE_T), foundetadr,(maxlien)]
&BEGIN_PROVIDER[real, foundetdmat,(maxlien)]
&BEGIN_PROVIDER[integer(C_SIZE_T), foundadd,(maxlien,3)]
&BEGIN_PROVIDER[integer(C_SIZE_T), foundaddh,(maxlien,3)]
&BEGIN_PROVIDER[integer, detfound]
    BEGIN_DOC
    ! provides all found determinants
    END_DOC
    detfound=0
    founddet=0
    foundetdmat=0d0
    founddetadr=0
    foundadd=0
    foundaddh=0
END_PROVIDER
