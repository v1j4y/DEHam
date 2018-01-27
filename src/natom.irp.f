BEGIN_PROVIDER [integer, natom]
&BEGIN_PROVIDER [integer, natrest]
&BEGIN_PROVIDER [integer, ial0]
&BEGIN_PROVIDER [logical*1, yham]
&BEGIN_PROVIDER [integer, nlientot]
&BEGIN_PROVIDER [real*8, xt,(maxlien)]
&BEGIN_PROVIDER [real*8 , xjz,(maxlien)]
&BEGIN_PROVIDER [real*8 , xv1]
&BEGIN_PROVIDER [integer, iliatom1,(maxlien)]
&BEGIN_PROVIDER [integer, iliatom2,(maxlien)]
&BEGIN_PROVIDER [integer, ivoisl1,(6,maxlien)]
&BEGIN_PROVIDER [real*8, Emin]
&BEGIN_PROVIDER [real*8, Emax]
&BEGIN_PROVIDER [integer, M0]
&BEGIN_PROVIDER [integer, nalpha]
&BEGIN_PROVIDER [integer, nbeta]
    BEGIN_DOC
    ! read data
    END_DOC

    implicit none
      logical*1 ysuiv
      logical*1 yec
      logical*1 yw
      logical y(natomax)
      logical ykl(maxlien)
      logical yplac2(maxlien)
      logical yplac(maxlien)
!   integer,allocatable :: l1(:),l2(:),ktyp(:)
    integer :: ityp
!C  real*8,allocatable :: xjjz(:),xjjxy(:)
!C  real*8,allocatable :: xtt(:)
    integer :: iprec,maxda2,niter,ngao,nvec,numero,nes4
    integer :: ilien(natomax,natomax)
    integer :: nlien,ilien2(natomax,natomax)
    real*8 xjxy(maxlien)
    integer :: nvois1(natomax),ikx(maxlien)
    integer :: iech(natomax,natomax)
    integer :: ltyp(maxlien),nvois(natomax),ivois1(6,natomax)
    integer :: nvoisl1(maxlien)
    integer :: i,ik,ikk,ikl,ik1,ik2,iiki,j,il,il2,il3,il4
    integer :: kko,kkio1,kkio2,kkio1ou,klo1,l,nlientot2
    integer :: klo2,kat1,kat2,ii,jclust,kl,jk,nplac
    integer :: jplac(2,maxlien),iplac(4,maxplac)
      real*8 xj1
      real*8 xj2
      real*8 xt1
      real*8 xt2
      real*8 xv2
!     real*8 xv1
      real*8 xv3
      real*8 xbJ
      real*8 xbT
      real*8 xspar
      real*8 xsperp
      real*8 xseuil
      real*8 xeneparJ
      real*8 xeneparT
      real*8 xeneperpJ
      real*8 xeneperpT
      real*8 xenediagJ
      real*8 xenediagT
!C  allocate (xtt(maxlien))
!C  allocate (xjjz(maxlien))
!C  allocate (xjjxy(maxlien))

!C   NAMELIST/hamilton/yham,FAM1
!C   NAMELIST/clust/l1,l2,ktyp
!C   NAMELIST/integ/ityp,xjjz,xjjxy,xtt
!C   NAMELIST/param/xj1,xj2,xt1,xt2,xv1,xv2,xv3,xbJ,xbT,xeneperpJ, &
!C   xeneparJ,xenediagJ,xeneperpT,xeneparT,xenediagT,xspar,xsperp
!C   NAMELIST/infdia/isz,iprec,maxda2,niter,ngao,nvec,numero,nes4, &
!C   xseuil,ysuiv,yec
!   NAMELIST/scsrev1/Emin,Emax,M0

    ilien2=0

      yw=.FALSE.
!-----mise a zero
      jclust=0
      natom=0
      nlientot=0
      do ikk=1,maxlien
	iliatom1(ikk)=0
	iliatom2(ikk)=0
      enddo
      do ik=1,natomax
        do j=1,natomax
           ilien(j,ik)=0
           iech(j,ik)=0
        enddo
      enddo
!------------------Lecture Hamiltonien

!      FAM1=.TRUE. 
       yham=.TRUE.
       write(6,*)'HAMILTONIEN t-J'
       write(6,*)'Le nombre de trou est : ',ntrou
       write(6,*)'Famille 1 : ',FAM1
!---------------------------------------------
      write(6,*)' '
      write(6,*)' '
      write(6,*)'LECTURE DES ATOMES, DES LIAISONS, DES INTEGRALES'
      write(6,*)' '
      write(6,*)' '
!-----------Lecture 1ier voisin
!     read(5,clust)
        jclust=jclust+1
	if(yw)write(6,*)' '
	write(6,*)'================ CLUSTER',jclust,'=================='
	write(6,*)' '
!------------------------
       do i=1,maxlien
	  if(l1(i).eq.0)EXIT
	  nlien=i
       enddo

        if(yw)write(6,*)(l1(i),i=1,maxlien)

       write(6,*)'Liaisons entre les atomes',nlien

       do i=1,natomax
	  y(i)=.false.
       enddo
       do i=1,nlien
	  ilien(l1(i),l2(i))=nlientot+i
	  ilien(l2(i),l1(i))=nlientot+i
	  if(l1(i).lt.l2(i))then
	  iliatom1(nlientot+i)=l1(i)
	  iliatom2(nlientot+i)=l2(i)
	  else
	  iliatom1(nlientot+i)=l2(i)
	  iliatom2(nlientot+i)=l1(i)
	  endif
	  iech(l1(i),l2(i))=(-1)**(l2(i)-l1(i)+1)
	  iech(l2(i),l1(i))=(-1)**(l2(i)-l1(i)+1)
	  y(l1(i))=.true.
          y(l2(i))=.true.
          ltyp(nlientot+i)=ktyp(i)
       write(6,*)'Les atomes ',l1(i),' et ',l2(i),' forment la liaison', &
      ilien(l1(i),l2(i)),'qui est de type',ltyp(i)
          if(yw)write(6,*)'iliatom1',iliatom1(i)
       enddo
          nlientot=nlientot+nlien

      do i=1,natomax
       if(y(i))natom=natom+1
      enddo
      write(6,*)'============================================='
      write(6,*)'Le nombre total d atomes est ',natom
      write(6,*)'============================================='

!----------------------------------------------------------------
!----------------Voisins---------------------
!       nvois(ik)=nbre de voisins de ik
!       (ivois(ii,ik),ii=1,nvois(ik))=No des voisins
!--------------------------------------------
      do ik=1,natom
	 nvois1(ik)=0
	 do ii=1,6
	    ivois1(ii,ik)=0
	 enddo
      enddo
      do ik=1,natom
	do ii=1,natom
	    if (ilien(ik,ii).ne.0) then
	       nvois1(ik)=nvois1(ik)+1
	       ivois1(nvois1(ik),ik)=ii
	    endif
	enddo
      enddo

!----------------------------------------------------------------
!      calcul des liaisons voisines de type 1
!-------------------------------------------------------
      do i=1,nlientot
	 nvoisl1(i)=0
	 do j=1,6
	    ivoisl1(j,i)=0
	 enddo
      enddo
      do i=1,natom
	 do j=i+1,natom
	    if(ilien(i,j).ne.0)then
	       l=ilien(i,j)
	    else
	      EXIT
	    endif
	    do ik=1,nvois1(i)
	       if(ivois1(ik,i).ne.j)then
		  nvoisl1(l)=nvoisl1(l)+1
		  ivoisl1(nvoisl1(l),l)=ilien(ivois1(ik,i),i)
	       endif
            enddo
	    do jk=1,nvois1(j)
	       if(ivois1(jk,j).ne.i)then
	         nvoisl1(l)=nvoisl1(l)+1
		 ivoisl1(nvoisl1(l),l)=ilien(ivois1(jk,j),j)
	       endif
            enddo
        enddo
      enddo


      ityp=3
!C    xjjz=(/.1000d0,-0.8d0,0.000d0/)
!C    xjjxy=(/.1000d0,-0.8d0,0.000d0/)
!C    xtt=(/-1.0000d0,0.d0,0.0d0/)

       write(6,*)'Nombre de J differents',ityp
       do ikl=1,nlientot
         Ykl(ikl)=.false.
       enddo
       do il=1,nlientot
          write(6,*)'type de liaison',il,ltyp(il)
       enddo
       do iiki=1,ityp
          write(6,*)'type de J',xjjz(iiki)
       enddo
       do il=1,nlientot
         if(.not.ykl(il))then
          xjz(il)=xjjz(ltyp(il))
          write(6,*)'Parametres : Jz',il,'=',xjz(il)
          xjxy(il)=xjjxy(ltyp(il))
          write(6,*)'Parametres : Jxy',il,'=',xjxy(il)
          xt(il)=xtt(ltyp(il))
          write(6,*)'Parametre : t',il,'=',xt(il)
          ykl(il)=.true.
         endif
       enddo
        xbJ=0.00d0
        xbT=0.0d0
        xv1=0.0d0
        xv2=0.0d0
        xv3=0.0d0
        xt1=-0.20d0
        xt2=0.0d0
        xj1=0.01d0
        xj2=0.d0
        xeneparJ=0.000d0
        xeneparT=0.000d0
        xeneperpJ=0.000d0
        xeneperpT=0.000d0
        xenediagJ=0.000d0
        xenediagT=0.000d0
        xspar=-0.00d0
        xsperp=-0.00d0
       write(6,*)'coucoudslect3'
      write(6,*)'coucou'
	write(6,*)'Parametres pour le t-J'
	 write(6,*)'xj1 = ',xj1
	 write(6,*)'xj2 = ',xj2
	 write(6,*)'xt1 = ',xt1
	 write(6,*)'xt2 = ',xt2
	 write(6,*)'xv1 = ',xv1
	 write(6,*)'xv2 = ',xv2
	 write(6,*)'xv3 = ',xv3
	 write(6,*)'xbj = ',xbj
	 write(6,*)'xbt = ',xbt
	 write(6,*)'xeneparJ = ',xeneparJ
	 write(6,*)'xeneperpJ = ',xeneperpJ
	 write(6,*)'xeneparT = ',xeneparT
	 write(6,*)'xeneperpT = ',xeneperpT
	 write(6,*)'xenediagJ = ',xenediagJ
	 write(6,*)'xenediagT = ',xenediagT
	 write(6,*)'xspar = ',xspar
	 write(6,*)'xsperp = ',xsperp
!===================================================================
!====================================================================
!       calcul des plaquettes pour un 2D t-J
!===================================================================
      do ikl=1,nlientot
         yplac(ikl)=.false.
	 IKX(ikl)=0
      enddo
      nplac=0
      do ikl=1,nlientot
           ik1=iliatom1(ikl)
           ik2=iliatom2(ikl)
	   do kkio1=1,nvois1(ik1)
	      do kkio2=1,nvois1(ik2)
	        kat1=ivois1(kkio1,ik1)
	        kat2=ivois1(kkio2,ik2)
	        if(kat1.ne.ik2.and.kat2.ne.ik1)then
	         if(ilien(kat1,kat2).ne.0)then
		   klo1=ilien2(kat1,ik2)
                   if (klo1 == 0) cycle
		   klo2=ilien2(kat2,ik1)
                   if (klo2 == 0) cycle
		   if(.not.yplac2(klo1))then
		     yplac2(klo1)=.true.
		     yplac2(klo2)=.true.
		     nplac=nplac+1
		     il2=ilien(ik1,kat1)
		     il3=ilien(ik2,kat2)
		     il4=ilien(kat1,kat2)
		     iplac(1,nplac)=ik1
		     iplac(2,nplac)=ik2
		     iplac(3,nplac)=kat1
		     iplac(4,nplac)=kat2
		     IKX(ikl)=IKX(ikl)+1
		     IKX(il2)=IKX(il2)+1
		     IKX(il3)=IKX(il3)+1
		     IKX(il4)=IKX(il4)+1
		     jplac(IKX(ikl),ikl)=nplac
		     jplac(IKX(il2),il2)=nplac
		     jplac(IKX(il3),il3)=nplac
		     jplac(IKX(il4),il4)=nplac
		   endif
	         endif
                endif
	      enddo
           enddo
      enddo
      write(6,*)'Le systeme comporte ',nplac,' plaquettes.'
      do kko=1,nplac
         write(6,*)'La plaquette ',kko,' est contituee des atomes',&
      iplac(1,kko),' ',iplac(2,kko),' ',iplac(3,kko),' et ',iplac(4,kko)
      enddo
!===================================================================
!   isz=0
    IPREC=8
    maxda2=20
    NITER=280
    ngao=100
    nvec=8
    numero=1
    NES4=0
    xseuil=1.0E-008
    ysuiv=.FALSE.
    yec=.TRUE.
      write(6,*)'Spin total',isz
      write(6,*)'Nombre de vecteurs demande',nvec
      write(6,*)'Nombre maximal d iterations de Davidson',niter
      write(6,*)'Vecteur calcule le plus bas',numero
      write(6,*)'Variable Nes4 (vecteurs d essai)',nes4
      write(6,*)'Nombre de determinants en donnees',ngao
      write(6,*)'Variable Ysuiv (suivre le vecteur initial)',ysuiv
      write(6,*)'Seuil au dela duquel seront ecrits les vecteurs',xseuil
      write(6,*)'Option d ecriture des determinants sur FIL2',yec
!     write(6,*)Emin
!     write(6,*)Emax
!     write(6,*)M0
      if(yham)then
      IAL0=(natom-ntrou)/2+mod(natom-ntrou,2)+isz
      else
      IAL0=(natom)/2+mod(natom,2)+isz
      endif
      write(6,*)'=======nombre de centres de spin alpha=====',ial0
      natrest=natom-ntrou

    !C calculating nalpha and nbeta
    if(mod(natom-ntrou+2*isz,2).eq.0)then
        nalpha=(natom-ntrou+2*isz)/2
        nbeta=(natom -ntrou-2*isz)/2
        if(((natom-ntrou)/2).eq.isz)then
            nbeta=0
        endif
    else
        nalpha=(natom-ntrou+2*isz+1)/2
        nbeta=(natom -ntrou-2*isz-1)/2
        if(((natom-ntrou+1)/2).eq.isz)then
            nbeta=0
        endif
    endif

END_PROVIDER
