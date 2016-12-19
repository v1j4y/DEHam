      SUBROUTINE COMBIN(NS1,N,M,NT,MD,ND)
      INTEGER NS1(MD,ND)
      INTEGER :: i,j,MMJ1,J1
      I=1
      DO J=1,M
       NS1(j,i)=J
      enddo
      J1=0
      do while (J1.LT.M)
      J1=0
      DO J=1,M
      IF(NS1(j,i).EQ.N-M+J) J1=J1+1
      enddo
      IF(J1.EQ.M) EXIT
      I=I+1
      DO J=1,M
      NS1(j,i)=NS1(j,I-1)
      enddo
      NS1(m-j1,I)=NS1(m-j1,I-1)+1
      IF(J1.EQ.0) CYCLE
      MMJ1=M-J1+1
      DO J=MMJ1,M
      NS1(j,I)=NS1(j-1,I)+1
      enddo
      enddo
      NT=I
      RETURN
      END
