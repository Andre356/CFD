      SUBROUTINE PRINTTEC
      include 'common.f'
      INTEGER ICELL(NCELL,4)
      OPEN(UNIT=26,FILE='varout.dat',FORM='FORMATTED')
      REWIND(26)
      WRITE(26,*) 'TITLE = "2D Jameson Output"'
      WRITE(26,*) 'VARIABLES = "X","Y","Rho","U","V","P"'
      WRITE(26,*) 'ZONE T="Zone 1" Nodes =',NPTI,', Elements=',NC,
     :  ', ZONETYPE=FETriangle, DATAPACKING=BLOCK'
      WRITE(26,*) 'VARLOCATION=([3,4,5,6]=CELLCENTERED)'
      WRITE(26,*) 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)'
C     Print out X co-ordinates in blocks of 10
      DO I=1,NPTI,10
         WRITE(26,901) (XP(J),J=I,MIN(I+9,NPTI))
      ENDDO
C     Print out Y co-ordinates in blocks of 10
      DO I=1,NPTI,10
         WRITE(26,901) (YP(J),J=I,MIN(I+9,NPTI))
      ENDDO
C     Print out Rho 
      DO I=1,NC,10
         WRITE(26,901) (W(J,1),J=I,MIN(I+9,NC))
      ENDDO
C     Print out U
      DO I=1,NC,10
         WRITE(26,901) (W(J,2)/W(J,1),J=I,MIN(I+9,NC))
      ENDDO
C     Print out V
      DO I=1,NC,10
         WRITE(26,901) (W(J,3)/W(J,1),J=I,MIN(I+9,NC))
      ENDDO
C     Print out Pressure
      DO I=1,NC,10
         WRITE(26,901) (P(J),J=I,MIN(I+9,NC))
      ENDDO
C     Zero ICELL - will be used to build connectivity data
      DO I = 1,NC
         ICELL(I,1) = 0
         ICELL(I,2) = 0
         ICELL(I,3) = 0
         ICELL(I,4) = 0
      ENDDO
C     Build cell list - can do something interesting here.
C
C     First we loop over the boundary elements - i.e. all the one sided elements
C    
C     If a given cell has no current points associated to it then both the points on
C     the current edge must be new. (This means its impossible to get to a state where ICELL(cell,1) = 1)
C
C     IF ICELl(cell,1) is not zero it can either be 2 - from above and hence one more point is needed or
C     3 and no more points are required hence only have to check for ICELL(cell,1) equal 2 now
C     If it is 2 then one of the current edge must contain the extra point find which edge it is and add it
C     into the list. The list is now complete as it contains 3 nodes and no more checks are required.

      DO I = 1,NBD
        K  = JLIST(I,2)
        IA = JLIST(I,3)
        IB = JLIST(I,4)
        IF (ICELL(K,1).eq.0) THEN
           ICELL(K,2) = IA
           ICELL(K,3) = IB
           ICELL(K,1) = 2
        ELSE IF (ICELL(K,1).eq.2) THEN
           IF(ICELL(K,2).NE.IA.AND.ICELL(K,3).NE.IA) THEN
              ICELL(K,4) = IA
              ICELL(K,1) = 3
           ELSE
              ICELL(K,4) = IB
              ICELL(K,1) = 3
           ENDIF
        ENDIF
      ENDDO
      DO I=NBD+1,NEDGE
        K  = JLIST(I,2)
        IA = JLIST(I,3)
        IB = JLIST(I,4)
        IP = JLIST(I,5)
C
        IF (ICELL(K,1).eq.0) THEN
           ICELL(K,2) = IA
           ICELL(K,3) = IB
           ICELL(K,1) = 2
        ELSE IF (ICELL(K,1).eq.2) THEN
           IF(ICELL(K,2).NE.IA.AND.ICELL(K,3).NE.IA) THEN
              ICELL(K,4) = IA
              ICELL(K,1) = ICELL(K,1) + 1
           ELSE
              ICELL(K,4) = IB
              ICELL(K,1) = ICELL(K,1) + 1
           ENDIF
        ENDIF
C
        IF (ICELL(IP,1).eq.0) THEN
           ICELL(IP,2) = IA
           ICELL(IP,3) = IB
           ICELL(IP,1) = 2
        ELSE IF (ICELL(IP,1).eq.2) THEN
           IF(ICELL(IP,2).NE.IA.AND.ICELL(IP,3).NE.IA) THEN
              ICELL(IP,4) = IA
              ICELL(IP,1) = ICELL(IP,1)+1
           ELSE
              ICELL(IP,4) = IB
              ICELL(IP,1) = ICELL(IP,1)+1
           ENDIF
        ENDIF
      ENDDO
C
      DO I=1,NC
         WRITE(26,902) ICELL(I,2),ICELL(I,3),ICELL(I,4)
      ENDDO
      CLOSE(26)
      RETURN
901   FORMAT(10E16.8)
902   FORMAT(3I8)
      END
