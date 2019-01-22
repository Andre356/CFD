      PROGRAM JAMESON ! UNSTRUCTURED GRID, EULER EQUATIONS, EXPLICIT METHOD
      include 'common.f'
      OPEN(UNIT=10,FILE='grid.dat'  ,FORM='FORMATTED') ! Input grid and flow params
      OPEN(UNIT=20,FILE='output.dat',FORM='FORMATTED') ! Nonvergence history, loads, vars
      CALL INIT                                        ! READ INPUT, INITIALISE VARIABLES
      CLOSE(10)                                        ! Done with input
      DO N=1,NCYC                                      ! START TIME STEPS (CYCLES)
         CALL BCON                                     ! APPLY BOUNDARY CONDITIONS
         CALL TIME                                     ! CALCULATE LOCAL TIME STEP
         DO NS=1,NRKS                                  ! START RUNGE-KUTTA STEPS
            IF(NS.NE.1) CALL BCON                      ! Apply boundary conditions
            CALL FLUX                                  ! CALCULATE FLUXES
            IF(NS.EQ.1) CALL DISS                      ! ADD DISSIPATION
            COEFF = ST(NS)                             ! R-K Coef for step NS
            DO I=1,NC                                  ! UPDATE SOLUTION, CALCULATE PRESSURE over all cells
               COEFF2 = COEFF*DT(I)
               W(I,1) = ABS( WS(I,1)-COEFF2*(Q(I,1)-D(I,1)) )
               W(I,2) =      WS(I,2)-COEFF2*(Q(I,2)-D(I,2))
               W(I,3) =      WS(I,3)-COEFF2*(Q(I,3)-D(I,3))
               W(I,4) = ABS( WS(I,4)-COEFF2*(Q(I,4)-D(I,4)) )
               U      = W(I,2)/W(I,1)
               V      = W(I,3)/W(I,1)
               P(I)   = GM1*(W(I,4)-0.5*W(I,1)*(U*U+V*V))
            ENDDO
         ENDDO                                         ! END OF R-K STEPS
         DO I=1,NC                                     ! STORE OLD SOLUTION
            WS(I,1) = W(I,1)
            WS(I,2) = W(I,2)
            WS(I,3) = W(I,3)
            WS(I,4) = W(I,4)
         ENDDO
         CX  = 0.0                                     ! CALCULATE INTEGRATED LOADS
         CY  = 0.0
         CM  = 0.0
         DO I=1,ICU-ICL
            K       = JLIST(I,2)
            CP(I)   = (P(K)-1.)/CPD
            IA      = JLIST(I,3)
            IB      = JLIST(I,4)
            XCM     = (XP(IB)+XP(IA))*0.5-0.25
            YCM     = (YP(IB)+YP(IA))*0.5
            DCX     =  CP(I)*(YP(IB)-YP(IA))
            DCY     = -CP(I)*(XP(IB)-XP(IA))
            CM      = CM-DCY*XCM+DCX*YCM
            CX      = CX+DCX
            CY      = CY+DCY
            XF(I)   = XCM+0.25
         ENDDO
         CL         = CY*CA-CX*SA
         CD         = CY*SA+CX*CA
         WRITE(20,*) N,CL,CD,CM                        ! CONVERGENCE HISTORY OF LOADS
      ENDDO                                            ! END OF TIME STEPS
      WRITE(20,*) '--- Flow Conditions ---'
      WRITE(20,*) FMACH,AL*180/PI,CFL
      WRITE(20,*) '--- Suface Loads ---'
      DO I = 1,ICU-ICL
         K = JLIST(I,2)
         WRITE(20,*) XF(I),(P(K)-1.)/CPD
      ENDDO
      WRITE(20,*) '--- Integral Loads ---'
      WRITE(20,*) 'CL=',CL,'CD=',CD,'CM=',CM
      WRITE(20,*) '--- Flow Variables ---'
      DO I=1,NC
         WRITE(20,*) W(I,1),W(I,2),W(I,3),W(I,4) 
      ENDDO
      CLOSE(20)
      Call printtec
      STOP
      END
C -------------------------------------------------------------------
      SUBROUTINE INIT ! READS FLOW PARAMS AND GRID, INITIALISES VARS
      include 'common.f'
      READ(10,*) FMACH,AL,CFL,NCYC ! Mach, Incidence, CFL, Iterations
      NRKS  = 4                    ! 4th order R-K
      ST(1) = 0.25   *ABS(CFL)     ! Coefs for the R-K steps
      ST(2) = 0.33333*ABS(CFL)
      ST(3) = 0.5    *ABS(CFL)
      ST(4) = 1.0    *ABS(CFL)
      RK2   = 0.5                  ! 1st dissipation parameter
      RK4   = 0.008                ! 2nd dissipation parameter
      GM    = 1.4                  ! Gamma of air
      GM1   = GM-1.                ! Gamma - 1
      PI    = 4.0*ATAN(1.)         ! 3.1415926535...
      AL    = AL*PI/180.           ! Alpha is in rads now more
      SA    = SIN(AL)              ! Sin of alpha
      CA    = COS(AL)              ! Cos of alpha
      R0    = 1.0                  ! Free stream density
      P0    = 1.0                  ! Free stream pressure
      C0    = SQRT(GM*P0/R0)       ! Free stream speed of sound
      U0    = FMACH*C0*COS(AL)     ! Free stream u-velocity scaled with Mach
      V0    = FMACH*C0*SIN(AL)     ! Free stream v-velocity scaled with Mach
      E0    = P0/(R0*GM1)+0.5*(U0*U0+V0*V0) ! Free stream energy
      CPD   = GM*FMACH*FMACH*0.5   ! Reference Dynamic Head
      READ(10,*) ICL,ICU,NBD,NPTI,NC,NEDGE ! ICL=1 ICU is the number of nodes around the aerofoil, NBD number of Far Field nodes
      DO I=1,NEDGE ! No of edges
         READ(10,*) (JLIST(I,K),K=1,5) ! K IA IB IP
      ENDDO
      DO I=1,NPTI ! No of nodes
         READ(10,*) XP(I),YP(I)
      ENDDO
      DO I=1,NC ! Inside the domain
         W(I,1)  = R0    ! W is the currnet solution
         W(I,2)  = R0*U0
         W(I,3)  = R0*V0
         W(I,4)  = R0*E0
         P(I)    = P0
         WS(I,1) = W(I,1) ! WS is a stored solution at the previous time step
         WS(I,2) = W(I,2)
         WS(I,3) = W(I,3)
         WS(I,4) = W(I,4)
      ENDDO
      DO I=1,NBD ! Bounbdary values
         BW(I,1) = R0
         BW(I,2) = R0*U0
         BW(I,3) = R0*V0
         BW(I,4) = R0*E0
         BW(I,5) = P0
      ENDDO
      RETURN
      END
C -------------------------------------------------------------------
      SUBROUTINE TIME ! CALCULATION OF THE LOCAL TIME STEP
      include 'common.f'
      VK2=0.15
      DO I=1,NC
         DT(I)=0.0
      ENDDO
      DO I=1,ICU-ICL ! Around the aerofoil
         K    =JLIST(I,2)
         IA   =JLIST(I,3)
         IB   =JLIST(I,4)
         C2   =GM*P(K)/W(K,1)
         DX   =XP(IB)-XP(IA)
         DY   =YP(IB)-YP(IA)
         DL2  =DX*DX+DY*DY
         T2   =SQRT(ABS(C2*DL2))
         U    =W(K,2)/W(K,1)
         V    =W(K,3)/W(K,1)
         T1   =ABS(DY*U-DX*V)
         DT(K)=DT(K)+T1+T2
      ENDDO
      DO I=ICU-ICL+1,NBD ! Around the boundary
         K       = JLIST(I,2)
         IA      = JLIST(I,3)
         IB      = JLIST(I,4)
         C2      = GM*BW(I,5)/BW(I,1)
         DX      = XP(IB)-XP(IA)
         DY      = YP(IB)-YP(IA)
         DL2     = DX*DX+DY*DY
         T2      = SQRT(ABS(C2*DL2))
         U       = BW(I,2)/BW(I,1)
         V       = BW(I,3)/BW(I,1)
         T1      = ABS(DY*U-DX*V)
         DT(K)   = DT(K)+T1+T2
      ENDDO
      DO I=NBD+1,NEDGE ! Inside the domain
         K       = JLIST(I,2)
         IA      = JLIST(I,3)
         IB      = JLIST(I,4)
         IP      = JLIST(I,5)
         C2K     = GM*(P(K)+P(IP))/(W(K,1)+W(IP,1))
         UK      = (W(K,2)/W(K,1)+W(IP,2)/W(IP,1))*0.5
         VK      = (W(K,3)/W(K,1)+W(IP,3)/W(IP,1))*0.5
         DX      = XP(IB)-XP(IA)
         DY      = YP(IB)-YP(IA)
         DL2     = DX*DX+DY*DY
         DTTT    = ABS(UK*DY - VK*DX)+SQRT(C2K*DL2)
         DT(K)   = DT(K)  + DTTT
         DT(IP)  = DT(IP) + DTTT
      ENDDO
      DO I=1,NC ! We store 1/Dt
         IF (DT(I).NE.0.0) DT(I) = 1./DT(I)
      ENDDO
      RETURN
      END
C -------------------------------------------------------------------
      SUBROUTINE BCON
      include 'common.f'
      DO I=ICU-ICL+1,NBD ! SETTING FARFIELD BOUNDARY CONDITION
        K       = JLIST(I,2)
        IA      = JLIST(I,3)
        IB      = JLIST(I,4)
        DX      = XP(IB)-XP(IA)
        DY      = YP(IB)-YP(IA)
        DH      = SQRT(DX*DX+DY*DY)
        STH     = DY/DH
        CTH     = DX/DH
        RE      = W(K,1)
        U       = W(K,2)/RE
        V       = W(K,3)/RE
        PE      = P(K)
        CE      = SQRT(ABS(GM*PE/RE))
        XMACH   = SQRT((U*U+V*V))/CE
        S0      = P0/R0**GM
        QN0     = U0*STH-V0*CTH
        QT0     = U0*CTH+V0*STH
        RI0     = QN0-2.*C0/GM1
        SE      = PE/RE**GM
        QNE     = U*STH-V*CTH
        QTE     = U*CTH+V*STH
        RIE     = QNE+2.*CE/GM1
        IF (XMACH.LT.1.0) THEN ! SUBSONIC FLOW
           QN     =     0.50*(RIE+RI0)
           C      = GM1*0.25*(RIE-RI0)
           IF (QN.LT.0.0) THEN
              QTT = QT0
              SB  = S0
           ELSE
              QTT = QTE
              SB  = SE
           ENDIF
        ELSE                   ! SUPERSONIC FLOW
           IF (QN0.LT.0.0) THEN
              QN  = QN0
              C   = C0
              QTT = QT0
              SB  = S0
           ELSE 
              QN  = QNE
              C   = CE
              QTT = QTE
              SB  = SE
           ENDIF
        ENDIF
        UB      = QN*STH+QTT*CTH ! CALCULATE PRIMARY VARIABLES
        VB      = QTT*STH-QN*CTH
        RB      = (C*C/GM/SB)**(1.0/GM1)
        PB      = C*C*RB/GM
        REB     = PB/GM1 + 0.5*RB*(UB*UB+VB*VB)
        BW(I,1) = RB ! CALCULATE CONSERVATIVE VARIABLES
        BW(I,2) = RB*UB
        BW(I,3) = RB*VB
        BW(I,4) = REB
        BW(I,5) = PB
      ENDDO
      RETURN
      END
C -------------------------------------------------------------------
      SUBROUTINE FLUX
      include 'common.f'
      DO I=1,NC
        Q(I,1) = 0.0
        Q(I,2) = 0.0
        Q(I,3) = 0.0
        Q(I,4) = 0.0
      ENDDO
      DO I=1,ICU-ICL ! FLUX CALCULATION (ONLY PRESSURE CONTRIBUTION) AIRFOIL SURFACE
         K       = JLIST(I,2)
         IA      = JLIST(I,3)
         IB      = JLIST(I,4)
         DX      = XP(IB)-XP(IA)
         DY      = YP(IB)-YP(IA)
         Q(K,2)  = + DY* P(K)
         Q(K,3)  = - DX* P(K)
      ENDDO
      DO I=ICU-ICL+1,NBD ! FLUX CALCULATION OUTER BOUNDARY
         K      = JLIST(I,2)
         IA     = JLIST(I,3)
         IB     = JLIST(I,4)
         U      = BW(I,2)/BW(I,1)
         V      = BW(I,3)/BW(I,1)
         DX     = XP(IB)-XP(IA)
         DY     = YP(IB)-YP(IA)
         QK     = DY*U-DX*V
         DP     = BW(I,5)
         F1     = QK*BW(I,1)
         Q(K,1) = Q(K,1)+F1
         F2     = QK*BW(I,2)+DP*DY
         Q(K,2) = Q(K,2)+F2
         F3     = QK*BW(I,3)-DP*DX
         Q(K,3) = Q(K,3)+F3
         F4     = QK*(BW(I,4)+DP)
         Q(K,4) = Q(K,4)+F4
      ENDDO
      DO I=NBD+1,NEDGE ! Flux calculation inside the domain
        K       = JLIST(I,2)
        IA      = JLIST(I,3)
        IB      = JLIST(I,4)
        IP      = JLIST(I,5)
        U       = .5*(W(K,2)/W(K,1)+W(IP,2)/W(IP,1))
        V       = .5*(W(K,3)/W(K,1)+W(IP,3)/W(IP,1))
        DX      = XP(IB)-XP(IA)
        DY      = YP(IB)-YP(IA)
        QK      = DY*U-DX*V
        DP      = P(K)+P(IP)
        F1      = QK*(W(K,1)+W(IP,1))*.5
        Q(K,1)  = Q(K,1)+F1
        Q(IP,1) = Q(IP,1)-F1
        F2      = (QK*(W(K,2)+W(IP,2))+DP*DY)*0.5
        Q(K,2)  = Q(K,2)+F2
        Q(IP,2) = Q(IP,2)-F2
        F3      = (QK*(W(K,3)+W(IP,3))-DP*DX)*0.5
        Q(K,3)  = Q(K,3)+F3
        Q(IP,3) = Q(IP,3)-F3
        F4      = QK*(W(K,4)+W(IP,4)+DP)*0.5
        Q(K,4)  = Q(K,4)+F4
        Q(IP,4) = Q(IP,4)-F4
      ENDDO
      RETURN
      END
C -------------------------------------------------------------------
      SUBROUTINE DISS ! CALCULATION OF THE DISSIPATION FUNCTION
      include 'common.f'
      DIMENSION DW2(NPOINT,4)
      ZNIL = 0.0
      DO I=1,NC
         DW2(I,1) = 0.0
         DW2(I,2) = 0.0
         DW2(I,3) = 0.0
         DW2(I,4) = 0.0
         D(I,1)   = 0.0
         D(I,2)   = 0.0
         D(I,3)   = 0.0
         D(I,4)   = 0.0
      ENDDO
      DO I=NBD+1,NEDGE
        K         = JLIST(I,2)
        IP        = JLIST(I,5)
        DELW1     = (W(IP,1)-W(K,1))
        DELW2     = (W(IP,2)-W(K,2))
        DELW3     = (W(IP,3)-W(K,3))
        DELW4     = (W(IP,4)-W(K,4)+P(IP)-P(K))
        DW2(K,1)  = DW2(K,1)+DELW1
        DW2(K,2)  = DW2(K,2)+DELW2
        DW2(K,3)  = DW2(K,3)+DELW3
        DW2(K,4)  = DW2(K,4)+DELW4
        DW2(IP,1) = DW2(IP,1)-DELW1
        DW2(IP,2) = DW2(IP,2)-DELW2
        DW2(IP,3) = DW2(IP,3)-DELW3
        DW2(IP,4) = DW2(IP,4)-DELW4
      ENDDO
      DO I=NBD+1,NEDGE
        K       = JLIST(I,2)
        IA      = JLIST(I,3)
        IB      = JLIST(I,4)
        IP      = JLIST(I,5)
        C2      = GM*(P(K)+P(IP))/(W(K,1)+W(IP,1))
        U       = (W(K,2)/W(K,1)+W(IP,2)/W(IP,1))*0.5
        V       = (W(K,3)/W(K,1)+W(IP,3)/W(IP,1))*0.5
        DX      = XP(IB)-XP(IA)
        DY      = YP(IB)-YP(IA)
        DL2     = DX*DX+DY*DY
        A       = ABS(DY*U-DX*V) + SQRT(C2*DL2)
        EPS2    = RK2*ABS((P(IP)-P(K))/(P(IP)+P(K)))
        EPS4    = DMAX1(ZNIL,(RK4-EPS2))
        DELW1   = W(IP,1)-W(K,1)
        DELW2   = W(IP,2)-W(K,2)
        DELW3   = W(IP,3)-W(K,3)
        DELW4   = W(IP,4)-W(K,4)+P(IP)-P(K)
        T1      = EPS2*DELW1
        T2      = EPS2*DELW2
        T3      = EPS2*DELW3
        T4      = EPS2*DELW4
        S1      = EPS4*(DW2(IP,1)-DW2(K,1))
        S2      = EPS4*(DW2(IP,2)-DW2(K,2))
        S3      = EPS4*(DW2(IP,3)-DW2(K,3))
        S4      = EPS4*(DW2(IP,4)-DW2(K,4))
        T1      = A*(T1-S1)
        T2      = A*(T2-S2)
        T3      = A*(T3-S3)
        T4      = A*(T4-S4)
        D(K,1)  = D(K,1)+T1
        D(K,2)  = D(K,2)+T2
        D(K,3)  = D(K,3)+T3
        D(K,4)  = D(K,4)+T4
        D(IP,1) = D(IP,1)-T1
        D(IP,2) = D(IP,2)-T2
        D(IP,3) = D(IP,3)-T3
        D(IP,4) = D(IP,4)-T4
      ENDDO
      RETURN
      END
