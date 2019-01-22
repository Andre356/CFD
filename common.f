      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NCELL=25000, NPOINT=40000, IPT=800)
      COMMON/GEO/ XP(NPOINT),      ! X-coords of nodes
     1            YP(NPOINT),      ! Y-coords of nodes
     1            JLIST(NPOINT,5), ! Connectivity array
     1            NC,              ! No of cells
     1            NEDGE,           ! No of edges
     1            NBD,             ! No of boundary nodes
     1            ICU,             ! Aerofoil poins
     1            ICL,             ! Aerofoil first point
     1            NPTI             ! No of points on mesh
      COMMON/VAR/ W (NCELL,4),     ! Variables new
     1            WS(NCELL,4),     ! Variables old
     1            P (NCELL),       ! Pressure
     1            Q (NCELL,4),     ! Inviscid flux at each cell and for each variable
     1            D (NCELL,4),     ! Dissipation flux for each cell and for each variable
     1            DT(NCELL),       ! Local time step for each cell , we store 1/DT
     1            CP(IPT),         ! surface Cp
     1            XF(IPT),         ! X-coords of surface mid-points between nodes where CP is computed
     1            BW(NPOINT,5)     ! Flow vars at boundary Rho, RhoU, RhoV, E, Pressure
      COMMON/EXT/ GM,              ! gamma = 1.4
     1            GM1,             ! gamma -1 = 0.4
     1            FMACH,           ! Free stream mach
     1            AL,              ! Angle of attack, read in degrees and converted to rads
     1            PI,              ! 3.1415926535....
     1            R0,              ! Free stream density
     1            P0,              ! Free stream pressure
     1            C0,              ! Free stream speed of sound
     1            U0,              ! Free stream u-velocity
     1            V0,              ! Free stream v-velocity
     1            E0,              ! Free stream energy
     1            CL,              ! Lift coef
     1            CD,              ! Drag coef
     1            CM,              ! Pitching moment coef about c/4
     1            CA,              ! Cos of Alpha
     1            SA,              ! Sin of Alpha
     1            CPD,             ! Reference Dynamic Head, used for CP calcs
     1            CFL,             ! CFL number
     1            RK2,RK4,         ! Dissipation coefs
     1            ST(4),           ! R-K coefs
     1            NRKS,            ! No of R-K steps 4
     1            N,               ! Time step counter
     1            NCYC             ! Max number of time steps
