1         ******************************************************
          *            GAMESS VERSION =  6 JUN 1999            *
          *             FROM IOWA STATE UNIVERSITY             *
          * M.W.SCHMIDT, K.K.BALDRIDGE, J.A.BOATZ, S.T.ELBERT, *
          *   M.S.GORDON, J.H.JENSEN, S.KOSEKI, N.MATSUNAGA,   *
          *          K.A.NGUYEN, S.J.SU, T.L.WINDUS,           *
          *       TOGETHER WITH M.DUPUIS, J.A.MONTGOMERY       *
          *         J.COMPUT.CHEM.  14, 1347-1363(1993)        *
          *******Intel x86 (Win32,Linux,OS/2,DOS) VERSION*******
          * PC GAMESS version 7.1 (Tornado), build number 4470 *
          *   Compiled on    Tuesday,   21-08-2007, 17:21:28   *
          *      Intel specific optimization, bug fixes,       *
          *    code changes, and additional functionality -    *
          *  copyright (c) 1994, 2007 by  Alex A. Granovsky,   *
          *        Laboratory of Chemical Cybernetics,         *
          *      Moscow State University, Moscow, Russia.      *
          *   Some parts of this program include code due to   *
          * work of Jim Kress, Peter Burger, and Robert Ponec. *
          ******************************************************
          *                PC GAMESS homepage:                 *
          * http://classic.chem.msu.su/gran/gamess/index.html  *
          *                      e-mail:                       *
          *               gran@classic.chem.msu.su             *
          *   This program may not be redistributed without    *
          * the specific, written permission of its developers.*
          ******************************************************


 AMD Opteron/ Win32  PC GAMESS version running under Windows NT
 Running on AMD CPU  :  CPU Generation 15, Family 15, Model  5, Stepping  8
 CPU Brand String    :  AMD Opteron(tm) Processor 144                   
 CPU Features        :  CMOV, MMX, SSE, SSE2, AMD64
 Data cache size     :  L1 64 KB, L2 1024 KB, L3     0 KB
 Operating System successfully passed SSE support test.


 EXECUTION OF GAMESS BEGUN 20:41:24 LT  27-AUG-2007

            ECHO OF THE FIRST FEW INPUT CARDS -
 INPUT CARD>! EXAM 12.                                                                      
 INPUT CARD>!   This job illustrates linear bends, for acetylene.                           
 INPUT CARD>!   The optimal RHF/STO-2G geometry is located.                                 
 INPUT CARD>!                                                                               
 INPUT CARD>!   At the input geometry,                                                      
 INPUT CARD>!   the FINAL E= -73.5036974734 after 7 iterations,                             
 INPUT CARD>!   and the RMS gradient is 0.1506891.                                          
 INPUT CARD>!                                                                               
 INPUT CARD>!   At the final geometry, 7 steps later,                                       
 INPUT CARD>!   the FINAL E= -73.6046483165, RMS gradient=0.0000028,                        
 INPUT CARD>!   R(CC)=1.1777007 and R(CH)=1.0749435.                                        
 INPUT CARD>!                                                                               
 INPUT CARD> $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE NZVAR=5 $END                                
 INPUT CARD> $SYSTEM TIMLIM=6 MEMORY=100000 $END                                            
 INPUT CARD> $BASIS  GBASIS=STO NGAUSS=2 $END                                               
 INPUT CARD> $GUESS  GUESS=HUCKEL $END                                                      
 INPUT CARD> $DATA                                                                          
 INPUT CARD>Acetylene geometry optimization in internal coordinates                         
 INPUT CARD>Dnh      4                                                                      
 INPUT CARD>                                                                                
 INPUT CARD>CARBON      6.0    0.0  0.0  0.70                                               
 INPUT CARD>HYDROGEN    1.0    0.0  0.0  1.78                                               
 INPUT CARD> $END                                                                           
 INPUT CARD> $ZMAT  IZMAT(1)=1,1,2,   1,1,3,   1,2,4,                                       
 INPUT CARD>                 5,1,2,4,    5,2,1,3  $END                                      
 INPUT CARD>------- XZ is 1st plane for both bends -------                                  
 INPUT CARD> $LIBE  APTS(1)=1.0,0.0,0.0,1.0,0.0,0.0 $END                                    
     200000 WORDS OF MEMORY AVAILABLE

     BASIS OPTIONS
     -------------
     GBASIS=STO          IGAUSS=       2      POLAR=NONE    
     NDFUNC=       0     NFFUNC=       0     DIFFSP=       F
     NPFUNC=       0      DIFFS=       F


     RUN TITLE
     ---------
 Acetylene geometry optimization in internal coordinates                         

 THE POINT GROUP OF THE MOLECULE IS DNH     
 THE ORDER OF THE PRINCIPAL AXIS IS     4

 ATOM      ATOMIC                      COORDINATES (BOHR)
           CHARGE         X                   Y                   Z
 CARBON      6.0     0.0000000000        0.0000000000       -1.3228081914
 CARBON      6.0     0.0000000000        0.0000000000        1.3228081914
 HYDROGEN    1.0     0.0000000000        0.0000000000       -3.3637122581
 HYDROGEN    1.0     0.0000000000        0.0000000000        3.3637122581

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    CARBON         CARBON         HYDROGEN       HYDROGEN  

  1  CARBON          0.0000000      1.4000000 *    1.0800000 *    2.4800000 *  
  2  CARBON          1.4000000 *    0.0000000      2.4800000 *    1.0800000 *  
  3  HYDROGEN        1.0800000 *    2.4800000 *    0.0000000      3.5600000    
  4  HYDROGEN        2.4800000 *    1.0800000 *    3.5600000      0.0000000    

  * ... LESS THAN  3.000


     ATOMIC BASIS SET
     ----------------
 THE CONTRACTED PRIMITIVE FUNCTIONS HAVE BEEN UNNORMALIZED
 THE CONTRACTED BASIS FUNCTIONS ARE NOW NORMALIZED TO UNITY

 SHELL TYPE PRIM    EXPONENT          CONTRACTION COEFFICIENTS

 CARBON    

   3   S    1      27.385033    3.669807 (  0.430128) 
   3   S    2       4.874522    1.587353 (  0.678914) 

   4   L    3       1.136748    0.038816 (  0.049472)     0.855856 (  0.511541) 
   4   L    4       0.288309    0.270261 (  0.963782)     0.184542 (  0.612820) 

 HYDROGEN  

   6   S    5       1.309756    0.375320 (  0.430128) 
   6   S    6       0.233136    0.162342 (  0.678914) 

 TOTAL NUMBER OF SHELLS              =    6
 TOTAL NUMBER OF BASIS FUNCTIONS     =   12
 NUMBER OF ELECTRONS                 =   14
 CHARGE OF MOLECULE                  =    0
 STATE MULTIPLICITY                  =    1
 NUMBER OF OCCUPIED ORBITALS (ALPHA) =    7
 NUMBER OF OCCUPIED ORBITALS (BETA ) =    7
 TOTAL NUMBER OF ATOMS               =    4
 THE NUCLEAR REPULSION ENERGY IS       22.1963425659

 THIS MOLECULE IS RECOGNIZED AS BEING LINEAR.

     $CONTRL OPTIONS
     ---------------
     SCFTYP=RHF          RUNTYP=OPTIMIZE     EXETYP=RUN     
     MPLEVL=       0     LOCAL =NONE         UNITS =ANGS    
     MULT  =       1     ICHARG=       0     MAXIT =      30
     NPRINT=       7     IREST =       0     COORD =UNIQUE  
     ECP   =NONE         NORMF =       0     NORMP =       0
     ITOL  =      20     ICUT  =       9     NZVAR =       5
     NOSYM =       0     INTTYP=POPLE        GEOM  =INPUT   
     PLTORB=       F     MOLPLT=       F     RPAC  =       F
     AIMPAC=       F     FRIEND=             CITYP =NONE    
     DFTTYP=NONE    

     $SYSTEM OPTIONS
     ---------------
     KDIAG =       0     MEMORY=   200000     TIMLIM=      360.0 SEC.
     COREFL=       F     PTIME =        F     XDR   =       F
     BALTYP=LOOP    

          ----------------
          PROPERTIES INPUT
          ----------------

     MOMENTS            FIELD           POTENTIAL          DENSITY
 IEMOM =       1   IEFLD =       0   IEPOT =       0   IEDEN =       0
 WHERE =COMASS     WHERE =NUCLEI     WHERE =NUCLEI     WHERE =NUCLEI  
 OUTPUT=BOTH       OUTPUT=BOTH       OUTPUT=BOTH       OUTPUT=BOTH    
 IEMINT=       0   IEFINT=       0                     IEDINT=       0
                                                       MORB  =       0

          EXTRAPOLATION IN EFFECT
          SOSCF IN EFFECT

          ----------------------
          INTEGRAL INPUT OPTIONS
          ----------------------
 NOPK  =       1 NORDER=       0 SCHWRZ=       T

 ATTENTION! AO INTEGRALS WILL BE PACKED.
 THRESHOLD FOR PACKING PKTHR =  0.10000000D-01

     -------------------------------
     INTEGRAL TRANSFORMATION OPTIONS
     -------------------------------
     NWORD  =       0     CUTTRF = 1.0E-09
     MPTRAN =       0     DIRTRF =       F
     AOINTS =DUP          IREST  =       0

   --- ENCODED Z MATRIX ---
 COORD    TYPE     I    J    K    L    M    N
     1       1     1    2
     2       1     1    3
     3       1     2    4
     4       5     1    2    4
     5       5     2    1    3

 THE DETERMINANT OF THE G MATRIX IS 10**(    -3)


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -              COORDINATE      COORDINATE
 NO.     TYPE      I    J    K    L    M    N        (BOHR,RAD)       (ANG,DEG)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           2.6456164       1.4000000
     2 STRETCH     1    3                           2.0409041       1.0800000
     3 STRETCH     2    4                           2.0409041       1.0800000
     4 LIN.BEND    1    2    4                      3.1415927     180.0000000
     5 LIN.BEND    1    2    4                      3.1415927     180.0000000
     6 LIN.BEND    2    1    3                      3.1415927     180.0000000
     7 LIN.BEND    2    1    3                      3.1415927     180.0000000

     ------------------------------------------
     THE POINT GROUP IS DNH, NAXIS= 4, ORDER=16
     ------------------------------------------

     DIMENSIONS OF THE SYMMETRY SUBSPACES ARE
 A1G =   4      A1U =   0      B1G =   0      B1U =   0      A2G =   0
 A2U =   4      B2G =   0      B2U =   0      EG  =   1      EU  =   1

 ..... DONE SETTING UP THE RUN .....

 CPU        TIME:   STEP =      0.02 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.01 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =    136.59%,  TOTAL =     273.13%


          -----------------------------
          STATIONARY POINT LOCATION RUN
          -----------------------------

 OBTAINING INITIAL HESSIAN, HESS=GUESS   
 DIAGONAL GUESS HESSIAN IN INTERNAL COORDS IS
     1  0.4727     2  0.3532     3  0.3532     4  0.2500     5  0.2500
     6  0.2500     7  0.2500

          PARAMETERS CONTROLLING GEOMETRY SEARCH ARE
          METHOD =QA                  UPHESS =BFGS    
          NNEG   =         0          NFRZ   =         0
          NSTEP  =        20          IFOLOW =         1
          HESS   =GUESS               RESTAR =         F
          IHREP  =         0          HSSEND =         F
          NPRT   =         0          NPUN   =         0
          OPTTOL = 1.000E-04          RMIN   = 1.500E-03
          RMAX   = 1.000E-01          RLIM   = 7.000E-02
          DXMAX  = 3.000E-01          PURIFY =         F
          MOVIE  =         F          TRUPD  =         T
          TRMAX  = 5.000E-01          TRMIN  = 5.000E-02
          ITBMAT =        10          STPT   =         F
          STSTEP = 1.000E-02          PROJCT=          T
          MAXDII =        20          NSKIP  =         2
1NSERCH=   0

 COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000   0.7000000000
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.7800000000

 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000  -0.7000000000
 CARBON      6.0   0.0000000000   0.0000000000   0.7000000000
 HYDROGEN    1.0   0.0000000000   0.0000000000  -1.7800000000
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.7800000000


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -              COORDINATE      COORDINATE
 NO.     TYPE      I    J    K    L    M    N        (BOHR,RAD)       (ANG,DEG)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           2.6456164       1.4000000
     2 STRETCH     1    3                           2.0409041       1.0800000
     3 STRETCH     2    4                           2.0409041       1.0800000
     4 LIN.BEND    1    2    4                      3.1415927     180.0000000
     5 LIN.BEND    1    2    4                      3.1415927     180.0000000
     6 LIN.BEND    2    1    3                      3.1415927     180.0000000
     7 LIN.BEND    2    1    3                      3.1415927     180.0000000

          ********************
          1 ELECTRON INTEGRALS
          ********************
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     250.25%

          -------------
          GUESS OPTIONS
          -------------
          GUESS =HUCKEL            NORB  =       0          NORDER=       0
          MIX   =       F          PRTMO =       F          SYMDEN=       F
          TOLZ  = 1.0E-08          TOLE  = 1.0E-05

 INITIAL GUESS ORBITALS GENERATED BY HUCKEL   ROUTINE.
 HUCKEL GUESS REQUIRES      4072 WORDS.

 SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW.   BOTH SET(S).
     7 ORBITALS ARE OCCUPIED (    2 CORE ORBITALS).
     3=A1G      4=A2U      5=A1G      6=EU       7=EU       8=EG       9=EG  
    10=A2U     11=A1G     12=A2U 
 ...... END OF INITIAL ORBITAL SELECTION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     227.08%

          --------------------
          2 ELECTRON INTEGRALS
          --------------------

 THE -PK- OPTION IS OFF, THE INTEGRALS ARE NOT IN SUPERMATRIX FORM.
 STORING    4998 INTEGRALS/RECORD ON DISK, USING 12 BYTES/INTEGRAL.
 TWO ELECTRON INTEGRAL EVALUATION REQUIRES   34429 WORDS OF MEMORY.
 SCHWARZ INEQUALITY OVERHEAD:        78 INTEGRALS, T=        0.00
 II,JST,KST,LST =  1  1  1  1 NREC =         1 INTLOC =    1
 II,JST,KST,LST =  2  1  1  1 NREC =         1 INTLOC =    1
 II,JST,KST,LST =  3  1  1  1 NREC =         1 INTLOC =    1
 II,JST,KST,LST =  4  1  1  1 NREC =         1 INTLOC =    4
 II,JST,KST,LST =  5  1  1  1 NREC =         1 INTLOC =  305
 II,JST,KST,LST =  6  1  1  1 NREC =         1 INTLOC =  305
 SCHWARZ INEQUALITY TEST SKIPPED         3 INTEGRAL BLOCKS.
 TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =                 634
          1 INTEGRAL RECORDS WERE STORED ON DISK FILE  8.
 ...... END OF TWO-ELECTRON INTEGRALS .....

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     181.42%

          -------------------
          RHF SCF CALCULATION
          -------------------

     NUCLEAR ENERGY =        22.1963425659
     MAXIT =   30     NPUNCH=    2
     EXTRAP=T  DAMP=F  SHIFT=F  RSTRCT=F  DIIS=F  DEM=F  SOSCF=T
     DENSITY CONV=  1.00E-05
     SOSCF WILL OPTIMIZE      35 ORBITAL ROTATIONS, SOGTOL=   0.250
     MEMORY REQUIRED FOR RHF STEP=      8813 WORDS.

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE     ORB. GRAD
   1  0  0   -73.449308983   -73.449308983   0.109062893   0.000000000
          ---------------START SECOND ORDER SCF---------------
   2  1  0   -73.503292358    -0.053983375   0.019983731   0.019966228
   3  2  0   -73.503682554    -0.000390196   0.005215029   0.002167368
   4  3  0   -73.503697280    -0.000014726   0.000326539   0.000370492
   5  4  0   -73.503697469    -0.000000189   0.000048837   0.000059138
   6  5  0   -73.503697472    -0.000000004   0.000002034   0.000002012
   7  6  0   -73.503697472     0.000000000   0.000000387   0.000000210

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -73.5036974723 AFTER   7 ITERATIONS

          ------------
          EIGENVECTORS
          ------------

                      1          2          3          4          5
                  -10.5911   -10.5897    -0.8499    -0.6785    -0.5528
                     A2U        A1G        A1G        A2U        A1G 
    1  C   1  S   0.694565   0.696395  -0.207795  -0.140682  -0.020381
    2  C   1  S   0.054894   0.035829   0.525338   0.359296   0.038641
    3  C   1  X   0.000000   0.000000   0.000000   0.000000   0.000000
    4  C   1  Y   0.000000   0.000000   0.000000   0.000000   0.000000
    5  C   1  Z   0.009198  -0.006212   0.089835  -0.233197   0.435859
    6  C   2  S  -0.694565   0.696395  -0.207795   0.140682  -0.020381
    7  C   2  S  -0.054894   0.035829   0.525338  -0.359296   0.038641
    8  C   2  X   0.000000   0.000000   0.000000   0.000000   0.000000
    9  C   2  Y   0.000000   0.000000   0.000000   0.000000   0.000000
   10  C   2  Z   0.009198   0.006212  -0.089835  -0.233197  -0.435859
   11  H   3  S  -0.009768  -0.010735   0.165634   0.352692  -0.348969
   12  H   4  S   0.009768  -0.010735   0.165634  -0.352692  -0.348969

                      6          7          8          9         10
                   -0.2516    -0.2516     0.3618     0.3618     0.4782
                     EU         EU         EG         EG         A2U 
    1  C   1  S   0.000000   0.000000   0.000000   0.000000  -0.227829
    2  C   1  S   0.000000   0.000000   0.000000   0.000000   1.062588
    3  C   1  X   0.641670   0.000000   0.797760   0.000000   0.000000
    4  C   1  Y   0.000000   0.641670   0.000000   0.797760   0.000000
    5  C   1  Z   0.000000   0.000000   0.000000   0.000000   0.091147
    6  C   2  S   0.000000   0.000000   0.000000   0.000000   0.227829
    7  C   2  S   0.000000   0.000000   0.000000   0.000000  -1.062588
    8  C   2  X   0.641670   0.000000  -0.797760   0.000000   0.000000
    9  C   2  Y   0.000000   0.641670   0.000000  -0.797760   0.000000
   10  C   2  Z   0.000000   0.000000   0.000000   0.000000   0.091147
   11  H   3  S   0.000000   0.000000   0.000000   0.000000  -0.700408
   12  H   4  S   0.000000   0.000000   0.000000   0.000000   0.700408

                     11         12
                    0.7498     1.3028
                     A1G        A2U 
    1  C   1  S   0.120331  -0.074450
    2  C   1  S  -0.666863   0.461626
    3  C   1  X   0.000000   0.000000
    4  C   1  Y   0.000000   0.000000
    5  C   1  Z   0.647908   1.219721
    6  C   2  S   0.120331   0.074450
    7  C   2  S  -0.666863  -0.461626
    8  C   2  X   0.000000   0.000000
    9  C   2  Y   0.000000   0.000000
   10  C   2  Z  -0.647908   1.219721
   11  H   3  S   0.897733   0.589874
   12  H   4  S   0.897733  -0.589874
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.02 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =    323.52%,  TOTAL =     212.54%


                         ------------------------------
                         properties for the RHF density
                         ------------------------------

          -----------------
          ENERGY COMPONENTS
          -----------------

         WAVEFUNCTION NORMALIZATION =       1.0000000000

                ONE ELECTRON ENERGY =    -143.8697110570
                TWO ELECTRON ENERGY =      48.1696710187
           NUCLEAR REPULSION ENERGY =      22.1963425659
                                      ------------------
                       TOTAL ENERGY =     -73.5036974723

 ELECTRON-ELECTRON POTENTIAL ENERGY =      48.1696710187
  NUCLEUS-ELECTRON POTENTIAL ENERGY =    -216.0280545139
   NUCLEUS-NUCLEUS POTENTIAL ENERGY =      22.1963425659
                                      ------------------
             TOTAL POTENTIAL ENERGY =    -145.6620409293
               TOTAL KINETIC ENERGY =      72.1583434570
                 VIRIAL RATIO (V/T) =       2.0186444692

  ...... PI ENERGY ANALYSIS ......

 ENERGY ANALYSIS:
            FOCK ENERGY=    -47.5303690870
          BARE H ENERGY=   -143.8697110570
     ELECTRONIC ENERGY =    -95.7000400720
         KINETIC ENERGY=     72.1583434570
          N-N REPULSION=     22.1963425659
           TOTAL ENERGY=    -73.5036975061
        SIGMA PART(1+2)=    -82.4325227136
               (K,V1,2)=     67.1785217856   -185.5195470405     35.9085025413
           PI PART(1+2)=    -13.2675173584
               (K,V1,2)=      4.9798216714    -30.5085074735     12.2611684437
  SIGMA SKELETON, ERROR=    -60.2361801476      0.0000000000
             MIXED PART= 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00
 ...... END OF PI ENERGY ANALYSIS ......

          ---------------------------------------
          MULLIKEN AND LOWDIN POPULATION ANALYSES
          ---------------------------------------

     MULLIKEN ATOMIC POPULATION IN EACH MOLECULAR ORBITAL

                      1          2          3          4          5

                  2.000000   2.000000   2.000000   2.000000   2.000000

    1             1.001016   1.001165   0.868223   0.554377   0.645686
    2             1.001016   1.001165   0.868223   0.554377   0.645686
    3            -0.001016  -0.001165   0.131777   0.445623   0.354314
    4            -0.001016  -0.001165   0.131777   0.445623   0.354314

                      6          7

                  2.000000   2.000000

    1             1.000000   1.000000
    2             1.000000   1.000000
    3             0.000000   0.000000
    4             0.000000   0.000000

               ----- POPULATIONS IN EACH AO -----
                             MULLIKEN      LOWDIN
              1  C   1  S     1.98375     1.97631
              2  C   1  S     1.17763     1.10745
              3  C   1  X     1.00000     1.00000
              4  C   1  Y     1.00000     1.00000
              5  C   1  Z     0.90909     0.95938
              6  C   2  S     1.98375     1.97631
              7  C   2  S     1.17763     1.10745
              8  C   2  X     1.00000     1.00000
              9  C   2  Y     1.00000     1.00000
             10  C   2  Z     0.90909     0.95938
             11  H   3  S     0.92953     0.95686
             12  H   4  S     0.92953     0.95686

          ----- MULLIKEN ATOMIC OVERLAP POPULATIONS -----
          (OFF-DIAGONAL ELEMENTS NEED TO BE MULTIPLIED BY 2)

             1           2           3           4

    1    4.9300951
    2    0.7585908   4.9300951
    3    0.3930396  -0.0112588   0.5476331
    4   -0.0112588   0.3930396   0.0001194   0.5476331

          TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS
       ATOM         MULL.POP.    CHARGE          LOW.POP.     CHARGE
    1 CARBON        6.070467   -0.070467         6.043141   -0.043141
    2 CARBON        6.070467   -0.070467         6.043141   -0.043141
    3 HYDROGEN      0.929533    0.070467         0.956859    0.043141
    4 HYDROGEN      0.929533    0.070467         0.956859    0.043141

          -------------------------------
          BOND ORDER AND VALENCE ANALYSIS     BOND ORDER THRESHOLD=0.050
          -------------------------------

                   BOND                       BOND                       BOND
  ATOM PAIR DIST  ORDER      ATOM PAIR DIST  ORDER      ATOM PAIR DIST  ORDER
    1   2  1.400  2.970        1   3  1.080  0.979        2   4  1.080  0.979

                       TOTAL       BONDED        FREE
      ATOM            VALENCE     VALENCE     VALENCE
    1 CARBON            3.963       3.963       0.000
    2 CARBON            3.963       3.963       0.000
    3 HYDROGEN          0.995       0.995       0.000
    4 HYDROGEN          0.995       0.995       0.000

          ---------------------
          ELECTROSTATIC MOMENTS
          ---------------------

 POINT   1           X           Y           Z (BOHR)    CHARGE
                 0.000000    0.000000    0.000000        0.00 (A.U.)
         DX          DY          DZ         /D/  (DEBYE)
     0.000000    0.000000    0.000000    0.000000
 ...... END OF PROPERTY EVALUATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     202.13%
 ......END OF NBO ANALYSIS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     194.36%

 BEGINNING ONE ELECTRON GRADIENT...
 ..... END OF 1-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     182.49%

          ----------------------
          GRADIENT OF THE ENERGY
          ----------------------
 SCHWARZ SCREENING SKIPPED          4 BLOCKS, COMPUTED        112 BLOCKS

 ...... END OF 2-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     155.31%

          NSERCH=  0     ENERGY=     -73.5036975

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  CARBON       6.0     0.0000000     0.0000000    -0.3975936
  2  CARBON       6.0     0.0000000     0.0000000     0.3975936
  3  HYDROGEN     1.0     0.0000000     0.0000000    -0.0010894
  4  HYDROGEN     1.0     0.0000000     0.0000000     0.0010894


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -            COORDINATE        GRADIENT
 NO.     TYPE      I    J    K    L    M    N       (ANG,DEG)     (H/B,H/RAD)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           1.4000000       0.3986830
     2 STRETCH     1    3                           1.0800000       0.0010894
     3 STRETCH     2    4                           1.0800000       0.0010894
     4 LIN.BEND    1    2    4                    180.0000000       0.0000000
     5 LIN.BEND    1    2    4                    180.0000000       0.0000000
     6 LIN.BEND    2    1    3                    180.0000000       0.0000000
     7 LIN.BEND    2    1    3                    180.0000000       0.0000000

          MAXIMUM GRADIENT =  0.3986830    RMS GRADIENT = 0.1506891
          FORCE CONSTANT MATRIX NOT UPDATED --- TAKING FIRST STEP
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.843360
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.85621714
          TRIM/QA STEP HAS LENGTH         =   0.300000
          RADIUS OF STEP TAKEN=   0.30000  CURRENT TRUST RADIUS=   0.30000
          TRANSFORMING DISPLACEMENT FROM INTERNALS TO CARTESIANS
          THE ROOT MEAN SQUARE ERROR IN ITERATION   1 IS   0.00000000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     150.19%

1NSERCH=   1

 COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000   0.6206241283
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.7001474454

 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000  -0.6206241283
 CARBON      6.0   0.0000000000   0.0000000000   0.6206241283
 HYDROGEN    1.0   0.0000000000   0.0000000000  -1.7001474454
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.7001474454


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -              COORDINATE      COORDINATE
 NO.     TYPE      I    J    K    L    M    N        (BOHR,RAD)       (ANG,DEG)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           2.3456191       1.2412483
     2 STRETCH     1    3                           2.0400033       1.0795233
     3 STRETCH     2    4                           2.0400033       1.0795233
     4 LIN.BEND    1    2    4                      3.1415927     180.0000000
     5 LIN.BEND    1    2    4                      3.1415927     180.0000000
     6 LIN.BEND    2    1    3                      3.1415927     180.0000000
     7 LIN.BEND    2    1    3                      3.1415927     180.0000000

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    CARBON         CARBON         HYDROGEN       HYDROGEN  

  1  CARBON          0.0000000      1.2412483 *    1.0795233 *    2.3207716 *  
  2  CARBON          1.2412483 *    0.0000000      2.3207716 *    1.0795233 *  
  3  HYDROGEN        1.0795233 *    2.3207716 *    0.0000000      3.4002949    
  4  HYDROGEN        2.3207716 *    1.0795233 *    3.4002949      0.0000000    

  * ... LESS THAN  3.000

 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     147.63%
 TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =                 636
          1 INTEGRAL RECORDS WERE STORED ON DISK FILE  8.
 ...... END OF TWO-ELECTRON INTEGRALS .....

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.01 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     120.42%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE     ORB. GRAD
          ---------------START SECOND ORDER SCF---------------
   1  0  0   -73.591291154   -73.591291154   0.042990734   0.039546796
   2  1  0   -73.593943016    -0.002651862   0.010290820   0.003869111
   3  2  0   -73.594016684    -0.000073668   0.000398647   0.000488951
   4  3  0   -73.594016993    -0.000000309   0.000026500   0.000054437
   5  4  0   -73.594016995    -0.000000002   0.000000599   0.000000537

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -73.5940169953 AFTER   5 ITERATIONS
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     111.36%
 ..... END OF 1-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     107.68%

 ...... END OF 2-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =      97.74%

          NSERCH=  1     ENERGY=     -73.5940170

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  CARBON       6.0     0.0000000     0.0000000    -0.1642137
  2  CARBON       6.0     0.0000000     0.0000000     0.1642137
  3  HYDROGEN     1.0     0.0000000     0.0000000    -0.0031594
  4  HYDROGEN     1.0     0.0000000     0.0000000     0.0031594


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -            COORDINATE        GRADIENT
 NO.     TYPE      I    J    K    L    M    N       (ANG,DEG)     (H/B,H/RAD)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           1.2412483       0.1673731
     2 STRETCH     1    3                           1.0795233       0.0031594
     3 STRETCH     2    4                           1.0795233       0.0031594
     4 LIN.BEND    1    2    4                    180.0000000       0.0000000
     5 LIN.BEND    1    2    4                    180.0000000       0.0000000
     6 LIN.BEND    2    1    3                    180.0000000       0.0000000
     7 LIN.BEND    2    1    3                    180.0000000       0.0000000

          MAXIMUM GRADIENT =  0.1673731    RMS GRADIENT = 0.0632836
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0903195230
          PREDICTED ENERGY CHANGE WAS  -0.0983326647 RATIO=  0.919
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.218226
          RADIUS OF STEP TAKEN=   0.21823  CURRENT TRUST RADIUS=   0.50000
          TRANSFORMING DISPLACEMENT FROM INTERNALS TO CARTESIANS
          THE ROOT MEAN SQUARE ERROR IN ITERATION   1 IS   0.00000000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =      96.07%

1NSERCH=   2

 COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000   0.5631166234
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6353157425

 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000  -0.5631166234
 CARBON      6.0   0.0000000000   0.0000000000   0.5631166234
 HYDROGEN    1.0   0.0000000000   0.0000000000  -1.6353157425
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6353157425


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -              COORDINATE      COORDINATE
 NO.     TYPE      I    J    K    L    M    N        (BOHR,RAD)       (ANG,DEG)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           2.1282722       1.1262332
     2 STRETCH     1    3                           2.0261625       1.0721991
     3 STRETCH     2    4                           2.0261625       1.0721991
     4 LIN.BEND    1    2    4                      3.1415927     180.0000000
     5 LIN.BEND    1    2    4                      3.1415927     180.0000000
     6 LIN.BEND    2    1    3                      3.1415927     180.0000000
     7 LIN.BEND    2    1    3                      3.1415927     180.0000000

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    CARBON         CARBON         HYDROGEN       HYDROGEN  

  1  CARBON          0.0000000      1.1262332 *    1.0721991 *    2.1984324 *  
  2  CARBON          1.1262332 *    0.0000000      2.1984324 *    1.0721991 *  
  3  HYDROGEN        1.0721991 *    2.1984324 *    0.0000000      3.2706315    
  4  HYDROGEN        2.1984324 *    1.0721991 *    3.2706315      0.0000000    

  * ... LESS THAN  3.000

 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.02 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =   2646.97%,  TOTAL =     126.57%
 TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =                 636
          1 INTEGRAL RECORDS WERE STORED ON DISK FILE  8.
 ...... END OF TWO-ELECTRON INTEGRALS .....

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.01 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     113.29%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE     ORB. GRAD
          ---------------START SECOND ORDER SCF---------------
   1  0  0   -73.593819887   -73.593819887   0.043424807   0.033632738
   2  1  0   -73.596122636    -0.002302749   0.009405529   0.004114442
   3  2  0   -73.596184365    -0.000061730   0.000365930   0.000386778
   4  3  0   -73.596184657    -0.000000292   0.000027542   0.000040713
   5  4  0   -73.596184659    -0.000000002   0.000000750   0.000000682

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -73.5961846590 AFTER   5 ITERATIONS
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     107.30%
 ..... END OF 1-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     104.72%

 ...... END OF 2-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =      97.33%

          NSERCH=  2     ENERGY=     -73.5961847

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  CARBON       6.0     0.0000000     0.0000000     0.1803572
  2  CARBON       6.0     0.0000000     0.0000000    -0.1803572
  3  HYDROGEN     1.0     0.0000000     0.0000000     0.0013949
  4  HYDROGEN     1.0     0.0000000     0.0000000    -0.0013949


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -            COORDINATE        GRADIENT
 NO.     TYPE      I    J    K    L    M    N       (ANG,DEG)     (H/B,H/RAD)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           1.1262332      -0.1817521
     2 STRETCH     1    3                           1.0721991      -0.0013949
     3 STRETCH     2    4                           1.0721991      -0.0013949
     4 LIN.BEND    1    2    4                    180.0000000       0.0000000
     5 LIN.BEND    1    2    4                    180.0000000       0.0000000
     6 LIN.BEND    2    1    3                    180.0000000       0.0000000
     7 LIN.BEND    2    1    3                    180.0000000       0.0000000

          MAXIMUM GRADIENT =  0.1817521    RMS GRADIENT = 0.0686999
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0021676637
          PREDICTED ENERGY CHANGE WAS  -0.0182327394 RATIO=  0.119
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.113318
          RADIUS OF STEP TAKEN=   0.11332  CURRENT TRUST RADIUS=   0.21823
          TRANSFORMING DISPLACEMENT FROM INTERNALS TO CARTESIANS
          THE ROOT MEAN SQUARE ERROR IN ITERATION   1 IS   0.00000000

 CPU        TIME:   STEP =      0.02 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =   1831.98%,  TOTAL =     120.07%

1NSERCH=   3

 COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000   0.5930532166
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6676032154

 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000  -0.5930532166
 CARBON      6.0   0.0000000000   0.0000000000   0.5930532166
 HYDROGEN    1.0   0.0000000000   0.0000000000  -1.6676032154
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6676032154


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -              COORDINATE      COORDINATE
 NO.     TYPE      I    J    K    L    M    N        (BOHR,RAD)       (ANG,DEG)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           2.2414162       1.1861064
     2 STRETCH     1    3                           2.0306051       1.0745500
     3 STRETCH     2    4                           2.0306051       1.0745500
     4 LIN.BEND    1    2    4                      3.1415927     180.0000000
     5 LIN.BEND    1    2    4                      3.1415927     180.0000000
     6 LIN.BEND    2    1    3                      3.1415927     180.0000000
     7 LIN.BEND    2    1    3                      3.1415927     180.0000000

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    CARBON         CARBON         HYDROGEN       HYDROGEN  

  1  CARBON          0.0000000      1.1861064 *    1.0745500 *    2.2606564 *  
  2  CARBON          1.1861064 *    0.0000000      2.2606564 *    1.0745500 *  
  3  HYDROGEN        1.0745500 *    2.2606564 *    0.0000000      3.3352064    
  4  HYDROGEN        2.2606564 *    1.0745500 *    3.3352064      0.0000000    

  * ... LESS THAN  3.000

 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     119.13%
 TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =                 636
          1 INTEGRAL RECORDS WERE STORED ON DISK FILE  8.
 ...... END OF TWO-ELECTRON INTEGRALS .....

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.01 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     109.55%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE     ORB. GRAD
          ---------------START SECOND ORDER SCF---------------
   1  0  0   -73.603706493   -73.603706493   0.020755917   0.019371414
   2  1  0   -73.604405342    -0.000698849   0.007717317   0.003088747
   3  2  0   -73.604443817    -0.000038475   0.000352023   0.000381058
   4  3  0   -73.604443960    -0.000000143   0.000033446   0.000041675
   5  4  0   -73.604443963    -0.000000002   0.000001371   0.000001579

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -73.6044439625 AFTER   5 ITERATIONS
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     105.06%
 ..... END OF 1-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     102.98%

 ...... END OF 2-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =      97.28%

          NSERCH=  3     ENERGY=     -73.6044440

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  CARBON       6.0     0.0000000     0.0000000    -0.0260843
  2  CARBON       6.0     0.0000000     0.0000000     0.0260843
  3  HYDROGEN     1.0     0.0000000     0.0000000     0.0006012
  4  HYDROGEN     1.0     0.0000000     0.0000000    -0.0006012


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -            COORDINATE        GRADIENT
 NO.     TYPE      I    J    K    L    M    N       (ANG,DEG)     (H/B,H/RAD)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           1.1861064       0.0254830
     2 STRETCH     1    3                           1.0745500      -0.0006012
     3 STRETCH     2    4                           1.0745500      -0.0006012
     4 LIN.BEND    1    2    4                    180.0000000       0.0000000
     5 LIN.BEND    1    2    4                    180.0000000       0.0000000
     6 LIN.BEND    2    1    3                    180.0000000       0.0000000
     7 LIN.BEND    2    1    3                    180.0000000       0.0000000

          MAXIMUM GRADIENT =  0.0254830    RMS GRADIENT = 0.0096370
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0082593035
          PREDICTED ENERGY CHANGE WAS  -0.0102882694 RATIO=  0.803
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.014045
          RADIUS OF STEP TAKEN=   0.01405  CURRENT TRUST RADIUS=   0.16026
          TRANSFORMING DISPLACEMENT FROM INTERNALS TO CARTESIANS
          THE ROOT MEAN SQUARE ERROR IN ITERATION   1 IS   0.00000000

 CPU        TIME:   STEP =      0.02 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =   1850.77%,  TOTAL =     115.53%

1NSERCH=   4

 COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000   0.5893759273
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6646841533

 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000  -0.5893759273
 CARBON      6.0   0.0000000000   0.0000000000   0.5893759273
 HYDROGEN    1.0   0.0000000000   0.0000000000  -1.6646841533
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6646841533


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -              COORDINATE      COORDINATE
 NO.     TYPE      I    J    K    L    M    N        (BOHR,RAD)       (ANG,DEG)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           2.2275180       1.1787519
     2 STRETCH     1    3                           2.0320379       1.0753082
     3 STRETCH     2    4                           2.0320379       1.0753082
     4 LIN.BEND    1    2    4                      3.1415927     180.0000000
     5 LIN.BEND    1    2    4                      3.1415927     180.0000000
     6 LIN.BEND    2    1    3                      3.1415927     180.0000000
     7 LIN.BEND    2    1    3                      3.1415927     180.0000000

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    CARBON         CARBON         HYDROGEN       HYDROGEN  

  1  CARBON          0.0000000      1.1787519 *    1.0753082 *    2.2540601 *  
  2  CARBON          1.1787519 *    0.0000000      2.2540601 *    1.0753082 *  
  3  HYDROGEN        1.0753082 *    2.2540601 *    0.0000000      3.3293683    
  4  HYDROGEN        2.2540601 *    1.0753082 *    3.3293683      0.0000000    

  * ... LESS THAN  3.000

 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     114.70%
 TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =                 636
          1 INTEGRAL RECORDS WERE STORED ON DISK FILE  8.
 ...... END OF TWO-ELECTRON INTEGRALS .....

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.01 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     107.26%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE     ORB. GRAD
          ---------------START SECOND ORDER SCF---------------
   1  0  0   -73.604636248   -73.604636248   0.002257750   0.002078926
   2  1  0   -73.604644517    -0.000008269   0.000716463   0.000294594
   3  2  0   -73.604644857    -0.000000340   0.000035764   0.000031201
   4  3  0   -73.604644858    -0.000000002   0.000003326   0.000003938
   5  4  0   -73.604644859     0.000000000   0.000000128   0.000000152

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -73.6046448585 AFTER   5 ITERATIONS
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     103.38%
 ..... END OF 1-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     101.77%

 ...... END OF 2-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.02 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =    350.24%,  TOTAL =     113.24%

          NSERCH=  4     ENERGY=     -73.6046449

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  CARBON       6.0     0.0000000     0.0000000    -0.0028922
  2  CARBON       6.0     0.0000000     0.0000000     0.0028922
  3  HYDROGEN     1.0     0.0000000     0.0000000    -0.0003410
  4  HYDROGEN     1.0     0.0000000     0.0000000     0.0003410


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -            COORDINATE        GRADIENT
 NO.     TYPE      I    J    K    L    M    N       (ANG,DEG)     (H/B,H/RAD)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           1.1787519       0.0032332
     2 STRETCH     1    3                           1.0753082       0.0003410
     3 STRETCH     2    4                           1.0753082       0.0003410
     4 LIN.BEND    1    2    4                    180.0000000       0.0000000
     5 LIN.BEND    1    2    4                    180.0000000       0.0000000
     6 LIN.BEND    2    1    3                    180.0000000       0.0000000
     7 LIN.BEND    2    1    3                    180.0000000       0.0000000

          MAXIMUM GRADIENT =  0.0032332    RMS GRADIENT = 0.0012356
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0002008960
          PREDICTED ENERGY CHANGE WAS  -0.0001779448 RATIO=  1.129
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.002622
          RADIUS OF STEP TAKEN=   0.00262  CURRENT TRUST RADIUS=   0.05000
          TRANSFORMING DISPLACEMENT FROM INTERNALS TO CARTESIANS
          THE ROOT MEAN SQUARE ERROR IN ITERATION   1 IS   0.00000000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     112.28%

1NSERCH=   5

 COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000   0.5888277365
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6635349260

 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000  -0.5888277365
 CARBON      6.0   0.0000000000   0.0000000000   0.5888277365
 HYDROGEN    1.0   0.0000000000   0.0000000000  -1.6635349260
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6635349260


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -              COORDINATE      COORDINATE
 NO.     TYPE      I    J    K    L    M    N        (BOHR,RAD)       (ANG,DEG)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           2.2254462       1.1776555
     2 STRETCH     1    3                           2.0309021       1.0747072
     3 STRETCH     2    4                           2.0309021       1.0747072
     4 LIN.BEND    1    2    4                      3.1415927     180.0000000
     5 LIN.BEND    1    2    4                      3.1415927     180.0000000
     6 LIN.BEND    2    1    3                      3.1415927     180.0000000
     7 LIN.BEND    2    1    3                      3.1415927     180.0000000

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    CARBON         CARBON         HYDROGEN       HYDROGEN  

  1  CARBON          0.0000000      1.1776555 *    1.0747072 *    2.2523627 *  
  2  CARBON          1.1776555 *    0.0000000      2.2523627 *    1.0747072 *  
  3  HYDROGEN        1.0747072 *    2.2523627 *    0.0000000      3.3270699    
  4  HYDROGEN        2.2523627 *    1.0747072 *    3.3270699      0.0000000    

  * ... LESS THAN  3.000

 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     111.59%
 TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =                 636
          1 INTEGRAL RECORDS WERE STORED ON DISK FILE  8.
 ...... END OF TWO-ELECTRON INTEGRALS .....

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.01 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     105.36%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE     ORB. GRAD
          ---------------START SECOND ORDER SCF---------------
   1  0  0   -73.604647846   -73.604647846   0.000514922   0.000442272
   2  1  0   -73.604648188    -0.000000342   0.000163224   0.000062191
   3  2  0   -73.604648204    -0.000000016   0.000004715   0.000006682
   4  3  0   -73.604648204     0.000000000   0.000000588   0.000000638
   5  4  0   -73.604648204     0.000000000   0.000000017   0.000000022

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -73.6046482043 AFTER   5 ITERATIONS
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     102.28%
 ..... END OF 1-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.88%

 ...... END OF 2-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.02 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =    355.09%,  TOTAL =     110.79%

          NSERCH=  5     ENERGY=     -73.6046482

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  CARBON       6.0     0.0000000     0.0000000    -0.0001163
  2  CARBON       6.0     0.0000000     0.0000000     0.0001163
  3  HYDROGEN     1.0     0.0000000     0.0000000     0.0002372
  4  HYDROGEN     1.0     0.0000000     0.0000000    -0.0002372


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -            COORDINATE        GRADIENT
 NO.     TYPE      I    J    K    L    M    N       (ANG,DEG)     (H/B,H/RAD)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           1.1776555      -0.0001209
     2 STRETCH     1    3                           1.0747072      -0.0002372
     3 STRETCH     2    4                           1.0747072      -0.0002372
     4 LIN.BEND    1    2    4                    180.0000000       0.0000000
     5 LIN.BEND    1    2    4                    180.0000000       0.0000000
     6 LIN.BEND    2    1    3                    180.0000000       0.0000000
     7 LIN.BEND    2    1    3                    180.0000000       0.0000000

          MAXIMUM GRADIENT =  0.0002372    RMS GRADIENT = 0.0001348
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0000033458
          PREDICTED ENERGY CHANGE WAS  -0.0000037367 RATIO=  0.895
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.000819
          RADIUS OF STEP TAKEN=   0.00082  CURRENT TRUST RADIUS=   0.05000
          TRANSFORMING DISPLACEMENT FROM INTERNALS TO CARTESIANS
          THE ROOT MEAN SQUARE ERROR IN ITERATION   1 IS   0.00000000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     109.91%

1NSERCH=   6

 COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000   0.5888372303
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6638504869

 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000  -0.5888372303
 CARBON      6.0   0.0000000000   0.0000000000   0.5888372303
 HYDROGEN    1.0   0.0000000000   0.0000000000  -1.6638504869
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6638504869


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -              COORDINATE      COORDINATE
 NO.     TYPE      I    J    K    L    M    N        (BOHR,RAD)       (ANG,DEG)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           2.2254820       1.1776745
     2 STRETCH     1    3                           2.0314805       1.0750133
     3 STRETCH     2    4                           2.0314805       1.0750133
     4 LIN.BEND    1    2    4                      3.1415927     180.0000000
     5 LIN.BEND    1    2    4                      3.1415927     180.0000000
     6 LIN.BEND    2    1    3                      3.1415927     180.0000000
     7 LIN.BEND    2    1    3                      3.1415927     180.0000000

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    CARBON         CARBON         HYDROGEN       HYDROGEN  

  1  CARBON          0.0000000      1.1776745 *    1.0750133 *    2.2526877 *  
  2  CARBON          1.1776745 *    0.0000000      2.2526877 *    1.0750133 *  
  3  HYDROGEN        1.0750133 *    2.2526877 *    0.0000000      3.3277010    
  4  HYDROGEN        2.2526877 *    1.0750133 *    3.3277010      0.0000000    

  * ... LESS THAN  3.000

 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     109.18%
 TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =                 636
          1 INTEGRAL RECORDS WERE STORED ON DISK FILE  8.
 ...... END OF TWO-ELECTRON INTEGRALS .....

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.01 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     104.01%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE     ORB. GRAD
          ---------------START SECOND ORDER SCF---------------
   1  0  0   -73.604648296   -73.604648296   0.000081772   0.000061778
   2  1  0   -73.604648304    -0.000000007   0.000024662   0.000008435
   3  2  0   -73.604648304     0.000000000   0.000000634   0.000000991
   4  3  0   -73.604648304     0.000000000   0.000000085   0.000000098

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -73.6046483041 AFTER   4 ITERATIONS
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     101.76%
 ..... END OF 1-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.59%

 ...... END OF 2-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.02 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =    350.51%,  TOTAL =     109.24%

          NSERCH=  6     ENERGY=     -73.6046483

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  CARBON       6.0     0.0000000     0.0000000     0.0001479
  2  CARBON       6.0     0.0000000     0.0000000    -0.0001479
  3  HYDROGEN     1.0     0.0000000     0.0000000    -0.0000707
  4  HYDROGEN     1.0     0.0000000     0.0000000     0.0000707


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -            COORDINATE        GRADIENT
 NO.     TYPE      I    J    K    L    M    N       (ANG,DEG)     (H/B,H/RAD)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           1.1776745      -0.0000772
     2 STRETCH     1    3                           1.0750133       0.0000707
     3 STRETCH     2    4                           1.0750133       0.0000707
     4 LIN.BEND    1    2    4                    180.0000000       0.0000000
     5 LIN.BEND    1    2    4                    180.0000000       0.0000000
     6 LIN.BEND    2    1    3                    180.0000000       0.0000000
     7 LIN.BEND    2    1    3                    180.0000000       0.0000000

          MAXIMUM GRADIENT =  0.0000772    RMS GRADIENT = 0.0000477
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0000000998
          PREDICTED ENERGY CHANGE WAS  -0.0000001394 RATIO=  0.716
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.000193
          RADIUS OF STEP TAKEN=   0.00019  CURRENT TRUST RADIUS=   0.05000
          TRANSFORMING DISPLACEMENT FROM INTERNALS TO CARTESIANS
          THE ROOT MEAN SQUARE ERROR IN ITERATION   1 IS   0.00000000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     108.41%

1NSERCH=   7

 COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000   0.5888503322
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6637938208

 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000  -0.5888503322
 CARBON      6.0   0.0000000000   0.0000000000   0.5888503322
 HYDROGEN    1.0   0.0000000000   0.0000000000  -1.6637938208
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6637938208


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -              COORDINATE      COORDINATE
 NO.     TYPE      I    J    K    L    M    N        (BOHR,RAD)       (ANG,DEG)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           2.2255316       1.1777007
     2 STRETCH     1    3                           2.0313486       1.0749435
     3 STRETCH     2    4                           2.0313486       1.0749435
     4 LIN.BEND    1    2    4                      3.1415927     180.0000000
     5 LIN.BEND    1    2    4                      3.1415927     180.0000000
     6 LIN.BEND    2    1    3                      3.1415927     180.0000000
     7 LIN.BEND    2    1    3                      3.1415927     180.0000000

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    CARBON         CARBON         HYDROGEN       HYDROGEN  

  1  CARBON          0.0000000      1.1777007 *    1.0749435 *    2.2526442 *  
  2  CARBON          1.1777007 *    0.0000000      2.2526442 *    1.0749435 *  
  3  HYDROGEN        1.0749435 *    2.2526442 *    0.0000000      3.3275876    
  4  HYDROGEN        2.2526442 *    1.0749435 *    3.3275876      0.0000000    

  * ... LESS THAN  3.000

 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     107.99%
 TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =                 636
          1 INTEGRAL RECORDS WERE STORED ON DISK FILE  8.
 ...... END OF TWO-ELECTRON INTEGRALS .....

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.01 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     103.47%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE     ORB. GRAD
          ---------------START SECOND ORDER SCF---------------
   1  0  0   -73.604648315   -73.604648315   0.000008422   0.000018109
   2  1  0   -73.604648315     0.000000000   0.000002077   0.000000942
   3  2  0   -73.604648315     0.000000000   0.000000236   0.000000196

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -73.6046483151 AFTER   3 ITERATIONS
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.77%
 ..... END OF 1-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =      99.70%

 ...... END OF 2-ELECTRON GRADIENT ......

 CPU        TIME:   STEP =      0.02 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =    348.87%,  TOTAL =     107.37%

          NSERCH=  7     ENERGY=     -73.6046483

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  CARBON       6.0     0.0000000     0.0000000    -0.0000075
  2  CARBON       6.0     0.0000000     0.0000000     0.0000075
  3  HYDROGEN     1.0     0.0000000     0.0000000     0.0000002
  4  HYDROGEN     1.0     0.0000000     0.0000000    -0.0000002


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -            COORDINATE        GRADIENT
 NO.     TYPE      I    J    K    L    M    N       (ANG,DEG)     (H/B,H/RAD)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           1.1777007       0.0000073
     2 STRETCH     1    3                           1.0749435      -0.0000002
     3 STRETCH     2    4                           1.0749435      -0.0000002
     4 LIN.BEND    1    2    4                    180.0000000       0.0000000
     5 LIN.BEND    1    2    4                    180.0000000       0.0000000
     6 LIN.BEND    2    1    3                    180.0000000       0.0000000
     7 LIN.BEND    2    1    3                    180.0000000       0.0000000

          MAXIMUM GRADIENT =  0.0000073    RMS GRADIENT = 0.0000028
1     ***** EQUILIBRIUM GEOMETRY LOCATED *****

 Acetylene geometry optimization in internal coordinates                         
 COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000   0.5888503322
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6637938208

 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 CARBON      6.0   0.0000000000   0.0000000000  -0.5888503322
 CARBON      6.0   0.0000000000   0.0000000000   0.5888503322
 HYDROGEN    1.0   0.0000000000   0.0000000000  -1.6637938208
 HYDROGEN    1.0   0.0000000000   0.0000000000   1.6637938208


                     --------------------
                     INTERNAL COORDINATES
                     --------------------

                          - - ATOMS - -              COORDINATE      COORDINATE
 NO.     TYPE      I    J    K    L    M    N        (BOHR,RAD)       (ANG,DEG)
 ------------------------------------------------------------------------------
     1 STRETCH     1    2                           2.2255316       1.1777007
     2 STRETCH     1    3                           2.0313486       1.0749435
     3 STRETCH     2    4                           2.0313486       1.0749435
     4 LIN.BEND    1    2    4                      3.1415927     180.0000000
     5 LIN.BEND    1    2    4                      3.1415927     180.0000000
     6 LIN.BEND    2    1    3                      3.1415927     180.0000000
     7 LIN.BEND    2    1    3                      3.1415927     180.0000000

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    CARBON         CARBON         HYDROGEN       HYDROGEN  

  1  CARBON          0.0000000      1.1777007 *    1.0749435 *    2.2526442 *  
  2  CARBON          1.1777007 *    0.0000000      2.2526442 *    1.0749435 *  
  3  HYDROGEN        1.0749435 *    2.2526442 *    0.0000000      3.3275876    
  4  HYDROGEN        2.2526442 *    1.0749435 *    3.3275876      0.0000000    

  * ... LESS THAN  3.000


          NUCLEAR ENERGY    =       25.0613094487
          ELECTRONIC ENERGY =      -98.6659577638
          TOTAL ENERGY      =      -73.6046483151

          ------------------
          MOLECULAR ORBITALS
          ------------------

                      1          2          3          4          5
                  -10.6116   -10.6091    -0.9304    -0.6757    -0.5880
                     A1G        A2U        A1G        A2U        A1G 
    1  C   1  S   0.695522   0.694586  -0.220393  -0.134164  -0.000042
    2  C   1  S   0.034586   0.066095   0.483769   0.324041  -0.065440
    3  C   1  X   0.000000   0.000000   0.000000   0.000000   0.000000
    4  C   1  Y   0.000000   0.000000   0.000000   0.000000   0.000000
    5  C   1  Z  -0.002464   0.013467   0.171676  -0.258635   0.440249
    6  C   2  S   0.695522  -0.694586  -0.220393   0.134164  -0.000042
    7  C   2  S   0.034586  -0.066095   0.483769  -0.324041  -0.065440
    8  C   2  X   0.000000   0.000000   0.000000   0.000000   0.000000
    9  C   2  Y   0.000000   0.000000   0.000000   0.000000   0.000000
   10  C   2  Z   0.002464   0.013467  -0.171676  -0.258635  -0.440249
   11  H   3  S  -0.009030  -0.010675   0.129637   0.356497  -0.354979
   12  H   4  S  -0.009030   0.010675   0.129637  -0.356497  -0.354979

                      6          7          8          9         10
                   -0.3377    -0.3377     0.4484     0.4484     0.5674
                     EU         EU         EG         EG         A2U 
    1  C   1  S   0.000000   0.000000   0.000000   0.000000  -0.220768
    2  C   1  S   0.000000   0.000000   0.000000   0.000000   1.156159
    3  C   1  X   0.616462   0.000000   0.854796   0.000000   0.000000
    4  C   1  Y   0.000000   0.616462   0.000000   0.854796   0.000000
    5  C   1  Z   0.000000   0.000000   0.000000   0.000000  -0.026792
    6  C   2  S   0.000000   0.000000   0.000000   0.000000   0.220768
    7  C   2  S   0.000000   0.000000   0.000000   0.000000  -1.156159
    8  C   2  X   0.616462   0.000000  -0.854796   0.000000   0.000000
    9  C   2  Y   0.000000   0.616462   0.000000  -0.854796   0.000000
   10  C   2  Z   0.000000   0.000000   0.000000   0.000000  -0.026792
   11  H   3  S   0.000000   0.000000   0.000000   0.000000  -0.788762
   12  H   4  S   0.000000   0.000000   0.000000   0.000000   0.788762

                     11         12
                    0.7503     1.6599
                     A1G        A2U 
    1  C   1  S   0.113967  -0.126627
    2  C   1  S  -0.682166   1.146495
    3  C   1  X   0.000000   0.000000
    4  C   1  Y   0.000000   0.000000
    5  C   1  Z   0.659269   1.491827
    6  C   2  S   0.113967   0.126627
    7  C   2  S  -0.682166  -1.146495
    8  C   2  X   0.000000   0.000000
    9  C   2  Y   0.000000   0.000000
   10  C   2  Z  -0.659269   1.491827
   11  H   3  S   0.906345   0.476361
   12  H   4  S   0.906345  -0.476361


                         ------------------------------
                         properties for the RHF density
                         ------------------------------

          -----------------
          ENERGY COMPONENTS
          -----------------

         WAVEFUNCTION NORMALIZATION =       1.0000000000

                ONE ELECTRON ENERGY =    -149.1514355874
                TWO ELECTRON ENERGY =      50.4854778235
           NUCLEAR REPULSION ENERGY =      25.0613094487
                                      ------------------
                       TOTAL ENERGY =     -73.6046483151

 ELECTRON-ELECTRON POTENTIAL ENERGY =      50.4854778235
  NUCLEUS-ELECTRON POTENTIAL ENERGY =    -221.6237407379
   NUCLEUS-NUCLEUS POTENTIAL ENERGY =      25.0613094487
                                      ------------------
             TOTAL POTENTIAL ENERGY =    -146.0769534657
               TOTAL KINETIC ENERGY =      72.4723051505
                 VIRIAL RATIO (V/T) =       2.0156244949

  ...... PI ENERGY ANALYSIS ......

 ENERGY ANALYSIS:
            FOCK ENERGY=    -48.1804800578
          BARE H ENERGY=   -149.1514355874
     ELECTRONIC ENERGY =    -98.6659578226
         KINETIC ENERGY=     72.4723051505
          N-N REPULSION=     25.0613094487
           TOTAL ENERGY=    -73.6046483738
        SIGMA PART(1+2)=    -84.4762118307
               (K,V1,2)=     67.6284747354   -189.7511460671     37.6464595011
           PI PART(1+2)=    -14.1897459919
               (K,V1,2)=      4.8438304152    -31.8725946708     12.8390182637
  SIGMA SKELETON, ERROR=    -59.4149023820      0.0000000000
             MIXED PART= 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00
 ...... END OF PI ENERGY ANALYSIS ......

          ---------------------------------------
          MULLIKEN AND LOWDIN POPULATION ANALYSES
          ---------------------------------------

     MULLIKEN ATOMIC POPULATION IN EACH MOLECULAR ORBITAL

                      1          2          3          4          5

                  2.000000   2.000000   2.000000   2.000000   2.000000

    1             1.001014   1.001094   0.915506   0.547299   0.604449
    2             1.001014   1.001094   0.915506   0.547299   0.604449
    3            -0.001014  -0.001094   0.084494   0.452701   0.395551
    4            -0.001014  -0.001094   0.084494   0.452701   0.395551

                      6          7

                  2.000000   2.000000

    1             1.000000   1.000000
    2             1.000000   1.000000
    3             0.000000   0.000000
    4             0.000000   0.000000

               ----- POPULATIONS IN EACH AO -----
                             MULLIKEN      LOWDIN
              1  C   1  S     1.98383     1.97896
              2  C   1  S     1.09313     1.02552
              3  C   1  X     1.00000     1.00000
              4  C   1  Y     1.00000     1.00000
              5  C   1  Z     0.99240     1.03833
              6  C   2  S     1.98383     1.97896
              7  C   2  S     1.09313     1.02552
              8  C   2  X     1.00000     1.00000
              9  C   2  Y     1.00000     1.00000
             10  C   2  Z     0.99240     1.03833
             11  H   3  S     0.93064     0.95719
             12  H   4  S     0.93064     0.95719

          ----- MULLIKEN ATOMIC OVERLAP POPULATIONS -----
          (OFF-DIAGONAL ELEMENTS NEED TO BE MULTIPLIED BY 2)

             1           2           3           4

    1    4.7775980
    2    0.9014781   4.7775980
    3    0.4025468  -0.0122602   0.5402023
    4   -0.0122602   0.4025468   0.0001485   0.5402023

          TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS
       ATOM         MULL.POP.    CHARGE          LOW.POP.     CHARGE
    1 CARBON        6.069363   -0.069363         6.042811   -0.042811
    2 CARBON        6.069363   -0.069363         6.042811   -0.042811
    3 HYDROGEN      0.930637    0.069363         0.957189    0.042811
    4 HYDROGEN      0.930637    0.069363         0.957189    0.042811

          -------------------------------
          BOND ORDER AND VALENCE ANALYSIS     BOND ORDER THRESHOLD=0.050
          -------------------------------

                   BOND                       BOND                       BOND
  ATOM PAIR DIST  ORDER      ATOM PAIR DIST  ORDER      ATOM PAIR DIST  ORDER
    1   2  1.178  2.997        1   3  1.075  0.993        2   4  1.075  0.993

                       TOTAL       BONDED        FREE
      ATOM            VALENCE     VALENCE     VALENCE
    1 CARBON            3.992       3.992       0.000
    2 CARBON            3.992       3.992       0.000
    3 HYDROGEN          0.995       0.995       0.000
    4 HYDROGEN          0.995       0.995       0.000

          ---------------------
          ELECTROSTATIC MOMENTS
          ---------------------

 POINT   1           X           Y           Z (BOHR)    CHARGE
                 0.000000    0.000000    0.000000        0.00 (A.U.)
         DX          DY          DZ         /D/  (DEBYE)
     0.000000    0.000000    0.000000    0.000000
 ...... END OF PROPERTY EVALUATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     106.01%
 ......END OF NBO ANALYSIS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     105.37%
  $VIB   
          IVIB=   0 IATOM=   0 ICOORD=   0 E=      -73.6046483151
 -8.020944297E-50 1.657114120E-49-7.508778270E-06-8.020944297E-50 1.657114120E-49
  7.508778270E-06-1.095505175E-50-3.634960006E-49 1.962522075E-07-1.095505175E-50
 -3.634960006E-49-1.962522075E-07
 -2.762709391E-34 2.762709391E-34 6.772625127E-15
 ......END OF GEOMETRY SEARCH......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     105.12%

                         I/O STATISTICS:
 DATA READ TOTAL =        3.126 MB,  DATA WRITTEN TOTAL =        1.042 MB

      197292 WORDS OF    DYNAMIC MEMORY USED
 EXECUTION OF GAMESS TERMINATED NORMALLY 20:41:25 LT  27-AUG-2007
