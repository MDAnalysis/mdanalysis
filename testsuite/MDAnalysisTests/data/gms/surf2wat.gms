1         ******************************************************
          *            GAMESS VERSION =  6 JUN 1999            *
          *             FROM IOWA STATE UNIVERSITY             *
          * M.W.SCHMIDT, K.K.BALDRIDGE, J.A.BOATZ, S.T.ELBERT, *
          *   M.S.GORDON, J.H.JENSEN, S.KOSEKI, N.MATSUNAGA,   *
          *          K.A.NGUYEN, S.J.SU, T.L.WINDUS,           *
          *       TOGETHER WITH M.DUPUIS, J.A.MONTGOMERY       *
          *         J.COMPUT.CHEM.  14, 1347-1363(1993)        *
          *******Intel x86 (Win32,Linux,OS/2,DOS) VERSION*******
          * PC GAMESS version 7.1 (Tornado), build number 4630 *
          *   Compiled on    Friday,    21-12-2007, 17:57:38   *
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


 Intel Core2/ Linux  PC GAMESS version running under Linux.
 Running on Intel CPU:  Brand ID  0, Family  6, Model 42, Stepping  7
 CPU Brand String    :  Intel(R) Core(TM) i7-2640M CPU @ 2.80GHz        
 CPU Features        :  CMOV, MMX, SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2, HTT, MWAIT, EM64T
 Data cache size     :  L1 32 KB, L2  256 KB, L3  4096 KB
 # of cores/package  :  8
 Operating System successfully passed SSE support test.


 PARALLEL VERSION (MPICH) RUNNING SEQUENTIALLY ON SINGLE NODE


 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 EXECUTION OF GAMESS BEGUN  9:54:58 LT  20-MAR-2015

            ECHO OF THE FIRST FEW INPUT CARDS -
 INPUT CARD> $CONTRL EXETYP=RUN RUNTYP=SURFACE SCFTYP=RHF COORD=CART UNITS=ANGS $END        
 INPUT CARD> $SURF IVEC1(1)=1,4 IGRP1(1)=4,5,6 ORIG1=0 DISP1=0.1 NDISP1=10 $END             
 INPUT CARD> $BASIS GBASIS=PM3 $END                                                         
 INPUT CARD>                                                                                
 INPUT CARD> $DATA                                                                          
 INPUT CARD>2wat.pdb                                                                        
 INPUT CARD>C1                                                                              
 INPUT CARD>O      8.0     -0.0000000000   -0.0000000000    0.0000000000                    
 INPUT CARD>H      1.0     -0.5390000000   -0.6090000000    0.5110000000                    
 INPUT CARD>H      1.0      0.1250000000    0.8110000000    0.4980000000                    
 INPUT CARD>O      8.0      0.2900000000    2.8640000000    0.1790000000                    
 INPUT CARD>H      1.0     -0.5190000000    2.3480000000    0.1430000000                    
 INPUT CARD>H      1.0      0.5640000000    2.9670000000    1.0930000000                    
 INPUT CARD> $END                                                                           
 INPUT CARD>                                                                                
 INPUT CARD>                                                                                
    2000000 WORDS OF MEMORY AVAILABLE

     BASIS OPTIONS
     -------------
     GBASIS=PM3          IGAUSS=       0      POLAR=NONE    
     NDFUNC=       0     NFFUNC=       0     DIFFSP=       F
     NPFUNC=       0      DIFFS=       F


     RUN TITLE
     ---------
 2wat.pdb                                                                        

 THE POINT GROUP OF THE MOLECULE IS C1      
 THE ORDER OF THE PRINCIPAL AXIS IS     0

 THE MOMENTS OF INERTIA ARE (AMU-ANGSTROM**2)
 IXX=       2.067   IYY=      75.433   IZZ=      75.833

 ATOM      ATOMIC                      COORDINATES (BOHR)
           CHARGE         X                   Y                   Z
 O           8.0    -2.7147198593       -0.0889313894        0.0863298571
 H           1.0    -3.9079005753        0.8531456570       -0.9052343501
 H           1.0    -1.1090827079        0.7545511549        0.0605105544
 O           8.0     2.7353529766       -0.1007096046        0.0319695322
 H           1.0     1.5947702889       -0.2482258057       -1.3715792702
 H           1.0     3.0947505344        1.6502684137        0.3388063934

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9604806 *    0.9598698 *    2.8842047 *  
  2  H               0.9604806 *    0.0000000      1.5676304 *    3.5859718    
  3  H               0.9598698 *    1.5676304 *    0.0000000      2.0841773 *  
  4  O               2.8842047 *    3.5859718      2.0841773 *    0.0000000    
  5  H               2.4089238 *    2.9798780 *    1.7038574 *    0.9602255 *  
  6  H               3.2118272      3.7872297      2.2792723 *    0.9597296 *  

                    H              H         

  1  O               2.4089238 *    3.2118272    
  2  H               2.9798780 *    3.7872297    
  3  H               1.7038574 *    2.2792723 *  
  4  O               0.9602255 *    0.9597296 *  
  5  H               0.0000000      1.5679764 *  
  6  H               1.5679764 *    0.0000000    

  * ... LESS THAN  3.000


 TOTAL NUMBER OF SHELLS              =    6
 TOTAL NUMBER OF BASIS FUNCTIONS     =   12
 NUMBER OF ELECTRONS                 =   16
 CHARGE OF MOLECULE                  =    0
 STATE MULTIPLICITY                  =    1
 NUMBER OF OCCUPIED ORBITALS (ALPHA) =    8
 NUMBER OF OCCUPIED ORBITALS (BETA ) =    8
 TOTAL NUMBER OF ATOMS               =    6
 THE NUCLEAR REPULSION ENERGY IS       37.2025308153

 THE PARAMETERS USED IN THIS CALCULATION ARE DESCRIBED IN:

  H: (PM3): J. J. P. STEWART, J. COMP. CHEM.     10, 209 (1989).                
  O: (PM3): J. J. P. STEWART, J. COMP. CHEM.     10, 209 (1989).                

 THERE ARE    2 HEAVY AND    4 LIGHT ATOMS,
 YIELDING A TOTAL OF       186 MOPAC 2E- INTEGRALS.

     $CONTRL OPTIONS
     ---------------
     SCFTYP=RHF          RUNTYP=SURFACE      EXETYP=RUN     
     MPLEVL=       0     LOCAL =NONE         UNITS =ANGS    
     MULT  =       1     ICHARG=       0     MAXIT =      30
     NPRINT=       7     IREST =       0     COORD =CART    
     ECP   =NONE         NORMF =       0     NORMP =       0
     ITOL  =      20     ICUT  =       9     NZVAR =       0
     NOSYM =       0     INTTYP=POPLE        GEOM  =INPUT   
     PLTORB=       F     MOLPLT=       F     RPAC  =       F
     AIMPAC=       F     FRIEND=             CITYP =NONE    
     DFTTYP=NONE    

     $SYSTEM OPTIONS
     ---------------
     KDIAG =       0     MEMORY=  2000000     TIMLIM=    36000.0 SEC.
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

     ------------------------------------------
     THE POINT GROUP IS C1 , NAXIS= 0, ORDER= 1
     ------------------------------------------

     DIMENSIONS OF THE SYMMETRY SUBSPACES ARE
 A   =  12

 ..... DONE SETTING UP THE RUN .....

 CPU        TIME:   STEP =      0.06 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.06 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =     95.86%,  TOTAL =     100.00%

     ---------------------------
     POTENTIAL SURFACE MAP INPUT
     ---------------------------
 COORD 1 LYING ALONG ATOM PAIR    1    4
 HAS ORIGIN= 0.000, DISPLACEMENT= 0.100  AND 10 STEPS.
 GROUP 1 CONTAINS    3 ATOMS:
     4    5    6

 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          -------------
          GUESS OPTIONS
          -------------
          GUESS =HUCKEL            NORB  =       0          NORDER=       0
          MIX   =       F          PRTMO =       F          SYMDEN=       F
          TOLZ  = 0.0E+00          TOLE  = 0.0E+00

 INITIAL GUESS ORBITALS GENERATED BY HUCKEL   ROUTINE.

 SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW.   BOTH SET(S).
     8 ORBITALS ARE OCCUPIED (    0 CORE ORBITALS).
     1=A        2=A        3=A        4=A        5=A        6=A        7=A   
     8=A        9=A       10=A       11=A       12=A   
 ...... END OF INITIAL ORBITAL SELECTION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.694356453   -23.694356453   0.186490456   0.000000000
   2  1  0   -23.861520279    -0.167163826   0.070924486   0.000000000
   3  2  0   -23.874677468    -0.013157189   0.026633612   0.000000000
   4  3  0   -23.876150396    -0.001472928   0.010257081   0.000000000
   5  0  0   -23.876359478    -0.000209083   0.006749072   0.000000000
   6  1  0   -23.876400523    -0.000041045   0.000281158   0.000000000
   7  2  0   -23.876400692    -0.000000169   0.000173226   0.000000000
   8  3  0   -23.876400755    -0.000000063   0.000106665   0.000000000
   9  4  0   -23.876400779    -0.000000024   0.000065743   0.000000000
  10  5  0   -23.876400788    -0.000000009   0.000040570   0.000000000
  11  6  0   -23.876400792    -0.000000004   0.000025066   0.000000000
  12  7  0   -23.876400793    -0.000000001   0.000015503   0.000000000
  13  8  0   -23.876400793    -0.000000001   0.000009597   0.000000000
  14  9  0   -23.876400794     0.000000000   0.000005946   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8764007937 AFTER  14 ITERATIONS

 HEAT OF FORMATION IS     -104.60522 KCAL/MOL
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ---- SURFACE MAPPING GEOMETRY ----
 COORD 1=   0.00000  COORD 2=   0.00000
 HAS ENERGY VALUE        -23.8764007937
 O     -1.43657  -0.04706   0.04568
 H     -2.06797   0.45147  -0.47903
 H     -0.58690   0.39929   0.03202
 O      1.44749  -0.05329   0.01692
 H      0.84392  -0.13136  -0.72581
 H      1.63767   0.87328   0.17929
 ----------------------------------

 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.877314219   -23.877314219   0.008375251   0.000000000
   2  1  0   -23.877503208    -0.000188990   0.002860388   0.000000000
   3  2  0   -23.877517479    -0.000014271   0.000963002   0.000000000
   4  3  0   -23.877519041    -0.000001562   0.000326006   0.000000000
   5  0  0   -23.877519256    -0.000000215   0.000189167   0.000000000
   6  1  0   -23.877519305    -0.000000049   0.000012747   0.000000000
   7  2  0   -23.877519305     0.000000000   0.000007128   0.000000000
   8  3  0   -23.877519305     0.000000000   0.000004098   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8775193054 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -105.30712 KCAL/MOL
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ---- SURFACE MAPPING GEOMETRY ----
 COORD 1=   0.10000  COORD 2=   0.00000
 HAS ENERGY VALUE        -23.8775193054
 O     -1.43657  -0.04706   0.04568
 H     -2.06797   0.45147  -0.47903
 H     -0.58690   0.39929   0.03202
 O      1.54748  -0.05351   0.01592
 H      0.94391  -0.13157  -0.72681
 H      1.73767   0.87307   0.17829
 ----------------------------------

 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.878001267   -23.878001267   0.007486083   0.000000000
   2  1  0   -23.878151862    -0.000150595   0.002460939   0.000000000
   3  2  0   -23.878162312    -0.000010450   0.000796017   0.000000000
   4  3  0   -23.878163337    -0.000001025   0.000258401   0.000000000
   5  0  0   -23.878163463    -0.000000126   0.000133270   0.000000000
   6  1  0   -23.878163489    -0.000000026   0.000009299   0.000000000
   7  2  0   -23.878163489     0.000000000   0.000005239   0.000000000
   8  3  0   -23.878163489     0.000000000   0.000003030   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8781634891 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -105.71136 KCAL/MOL
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ---- SURFACE MAPPING GEOMETRY ----
 COORD 1=   0.20000  COORD 2=   0.00000
 HAS ENERGY VALUE        -23.8781634891
 O     -1.43657  -0.04706   0.04568
 H     -2.06797   0.45147  -0.47903
 H     -0.58690   0.39929   0.03202
 O      1.64748  -0.05373   0.01492
 H      1.04391  -0.13179  -0.72780
 H      1.83766   0.87285   0.17729
 ----------------------------------

 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.878369863   -23.878369863   0.006645530   0.000000000
   2  1  0   -23.878489560    -0.000119697   0.002105969   0.000000000
   3  2  0   -23.878497217    -0.000007657   0.000655579   0.000000000
   4  3  0   -23.878497893    -0.000000675   0.000204432   0.000000000
   5  0  0   -23.878497967    -0.000000074   0.000098989   0.000000000
   6  1  0   -23.878497980    -0.000000013   0.000006577   0.000000000
   7  2  0   -23.878497980     0.000000000   0.000003735   0.000000000
   8  3  0   -23.878497980     0.000000000   0.000002174   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8784979803 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -105.92126 KCAL/MOL
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ---- SURFACE MAPPING GEOMETRY ----
 COORD 1=   0.30000  COORD 2=   0.00000
 HAS ENERGY VALUE        -23.8784979803
 O     -1.43657  -0.04706   0.04568
 H     -2.06797   0.45147  -0.47903
 H     -0.58690   0.39929   0.03202
 O      1.74747  -0.05394   0.01393
 H      1.14390  -0.13200  -0.72880
 H      1.93766   0.87264   0.17630
 ----------------------------------

 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.878582830   -23.878582830   0.005948268   0.000000000
   2  1  0   -23.878677703    -0.000094873   0.001793754   0.000000000
   3  2  0   -23.878683317    -0.000005614   0.000538166   0.000000000
   4  3  0   -23.878683765    -0.000000448   0.000161465   0.000000000
   5  4  0   -23.878683808    -0.000000044   0.000048925   0.000000000
   6  5  0   -23.878683814    -0.000000006   0.000016073   0.000000000
   7  6  0   -23.878683815    -0.000000001   0.000009209   0.000000000
   8  7  0   -23.878683815     0.000000000   0.000005411   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8786838155 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -106.03788 KCAL/MOL
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.01 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =    100.00%,  TOTAL =     100.00%
 ---- SURFACE MAPPING GEOMETRY ----
 COORD 1=   0.40000  COORD 2=   0.00000
 HAS ENERGY VALUE        -23.8786838155
 O     -1.43657  -0.04706   0.04568
 H     -2.06797   0.45147  -0.47903
 H     -0.58690   0.39929   0.03202
 O      1.84747  -0.05416   0.01293
 H      1.24390  -0.13222  -0.72980
 H      2.03765   0.87242   0.17530
 ----------------------------------

 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.878741417   -23.878741417   0.005527910   0.000000000
   2  1  0   -23.878816405    -0.000074988   0.001521838   0.000000000
   3  2  0   -23.878820523    -0.000004118   0.000440622   0.000000000
   4  3  0   -23.878820822    -0.000000299   0.000127380   0.000000000
   5  4  0   -23.878820848    -0.000000026   0.000037122   0.000000000
   6  5  0   -23.878820851    -0.000000003   0.000011007   0.000000000
   7  6  0   -23.878820852    -0.000000001   0.000006303   0.000000000
   8  7  0   -23.878820852     0.000000000   0.000003706   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8788208518 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -106.12387 KCAL/MOL
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ---- SURFACE MAPPING GEOMETRY ----
 COORD 1=   0.50000  COORD 2=   0.00000
 HAS ENERGY VALUE        -23.8788208518
 O     -1.43657  -0.04706   0.04568
 H     -2.06797   0.45147  -0.47903
 H     -0.58690   0.39929   0.03202
 O      1.94746  -0.05437   0.01193
 H      1.34389  -0.13244  -0.73080
 H      2.13765   0.87220   0.17430
 ----------------------------------

 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.878872031   -23.878872031   0.005100350   0.000000000
   2  1  0   -23.878931120    -0.000059089   0.001286481   0.000000000
   3  2  0   -23.878934140    -0.000003020   0.000359862   0.000000000
   4  3  0   -23.878934340    -0.000000200   0.000100365   0.000000000
   5  4  0   -23.878934356    -0.000000016   0.000028168   0.000000000
   6  5  0   -23.878934358    -0.000000002   0.000008005   0.000000000
   7  6  0   -23.878934358     0.000000000   0.000004273   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8789343580 AFTER   7 ITERATIONS

 HEAT OF FORMATION IS     -106.19510 KCAL/MOL
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ---- SURFACE MAPPING GEOMETRY ----
 COORD 1=   0.60000  COORD 2=   0.00000
 HAS ENERGY VALUE        -23.8789343580
 O     -1.43657  -0.04706   0.04568
 H     -2.06797   0.45147  -0.47903
 H     -0.58690   0.39929   0.03202
 O      2.04746  -0.05459   0.01093
 H      1.44388  -0.13265  -0.73179
 H      2.23764   0.87199   0.17330
 ----------------------------------

 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.878961563   -23.878961563   0.004676991   0.000000000
   2  1  0   -23.879008001    -0.000046438   0.001084511   0.000000000
   3  2  0   -23.879010216    -0.000002215   0.000293408   0.000000000
   4  3  0   -23.879010351    -0.000000135   0.000079043   0.000000000
   5  4  0   -23.879010360    -0.000000010   0.000021391   0.000000000
   6  5  0   -23.879010361    -0.000000001   0.000005849   0.000000000
   7  6  0   -23.879010361     0.000000000   0.000003000   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8790103613 AFTER   7 ITERATIONS

 HEAT OF FORMATION IS     -106.24279 KCAL/MOL
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ---- SURFACE MAPPING GEOMETRY ----
 COORD 1=   0.70000  COORD 2=   0.00000
 HAS ENERGY VALUE        -23.8790103613
 O     -1.43657  -0.04706   0.04568
 H     -2.06797   0.45147  -0.47903
 H     -0.58690   0.39929   0.03202
 O      2.14745  -0.05481   0.00994
 H      1.54388  -0.13287  -0.73279
 H      2.33764   0.87177   0.17231
 ----------------------------------

 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.878998279   -23.878998279   0.004265681   0.000000000
   2  1  0   -23.879034671    -0.000036393   0.000911894   0.000000000
   3  2  0   -23.879036295    -0.000001624   0.000238845   0.000000000
   4  3  0   -23.879036387    -0.000000091   0.000062220   0.000000000
   5  4  0   -23.879036393    -0.000000006   0.000016254   0.000000000
   6  5  0   -23.879036393    -0.000000001   0.000004280   0.000000000
   7  6  0   -23.879036393     0.000000000   0.000002262   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8790363931 AFTER   7 ITERATIONS

 HEAT OF FORMATION IS     -106.25913 KCAL/MOL
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ---- SURFACE MAPPING GEOMETRY ----
 COORD 1=   0.80000  COORD 2=   0.00000
 HAS ENERGY VALUE        -23.8790363931
 O     -1.43657  -0.04706   0.04568
 H     -2.06797   0.45147  -0.47903
 H     -0.58690   0.39929   0.03202
 O      2.24744  -0.05502   0.00894
 H      1.64387  -0.13308  -0.73379
 H      2.43763   0.87156   0.17131
 ----------------------------------

 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.878988307   -23.878988307   0.003872333   0.000000000
   2  1  0   -23.879016755    -0.000028448   0.000765221   0.000000000
   3  2  0   -23.879017946    -0.000001191   0.000194224   0.000000000
   4  3  0   -23.879018008    -0.000000062   0.000048977   0.000000000
   5  4  0   -23.879018012    -0.000000004   0.000012365   0.000000000
   6  5  0   -23.879018012     0.000000000   0.000003138   0.000000000
   7  6  0   -23.879018012     0.000000000   0.000001812   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8790180122 AFTER   7 ITERATIONS

 HEAT OF FORMATION IS     -106.24760 KCAL/MOL
 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ---- SURFACE MAPPING GEOMETRY ----
 COORD 1=   0.90000  COORD 2=   0.00000
 HAS ENERGY VALUE        -23.8790180122
 O     -1.43657  -0.04706   0.04568
 H     -2.06797   0.45147  -0.47903
 H     -0.58690   0.39929   0.03202
 O      2.34744  -0.05524   0.00794
 H      1.74387  -0.13330  -0.73478
 H      2.53762   0.87134   0.17031
 ----------------------------------


 OVERALL RESULTS OF THE POTENTIAL SURFACE SCAN

 ------------------------------------
   ICOORD1,     |
    COORD1      |       ENERGY
 ---------------+--------------------
  1 (   0.00000)|      -23.8764007937  
  2 (   0.10000)|      -23.8775193054  
  3 (   0.20000)|      -23.8781634891  
  4 (   0.30000)|      -23.8784979803  
  5 (   0.40000)|      -23.8786838155  
  6 (   0.50000)|      -23.8788208518  
  7 (   0.60000)|      -23.8789343580  
  8 (   0.70000)|      -23.8790103613  
  9 (   0.80000)|      -23.8790363931  
 10 (   0.90000)|      -23.8790180122  


 ENERGY DELTA MAP(S) (W.R. TO THE LOWEST FOUND)

 ------------------------------------
   ICOORD1,     |         ENERGY
    COORD1      |         DELTA
 ---------------+--------------------
  1 (   0.00000)|        0.0026355994  
  2 (   0.10000)|        0.0015170877  
  3 (   0.20000)|        0.0008729040  
  4 (   0.30000)|        0.0005384128  
  5 (   0.40000)|        0.0003525776  
  6 (   0.50000)|        0.0002155413  
  7 (   0.60000)|        0.0001020351  
  8 (   0.70000)|        0.0000260318  
  9 (   0.80000)|        0.0000000000  
 10 (   0.90000)|        0.0000183809  

 ... DONE WITH POTENTIAL SURFACE SCAN ...

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
      197280 WORDS OF    DYNAMIC MEMORY USED

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 EXECUTION OF GAMESS TERMINATED NORMALLY  9:54:59 LT  20-MAR-2015
