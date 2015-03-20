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

 EXECUTION OF GAMESS BEGUN 10:03:09 LT  20-MAR-2015

            ECHO OF THE FIRST FEW INPUT CARDS -
 INPUT CARD> $CONTRL EXETYP=RUN RUNTYP=OPTIMIZE SCFTYP=RHF COORD=CART UNITS=ANGS $END       
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
     SCFTYP=RHF          RUNTYP=OPTIMIZE     EXETYP=RUN     
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

 CPU        TIME:   STEP =      0.05 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.05 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =    100.00%,  TOTAL =     100.00%


          -----------------------------
          STATIONARY POINT LOCATION RUN
          -----------------------------

 OBTAINING INITIAL HESSIAN, HESS=GUESS   
 DIAGONAL GUESS HESSIAN IN CARTESIAN COORDS IS H(I,I)=  0.3333

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


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.4365679876  -0.0470604680   0.0456837963
 H           1.0  -2.0679720768   0.4514652720  -0.4790294233
 H           1.0  -0.5869013365   0.3992913046   0.0320208087
 O           8.0   1.4474865638  -0.0532932315   0.0169175491
 H           1.0   0.8439161546  -0.1313554490  -0.7258085453
 H           1.0   1.6376715749   0.8732844997   0.1792886353

          ********************
          1 ELECTRON INTEGRALS
          ********************

          **************************
          1 AND 2 ELECTRON INTEGRALS
          **************************

 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
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

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          -------------------
          RHF SCF CALCULATION
          -------------------

     NUCLEAR ENERGY =        22.4188147370
     MAXIT =   30     NPUNCH=    2
     EXTRAP=T  DAMP=F  SHIFT=F  RSTRCT=F  DIIS=F  DEM=F  SOSCF=F
     DENSITY CONV=  1.00E-05
     MEMORY REQUIRED FOR RHF STEP=      1458 WORDS.

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

          ------------
          EIGENVECTORS
          ------------

                      1          2          3          4          5
                   -1.3924    -1.3208    -0.6731    -0.6398    -0.5489
                     A          A          A          A          A   
    1  O   1  S   0.565432   0.673197  -0.052746   0.024919  -0.250005
    2  O   1  X   0.034724  -0.002370   0.488753   0.509330  -0.013906
    3  O   1  Y   0.061067   0.062490   0.008203  -0.026485   0.452473
    4  O   1  Z  -0.034009  -0.037678   0.119639   0.206876  -0.335798
    5  H   2  S   0.209406   0.250961  -0.281798  -0.335622   0.284798
    6  H   3  S   0.263406   0.205599   0.332028   0.309952   0.151816
    7  O   4  S   0.646876  -0.592512  -0.081515  -0.090060  -0.251302
    8  O   4  X  -0.046092   0.009696  -0.359399   0.169480  -0.255535
    9  O   4  Y   0.060590  -0.051661  -0.287949   0.423882   0.476973
   10  O   4  Z  -0.044103   0.031087  -0.326357   0.287461  -0.247918
   11  H   5  S   0.271096  -0.188685   0.373537  -0.253475   0.202459
   12  H   6  S   0.247546  -0.214060  -0.299232   0.352643   0.239263

                      6          7          8          9         10
                   -0.5208    -0.4589    -0.4543     0.1231     0.1512
                     A          A          A          A          A   
    1  O   1  S  -0.211049  -0.025376   0.005727  -0.215014  -0.263579
    2  O   1  X   0.184088   0.084550  -0.287241  -0.113226  -0.031207
    3  O   1  Y   0.545101   0.223436   0.480056  -0.283144  -0.357496
    4  O   1  Z  -0.224595   0.073426   0.812226   0.162895   0.204533
    5  H   2  S   0.175726   0.006355  -0.000936   0.309416   0.454768
    6  H   3  S   0.252784   0.054063  -0.011186   0.327893   0.424738
    7  O   4  S   0.189526   0.005109  -0.008523  -0.279547   0.195738
    8  O   4  X   0.104162   0.768922  -0.156571   0.201972  -0.050230
    9  O   4  Y  -0.396731  -0.071847   0.034136  -0.327727   0.260377
   10  O   4  Z   0.437573  -0.580807   0.032252   0.231766  -0.147448
   11  H   5  S  -0.227588  -0.011792   0.019216   0.426408  -0.321116
   12  H   6  S  -0.178273  -0.004990   0.002533   0.411228  -0.373109

                     11         12
                    0.1851     0.1927
                     A          A   
    1  O   1  S   0.016187  -0.000894
    2  O   1  X   0.133986   0.587374
    3  O   1  Y   0.027757  -0.024654
    4  O   1  Z   0.044934   0.204007
    5  H   2  S   0.090869   0.527519
    6  H   3  S  -0.174879  -0.523010
    7  O   4  S  -0.006923  -0.006870
    8  O   4  X   0.306135  -0.104655
    9  O   4  Y   0.391217  -0.107665
   10  O   4  Z   0.369873  -0.066134
   11  H   5  S   0.549863  -0.080888
   12  H   6  S  -0.504166   0.165652

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

                    -------------
                    MOPAC CHARGES
                    -------------
         ATOM NO.   TYPE          CHARGE        ATOM  ELECTRON DENSITY
           1         O             -0.3718          6.3718
           2         H              0.1782          0.8218
           3         H              0.1841          0.8159
           4         O             -0.3580          6.3580
           5         H              0.1877          0.8123
           6         H              0.1799          0.8201

          ---------------------
          ELECTROSTATIC MOMENTS
          ---------------------

 POINT   1           X           Y           Z (ANGS)      CHARGE
                 0.000000    0.000000    0.000000        0.00 (A.U.)
         DX          DY          DZ         /D/  (DEBYE)
    -0.175042    2.781578   -1.756571    3.294444
 ...... END OF PROPERTY EVALUATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH=  0     ENERGY=     -23.8764008

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0     0.0083916    -0.0047876     0.0044140
  2  H            1.0    -0.0095706     0.0037485    -0.0059210
  3  H            1.0     0.0085875     0.0007577    -0.0009514
  4  O            8.0     0.0002623    -0.0064840     0.0013898
  5  H            1.0    -0.0100330    -0.0021466    -0.0034411
  6  H            1.0     0.0023621     0.0089119     0.0045098

          MAXIMUM GRADIENT =  0.0100330    RMS GRADIENT = 0.0057506
          FORCE CONSTANT MATRIX NOT UPDATED --- TAKING FIRST STEP
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.073187
          RADIUS OF STEP TAKEN=   0.07319  CURRENT TRUST RADIUS=   0.30000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=   1


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.4498873518  -0.0394608355   0.0386774279
 H           1.0  -2.0527791822   0.4455152562  -0.4696313038
 H           1.0  -0.6005344657   0.3980892682   0.0335304926
 O           8.0   1.4470695935  -0.0430016610   0.0147123241
 H           1.0   0.8598430316  -0.1279477011  -0.7203454728
 H           1.0   1.6339212670   0.8591376008   0.1721293527

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9257743 *    0.9554459 *    2.8970582 *  
  2  H               0.9257743 *    0.0000000      1.5376722 *    3.5668164    
  3  H               0.9554459 *    1.5376722 *    0.0000000      2.0946593 *  
  4  O               2.8970582 *    3.5668164      2.0946593 *    0.0000000    
  5  H               2.4328584 *    2.9791082 *    1.7256148 *    0.9446485 *  
  6  H               3.2148350      3.7649303      2.2857313 *    0.9346384 *  

                    H              H         

  1  O               2.4328584 *    3.2148350    
  2  H               2.9791082 *    3.7649303    
  3  H               1.7256148 *    2.2857313 *  
  4  O               0.9446485 *    0.9346384 *  
  5  H               0.0000000      1.5394953 *  
  6  H               1.5394953 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.0 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.875627379   -23.875627379   0.010862610   0.000000000
   2  1  0   -23.876011729    -0.000384349   0.003133233   0.000000000
   3  2  0   -23.876051737    -0.000040008   0.001024009   0.000000000
   4  3  0   -23.876058784    -0.000007048   0.000610062   0.000000000
   5  0  0   -23.876060454    -0.000001669   0.000823927   0.000000000
   6  1  0   -23.876061090    -0.000000636   0.000039773   0.000000000
   7  2  0   -23.876061094    -0.000000004   0.000024248   0.000000000
   8  3  0   -23.876061096    -0.000000001   0.000014888   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8760610956 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -104.39205 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH=  1     ENERGY=     -23.8760611

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0    -0.0119838     0.0157341    -0.0151997
  2  H            1.0     0.0162758    -0.0150997     0.0145734
  3  H            1.0     0.0022563    -0.0006666    -0.0013045
  4  O            8.0    -0.0042971     0.0170913    -0.0063121
  5  H            1.0     0.0014956     0.0000622     0.0091091
  6  H            1.0    -0.0037468    -0.0171213    -0.0008661

          MAXIMUM GRADIENT =  0.0171213    RMS GRADIENT = 0.0107452
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS   0.0003396981
          PREDICTED ENERGY CHANGE WAS  -0.0008927936 RATIO= -0.380
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.044064
          RADIUS OF STEP TAKEN=   0.04406  CURRENT TRUST RADIUS=   0.05000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=   2


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.4452942995  -0.0478810278   0.0468929045
 H           1.0  -2.0597187113   0.4539214269  -0.4768618084
 H           1.0  -0.6053001327   0.3982329075   0.0347461674
 O           8.0   1.4497702787  -0.0516563153   0.0182949869
 H           1.0   0.8627175934  -0.1271661217  -0.7249632548
 H           1.0   1.6354581638   0.8668810579   0.1709638250

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9506009 *    0.9511863 *    2.8952083 *  
  2  H               0.9506009 *    0.0000000      1.5427823 *    3.5801260    
  3  H               0.9511863 *    1.5427823 *    0.0000000      2.1038026 *  
  4  O               2.8952083 *    3.5801260      2.1038026 *    0.0000000    
  5  H               2.4349470 *    2.9899584 *    1.7344390 *    0.9501397 *  
  6  H               3.2160875      3.7741947      2.2932911 *    0.9494728 *  

                    H              H         

  1  O               2.4349470 *    3.2160875    
  2  H               2.9899584 *    3.7741947    
  3  H               1.7344390 *    2.2932911 *  
  4  O               0.9501397 *    0.9494728 *  
  5  H               0.0000000      1.5452971 *  
  6  H               1.5452971 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.876841253   -23.876841253   0.008794520   0.000000000
   2  1  0   -23.877033553    -0.000192300   0.002746584   0.000000000
   3  2  0   -23.877057625    -0.000024072   0.001229543   0.000000000
   4  3  0   -23.877062210    -0.000004586   0.000595769   0.000000000
   5  0  0   -23.877063275    -0.000001065   0.000581310   0.000000000
   6  1  0   -23.877063643    -0.000000368   0.000029755   0.000000000
   7  2  0   -23.877063645    -0.000000002   0.000018313   0.000000000
   8  3  0   -23.877063645    -0.000000001   0.000011247   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8770636454 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -105.02118 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH=  2     ENERGY=     -23.8770636

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0     0.0084281     0.0022000     0.0001340
  2  H            1.0    -0.0002561    -0.0004483    -0.0002146
  3  H            1.0    -0.0016725    -0.0017222    -0.0016792
  4  O            8.0    -0.0036025     0.0014160    -0.0035770
  5  H            1.0    -0.0017632     0.0003042     0.0043024
  6  H            1.0    -0.0011337    -0.0017496     0.0010344

          MAXIMUM GRADIENT =  0.0084281    RMS GRADIENT = 0.0027858
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0010025498
          PREDICTED ENERGY CHANGE WAS  -0.0009706231 RATIO=  1.033
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.041063
          RADIUS OF STEP TAKEN=   0.04106  CURRENT TRUST RADIUS=   0.08813

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=   3


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.4589787248  -0.0533815399   0.0484848975
 H           1.0  -2.0601881398   0.4566086378  -0.4778784184
 H           1.0  -0.6047789417   0.4010622658   0.0379905198
 O           8.0   1.4565032068  -0.0557288972   0.0251089246
 H           1.0   0.8675549598  -0.1272208088  -0.7330541637
 H           1.0   1.6375205322   0.8709922700   0.1684210609

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9479457 *    0.9676190 *    2.9155766 *  
  2  H               0.9479457 *    0.0000000      1.5451285 *    3.5892345    
  3  H               0.9676190 *    1.5451285 *    0.0000000      2.1113285 *  
  4  O               2.9155766 *    3.5892345      2.1113285 *    0.0000000    
  5  H               2.4554052 *    2.9962729 *    1.7439495 *    0.9626954 *  
  6  H               3.2337531      3.7765680      2.2947229 *    0.9550486 *  

                    H              H         

  1  O               2.4554052 *    3.2337531    
  2  H               2.9962729 *    3.7765680    
  3  H               1.7439495 *    2.2947229 *  
  4  O               0.9626954 *    0.9550486 *  
  5  H               0.0000000      1.5498174 *  
  6  H               1.5498174 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.876907232   -23.876907232   0.005938072   0.000000000
   2  1  0   -23.877038157    -0.000130924   0.001877192   0.000000000
   3  2  0   -23.877057982    -0.000019825   0.000997760   0.000000000
   4  3  0   -23.877062359    -0.000004378   0.000563322   0.000000000
   5  0  0   -23.877063497    -0.000001138   0.000742170   0.000000000
   6  1  0   -23.877063959    -0.000000461   0.000029328   0.000000000
   7  2  0   -23.877063961    -0.000000002   0.000018150   0.000000000
   8  3  0   -23.877063962    -0.000000001   0.000011144   0.000000000
   9  4  0   -23.877063962     0.000000000   0.000006834   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8770639618 AFTER   9 ITERATIONS

 HEAT OF FORMATION IS     -105.02137 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH=  3     ENERGY=     -23.8770640

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0    -0.0077909    -0.0049245    -0.0004489
  2  H            1.0     0.0020719    -0.0013012     0.0010985
  3  H            1.0     0.0114128     0.0064243    -0.0020890
  4  O            8.0     0.0034822    -0.0042274     0.0060462
  5  H            1.0    -0.0086076     0.0000982    -0.0057819
  6  H            1.0    -0.0005683     0.0039305     0.0011751

          MAXIMUM GRADIENT =  0.0114128    RMS GRADIENT = 0.0050523
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0000003164
          PREDICTED ENERGY CHANGE WAS  -0.0002354534 RATIO=  0.001
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.030998
          RADIUS OF STEP TAKEN=   0.03100  CURRENT TRUST RADIUS=   0.05000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=   4


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.4605673515  -0.0524761897   0.0498093667
 H           1.0  -2.0616664724   0.4592764889  -0.4792128462
 H           1.0  -0.6147235680   0.3973696613   0.0415745242
 O           8.0   1.4574755890  -0.0544744622   0.0240163055
 H           1.0   0.8781465224  -0.1270858657  -0.7327440168
 H           1.0   1.6389681730   0.8697222950   0.1656294873

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9503028 *    0.9580610 *    2.9181576 *  
  2  H               0.9503028 *    0.0000000      1.5390568 *    3.5918714    
  3  H               0.9580610 *    1.5390568 *    0.0000000      2.1209622 *  
  4  O               2.9181576 *    3.5918714      2.1209622 *    0.0000000    
  5  H               2.4672939 *    3.0084214      1.7616140 *    0.9558142 *  
  6  H               3.2358901      3.7787543      2.3059994 *    0.9524356 *  

                    H              H         

  1  O               2.4672939 *    3.2358901    
  2  H               3.0084214      3.7787543    
  3  H               1.7616140 *    2.3059994 *  
  4  O               0.9558142 *    0.9524356 *  
  5  H               0.0000000      1.5425793 *  
  6  H               1.5425793 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.877356481   -23.877356481   0.002953084   0.000000000
   2  1  0   -23.877393544    -0.000037063   0.001308427   0.000000000
   3  2  0   -23.877398813    -0.000005270   0.000785206   0.000000000
   4  3  0   -23.877400132    -0.000001319   0.000476590   0.000000000
   5  0  0   -23.877400555    -0.000000422   0.000731808   0.000000000
   6  1  0   -23.877400785    -0.000000231   0.000012438   0.000000000
   7  2  0   -23.877400786    -0.000000001   0.000006771   0.000000000
   8  3  0   -23.877400786     0.000000000   0.000004075   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8774007859 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -105.23274 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH=  4     ENERGY=     -23.8774008

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0     0.0012419    -0.0021962     0.0008557
  2  H            1.0     0.0010727     0.0000012    -0.0000635
  3  H            1.0     0.0029908     0.0025101    -0.0019400
  4  O            8.0    -0.0001632    -0.0021966     0.0014211
  5  H            1.0    -0.0038153     0.0007746    -0.0007190
  6  H            1.0    -0.0013269     0.0011069     0.0004457

          MAXIMUM GRADIENT =  0.0038153    RMS GRADIENT = 0.0017158
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0003368241
          PREDICTED ENERGY CHANGE WAS  -0.0002418355 RATIO=  1.393
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.058794
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.01521076
          TRIM/QA STEP HAS LENGTH         =   0.050000
          RADIUS OF STEP TAKEN=   0.05000  CURRENT TRUST RADIUS=   0.05000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=   5


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.4702749566  -0.0516974537   0.0508643496
 H           1.0  -2.0661540515   0.4634204943  -0.4815347716
 H           1.0  -0.6263012433   0.3917291488   0.0492904303
 O           8.0   1.4622791311  -0.0528305362   0.0246484682
 H           1.0   0.8942077328  -0.1285963143  -0.7360648222
 H           1.0   1.6438762799   0.8703065888   0.1618691664

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9507193 *    0.9533736 *    2.9326715 *  
  2  H               0.9507193 *    0.0000000      1.5362588 *    3.6017464    
  3  H               0.9533736 *    1.5362588 *    0.0000000      2.1355113 *  
  4  O               2.9326715 *    3.6017464      2.1355113 *    0.0000000    
  5  H               2.4931806 *    3.0296883      1.7887060 *    0.9524339 *  
  6  H               3.2496694      3.7873276      2.3228035 *    0.9507835 *  

                    H              H         

  1  O               2.4931806 *    3.2496694    
  2  H               3.0296883      3.7873276    
  3  H               1.7887060 *    2.3228035 *  
  4  O               0.9524339 *    0.9507835 *  
  5  H               0.0000000      1.5382118 *  
  6  H               1.5382118 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.877556327   -23.877556327   0.003510802   0.000000000
   2  1  0   -23.877609463    -0.000053136   0.002074152   0.000000000
   3  2  0   -23.877622424    -0.000012961   0.001277682   0.000000000
   4  3  0   -23.877626370    -0.000003946   0.000785529   0.000000000
   5  0  0   -23.877627671    -0.000001301   0.001232363   0.000000000
   6  1  0   -23.877628369    -0.000000698   0.000013491   0.000000000
   7  2  0   -23.877628370    -0.000000001   0.000007958   0.000000000
   8  3  0   -23.877628370     0.000000000   0.000004731   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8776283702 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -105.37556 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH=  5     ENERGY=     -23.8776284

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0     0.0049384    -0.0002989     0.0010751
  2  H            1.0     0.0009426     0.0001624    -0.0002091
  3  H            1.0    -0.0013982     0.0006273    -0.0015799
  4  O            8.0    -0.0018377    -0.0008665    -0.0007548
  5  H            1.0    -0.0008529     0.0010425     0.0015083
  6  H            1.0    -0.0017921    -0.0006668    -0.0000397

          MAXIMUM GRADIENT =  0.0049384    RMS GRADIENT = 0.0015614
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0002275843
          PREDICTED ENERGY CHANGE WAS  -0.0001629883 RATIO=  1.396
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.120971
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.04519694
          TRIM/QA STEP HAS LENGTH         =   0.050000
          RADIUS OF STEP TAKEN=   0.05000  CURRENT TRUST RADIUS=   0.05000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=   6


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.4822762011  -0.0509102128   0.0511442997
 H           1.0  -2.0721237046   0.4674849859  -0.4840603439
 H           1.0  -0.6348627974   0.3855405558   0.0575371500
 O           8.0   1.4680191303  -0.0506125896   0.0253915808
 H           1.0   0.9084612419  -0.1311833188  -0.7399450804
 H           1.0   1.6504152234   0.8720125072   0.1590052147

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9503145 *    0.9532259 *    2.9504077 *  
  2  H               0.9503145 *    0.0000000      1.5381033 *    3.6139421    
  3  H               0.9532259 *    1.5381033 *    0.0000000      2.1478769 *  
  4  O               2.9504077 *    3.6139421      2.1478769 *    0.0000000    
  5  H               2.5195023 *    3.0508634      1.8124102 *    0.9514919 *  
  6  H               3.2675948      3.7992725      2.3386848 *    0.9499253 *  

                    H              H         

  1  O               2.5195023 *    3.2675948    
  2  H               3.0508634      3.7992725    
  3  H               1.8124102 *    2.3386848 *  
  4  O               0.9514919 *    0.9499253 *  
  5  H               0.0000000      1.5378586 *  
  6  H               1.5378586 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.877755485   -23.877755485   0.003628161   0.000000000
   2  1  0   -23.877804549    -0.000049065   0.001783440   0.000000000
   3  2  0   -23.877816957    -0.000012408   0.001105528   0.000000000
   4  3  0   -23.877820620    -0.000003663   0.000681670   0.000000000
   5  0  0   -23.877821777    -0.000001157   0.001076529   0.000000000
   6  1  0   -23.877822365    -0.000000589   0.000009539   0.000000000
   7  2  0   -23.877822366     0.000000000   0.000005521   0.000000000
   8  3  0   -23.877822366     0.000000000   0.000003337   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8778223656 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -105.49729 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.01 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =    100.00%,  TOTAL =     100.00%

          NSERCH=  6     ENERGY=     -23.8778224

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0     0.0043588     0.0003039     0.0006993
  2  H            1.0     0.0007843    -0.0002474     0.0000000
  3  H            1.0    -0.0013926     0.0005539    -0.0010919
  4  O            8.0    -0.0022174     0.0000890    -0.0013921
  5  H            1.0     0.0002914     0.0008040     0.0018296
  6  H            1.0    -0.0018245    -0.0015034    -0.0000449

          MAXIMUM GRADIENT =  0.0043588    RMS GRADIENT = 0.0014949
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0001939954
          PREDICTED ENERGY CHANGE WAS  -0.0001652132 RATIO=  1.174
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.222742
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.03375438
          TRIM/QA STEP HAS LENGTH         =   0.070711
          RADIUS OF STEP TAKEN=   0.07071  CURRENT TRUST RADIUS=   0.07071

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=   7


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.4998851564  -0.0495776651   0.0509417701
 H           1.0  -2.0804570383   0.4734951694  -0.4874496457
 H           1.0  -0.6466441249   0.3754351510   0.0688380446
 O           8.0   1.4776142921  -0.0480328385   0.0273724751
 H           1.0   0.9259478272  -0.1350694446  -0.7458588352
 H           1.0   1.6610570928   0.8760815554   0.1552290119

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9489648 *    0.9534026 *    2.9775931 *  
  2  H               0.9489648 *    0.0000000      1.5410682 *    3.6327544    
  3  H               0.9534026 *    1.5410682 *    0.0000000      2.1664529 *  
  4  O               2.9775931 *    3.6327544      2.1664529 *    0.0000000    
  5  H               2.5547731 *    3.0782457      1.8432014 *    0.9538333 *  
  6  H               3.2953417      3.8175961      2.3629632 *    0.9507818 *  

                    H              H         

  1  O               2.5547731 *    3.2953417    
  2  H               3.0782457      3.8175961    
  3  H               1.8432014 *    2.3629632 *  
  4  O               0.9538333 *    0.9507818 *  
  5  H               0.0000000      1.5410293 *  
  6  H               1.5410293 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.01 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =    100.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.877923485   -23.877923485   0.005529385   0.000000000
   2  1  0   -23.878019789    -0.000096303   0.002669136   0.000000000
   3  2  0   -23.878043247    -0.000023459   0.001318467   0.000000000
   4  3  0   -23.878049855    -0.000006608   0.000784755   0.000000000
   5  0  0   -23.878051847    -0.000001991   0.001253209   0.000000000
   6  1  0   -23.878052802    -0.000000955   0.000011471   0.000000000
   7  2  0   -23.878052803    -0.000000001   0.000007259   0.000000000
   8  3  0   -23.878052803     0.000000000   0.000004468   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8780528030 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -105.64190 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH=  7     ENERGY=     -23.8780528

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0     0.0028491     0.0015379    -0.0003556
  2  H            1.0     0.0009332    -0.0013009     0.0007314
  3  H            1.0    -0.0007940     0.0004502    -0.0004503
  4  O            8.0    -0.0010393    -0.0003232     0.0003747
  5  H            1.0    -0.0005894     0.0001176    -0.0006500
  6  H            1.0    -0.0013596    -0.0004815     0.0003498

          MAXIMUM GRADIENT =  0.0028491    RMS GRADIENT = 0.0010289
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0002304374
          PREDICTED ENERGY CHANGE WAS  -0.0002205012 RATIO=  1.045
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.266640
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.00953245
          TRIM/QA STEP HAS LENGTH         =   0.141421
          RADIUS OF STEP TAKEN=   0.14142  CURRENT TRUST RADIUS=   0.14142

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=   8


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.5371062483  -0.0519435861   0.0526896616
 H           1.0  -2.0971577006   0.4885589276  -0.4958290689
 H           1.0  -0.6685940736   0.3562991091   0.0899754707
 O           8.0   1.4959785365  -0.0419141374   0.0280895008
 H           1.0   0.9620530473  -0.1412359650  -0.7522698674
 H           1.0   1.6824593312   0.8825675793   0.1464171239

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9521940 *    0.9603987 *    3.0332011    
  2  H               0.9521940 *    0.0000000      1.5496624 *    3.6696758    
  3  H               0.9603987 *    1.5496624 *    0.0000000      2.2017671 *  
  4  O               3.0332011      3.6696758      2.2017671 *    0.0000000    
  5  H               2.6271144 *    3.1338752      1.9015595 *    0.9507376 *  
  6  H               3.3537589      3.8539886      2.4098955 *    0.9504962 *  

                    H              H         

  1  O               2.6271144 *    3.3537589    
  2  H               3.1338752      3.8539886    
  3  H               1.9015595 *    2.4098955 *  
  4  O               0.9507376 *    0.9504962 *  
  5  H               0.0000000      1.5410377 *  
  6  H               1.5410377 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.877842039   -23.877842039   0.011988989   0.000000000
   2  1  0   -23.878211753    -0.000369714   0.005804729   0.000000000
   3  2  0   -23.878299996    -0.000088242   0.002866852   0.000000000
   4  3  0   -23.878324215    -0.000024220   0.001431834   0.000000000
   5  0  0   -23.878331304    -0.000007089   0.002193138   0.000000000
   6  1  0   -23.878334561    -0.000003257   0.000023661   0.000000000
   7  2  0   -23.878334564    -0.000000003   0.000014458   0.000000000
   8  3  0   -23.878334564    -0.000000001   0.000008810   0.000000000
   9  4  0   -23.878334565     0.000000000   0.000005375   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8783345645 AFTER   9 ITERATIONS

 HEAT OF FORMATION IS     -105.81871 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH=  8     ENERGY=     -23.8783346

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0    -0.0024133    -0.0030235     0.0011638
  2  H            1.0    -0.0018382     0.0004347    -0.0013592
  3  H            1.0     0.0063436     0.0032571     0.0003903
  4  O            8.0    -0.0026916     0.0001335    -0.0021737
  5  H            1.0     0.0016078    -0.0001017     0.0013259
  6  H            1.0    -0.0010083    -0.0007000     0.0006529

          MAXIMUM GRADIENT =  0.0063436    RMS GRADIENT = 0.0022486
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0002817616
          PREDICTED ENERGY CHANGE WAS  -0.0003091101 RATIO=  0.912
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.296137
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.00030070
          TRIM/QA STEP HAS LENGTH         =   0.282843
          RADIUS OF STEP TAKEN=   0.28284  CURRENT TRUST RADIUS=   0.28284

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=   9


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.6062370332  -0.0517484086   0.0542735814
 H           1.0  -2.1246290087   0.5170317666  -0.5090905978
 H           1.0  -0.7247046699   0.3143259787   0.1301464951
 O           8.0   1.5383621154  -0.0316553715   0.0360882907
 H           1.0   1.0296429321  -0.1524127954  -0.7695439392
 H           1.0   1.7251985568   0.8967907580   0.1271989906

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9537402 *    0.9575314 *    3.1447159    
  2  H               0.9537402 *    0.0000000      1.5522570 *    3.7437657    
  3  H               0.9575314 *    1.5522570 *    0.0000000      2.2912925 *  
  4  O               3.1447159      3.7437657      2.2912925 *    0.0000000    
  5  H               2.7634529 *    3.2350307      2.0260858 *    0.9604274 *  
  6  H               3.4646079      3.9204915      2.5181937 *    0.9514311 *  

                    H              H         

  1  O               2.7634529 *    3.4646079    
  2  H               3.2350307      3.9204915    
  3  H               2.0260858 *    2.5181937 *  
  4  O               0.9604274 *    0.9514311 *  
  5  H               0.0000000      1.5455658 *  
  6  H               1.5455658 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.876939928   -23.876939928   0.022596203   0.000000000
   2  1  0   -23.878334167    -0.001394240   0.011174166   0.000000000
   3  2  0   -23.878669355    -0.000335188   0.005581393   0.000000000
   4  3  0   -23.878760726    -0.000091372   0.002804138   0.000000000
   5  0  0   -23.878787183    -0.000026457   0.004180278   0.000000000
   6  1  0   -23.878799152    -0.000011969   0.000045096   0.000000000
   7  2  0   -23.878799162    -0.000000010   0.000025744   0.000000000
   8  3  0   -23.878799164    -0.000000002   0.000016094   0.000000000
   9  4  0   -23.878799165    -0.000000001   0.000009941   0.000000000
  10  5  0   -23.878799165     0.000000000   0.000006120   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8787991651 AFTER  10 ITERATIONS

 HEAT OF FORMATION IS     -106.11026 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH=  9     ENERGY=     -23.8787992

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0    -0.0006855    -0.0021978     0.0016687
  2  H            1.0    -0.0033266     0.0011515    -0.0023922
  3  H            1.0     0.0060142     0.0013194     0.0006145
  4  O            8.0     0.0027390     0.0001172     0.0063270
  5  H            1.0    -0.0042813    -0.0010409    -0.0069189
  6  H            1.0    -0.0004598     0.0006505     0.0007009

          MAXIMUM GRADIENT =  0.0069189    RMS GRADIENT = 0.0031662
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0004646006
          PREDICTED ENERGY CHANGE WAS  -0.0003302379 RATIO=  1.407
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.921937
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.00389945
          TRIM/QA STEP HAS LENGTH         =   0.282843
          RADIUS OF STEP TAKEN=   0.28284  CURRENT TRUST RADIUS=   0.28284

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=  10


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.6769286713  -0.0516281014   0.0542885863
 H           1.0  -2.1475162113   0.5440722206  -0.5195369834
 H           1.0  -0.7846626036   0.2741034453   0.1687146539
 O           8.0   1.5782033123  -0.0262487374   0.0380666042
 H           1.0   1.0999771736  -0.1622670401  -0.7803709218
 H           1.0   1.7685598929   0.9143001406   0.1079108815

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9516235 *    0.9567304 *    3.2552713    
  2  H               0.9516235 *    0.0000000      1.5504655 *    3.8101409    
  3  H               0.9567304 *    1.5504655 *    0.0000000      2.3854592 *  
  4  O               3.2552713      3.8101409      2.3854592 *    0.0000000    
  5  H               2.9017415 *    3.3336411      2.1547737 *    0.9576227 *  
  6  H               3.5787266      3.9832664      2.6329630 *    0.9621570 *  

                    H              H         

  1  O               2.9017415 *    3.5787266    
  2  H               3.3336411      3.9832664    
  3  H               2.1547737 *    2.6329630 *  
  4  O               0.9576227 *    0.9621570 *  
  5  H               0.0000000      1.5475931 *  
  6  H               1.5475931 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.877534280   -23.877534280   0.023830255   0.000000000
   2  1  0   -23.878858479    -0.001324199   0.011667966   0.000000000
   3  2  0   -23.879174967    -0.000316489   0.005790109   0.000000000
   4  3  0   -23.879260199    -0.000085232   0.002894743   0.000000000
   5  0  0   -23.879284495    -0.000024296   0.003947319   0.000000000
   6  1  0   -23.879295198    -0.000010703   0.000044472   0.000000000
   7  2  0   -23.879295207    -0.000000009   0.000027123   0.000000000
   8  3  0   -23.879295209    -0.000000002   0.000016591   0.000000000
   9  4  0   -23.879295210    -0.000000001   0.000010181   0.000000000
  10  5  0   -23.879295210     0.000000000   0.000006263   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8792952101 AFTER  10 ITERATIONS

 HEAT OF FORMATION IS     -106.42154 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH= 10     ENERGY=     -23.8792952

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0    -0.0021494    -0.0005758     0.0001892
  2  H            1.0    -0.0021684    -0.0000702    -0.0011788
  3  H            1.0     0.0062617     0.0006740     0.0004125
  4  O            8.0     0.0002426    -0.0107040     0.0041503
  5  H            1.0    -0.0038816    -0.0003649    -0.0047375
  6  H            1.0     0.0016950     0.0110410     0.0011644

          MAXIMUM GRADIENT =  0.0110410    RMS GRADIENT = 0.0043878
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0004960450
          PREDICTED ENERGY CHANGE WAS  -0.0004632786 RATIO=  1.071
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   2.083707
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.00252586
          TRIM/QA STEP HAS LENGTH         =   0.500000
          RADIUS OF STEP TAKEN=   0.50000  CURRENT TRUST RADIUS=   0.50000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=  11


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.8002993624  -0.0551708013   0.0561467475
 H           1.0  -2.1865762503   0.5936121902  -0.5384460410
 H           1.0  -0.8948722046   0.2056293416   0.2355447042
 O           8.0   1.6517138822  -0.0067055546   0.0413984997
 H           1.0   1.2244441662  -0.1798165824  -0.7988361630
 H           1.0   1.8432226614   0.9347833342   0.0732650733

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9610774 *    0.9591656 *    3.4523849    
  2  H               0.9610774 *    0.0000000      1.5550214 *    3.9279858    
  3  H               0.9591656 *    1.5550214 *    0.0000000      2.5627875 *  
  4  O               3.4523849      3.9279858      2.5627875 *    0.0000000    
  5  H               3.1457282      3.5072860      2.3895636 *    0.9583951 *  
  6  H               3.7756529      4.0902161      2.8381620 *    0.9612973 *  

                    H              H         

  1  O               3.1457282      3.7756529    
  2  H               3.5072860      4.0902161    
  3  H               2.3895636 *    2.8381620 *  
  4  O               0.9583951 *    0.9612973 *  
  5  H               0.0000000      1.5445972 *  
  6  H               1.5445972 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.874793817   -23.874793817   0.042712620   0.000000000
   2  1  0   -23.878646533    -0.003852716   0.021189464   0.000000000
   3  2  0   -23.879592893    -0.000946360   0.010580689   0.000000000
   4  3  0   -23.879849278    -0.000256384   0.005305141   0.000000000
   5  0  0   -23.879922306    -0.000073028   0.007105376   0.000000000
   6  1  0   -23.879954283    -0.000031977   0.000074547   0.000000000
   7  2  0   -23.879954303    -0.000000020   0.000041110   0.000000000
   8  3  0   -23.879954307    -0.000000004   0.000025565   0.000000000
   9  4  0   -23.879954308    -0.000000001   0.000015811   0.000000000
  10  5  0   -23.879954309     0.000000000   0.000009766   0.000000000
  11  6  0   -23.879954309     0.000000000   0.000006034   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8799543088 AFTER  11 ITERATIONS

 HEAT OF FORMATION IS     -106.83515 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH= 11     ENERGY=     -23.8799543

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0    -0.0020786    -0.0078773     0.0051858
  2  H            1.0    -0.0056862     0.0065929    -0.0071237
  3  H            1.0     0.0080698     0.0014276     0.0014963
  4  O            8.0     0.0022288    -0.0098848     0.0060496
  5  H            1.0    -0.0039512    -0.0004861    -0.0058867
  6  H            1.0     0.0014174     0.0102277     0.0002788

          MAXIMUM GRADIENT =  0.0102277    RMS GRADIENT = 0.0057042
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0006590987
          PREDICTED ENERGY CHANGE WAS  -0.0008751211 RATIO=  0.753
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.172839
          RADIUS OF STEP TAKEN=   0.17284  CURRENT TRUST RADIUS=   0.50000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=  12


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.8376232584  -0.0553221389   0.0529408142
 H           1.0  -2.1875114864   0.6045160017  -0.5349269479
 H           1.0  -0.9485030776   0.1874224428   0.2498355691
 O           8.0   1.6785905145   0.0075786884   0.0327773609
 H           1.0   1.2646130440  -0.1839683354  -0.7933400016
 H           1.0   1.8680671564   0.9321052690   0.0617860261

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9504718 *    0.9424580 *    3.5168341    
  2  H               0.9504718 *    0.0000000      1.5247823 *    3.9528935    
  3  H               0.9424580 *    1.5247823 *    0.0000000      2.6421731 *  
  4  O               3.5168341      3.9528935      2.6421731 *    0.0000000    
  5  H               3.2181689      3.5504434      2.4746776 *    0.9436829 *  
  6  H               3.8350010      4.1123107      2.9194148 *    0.9441887 *  

                    H              H         

  1  O               3.2181689      3.8350010    
  2  H               3.5504434      4.1123107    
  3  H               2.4746776 *    2.9194148 *  
  4  O               0.9436829 *    0.9441887 *  
  5  H               0.0000000      1.5300385 *  
  6  H               1.5300385 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.879491853   -23.879491853   0.012476198   0.000000000
   2  1  0   -23.880017139    -0.000525286   0.006170477   0.000000000
   3  2  0   -23.880119245    -0.000102106   0.003159555   0.000000000
   4  3  0   -23.880143788    -0.000024543   0.001665093   0.000000000
   5  0  0   -23.880150095    -0.000006307   0.001899396   0.000000000
   6  1  0   -23.880152426    -0.000002331   0.000042596   0.000000000
   7  2  0   -23.880152431    -0.000000005   0.000026003   0.000000000
   8  3  0   -23.880152433    -0.000000002   0.000015929   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8801524330 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -106.95948 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH= 12     ENERGY=     -23.8801524

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0     0.0083156     0.0018523     0.0027503
  2  H            1.0     0.0016158    -0.0004326     0.0005708
  3  H            1.0    -0.0101098    -0.0012362    -0.0036225
  4  O            8.0    -0.0008475     0.0049478    -0.0073494
  5  H            1.0     0.0029599     0.0019297     0.0075027
  6  H            1.0    -0.0019340    -0.0070610     0.0001482

          MAXIMUM GRADIENT =  0.0101098    RMS GRADIENT = 0.0047243
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0001981242
          PREDICTED ENERGY CHANGE WAS  -0.0006056507 RATIO=  0.327
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.140082
          RADIUS OF STEP TAKEN=   0.14008  CURRENT TRUST RADIUS=   0.17284

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=  13


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.8731999158  -0.0524848480   0.0498061473
 H           1.0  -2.2032521786   0.6170359254  -0.5406252899
 H           1.0  -0.9708077291   0.1658247535   0.2724397718
 O           8.0   1.6956382008   0.0120457938   0.0371019514
 H           1.0   1.2986097966  -0.1914195769  -0.8027963865
 H           1.0   1.8906447185   0.9413298797   0.0531466267

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9517362 *    0.9547442 *    3.5694441    
  2  H               0.9517362 *    0.0000000      1.5438865 *    3.9876218    
  3  H               0.9547442 *    1.5438865 *    0.0000000      2.6812247 *  
  4  O               3.5694441      3.9876218      2.6812247 *    0.0000000    
  5  H               3.2873410      3.6035221      2.5365355 *    0.9510305 *  
  6  H               3.8928403      4.1494245      2.9727778 *    0.9496599 *  

                    H              H         

  1  O               3.2873410      3.8928403    
  2  H               3.6035221      4.1494245    
  3  H               2.5365355 *    2.9727778 *  
  4  O               0.9510305 *    0.9496599 *  
  5  H               0.0000000      1.5382669 *  
  6  H               1.5382669 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.879878718   -23.879878718   0.013489785   0.000000000
   2  1  0   -23.880272684    -0.000393966   0.006370669   0.000000000
   3  2  0   -23.880367845    -0.000095161   0.003096045   0.000000000
   4  3  0   -23.880394704    -0.000026859   0.001784572   0.000000000
   5  0  0   -23.880402847    -0.000008143   0.002814725   0.000000000
   6  1  0   -23.880406787    -0.000003940   0.000030737   0.000000000
   7  2  0   -23.880406791    -0.000000004   0.000019004   0.000000000
   8  3  0   -23.880406791    -0.000000001   0.000011655   0.000000000
   9  4  0   -23.880406792     0.000000000   0.000007138   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8804067917 AFTER   9 ITERATIONS

 HEAT OF FORMATION IS     -107.11909 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH= 13     ENERGY=     -23.8804068

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0    -0.0021755    -0.0007039     0.0002128
  2  H            1.0    -0.0011866     0.0003956    -0.0011316
  3  H            1.0     0.0028886     0.0005164     0.0006986
  4  O            8.0     0.0015703     0.0008909    -0.0008282
  5  H            1.0    -0.0004522     0.0000331     0.0004824
  6  H            1.0    -0.0006445    -0.0011322     0.0005659

          MAXIMUM GRADIENT =  0.0028886    RMS GRADIENT = 0.0011455
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0002543587
          PREDICTED ENERGY CHANGE WAS  -0.0002663625 RATIO=  0.955
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.013011
          RADIUS OF STEP TAKEN=   0.01301  CURRENT TRUST RADIUS=   0.28016

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=  14


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.8744697413  -0.0532960445   0.0496344235
 H           1.0  -2.2021903946   0.6179228095  -0.5395739775
 H           1.0  -0.9748519877   0.1652965115   0.2729047559
 O           8.0   1.6946109487   0.0117519325   0.0388575823
 H           1.0   1.3020211963  -0.1916781464  -0.8043401072
 H           1.0   1.8925128710   0.9423348650   0.0515901438

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9513685 *    0.9523363 *    3.5696897    
  2  H               0.9513685 *    0.0000000      1.5399194 *    3.9858609    
  3  H               0.9523363 *    1.5399194 *    0.0000000      2.6840988 *  
  4  O               3.5696897      3.9858609      2.6840988 *    0.0000000    
  5  H               3.2921903      3.6062520      2.5440202 *    0.9520992 *  
  6  H               3.8963371      4.1498570      2.9790183 *    0.9514788 *  

                    H              H         

  1  O               3.2921903      3.8963371    
  2  H               3.6062520      4.1498570    
  3  H               2.5440202 *    2.9790183 *  
  4  O               0.9520992 *    0.9514788 *  
  5  H               0.0000000      1.5385976 *  
  6  H               1.5385976 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.880417998   -23.880417998   0.001794460   0.000000000
   2  1  0   -23.880426632    -0.000008635   0.001124475   0.000000000
   3  2  0   -23.880429327    -0.000002695   0.000699199   0.000000000
   4  3  0   -23.880430287    -0.000000960   0.000433538   0.000000000
   5  0  0   -23.880430641    -0.000000355   0.000703541   0.000000000
   6  1  0   -23.880430856    -0.000000215   0.000004123   0.000000000
   7  2  0   -23.880430856     0.000000000   0.000001957   0.000000000
   8  3  0   -23.880430857     0.000000000   0.000001132   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8804308565 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -107.13419 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH= 14     ENERGY=     -23.8804309

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0    -0.0002649    -0.0002073     0.0006581
  2  H            1.0    -0.0005731     0.0001457    -0.0007096
  3  H            1.0     0.0003512     0.0002733    -0.0001633
  4  O            8.0     0.0017081    -0.0008844     0.0003426
  5  H            1.0    -0.0008634     0.0000162    -0.0005194
  6  H            1.0    -0.0003579     0.0006566     0.0003917

          MAXIMUM GRADIENT =  0.0017081    RMS GRADIENT = 0.0006311
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0000240648
          PREDICTED ENERGY CHANGE WAS  -0.0000182165 RATIO=  1.321
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.053541
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.00021738
          TRIM/QA STEP HAS LENGTH         =   0.050000
          RADIUS OF STEP TAKEN=   0.05000  CURRENT TRUST RADIUS=   0.05000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=  15


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.8853163180  -0.0541971474   0.0480575656
 H           1.0  -2.2036689848   0.6225780192  -0.5389337200
 H           1.0  -0.9874274081   0.1588682095   0.2791469414
 O           8.0   1.6970495054   0.0144665896   0.0407476824
 H           1.0   1.3155500244  -0.1937298425  -0.8068054120
 H           1.0   1.9014460735   0.9443460992   0.0468597632

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9507533 *    0.9513168 *    3.5830313    
  2  H               0.9507533 *    0.0000000      1.5373764 *    3.9901673    
  3  H               0.9513168 *    1.5373764 *    0.0000000      2.6989076 *  
  4  O               3.5830313      3.9901673      2.6989076 *    0.0000000    
  5  H               3.3159924      3.6225704      2.5704713 *    0.9524882 *  
  6  H               3.9162047      4.1591656      3.0027525      0.9520983 *  

                    H              H         

  1  O               3.3159924      3.9162047    
  2  H               3.6225704      4.1591656    
  3  H               2.5704713 *    3.0027525    
  4  O               0.9524882 *    0.9520983 *  
  5  H               0.0000000      1.5385822 *  
  6  H               1.5385822 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.880359361   -23.880359361   0.004248007   0.000000000
   2  1  0   -23.880425727    -0.000066367   0.002664639   0.000000000
   3  2  0   -23.880445943    -0.000020216   0.001659515   0.000000000
   4  3  0   -23.880452600    -0.000006657   0.001030421   0.000000000
   5  0  0   -23.880454894    -0.000002294   0.001677931   0.000000000
   6  1  0   -23.880456188    -0.000001294   0.000005078   0.000000000
   7  2  0   -23.880456188     0.000000000   0.000002695   0.000000000
   8  3  0   -23.880456188     0.000000000   0.000001636   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8804561884 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -107.15009 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH= 15     ENERGY=     -23.8804562

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0     0.0003108     0.0002469     0.0005872
  2  H            1.0    -0.0000483    -0.0002783    -0.0002062
  3  H            1.0    -0.0008187     0.0002496    -0.0005744
  4  O            8.0     0.0017946    -0.0014889     0.0008196
  5  H            1.0    -0.0009717     0.0000100    -0.0009250
  6  H            1.0    -0.0002667     0.0012608     0.0002989

          MAXIMUM GRADIENT =  0.0017946    RMS GRADIENT = 0.0007939
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0000253319
          PREDICTED ENERGY CHANGE WAS  -0.0000152476 RATIO=  1.661
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.173319
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.00486041
          TRIM/QA STEP HAS LENGTH         =   0.050000
          RADIUS OF STEP TAKEN=   0.05000  CURRENT TRUST RADIUS=   0.05000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=  16


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.8935489864  -0.0561621523   0.0444388684
 H           1.0  -2.2031643477   0.6275011617  -0.5364261296
 H           1.0  -0.9984381579   0.1525504308   0.2855228281
 O           8.0   1.6923528760   0.0188778015   0.0441216659
 H           1.0   1.3287967378  -0.1961500095  -0.8094626742
 H           1.0   1.9116347706   0.9457146954   0.0408782620

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9490318 *    0.9502136 *    3.5866869    
  2  H               0.9490318 *    0.0000000      1.5338003 *    3.9852870    
  3  H               0.9502136 *    1.5338003 *    0.0000000      2.7049028 *  
  4  O               3.5866869      3.9852870      2.7049028 *    0.0000000    
  5  H               3.3365036      3.6369904      2.5954975 *    0.9523740 *  
  6  H               3.9348689      4.1672667      3.0261337      0.9524294 *  

                    H              H         

  1  O               3.3365036      3.9348689    
  2  H               3.6369904      4.1672667    
  3  H               2.5954975 *    3.0261337    
  4  O               0.9523740 *    0.9524294 *  
  5  H               0.0000000      1.5383871 *  
  6  H               1.5383871 *    0.0000000    

  * ... LESS THAN  3.000


 MOPAC ONE- AND TWO-ELECTRON INTEGRALS REQUIRE       342 WORDS.
 ...... END OF ONE- AND TWO-ELECTRON INTEGRALS......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
 ...... END OF ONE-ELECTRON INTEGRALS ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

 ITER EX DEM  TOTAL ENERGY      E CHANGE  DENSITY CHANGE    DIIS ERROR
   1  0  0   -23.880293344   -23.880293344   0.006947420   0.000000000
   2  1  0   -23.880429776    -0.000136432   0.004342130   0.000000000
   3  2  0   -23.880476766    -0.000046990   0.002706033   0.000000000
   4  3  0   -23.880493720    -0.000016954   0.001683992   0.000000000
   5  0  0   -23.880499983    -0.000006263   0.002763375   0.000000000
   6  1  0   -23.880503770    -0.000003788   0.000006994   0.000000000
   7  2  0   -23.880503770     0.000000000   0.000001980   0.000000000
   8  3  0   -23.880503770     0.000000000   0.000001187   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8805037704 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -107.17995 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.01 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =    100.00%,  TOTAL =     100.00%

          NSERCH= 16     ENERGY=     -23.8805038

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0     0.0006123     0.0014951    -0.0001013
  2  H            1.0     0.0009083    -0.0015356     0.0010208
  3  H            1.0    -0.0021291     0.0002657    -0.0010996
  4  O            8.0     0.0016964    -0.0018829     0.0008235
  5  H            1.0    -0.0008828     0.0000809    -0.0008726
  6  H            1.0    -0.0002052     0.0015767     0.0002291

          MAXIMUM GRADIENT =  0.0021291    RMS GRADIENT = 0.0011520
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0000475820
          PREDICTED ENERGY CHANGE WAS  -0.0000332442 RATIO=  1.431
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.378281
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.01588619
          TRIM/QA STEP HAS LENGTH         =   0.050000
          RADIUS OF STEP TAKEN=   0.05000  CURRENT TRUST RADIUS=   0.05000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=  17


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.8944165264  -0.0601176700   0.0389862515
 H           1.0  -2.1994210807   0.6322375538  -0.5321873535
 H           1.0  -1.0027952430   0.1482190243   0.2904433151
 O           8.0   1.6761571400   0.0249253609   0.0494545816
 H           1.0   1.3375060115  -0.1984819254  -0.8120250984
 H           1.0   1.9206025910   0.9455495840   0.0344011246

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9479572 *    0.9495385 *    3.5716016    
  2  H               0.9479572 *    0.0000000      1.5306562 *    3.9657587    
  3  H               0.9495385 *    1.5306562 *    0.0000000      2.6925941 *  
  4  O               3.5716016      3.9657587      2.6925941 *    0.0000000    
  5  H               3.3449497      3.6439343      2.6101050 *    0.9522304 *  
  6  H               3.9453464      4.1705853      3.0409782      0.9526432 *  

                    H              H         

  1  O               3.3449497      3.9453464    
  2  H               3.6439343      4.1705853    
  3  H               2.6101050 *    3.0409782    
  4  O               0.9522304 *    0.9526432 *  
  5  H               0.0000000      1.5379360 *  
  6  H               1.5379360 *    0.0000000    

  * ... LESS THAN  3.000


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
   1  0  0   -23.880193323   -23.880193323   0.009850865   0.000000000
   2  1  0   -23.880439541    -0.000246218   0.006182971   0.000000000
   3  2  0   -23.880531333    -0.000091791   0.003868525   0.000000000
   4  3  0   -23.880566094    -0.000034761   0.002415683   0.000000000
   5  0  0   -23.880579350    -0.000013257   0.003989060   0.000000000
   6  1  0   -23.880587595    -0.000008244   0.000006660   0.000000000
   7  2  0   -23.880587595     0.000000000   0.000001824   0.000000000
   8  3  0   -23.880587595     0.000000000   0.000000819   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8805875949 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -107.23255 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH= 17     ENERGY=     -23.8805876

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0     0.0007305     0.0021614    -0.0004297
  2  H            1.0     0.0016568    -0.0023378     0.0018205
  3  H            1.0    -0.0030173     0.0004100    -0.0015698
  4  O            8.0     0.0015652    -0.0022072     0.0008446
  5  H            1.0    -0.0007782     0.0002144    -0.0007769
  6  H            1.0    -0.0001569     0.0017592     0.0001112

          MAXIMUM GRADIENT =  0.0030173    RMS GRADIENT = 0.0015106
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0000838246
          PREDICTED ENERGY CHANGE WAS  -0.0000637189 RATIO=  1.316
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.853368
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.02272064
          TRIM/QA STEP HAS LENGTH         =   0.070711
          RADIUS OF STEP TAKEN=   0.07071  CURRENT TRUST RADIUS=   0.07071

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=  18


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.8907176209  -0.0660437166   0.0309781279
 H           1.0  -2.1935886149   0.6383575935  -0.5265490841
 H           1.0  -1.0027080842   0.1434167902   0.2962963275
 O           8.0   1.6488260848   0.0341771091   0.0565712697
 H           1.0   1.3449935417  -0.2015832151  -0.8146225834
 H           1.0   1.9308275860   0.9440073665   0.0263987631

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9480235 *    0.9501728 *    3.5410548    
  2  H               0.9480235 *    0.0000000      1.5297835 *    3.9330922    
  3  H               0.9501728 *    1.5297835 *    0.0000000      2.6645891 *  
  4  O               3.5410548      3.9330922      2.6645891 *    0.0000000    
  5  H               3.3471239      3.6482942      2.6200895 *    0.9523003 *  
  6  H               3.9527752      4.1725271      3.0527728      0.9530091 *  

                    H              H         

  1  O               3.3471239      3.9527752    
  2  H               3.6482942      4.1725271    
  3  H               2.6200895 *    3.0527728    
  4  O               0.9523003 *    0.9530091 *  
  5  H               0.0000000      1.5371715 *  
  6  H               1.5371715 *    0.0000000    

  * ... LESS THAN  3.000


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
   1  0  0   -23.879940824   -23.879940824   0.014721774   0.000000000
   2  1  0   -23.880433463    -0.000492640   0.009195011   0.000000000
   3  2  0   -23.880619111    -0.000185647   0.005729426   0.000000000
   4  3  0   -23.880689879    -0.000070768   0.003565129   0.000000000
   5  0  0   -23.880716972    -0.000027093   0.005849813   0.000000000
   6  1  0   -23.880733863    -0.000016890   0.000012870   0.000000000
   7  2  0   -23.880733863     0.000000000   0.000003568   0.000000000
   8  3  0   -23.880733863     0.000000000   0.000001386   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8807338632 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -107.32434 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH= 18     ENERGY=     -23.8807339

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0     0.0000677     0.0017893    -0.0004122
  2  H            1.0     0.0018302    -0.0022862     0.0018003
  3  H            1.0    -0.0025289     0.0007399    -0.0015710
  4  O            8.0     0.0013835    -0.0026948     0.0011869
  5  H            1.0    -0.0007040     0.0004017    -0.0008805
  6  H            1.0    -0.0000485     0.0020500    -0.0001235

          MAXIMUM GRADIENT =  0.0026948    RMS GRADIENT = 0.0015063
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0001462683
          PREDICTED ENERGY CHANGE WAS  -0.0001317622 RATIO=  1.110
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   1.412133
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.01913620
          TRIM/QA STEP HAS LENGTH         =   0.100000
          RADIUS OF STEP TAKEN=   0.10000  CURRENT TRUST RADIUS=   0.10000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=  19


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.8830260330  -0.0740056275   0.0195955675
 H           1.0  -2.1867514612   0.6473512949  -0.5201377543
 H           1.0  -0.9977690039   0.1358950362   0.3050264949
 O           8.0   1.6091316719   0.0495213871   0.0640211359
 H           1.0   1.3522770230  -0.2065115371  -0.8161922938
 H           1.0   1.9437706956   0.9400813741   0.0167596705

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9507455 *    0.9535246 *    3.4946242    
  2  H               0.9507455 *    0.0000000      1.5349797 *    3.8868202    
  3  H               0.9535246 *    1.5349797 *    0.0000000      2.6194418 *  
  4  O               3.4946242      3.8868202      2.6194418 *    0.0000000    
  5  H               3.3441419      3.6525953      2.6262312 *    0.9519994 *  
  6  H               3.9588829      4.1755434      3.0630817      0.9525303 *  

                    H              H         

  1  O               3.3441419      3.9588829    
  2  H               3.6525953      4.1755434    
  3  H               2.6262312 *    3.0630817    
  4  O               0.9519994 *    0.9525303 *  
  5  H               0.0000000      1.5356917 *  
  6  H               1.5356917 *    0.0000000    

  * ... LESS THAN  3.000


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
   1  0  0   -23.879455647   -23.879455647   0.020663057   0.000000000
   2  1  0   -23.880393964    -0.000938317   0.012872778   0.000000000
   3  2  0   -23.880746416    -0.000352451   0.007997885   0.000000000
   4  3  0   -23.880880652    -0.000134237   0.004962859   0.000000000
   5  0  0   -23.880932000    -0.000051347   0.008097978   0.000000000
   6  1  0   -23.880963951    -0.000031951   0.000019103   0.000000000
   7  2  0   -23.880963951    -0.000000001   0.000005893   0.000000000
   8  3  0   -23.880963952     0.000000000   0.000003001   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8809639517 AFTER   8 ITERATIONS

 HEAT OF FORMATION IS     -107.46873 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH= 19     ENERGY=     -23.8809640

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0    -0.0019807    -0.0008564     0.0002128
  2  H            1.0     0.0005178    -0.0001629    -0.0000583
  3  H            1.0     0.0008741     0.0012688    -0.0003432
  4  O            8.0     0.0012002    -0.0024474     0.0011796
  5  H            1.0    -0.0004640     0.0007069    -0.0006275
  6  H            1.0    -0.0001473     0.0014909    -0.0003635

          MAXIMUM GRADIENT =  0.0024474    RMS GRADIENT = 0.0010488
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0002300884
          PREDICTED ENERGY CHANGE WAS  -0.0002260194 RATIO=  1.018
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.827629
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.00669635
          TRIM/QA STEP HAS LENGTH         =   0.200000
          RADIUS OF STEP TAKEN=   0.20000  CURRENT TRUST RADIUS=   0.20000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1NSERCH=  20


 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.8638213181  -0.0883420931  -0.0036440133
 H           1.0  -2.1709738331   0.6622474588  -0.5045723887
 H           1.0  -0.9887973432   0.1231390636   0.3190769090
 O           8.0   1.5283042472   0.0793842627   0.0795839689
 H           1.0   1.3639938752  -0.2172906293  -0.8203983359
 H           1.0   1.9689272646   0.9331938649  -0.0009733193

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9532348 *    0.9563159 *    3.3972893    
  2  H               0.9532348 *    0.0000000      1.5383684 *    3.7902014    
  3  H               0.9563159 *    1.5383684 *    0.0000000      2.5288479 *  
  4  O               3.3972893      3.7902014      2.5288479 *    0.0000000    
  5  H               3.3320424      3.6564094      2.6362706 *    0.9617599 *  
  6  H               3.9665482      4.1792111      3.0833027      0.9641727 *  

                    H              H         

  1  O               3.3320424      3.9665482    
  2  H               3.6564094      4.1792111    
  3  H               2.6362706 *    3.0833027    
  4  O               0.9617599 *    0.9641727 *  
  5  H               0.0000000      1.5365599 *  
  6  H               1.5365599 *    0.0000000    

  * ... LESS THAN  3.000


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
   1  0  0   -23.875211834   -23.875211834   0.040452909   0.000000000
   2  1  0   -23.878836056    -0.003624222   0.024881683   0.000000000
   3  2  0   -23.880185640    -0.001349584   0.015302214   0.000000000
   4  3  0   -23.880696877    -0.000511238   0.009423985   0.000000000
   5  0  0   -23.880891626    -0.000194748   0.015186532   0.000000000
   6  1  0   -23.881012294    -0.000120669   0.000044631   0.000000000
   7  2  0   -23.881012302    -0.000000007   0.000021665   0.000000000
   8  3  0   -23.881012303    -0.000000002   0.000010563   0.000000000
   9  4  0   -23.881012304     0.000000000   0.000005821   0.000000000
  10  5  0   -23.881012304     0.000000000   0.000003562   0.000000000

          -----------------
          DENSITY CONVERGED
          -----------------
     TIME TO FORM FOCK OPERATORS=       0.0 SECONDS (       0.0 SEC/ITER)
     TIME TO SOLVE SCF EQUATIONS=       0.0 SECONDS (       0.0 SEC/ITER)

 FINAL ENERGY IS      -23.8810123038 AFTER  10 ITERATIONS

 HEAT OF FORMATION IS     -107.49907 KCAL/MOL

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 ...... END OF RHF CALCULATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

          NSERCH= 20     ENERGY=     -23.8810123

                                 -----------------------
                                 GRADIENT (HARTREE/BOHR)
                                 -----------------------
        ATOM     ZNUC       DE/DX         DE/DY         DE/DZ
 --------------------------------------------------------------
  1  O            8.0    -0.0034574    -0.0034787     0.0005682
  2  H            1.0    -0.0004870     0.0019002    -0.0014989
  3  H            1.0     0.0036269     0.0018308     0.0007398
  4  O            8.0    -0.0037164    -0.0106441     0.0134582
  5  H            1.0    -0.0010457    -0.0005961    -0.0100248
  6  H            1.0     0.0050796     0.0109880    -0.0032424

          MAXIMUM GRADIENT =  0.0134582    RMS GRADIENT = 0.0058441
          HESSIAN UPDATED USING THE BFGS FORMULA
             ACTUAL ENERGY CHANGE WAS  -0.0000483522
          PREDICTED ENERGY CHANGE WAS  -0.0003311356 RATIO=  0.146
          MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP
               NR STEP HAS LENGTH         =   0.452739
          TRIM/QA LAMBDA FOR NON-TS MODES =  -0.00173260
          TRIM/QA STEP HAS LENGTH         =   0.200000
          RADIUS OF STEP TAKEN=   0.20000  CURRENT TRUST RADIUS=   0.20000

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%

1     ***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****
     UPDATED HESSIAN, GEOMETRY, AND VECTORS WILL BE PUNCHED FOR RESTART

 2wat.pdb                                                                        
 **** THE GEOMETRY SEARCH IS NOT CONVERGED! ****

 THE NEXT PREDICTED SET OF COORDINATES FOLLOWS.  THEIR
 ENERGY AND GRADIENT IS UNKNOWN.  YOU MAY PREFER TO RESTART
 WITH SOME OTHER COORDINATES THAN THESE.

 YOU SHOULD RESTART "OPTIMIZE" RUNS WITH THE COORDINATES
 WHOSE ENERGY IS LOWEST.  RESTART "SADPOINT" RUNS WITH THE
 COORDINATES WHOSE RMS GRADIENT IS SMALLEST.  THESE ARE NOT
 ALWAYS THE LAST POINT COMPUTED!

 COORDINATES OF ALL ATOMS ARE (ANGS)
   ATOM   CHARGE       X              Y              Z
 ------------------------------------------------------------
 O           8.0  -1.8144054257  -0.0904653065  -0.0104107266
 H           1.0  -2.1602796712   0.6549843798  -0.5006942985
 H           1.0  -0.9366458862   0.1352060643   0.3058787705
 O           8.0   1.4815194448   0.0927404830   0.0705822198
 H           1.0   1.3183171465  -0.2176910686  -0.8070915818
 H           1.0   1.9491272843   0.9175573756   0.0108084374

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                    O              H              H              O         

  1  O               0.0000000      0.9569233 *    0.9599106 *    3.3020062    
  2  H               0.9569233 *    0.0000000      1.5549949 *    3.7289644    
  3  H               0.9599106 *    1.5549949 *    0.0000000      2.4299571 *  
  4  O               3.3020062      3.7289644      2.4299571 *    0.0000000    
  5  H               3.2349401      3.5994552      2.5393105 *    0.9451529 *  
  6  H               3.8962467      4.1494343      3.0044678      0.9500278 *  

                    H              H         

  1  O               3.2349401      3.8962467    
  2  H               3.5994552      4.1494343    
  3  H               2.5393105 *    3.0044678    
  4  O               0.9451529 *    0.9500278 *  
  5  H               0.0000000      1.5348195 *  
  6  H               1.5348195 *    0.0000000    

  * ... LESS THAN  3.000


          NUCLEAR ENERGY    =       20.6367543336
          ELECTRONIC ENERGY =        0.0000000000
          TOTAL ENERGY      =        0.0000000000

          ------------------
          MOLECULAR ORBITALS
          ------------------

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!


                      1          2          3          4          5
                   -1.3640    -1.3369    -0.6516    -0.6363    -0.5467
                     A          A          A          A          A   
    1  O   1  S   0.462368   0.747373   0.017730   0.012785   0.013930
    2  O   1  X   0.034217   0.039128  -0.200635   0.552012  -0.025009
    3  O   1  Y   0.048219   0.075580   0.071865  -0.264334  -0.019735
    4  O   1  Z  -0.007816  -0.016133  -0.122739   0.393480   0.011576
    5  H   2  S   0.175740   0.280891   0.137743  -0.433842  -0.010172
    6  H   3  S   0.193309   0.264341  -0.159679   0.418904  -0.009316
    7  O   4  S   0.741781  -0.472581   0.006520  -0.015801   0.334995
    8  O   4  X   0.014533  -0.020825   0.301962   0.054881  -0.200505
    9  O   4  Y   0.042789  -0.027845   0.539215   0.202415  -0.395391
   10  O   4  Z  -0.076153   0.047055   0.385661   0.143034   0.700678
   11  H   5  S   0.282267  -0.167465  -0.426120  -0.144560  -0.314241
   12  H   6  S   0.276858  -0.174224   0.427133   0.145366  -0.316247

                      6          7          8          9         10
                   -0.5284    -0.4624    -0.4464     0.1336     0.1590
                     A          A          A          A          A   
    1  O   1  S  -0.331125  -0.020842   0.000638  -0.131042  -0.315095
    2  O   1  X   0.414943   0.088511  -0.401716  -0.145530  -0.196165
    3  O   1  Y   0.703324   0.079758   0.387879  -0.165507  -0.445306
    4  O   1  Z  -0.128197   0.018833   0.829351   0.012718   0.111812
    5  H   2  S   0.309775   0.011741   0.000105   0.170141   0.549994
    6  H   3  S   0.309119   0.034972  -0.000306   0.228209   0.458053
    7  O   4  S  -0.003783   0.001479  -0.001130  -0.316129   0.115967
    8  O   4  X  -0.099250   0.879684  -0.009547  -0.113783   0.076061
    9  O   4  Y   0.053244  -0.458042   0.005675  -0.251181   0.100169
   10  O   4  Z   0.048199  -0.001839  -0.014713   0.423669  -0.162087
   11  H   5  S  -0.022372  -0.001447   0.003035   0.489873  -0.197673
   12  H   6  S   0.006316   0.000256   0.000101   0.507562  -0.210987

                     11         12
                    0.1811     0.2025
                     A          A   
    1  O   1  S   0.003958  -0.020651
    2  O   1  X   0.017481  -0.500347
    3  O   1  Y  -0.000679   0.198416
    4  O   1  Z   0.012609  -0.334990
    5  H   2  S   0.010096  -0.505251
    6  H   3  S  -0.028184   0.576391
    7  O   4  S  -0.007790   0.021224
    8  O   4  X   0.246379   0.041704
    9  O   4  Y   0.471395   0.036435
   10  O   4  Z   0.350243  -0.029379
   11  H   5  S   0.557113  -0.033468
   12  H   6  S  -0.531603  -0.068065

                    -------------
                    MOPAC CHARGES
                    -------------
         ATOM NO.   TYPE          CHARGE        ATOM  ELECTRON DENSITY
           1         O             -0.3654          6.3654
           2         H              0.1736          0.8264
           3         H              0.1898          0.8102
           4         O             -0.3579          6.3579
           5         H              0.1811          0.8189
           6         H              0.1787          0.8213

          ---------------------
          ELECTROSTATIC MOMENTS
          ---------------------

 POINT   1           X           Y           Z (ANGS)      CHARGE
                 0.000000    0.000000    0.000000        0.00 (A.U.)
         DX          DY          DZ         /D/  (DEBYE)
     1.369006    2.317321   -1.747454    3.209011
 ...... END OF PROPERTY EVALUATION ......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
  $VIB   
          IVIB=   0 IATOM=   0 ICOORD=   0 E=      -23.8810123038
 -3.457428572E-03-3.478746224E-03 5.681868924E-04-4.870481486E-04 1.900225172E-03
 -1.498930161E-03 3.626880276E-03 1.830770184E-03 7.397603631E-04-3.716354712E-03
 -1.064413610E-02 1.345820512E-02-1.045665836E-03-5.961090794E-04-1.002482699E-02
  5.079616992E-03 1.098799605E-02-3.242395226E-03
  1.369005732E+00 2.317321330E+00-1.747454003E+00
 ......END OF GEOMETRY SEARCH......

 CPU        TIME:   STEP =      0.00 ,  TOTAL =        0.2 SECONDS (    0.0 MIN)
 WALL CLOCK TIME:   STEP =      0.00 ,  TOTAL =        0.1 SECONDS (    0.0 MIN)
 CPU UTILIZATION:   STEP =      0.00%,  TOTAL =     100.00%
      197280 WORDS OF    DYNAMIC MEMORY USED

 WARNING! YOU ARE USING OUTDATED VERSION OF THE PC GAMESS!
 PLEASE CHECK PC GAMESS HOMEPAGE FOR INFORMATION ON UPDATES!

 EXECUTION OF GAMESS TERMINATED NORMALLY 10:03:09 LT  20-MAR-2015
