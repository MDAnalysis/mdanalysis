"""
Contains all of the elements in the periodic table and dictionaries that
store all of the data found in the periodic table for that element, indexed
by the element's symbol. For consistency with AMBER, a fictitious element
'EP' is added to all of the arrays that is just an Extra Point, with no mass
or any other meaningful attribute. It's just a container to put an extra
charge
"""
# Data descriptions:
#
#  KNOWN_ELEMENTS: number of known elements
#  Element:        array whose indices are the atomic number corresponding to
#                  the element
#  AtomicNum:      dictionary matching chemical symbol to atomic number
#  Mass:           dictionary matching chemical symbol to atomic mass
#  Name:           dicionary matching chemical symbol to their full name
#  OriginName:     dictionary matching chemical symbol to the name from which
#                  the symbol comes
#  Phase:          Lists standard phase that the element appears in

KNOWN_ELEMENTS = 118

AtomicNum = { 'H'  :  1, 'He' :  2, 'Li' :  3, 'Be' :  4, 'B'  :  5, 'C'  :  6,
              'N'  :  7, 'O'  :  8, 'F'  :  9, 'Ne' : 10, 'Na' : 11, 'Mg' : 12,
              'Al' : 13, 'Si' : 14, 'P'  : 15, 'S'  : 16, 'Cl' : 17, 'Ar' : 18,
              'K'  : 19, 'Ca' : 20, 'Sc' : 21, 'Ti' : 22, 'V'  : 23, 'Cr' : 24,
              'Mn' : 25, 'Fe' : 26, 'Co' : 27, 'Ni' : 28, 'Cu' : 29, 'Zn' : 30,
              'Ga' : 31, 'Ge' : 32, 'As' : 33, 'Se' : 34, 'Br' : 35, 'Kr' : 36,
              'Rb' : 37, 'Sr' : 38, 'Y'  : 39, 'Zr' : 40, 'Nb' : 41, 'Mo' : 42,
              'Tc' : 43, 'Ru' : 44, 'Rh' : 45, 'Pd' : 46, 'Ag' : 47, 'Cd' : 48,
              'In' : 49, 'Sn' : 50, 'Sb' : 51, 'Te' : 52, 'I'  : 53, 'Xe' : 54,
              'Cs' : 55, 'Ba' : 56, 'La' : 57, 'Ce' : 58, 'Pr' : 59, 'Nd' : 60,
              'Pm' : 61, 'Sm' : 62, 'Eu' : 63, 'Gd' : 64, 'Tb' : 65, 'Dy' : 66,
              'Ho' : 67, 'Er' : 68, 'Tm' : 69, 'Yb' : 70, 'Lu' : 71, 'Hf' : 72,
              'Ta' : 73, 'W'  : 74, 'Re' : 75, 'Os' : 76, 'Ir' : 77, 'Pt' : 78,
              'Au' : 79, 'Hg' : 80, 'Tl' : 81, 'Pb' : 82, 'Bi' : 83, 'Po' : 84,
              'At' : 85, 'Rn' : 86, 'Fr' : 87, 'Ra' : 88, 'Ac' : 89, 'Th' : 90,
              'Pa' : 91, 'U'  : 92, 'Np' : 93, 'Pu' : 94, 'Am' : 95, 'Cm' : 96,
              'Bk' : 97, 'Cf' : 98, 'Es' : 99, 'Fm' :100, 'Md' :101, 'No' :102,
              'Lr' :103, 'Rf' :104, 'Db' :105, 'Sg' :106, 'Bh' :107, 'Hs' :108,
              'Mt' :109, 'Ds' :110, 'Rg' :111, 'Cn' :112, 'Nh' :113, 'Fl' :114,
              'Mc' :115, 'Lv' :116, 'Ts' :117, 'Og' :118, 'EP' : 0 , 'LP' :  0,
              'Lp' :  0, 'Ep' :  0, 'D'  :  1, 'T'  :  1,}

Element = [ 'EP',
            'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne','Na','Mg',
            'Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca','Sc','Ti','V' ,'Cr',
            'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
            'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
            'In','Sn','Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
            'Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
            'At','Rn','Fr','Ra','Ac','Th','Pa','U' ,'Np','Pu','Am','Cm',
            'Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs',
            'Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' ]

Mass = { 'H'  :   1.0079 , 'He' :   4.0026 , 'Li' :   6.941  ,
         'Be' :   9.0122 , 'B'  :  10.811  , 'C'  :  12.0107 ,
         'N'  :  14.0067 , 'O'  :  15.9994 , 'F'  :  18.9984 ,
         'Ne' :  20.1797 , 'Na' :  22.9898 , 'Mg' :  24.3050 ,
         'Al' :  26.9815 , 'Si' :  28.0855 , 'P'  :  30.9738 ,
         'S'  :  32.065  , 'Cl' :  35.453  , 'Ar' :  39.948  ,
         'K'  :  39.0983 , 'Ca' :  40.078  , 'Sc' :  44.9559 ,
         'Ti' :  47.867  , 'V'  :  50.9415 , 'Cr' :  51.9961 ,
         'Mn' :  54.9380 , 'Fe' :  55.845  , 'Co' :  58.9331 ,
         'Ni' :  58.6934 , 'Cu' :  63.546  , 'Zn' :  65.409  ,
         'Ga' :  69.723  , 'Ge' :  72.64   , 'As' :  74.9216 ,
         'Se' :  78.96   , 'Br' :  79.904  , 'Kr' :  83.798  ,
         'Rb' :  85.4678 , 'Sr' :  87.62   , 'Y'  :  88.9059 ,
         'Zr' :  91.224  , 'Nb' :  92.9064 , 'Mo' :  95.94   ,
         'Tc' :  98.     , 'Ru' : 101.07   , 'Rh' : 102.9055 ,
         'Pd' : 106.42   , 'Ag' : 107.8682 , 'Cd' : 112.411  ,
         'In' : 114.818  , 'Sn' : 118.710  , 'Sb' : 121.760  ,
         'Te' : 127.60   , 'I'  : 126.9045 , 'Xe' : 131.293  ,
         'Cs' : 132.9055 , 'Ba' : 137.327  , 'La' : 138.9055 ,
         'Ce' : 140.116  , 'Pr' : 140.9077 , 'Nd' : 144.242  ,
         'Pm' : 145.     , 'Sm' : 150.36   , 'Eu' : 151.964  ,
         'Gd' : 157.25   , 'Tb' : 158.9254 , 'Dy' : 162.500  ,
         'Ho' : 164.9303 , 'Er' : 167.259  , 'Tm' : 168.9342 ,
         'Yb' : 173.04   , 'Lu' : 174.967  , 'Hf' : 178.49   ,
         'Ta' : 180.9479 , 'W'  : 183.84   , 'Re' : 186.207  ,
         'Os' : 190.23   , 'Ir' : 192.217  , 'Pt' : 195.084  ,
         'Au' : 196.9666 , 'Hg' : 200.59   , 'Tl' : 204.3833 ,
         'Pb' : 207.2    , 'Bi' : 208.9804 , 'Po' : 209.     ,
         'At' : 210.     , 'Rn' : 222.     , 'Fr' : 223.     ,
         'Ra' : 226.     , 'Ac' : 227.     , 'Th' : 232.0381 ,
         'Pa' : 231.0359 , 'U'  : 238.0289 , 'Np' : 237.     ,
         'Pu' : 244.     , 'Am' : 243.     , 'Cm' : 247.     ,
         'Bk' : 247.     , 'Cf' : 251.     , 'Es' : 252.     ,
         'Fm' : 257.     , 'Md' : 258.     , 'No' : 259.     ,
         'Lr' : 262.     , 'Rf' : 261.     , 'Db' : 262.     ,
         'Sg' : 266.     , 'Bh' : 264.     , 'Hs' : 277.     ,
         'Mt' : 268.     , 'Ds' : 281.     , 'Rg' : 272.     ,
         'Cn' : 285.     , 'Nh' : 286.     , 'Fl' : 289.     ,
         'Mc' : 289.     , 'Lv' : 293.     , 'Ts' : 294.     ,
         'Og' : 294.     , 'EP' : 0.000000 , 'D'  :   2.0141 ,
         'T'  :  3.0160}

Name = { 'H'  : 'Hydrogen'     ,'He' : 'Helium'       ,'Li' : 'Lithium'      ,
         'Be' : 'Beryllium'    ,'B'  : 'Boron'        ,'C'  : 'Carbon'       ,
         'N'  : 'Nitrogen'     ,'O'  : 'Oxygen'       ,'F'  : 'Fluorine'     ,
         'Ne' : 'Neon'         ,'Na' : 'Sodium'       ,'Mg' : 'Magnesium'    ,
         'Al' : 'Aluminum'     ,'Si' : 'Silicon'      ,'P'  : 'Phosphorus'   ,
         'S'  : 'Sulfur'       ,'Cl' : 'Chlorine'     ,'Ar' : 'Argon'        ,
         'K'  : 'Potassium'    ,'Ca' : 'Calcium'      ,'Sc' : 'Scandium'     ,
         'Ti' : 'Titanium'     ,'V'  : 'Vanadium'     ,'Cr' : 'Chromium'     ,
         'Mn' : 'Manganese'    ,'Fe' : 'Iron'         ,'Co' : 'Cobalt'       ,
         'Ni' : 'Nickel'       ,'Cu' : 'Copper'       ,'Zn' : 'Zinc'         ,
         'Ga' : 'Gallium'      ,'Ge' : 'Germanium'    ,'As' : 'Arsenic'      ,
         'Se' : 'Selenium'     ,'Br' : 'Bromine'      ,'Kr' : 'Krypton'      ,
         'Rb' : 'Rubidium'     ,'Sr' : 'Strontium'    ,'Y'  : 'Yttrium'      ,
         'Zr' : 'Zirconium'    ,'Nb' : 'Niobium'      ,'Mo' : 'Molybdenum'   ,
         'Tc' : 'Technetium'   ,'Ru' : 'Ruthenium'    ,'Rh' : 'Rhodium'      ,
         'Pd' : 'Palladium'    ,'Ag' : 'Silver'       ,'Cd' : 'Cadmium'      ,
         'In' : 'Indium'       ,'Sn' : 'Tin'          ,'Sb' : 'Antimony'     ,
         'Te' : 'Tellurium'    ,'I'  : 'Iodine'       ,'Xe' : 'Xenon'        ,
         'Cs' : 'Cesium'       ,'Ba' : 'Barium'       ,'La' : 'Lanthanum'    ,
         'Ce' : 'Cerium'       ,'Pr' : 'Praseodymium' ,'Nd' : 'Neodymium'    ,
         'Pm' : 'Promethium'   ,'Sm' : 'Samarium'     ,'Eu' : 'Europium'     ,
         'Gd' : 'Gadolinium'   ,'Tb' : 'Terbium'      ,'Dy' : 'Dysprosium'   ,
         'Ho' : 'Holmium'      ,'Er' : 'Erbium'       ,'Tm' : 'Thulium'      ,
         'Yb' : 'Ytterbium'    ,'Lu' : 'Lutetium'     ,'Hf' : 'Hafnium'      ,
         'Ta' : 'Tantalum'     ,'W'  : 'Tungsten'     ,'Re' : 'Rhenium'      ,
         'Os' : 'Osmium'       ,'Ir' : 'Iridium'      ,'Pt' : 'Platinum'     ,
         'Au' : 'Gold'         ,'Hg' : 'Mercury'      ,'Tl' : 'Thallium'     ,
         'Pb' : 'Lead'         ,'Bi' : 'Bismuth'      ,'Po' : 'Polonium'     ,
         'At' : 'Astatine'     ,'Rn' : 'Radon'        ,'Fr' : 'Francium'     ,
         'Ra' : 'Radium'       ,'Ac' : 'Actinium'     ,'Th' : 'Thorium'      ,
         'Pa' : 'Proactinium'  ,'U'  : 'Uranium'      ,'Np' : 'Neptunium'    ,
         'Pu' : 'Plutonium'    ,'Am' : 'Americium'    ,'Cm' : 'Curium'       ,
         'Bk' : 'Berkelium'    ,'Cf' : 'Californium'  ,'Es' : 'Einsteinium'  ,
         'Fm' : 'Fermium'      ,'Md' : 'Mendelevium'  ,'No' : 'Nobelium'     ,
         'Lr' : 'Lawrencium'   ,'Rf' : 'Rutherfordium','Db' : 'Dubnium'      ,
         'Sg' : 'Seaborgium'   ,'Bh' : 'Bohrium'      ,'Hs' : 'Hassium'      ,
         'Mt' : 'Meitnerium'   ,'Ds' : 'Darmstadtium' ,'Rg' : 'Roentgenium'  ,
         'Cn' : 'Copernicium'  ,'Nh' : 'Nihonium'     ,'Fl' : 'Flerovium'    ,
         'Mc' : 'Moscovium'    ,'Lv' : 'Livermorium'  ,'Ts' : 'Tennessine'   ,
         'Og' : 'Oganesson'    ,'EP' : 'Extra Point'  ,'LP' : 'Extra Point'  ,
         'Ep' : 'Extra Point'  ,'Lp' : 'Extra Point'  ,'D'  : 'Deuterium'    ,
         'T'  : 'Tritium'}

OriginName = {
        'H'  : 'Hydrogen'   ,'He' : 'Helium'       ,'Li' : 'Lithium'      ,
        'Be' : 'Beryllium'  ,'B'  : 'Boron'        ,'C'  : 'Carbon'       ,
        'N'  : 'Nitrogen'   ,'O'  : 'Oxygen'       ,'F'  : 'Fluorine'     ,
        'Ne' : 'Neon'       ,'Na' : 'Natrium'      ,'Mg' : 'Magnesium'    ,
        'Al' : 'Aluminum'   ,'Si' : 'Silicon'      ,'P'  : 'Phosphorus'   ,
        'S'  : 'Sulfur'     ,'Cl' : 'Chlorine'     ,'Ar' : 'Argon'        ,
        'K'  : 'Kalium'     ,'Ca' : 'Calcium'      ,'Sc' : 'Scandium'     ,
        'Ti' : 'Titanium'   ,'V'  : 'Vanadium'     ,'Cr' : 'Chromium'     ,
        'Mn' : 'Manganese'  ,'Fe' : 'Ferrum'       ,'Co' : 'Cobalt'       ,
        'Ni' : 'Nickel'     ,'Cu' : 'Cuprum'       ,'Zn' : 'Zinc'         ,
        'Ga' : 'Gallium'    ,'Ge' : 'Germanium'    ,'As' : 'Arsenic'      ,
        'Se' : 'Selenium'   ,'Br' : 'Bromine'      ,'Kr' : 'Krypton'      ,
        'Rb' : 'Rubidium'   ,'Sr' : 'Strontium'    ,'Y'  : 'Yttrium'      ,
        'Zr' : 'Zirconium'  ,'Nb' : 'Niobium'      ,'Mo' : 'Molybdenum'   ,
        'Tc' : 'Technetium' ,'Ru' : 'Ruthenium'    ,'Rh' : 'Rhodium'      ,
        'Pd' : 'Palladium'  ,'Ag' : 'Argentum'     ,'Cd' : 'Cadmium'      ,
        'In' : 'Indium'     ,'Sn' : 'Stannum'      ,'Sb' : 'Stibium'      ,
        'Te' : 'Tellurium'  ,'I'  : 'Iodine'       ,'Xe' : 'Xenon'        ,
        'Cs' : 'Cesium'     ,'Ba' : 'Barium'       ,'La' : 'Lanthanum'    ,
        'Ce' : 'Cerium'     ,'Pr' : 'Praseodymium' ,'Nd' : 'Neodymium'    ,
        'Pm' : 'Promethium' ,'Sm' : 'Samarium'     ,'Eu' : 'Europium'     ,
        'Gd' : 'Gadolinium' ,'Tb' : 'Terbium'      ,'Dy' : 'Dysprosium'   ,
        'Ho' : 'Holmium'    ,'Er' : 'Erbium'       ,'Tm' : 'Thulium'      ,
        'Yb' : 'Ytterbium'  ,'Lu' : 'Lutetium'     ,'Hf' : 'Hafnium'      ,
        'Ta' : 'Tantalum'   ,'W'  : 'Wolfram'      ,'Re' : 'Rhenium'      ,
        'Os' : 'Osmium'     ,'Ir' : 'Iridium'      ,'Pt' : 'Platinum'     ,
        'Au' : 'Aurum'      ,'Hg' : 'Hydrargyrum'  ,'Tl' : 'Thallium'     ,
        'Pb' : 'Plumbum'    ,'Bi' : 'Bismuth'      ,'Po' : 'Polonium'     ,
        'At' : 'Astatine'   ,'Rn' : 'Radon'        ,'Fr' : 'Francium'     ,
        'Ra' : 'Radium'     ,'Ac' : 'Actinium'     ,'Th' : 'Thorium'      ,
        'Pa' : 'Proactinium','U'  : 'Uranium'      ,'Np' : 'Neptunium'    ,
        'Pu' : 'Plutonium'  ,'Am' : 'Americium'    ,'Cm' : 'Curium'       ,
        'Bk' : 'Berkelium'  ,'Cf' : 'Californium'  ,'Es' : 'Einsteinium'  ,
        'Fm' : 'Fermium'    ,'Md' : 'Mendelevium'  ,'No' : 'Nobelium'     ,
        'Lr' : 'Lawrencium' ,'Rf' : 'Rutherfordium','Db' : 'Dubnium'      ,
        'Sg' : 'Seaborgium' ,'Bh' : 'Bohrium'      ,'Hs' : 'Hassium'      ,
        'Mt' : 'Meitnerium' ,'Ds' : 'Darmstadtium' ,'Rg' : 'Roentgenium'  ,
        'Cn' : 'Copernicium','Nh' : 'Nihonium'     ,'Fl' : 'Flerovium'    ,
        'Mc' : 'Moscovium'  ,'Lv' : 'Livermorium'  ,'Ts' : 'Tennessine'   ,
        'Og' : 'Oganesson'  ,'EP' : 'Extra Point'  ,'LP' : 'Extra Point'  ,
        'Ep' : 'Extra Point','Lp' : 'Extra Point'
}

Phase = { 'H'  : 'Gas'          ,'He' : 'Gas'          ,'Li' : 'Solid'        ,
          'Be' : 'Solid'        ,'B'  : 'Solid'        ,'C'  : 'Solid'        ,
          'N'  : 'Gas'          ,'O'  : 'Gas'          ,'F'  : 'Gas'          ,
          'Ne' : 'Gas'          ,'Na' : 'Solid'        ,'Mg' : 'Solid'        ,
          'Al' : 'Solid'        ,'Si' : 'Solid'        ,'P'  : 'Solid'        ,
          'S'  : 'Solid'        ,'Cl' : 'Gas'          ,'Ar' : 'Gas'          ,
          'K'  : 'Solid'        ,'Ca' : 'Solid'        ,'Sc' : 'Solid'        ,
          'Ti' : 'Solid'        ,'V'  : 'Solid'        ,'Cr' : 'Solid'        ,
          'Mn' : 'Solid'        ,'Fe' : 'Solid'        ,'Co' : 'Solid'        ,
          'Ni' : 'Solid'        ,'Cu' : 'Solid'        ,'Zn' : 'Solid'        ,
          'Ga' : 'Solid'        ,'Ge' : 'Solid'        ,'As' : 'Solid'        ,
          'Se' : 'Solid'        ,'Br' : 'Liquid'       ,'Kr' : 'Gas'          ,
          'Rb' : 'Solid'        ,'Sr' : 'Solid'        ,'Y'  : 'Solid'        ,
          'Zr' : 'Solid'        ,'Nb' : 'Solid'        ,'Mo' : 'Solid'        ,
          'Tc' : 'Solid'        ,'Ru' : 'Solid'        ,'Rh' : 'Solid'        ,
          'Pd' : 'Solid'        ,'Ag' : 'Solid'        ,'Cd' : 'Solid'        ,
          'In' : 'Solid'        ,'Sn' : 'Solid'        ,'Sb' : 'Solid'        ,
          'Te' : 'Solid'        ,'I'  : 'Solid'        ,'Xe' : 'Gas'          ,
          'Cs' : 'Solid'        ,'Ba' : 'Solid'        ,'La' : 'Solid'        ,
          'Ce' : 'Solid'        ,'Pr' : 'Solid'        ,'Nd' : 'Solid'        ,
          'Pm' : 'Solid'        ,'Sm' : 'Solid'        ,'Eu' : 'Solid'        ,
          'Gd' : 'Solid'        ,'Tb' : 'Solid'        ,'Dy' : 'Solid'        ,
          'Ho' : 'Solid'        ,'Er' : 'Solid'        ,'Tm' : 'Solid'        ,
          'Yb' : 'Solid'        ,'Lu' : 'Solid'        ,'Hf' : 'Solid'        ,
          'Ta' : 'Solid'        ,'W'  : 'Solid'        ,'Re' : 'Solid'        ,
          'Os' : 'Solid'        ,'Ir' : 'Solid'        ,'Pt' : 'Solid'        ,
          'Au' : 'Solid'        ,'Hg' : 'Liquid'       ,'Tl' : 'Solid'        ,
          'Pb' : 'Solid'        ,'Bi' : 'Solid'        ,'Po' : 'Solid'        ,
          'At' : 'Solid'        ,'Rn' : 'Solid'        ,'Fr' : 'Solid'        ,
          'Ra' : 'Gas'          ,'Ac' : 'Solid'        ,'Th' : 'Solid'        ,
          'Pa' : 'Solid'        ,'U'  : 'Solid'        ,'Np' : 'Solid'        ,
          'Pu' : 'Solid'        ,'Am' : 'Solid'        ,'Cm' : 'Solid'        ,
          'Bk' : 'Solid'        ,'Cf' : 'Solid'        ,'Es' : 'Solid'        ,
          'Fm' : 'Solid'        ,'Md' : 'Solid'        ,'No' : 'Solid'        ,
          'Lr' : 'Solid'        ,'Rf' : 'Unknown'      ,'Db' : 'Unknown'      ,
          'Sg' : 'Unknown'      ,'Bh' : 'Unknown'      ,'Hs' : 'Unknown'      ,
          'Mt' : 'Unknown'      ,'Ds' : 'Unknown'      ,'Rg' : 'Unknown'      ,
          'Cn' : 'Unknown'      ,'Nh' : 'Unknown'      ,'Fl' : 'Unknown'      ,
          'Mc' : 'Unknown'      ,'Lv' : 'Unknown'      ,'Ts' : 'Unknown'      ,
          'Og' : 'Unknown'      ,'EP' : 'N/A'          ,'Ep' : 'N/A'          ,
          'LP' : 'N/A'          ,'Lp' : 'N/A'
}

_sorted_masses = sorted([(k, v) for k, v in Mass.items() if k not in ("D", "T")], key=lambda x: x[1])

def element_by_mass(mass):
    """
    Determine the element that has a mass closest to the input mass

    Parameters
    ----------
    mass : float
        The atomic mass to compare to

    Returns
    -------
    element: str
        The returned string is the name of the element whose atomic mass is
        closest to the input mass

    Notes
    -----
    This actually fails (i.e., produces poor predictions) for some cases when
    masses have changed -- particularly in the case of Hydrogen mass
    repartitioning.
    """
    diff = mass + 1
    best_guess = 'EP'

    for element, element_mass in _sorted_masses:
        d = abs(element_mass - mass)
        if d < diff:
            best_guess = element
            diff = d
        else:
            break

    return best_guess

def element_by_name(name):
    """
    Determine the element based on the name of an atom. This is very naive. It
    first tries to match the first letter of the element. If that doesn't work,
    it tries to match the first *two* letters. If that still doesn't work, it
    defaults to an extra point

    Parameters
    ----------
    name : str
        Name of the atom to determine an element for

    Returns
    -------
    element : str
        The name of the best-matching element

    Notes
    -----
    This should be a last-case scenario for guessing element information. For
    instance, Ca will never be matched, since calcium atoms will be tagged as
    carbon before the second letter is tried. This is usually OK for
    biomolecules, but you are better off using the mass or, even better, an
    appropriate representation of the atomic number to begin with
    """
    name = name.strip()
    if len(name) == 0:
        return Element[0]
    try:
        atomic_number = AtomicNum[name[0].upper()]
    except KeyError:
        sym = name[:2]
        try:
            sym = '%s%s' % (sym[0].upper(), sym[1].lower())
            atomic_number = AtomicNum[sym]
        except (KeyError, IndexError):
            atomic_number = 0 # give up

    return Element[atomic_number]

# Add some mass aliases here. We need to do it *after* _sorted_masses is created above, since
# the element_by_mass routine which uses the _sorted_masses assumes that the masses are all
# monotonically strictly increasing
Mass.update(dict(Ep=0.0, LP=0.0, Lp=0.0))
