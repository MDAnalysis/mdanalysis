
exe_err = ('HOLE binary {name} not found. {name} must be on the '
           'PATH, or the path must provided with the keyword '
           'argument: {kw}')

IGNORE_RESIDUES = ["SOL", "WAT", "TIP", "HOH", "K  ", "NA ", "CL "]


#: Built-in HOLE radii (based on ``simple.rad`` from the HOLE_ distribution):
#: van der Waals radii are AMBER united atom from Weiner et al. (1984), JACS, vol 106 pp765-768.
#: *Simple* - Only use one value for each element C O H etc.
#: Added radii for K+, NA+, CL- (Pauling hydration radius from Hille 2002).
#: The data file can be written with the convenience function :func:`write_simplerad2`.
SIMPLE2_RAD = """
remark: Time-stamp: <2005-11-21 13:57:55 oliver> [OB]
remark: van der Waals radii: AMBER united atom
remark: from Weiner et al. (1984), JACS, vol 106 pp765-768
remark: Simple - Only use one value for each element C O H etc.
remark: van der Waals radii
remark: general last
VDWR C??? ??? 1.85
VDWR O??? ??? 1.65
VDWR S??? ??? 2.00
VDWR N??? ??? 1.75
VDWR H??? ??? 1.00
VDWR H?   ??? 1.00
VDWR P??? ??? 2.10
remark: ASN, GLN polar H (odd names for these atoms in xplor)
VDWR E2?  GLN 1.00
VDWR D2?  ASN 1.00
remark: amber lone pairs on sulphurs
VDWR LP?? ??? 0.00
remark: for some funny reason it wants radius for K even though
remark: it is on the exclude list
remark: Use Pauling hydration radius (Hille 2001) [OB]
VDWR K?   ??? 1.33
VDWR NA?  ??? 0.95
VDWR CL?  ??? 1.81
remark: funny hydrogens in gA structure [OB]
VDWR 1H?  ??? 1.00
VDWR 2H?  ??? 1.00
VDWR 3H?  ??? 1.00
remark: need bond rad for molqpt option
BOND C??? 0.85
BOND N??? 0.75
BOND O??? 0.7
BOND S??? 1.1
BOND H??? 0.5
BOND P??? 1.0
BOND ???? 0.85
"""

hole_input = """
! Input file for Oliver Smart's HOLE program
! written by MDAnalysis.analysis.hole.HOLE
! filename = {filename}
COORD  {coordinates}
RADIUS {radius}
SPHPDB {sphpdb}
SAMPLE {sample:f}
ENDRAD {end_radius:f}
IGNORE {ignore}
SHORTO {output_level:d}
"""

hole_lines = {
    'cpoint': 'CPOINT {} {} {}\n',
    'cvect': 'CVECT {} {} {}\n',
    'random_seed': 'RASEED {}\n',
    'dcd': '\nCHARMD {dcd}\nCHARMS {iniskip} {step}\n',
}

vmd_script_array = """\
set no_water_color {no_water_color}
set one_water_color {one_water_color}
set double_water_color {double_water_color}
array set triangles {{}}
"""

vmd_script_function = """
global vmd_frame;
trace add variable vmd_frame([molinfo top]) write drawFrame

proc drawFrame { name element op } {
            global vmd_frame triangles no_water_color one_water_color double_water_color;
            set frame $vmd_frame([molinfo top])
            draw delete all;

            draw color $no_water_color;
            foreach shape [lindex $triangles($frame) 0] {
                draw trinorm {*}$shape
            }
            draw color $one_water_color;
            foreach shape [lindex $triangles($frame) 1] {
                draw trinorm {*}$shape
            }
            draw color $double_water_color;
            foreach shape [lindex $triangles($frame) 2] {
                draw trinorm {*}$shape
            }
        }

drawFrame 0 0 0 
"""

