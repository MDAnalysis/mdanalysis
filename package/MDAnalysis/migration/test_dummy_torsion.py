# Torsions now called dihedrals throughout package
# (class names)

#Torsion -> Dihedral
class NewDihedral(Torsion, NewParent):
    __metaclass__ = MyMeta

class NewDihedral(NewParent, Torsion):
    __metaclass__ = MyMeta

class NewDihedral(Torsion):
    pass

class NewDihedral(NewParent, Torsion): pass

#Improper_Torsion -> ImproperDihedral
class NewImproperdihedral(Improper_Torsion):
    pass

