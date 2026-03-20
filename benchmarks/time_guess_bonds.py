class TimeGuessBonds:
    def setup(self):
        import MDAnalysis as mda
        from MDAnalysis.tests.datafiles import PSF, DCD
        
        self.u = mda.Universe(PSF, DCD)

        # Provide simple vdW radii (dummy but valid)
        self.vdwradii = {atom.type: 1.5 + (i % 3)*0.1 for i, atom in enumerate(self.u.atoms)}

    def time_guess_bonds(self):
        self.u.atoms.guess_bonds(vdwradii=self.vdwradii)
