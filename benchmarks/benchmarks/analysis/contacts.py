import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysisTests.datafiles import PSF, DCD


class ContactsBench:

    params = [4.5, 6.0]
    param_names = ["radius"]

    def setup(self, radius):
        self.u = mda.Universe(PSF, DCD)
        self.sel1 = "protein"
        self.sel2 = "name CA"

        g1 = self.u.select_atoms(self.sel1)
        g2 = self.u.select_atoms(self.sel2)

        self.analysis = contacts.Contacts(
            self.u,
            select=(self.sel1, self.sel2),
            refgroup=(g1, g2),
            radius=radius,
        )

    def time_contacts_run(self, radius):
        self.analysis.run()