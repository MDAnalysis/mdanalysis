import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysisTests.datafiles import PSF, DCD


class ContactsBench(object):
    """
    Benchmarks for MDAnalysis.analysis.contacts.Contacts
    """

    # Parameter combinations tested in the benchmark.
    # radius : cutoff distance used to define a contact
    # method : algorithm used to compute contacts
    # pbc    : whether periodic boundary conditions are applied
    params = [
        [4.5, 6.0],
        ["hard_cut", "soft_cut", "radius_cut"],
        [True, False],
    ]

    # Names corresponding to the parameters above
    param_names = ["radius", "method", "pbc"]

    def setup(self, radius, method, pbc):
        """
        Prepare the Universe and contact analysis object
        for benchmarking.
        """

        # Load test trajectory
        self.u = mda.Universe(PSF, DCD)

        # Define atom selections
        self.sel1 = "protein"
        self.sel2 = "name CA"

        # Create atom groups from the selections
        g1 = self.u.select_atoms(self.sel1)
        g2 = self.u.select_atoms(self.sel2)

        # Initialize the Contacts analysis
        # select   : atom selection strings
        # refgroup : reference atom groups used for contacts
        # radius   : contact cutoff distance
        # method   : contact calculation method
        # pbc      : periodic boundary conditions flag
        self.analysis = contacts.Contacts(
            self.u,
            select=(self.sel1, self.sel2),
            refgroup=(g1, g2),
            radius=radius,
            method=method,
            pbc=pbc,
        )

    def time_contacts_run(self, radius, method, pbc):
        """
        Benchmark execution of Contacts.run()
        over the full trajectory.
        """

        self.analysis.run()
