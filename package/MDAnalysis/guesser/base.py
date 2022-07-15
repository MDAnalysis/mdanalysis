class GuesserBase:
    context = "base"
    _guess = {}

    def is_guessed(self, to_guess):
        """check that the passed atrributes in the to_guess list can be guessed by the class

        Parameters
        ----------
        to_guess: list of atrributes to be guessed then added to the universe
        Returns
        -------
        True or False
        """
        for a in to_guess:
            if a.lower() not in self._guess:
                raise ValueError('{0} guesser can not guess the {1} atrribute'
                                 .format(self.context, a))
        return True

    def guessTopologyAttribute(self, to_guess):
        """map the attribute to be guessed with the apporpiate guessing method

        Parameters
        ----------
        to_guess: an atrribute to be guessed then added to the universe

        Returns
        -------
        values: list of guessed values
        """
        values = self._guess[to_guess]()
        return values

    def setAtoms(self, atoms):
        """get the AtomGroup of the universe to use its attributes
           inside the guessing methods

        Parameters
        ----------
        atoms: AtomGroup of the universe

        Returns
        -------
        values: list of guessed values
        """
        self._atoms = atoms
