class GuesserBase:
    context = "base"
    _guess = {}
<<<<<<< Updated upstream
=======
    #give a rank to each atrribute based on its dependency on other attributes
    #to be guessed, so that the attribute with lesser dependcy will be guessed 
    #first
    _rank = {}
    def __init__(self, atoms):
        self._atoms = atoms
>>>>>>> Stashed changes

    def is_guessed(self, to_guess):
        """check that the passed atrributes in the to_guess
        list can be guessed by the class

        Parameters
        ----------
        to_guess: list of atrributes to be guessed then added to the universe
        Returns
        -------
        True or False
        """
<<<<<<< Updated upstream
        for a in to_guess:
            if a.lower() not in self._guess:
                raise ValueError('{0} guesser can not guess the {1} atrribute'
                                 .format(self.context, a))
        return True

    def guessTopologyAttribute(self, to_guess):
=======
        for att in to_guess:
            if att.lower() not in self._guess:
                raise ValueError('{0} guesser can not guess the {1} atrribute'
                                 .format(self.context, att))
        return True

    def guess_topologyAttr(self, guess):
>>>>>>> Stashed changes
        """map the attribute to be guessed with the apporpiate guessing method

        Parameters
        ----------
<<<<<<< Updated upstream
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
=======
        guess: an atrribute to be guessed t
        
        Returns
        -------
        guessed values
        """
        return self._guess[guess]()
    
    def rank_attributes(self, attrs):
       """give a rank to each atrribute based on 
          its dependency on other attributes to be guessed,
          so that the attribute with lesser dependcy will 
          be guessed first
          
       Parameters
       ----------
       attrs: attributes list to be sorted
       
       Returns
       -------
       sorted attributes list
       """ 
       to_rank = {a: self._rank[a] for a in attrs}
       ranked_attrs = sorted(to_rank, key= to_rank.get)
       return ranked_attrs
       
       
       
       
       
       
>>>>>>> Stashed changes
