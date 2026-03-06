"""
List of topology file objects that can be edited in ParmEd. They can be indexed
via either the name of the topology file or by the order in which they were
loaded.
"""
from ..structure import Structure
from ..formats.registry import load_file
from ..exceptions import FormatNotFound
from .exceptions import DuplicateParm, ParmIndexError, ParmError

class ParmList:
    """
    List of Structure objects index-able via either file name or file index
    (based on order added)
    """
    def __init__(self):
        """ Must be instantiated by itself """
        self._parm_names = []
        self._parm_instances = []
        self.parm = None # The active parm instance
        self.current_index = 0

    def set_new_active(self, idx):
        """ Sets a new active parm """
        self.current_index = self.index(idx)
        self.parm = self[self.current_index]

    def add_parm(self, parm, name=None):
        """ Add a parm to the list """
        # Make sure this parm is not part of the list already
        if isinstance(parm, str):
            name = parm
            if name in self._parm_names:
                raise DuplicateParm(f'{name} already in ParmList')
            try:
                parm = load_file(name, structure=True)
            except FormatNotFound:
                raise ParmError(f'Could not determine file type of {name}')
            if not isinstance(parm, Structure):
                raise ParmError('Added parm must be Structure or a subclass')
        elif not isinstance(parm, Structure):
            raise ParmError('Added parm must be Structure or a subclass')
        else:
            name = name or str(parm)
            if name in self._parm_names:
                raise DuplicateParm(f'{parm} already in ParmList')
        # Otherwise, add in the new parm's name
        self._parm_names.append(name)
        self._parm_instances.append(parm)
        # A newly added topology file is the currently active parm
        self.current_index = len(self._parm_instances) - 1
        self.parm = parm

    def index(self, idx):
        """ See what the index of the requested parm is (can be int or str) """
        if isinstance(idx, int):
            if idx < -len(self._parm_instances) or idx >= len(self._parm_instances):
                raise ParmIndexError('index out of range for ParmList')
            return idx
        else:
            try:
                return self._parm_names.index(str(idx))
            except ValueError:
                raise ParmIndexError(f'{idx} prmtop not in ParmList')

    def __getitem__(self, idx):
        """ Allow an integer index or string index to identify a parm """
        return self._parm_instances[self.index(idx)]

    def __contains__(self, name):
        """ See if a given parm name is in the list """
        if isinstance(name, int):
            return name >= 0 and name < len(self)
        else:
            return str(name) in self._parm_names
   
    def __iter__(self):
        """ Iterate through the prmtop names """
        return iter(self._parm_instances)

    def __len__(self):
        return len(self._parm_instances)

    def empty(self):
        """ Returns True if the list is empty; False otherwise """
        return len(self._parm_instances) == 0
