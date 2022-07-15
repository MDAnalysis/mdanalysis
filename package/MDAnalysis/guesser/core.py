
from .base import GuesserBase
from .DefaultGuesser import DefaultGuesser

# guessers dictionary mimic registaration by metaclass
GUESSERS = {'DEFAULT': DefaultGuesser}


def get_guesser(atoms, context):
    """get an appropiate guesser to the universe and pass
       the atomGroup of the universe to the guesser

    Parameters
    ----------
    atoms: AtomGroup of the universe
    context: string or Guesser class
    Returns
    -------
    TGuesser class
    """
    if isinstance(context, GuesserBase):
        context.setAtoms(atoms)
        return context
    try:
        guesser = GUESSERS[context.upper()]
        guesser.setAtoms(atoms)
    except KeyError:
        raise TypeError("Unidentified guesser type {0}".format(context))
    return guesser
