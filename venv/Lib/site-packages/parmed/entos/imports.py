""" Imports from the Entos ecosystem """

__all__ = ['sierra', 'Molecule', 'QMMMInput', 'HAS_MISTRAL', 'HAS_ENTOS', 'constants']

try:
    import sierra
    from sierra import constants
    from sierra.inputs import Molecule
    try:
        from sierra.inputs import QMMMInput, QMMMSystem
    except ImportError:
        HAS_MISTRAL = False
        QMMMInput = QMMMSystem = None
    else:
        HAS_MISTRAL = True
except ImportError:
    sierra = Molecule = constants = QMMMInput = QMMMSystem = None
    HAS_ENTOS = HAS_MISTRAL = False
else:
    HAS_ENTOS = True
