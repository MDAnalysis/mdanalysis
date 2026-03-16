import pytest
from MDAnalysis.lib.transformations import rotation_matrix


def test_rotation_matrix_invalid_axis():
    """rotation_matrix should raise ValueError for invalid axis"""
    with pytest.raises(ValueError):
        rotation_matrix(90, "invalid_axis")


def test_rotation_matrix_invalid_angle():
    """rotation_matrix should raise TypeError for invalid angle"""
    with pytest.raises(TypeError):
        rotation_matrix("ninety", [1, 0, 0])