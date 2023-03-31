# Copyright (c) 2016-2023 The Regents of the University of Michigan
# Part of GSD, released under the BSD 2-Clause License.

"""Test the gsd.hoomd API."""

import gsd.fl
import gsd.hoomd
import numpy
import pickle
import pytest


def test_create(tmp_path):
    """Test that gsd files can be created."""
    with gsd.hoomd.open(name=tmp_path / "test_create.gsd", mode='wb') as hf:
        assert hf.file.schema == 'hoomd'
        assert hf.file.schema_version >= (1, 0)


def test_append(tmp_path, open_mode):
    """Test that gsd files can be appended to."""
    frame = gsd.hoomd.Frame()
    frame.particles.N = 10

    with gsd.hoomd.open(name=tmp_path / "test_append.gsd",
                        mode=open_mode.write) as hf:
        for i in range(5):
            frame.configuration.step = i + 1
            hf.append(frame)

    with gsd.hoomd.open(name=tmp_path / "test_append.gsd",
                        mode=open_mode.read) as hf:
        assert len(hf) == 5


def create_frame(i):
    """Helper function to create frame objects."""
    frame = gsd.hoomd.Frame()
    frame.configuration.step = i + 1
    return frame


def test_extend(tmp_path, open_mode):
    """Test that the extend method works."""
    frame = gsd.hoomd.Frame()
    frame.particles.N = 10

    with gsd.hoomd.open(name=tmp_path / "test_extend.gsd",
                        mode=open_mode.write) as hf:
        hf.extend((create_frame(i) for i in range(5)))

    with gsd.hoomd.open(name=tmp_path / "test_extend.gsd",
                        mode=open_mode.read) as hf:
        assert len(hf) == 5


def test_defaults(tmp_path, open_mode):
    """Test that the property defaults are properly set."""
    frame = gsd.hoomd.Frame()
    frame.particles.N = 2
    frame.bonds.N = 3
    frame.angles.N = 4
    frame.dihedrals.N = 5
    frame.impropers.N = 6
    frame.constraints.N = 4
    frame.pairs.N = 7

    with gsd.hoomd.open(name=tmp_path / "test_defaults.gsd",
                        mode=open_mode.write) as hf:
        hf.append(frame)

    with gsd.hoomd.open(name=tmp_path / "test_defaults.gsd",
                        mode=open_mode.read) as hf:
        s = hf[0]

        assert s.configuration.step == 0
        assert s.configuration.dimensions == 3
        numpy.testing.assert_array_equal(
            s.configuration.box,
            numpy.array([1, 1, 1, 0, 0, 0], dtype=numpy.float32))
        assert s.particles.N == 2
        assert s.particles.types == ['A']
        assert s.particles.type_shapes == [{}]
        numpy.testing.assert_array_equal(
            s.particles.typeid, numpy.array([0, 0], dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.particles.mass, numpy.array([1, 1], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.diameter, numpy.array([1, 1], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.body, numpy.array([-1, -1], dtype=numpy.int32))
        numpy.testing.assert_array_equal(
            s.particles.charge, numpy.array([0, 0], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.moment_inertia,
            numpy.array([[0, 0, 0], [0, 0, 0]], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.position,
            numpy.array([[0, 0, 0], [0, 0, 0]], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.orientation,
            numpy.array([[1, 0, 0, 0], [1, 0, 0, 0]], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.velocity,
            numpy.array([[0, 0, 0], [0, 0, 0]], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.angmom,
            numpy.array([[0, 0, 0, 0], [0, 0, 0, 0]], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.image,
            numpy.array([[0, 0, 0], [0, 0, 0]], dtype=numpy.int32))

        assert s.bonds.N == 3
        assert s.bonds.types == []
        numpy.testing.assert_array_equal(
            s.bonds.typeid, numpy.array([0, 0, 0], dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.bonds.group,
            numpy.array([[0, 0], [0, 0], [0, 0]], dtype=numpy.uint32))

        assert s.angles.N == 4
        assert s.angles.types == []
        numpy.testing.assert_array_equal(
            s.angles.typeid, numpy.array([0, 0, 0, 0], dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.angles.group,
            numpy.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
                        dtype=numpy.uint32))

        assert s.dihedrals.N == 5
        assert s.dihedrals.types == []
        numpy.testing.assert_array_equal(
            s.dihedrals.typeid, numpy.array([0, 0, 0, 0, 0],
                                            dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.dihedrals.group,
            numpy.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                         [0, 0, 0, 0]],
                        dtype=numpy.uint32))

        assert s.impropers.N == 6
        assert s.impropers.types == []
        numpy.testing.assert_array_equal(
            s.impropers.typeid,
            numpy.array([0, 0, 0, 0, 0, 0], dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.impropers.group,
            numpy.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                         [0, 0, 0, 0], [0, 0, 0, 0]],
                        dtype=numpy.uint32))

        assert s.constraints.N == 4
        numpy.testing.assert_array_equal(
            s.constraints.value, numpy.array([0, 0, 0, 0], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.constraints.group,
            numpy.array([[0, 0], [0, 0], [0, 0], [0, 0]], dtype=numpy.uint32))

        assert s.pairs.N == 7
        assert s.pairs.types == []
        numpy.testing.assert_array_equal(
            s.pairs.typeid, numpy.array([0] * 7, dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.pairs.group, numpy.array([[0, 0]] * 7, dtype=numpy.uint32))

        assert len(s.state) == 0


def make_nondefault_frame():
    """Make a frame with all non-default values."""
    frame0 = gsd.hoomd.Frame()
    frame0.configuration.step = 10000
    frame0.configuration.dimensions = 2
    frame0.configuration.box = [4, 5, 6, 1.0, 0.5, 0.25]
    frame0.particles.N = 2
    frame0.particles.types = ['A', 'B', 'C']
    frame0.particles.type_shapes = [
        {
            "type": "Sphere",
            "diameter": 2.0
        },
        {
            "type": "Sphere",
            "diameter": 3.0
        },
        {
            "type": "Sphere",
            "diameter": 4.0
        },
    ]
    frame0.particles.typeid = [1, 2]
    frame0.particles.mass = [2, 3]
    frame0.particles.diameter = [3, 4]
    frame0.particles.body = [10, 20]
    frame0.particles.charge = [0.5, 0.25]
    frame0.particles.moment_inertia = [[1, 2, 3], [3, 2, 1]]
    frame0.particles.position = [[0.1, 0.2, 0.3], [-1.0, -2.0, -3.0]]
    frame0.particles.orientation = [[1, 0.1, 0.2, 0.3], [0, -1.0, -2.0, -3.0]]
    frame0.particles.velocity = [[1.1, 2.2, 3.3], [-3.3, -2.2, -1.1]]
    frame0.particles.angmom = [[1, 1.1, 2.2, 3.3], [-1, -3.3, -2.2, -1.1]]
    frame0.particles.image = [[10, 20, 30], [5, 6, 7]]

    frame0.bonds.N = 1
    frame0.bonds.types = ['bondA', 'bondB']
    frame0.bonds.typeid = [1]
    frame0.bonds.group = [[0, 1]]

    frame0.angles.N = 1
    frame0.angles.typeid = [2]
    frame0.angles.types = ['angleA', 'angleB']
    frame0.angles.group = [[0, 1, 0]]

    frame0.dihedrals.N = 1
    frame0.dihedrals.typeid = [3]
    frame0.dihedrals.types = ['dihedralA', 'dihedralB']
    frame0.dihedrals.group = [[0, 1, 1, 0]]

    frame0.impropers.N = 1
    frame0.impropers.typeid = [4]
    frame0.impropers.types = ['improperA', 'improperB']
    frame0.impropers.group = [[1, 0, 0, 1]]

    frame0.constraints.N = 1
    frame0.constraints.value = [1.1]
    frame0.constraints.group = [[0, 1]]

    frame0.pairs.N = 1
    frame0.pairs.types = ['pairA', 'pairB']
    frame0.pairs.typeid = [1]
    frame0.pairs.group = [[0, 3]]

    frame0.log['value'] = [1, 2, 4, 10, 12, 18, 22]
    return frame0


def assert_frames_equal(s, frame0, check_position=True, check_step=True):
    """Assert that two frames are equal."""
    if check_step:
        assert s.configuration.step == frame0.configuration.step

    assert s.configuration.dimensions == frame0.configuration.dimensions
    numpy.testing.assert_array_equal(s.configuration.box,
                                     frame0.configuration.box)
    assert s.particles.N == frame0.particles.N
    assert s.particles.types == frame0.particles.types
    assert s.particles.type_shapes == frame0.particles.type_shapes
    numpy.testing.assert_array_equal(s.particles.typeid,
                                     frame0.particles.typeid)
    numpy.testing.assert_array_equal(s.particles.mass, frame0.particles.mass)
    numpy.testing.assert_array_equal(s.particles.diameter,
                                     frame0.particles.diameter)
    numpy.testing.assert_array_equal(s.particles.body, frame0.particles.body)
    numpy.testing.assert_array_equal(s.particles.charge,
                                     frame0.particles.charge)
    numpy.testing.assert_array_equal(s.particles.moment_inertia,
                                     frame0.particles.moment_inertia)
    if check_position:
        numpy.testing.assert_array_equal(s.particles.position,
                                         frame0.particles.position)
    numpy.testing.assert_array_equal(s.particles.orientation,
                                     frame0.particles.orientation)
    numpy.testing.assert_array_equal(s.particles.velocity,
                                     frame0.particles.velocity)
    numpy.testing.assert_array_equal(s.particles.angmom,
                                     frame0.particles.angmom)
    numpy.testing.assert_array_equal(s.particles.image, frame0.particles.image)

    assert s.bonds.N == frame0.bonds.N
    assert s.bonds.types == frame0.bonds.types
    numpy.testing.assert_array_equal(s.bonds.typeid, frame0.bonds.typeid)
    numpy.testing.assert_array_equal(s.bonds.group, frame0.bonds.group)

    assert s.angles.N == frame0.angles.N
    assert s.angles.types == frame0.angles.types
    numpy.testing.assert_array_equal(s.angles.typeid, frame0.angles.typeid)
    numpy.testing.assert_array_equal(s.angles.group, frame0.angles.group)

    assert s.dihedrals.N == frame0.dihedrals.N
    assert s.dihedrals.types == frame0.dihedrals.types
    numpy.testing.assert_array_equal(s.dihedrals.typeid,
                                     frame0.dihedrals.typeid)
    numpy.testing.assert_array_equal(s.dihedrals.group, frame0.dihedrals.group)

    assert s.impropers.N == frame0.impropers.N
    assert s.impropers.types == frame0.impropers.types
    numpy.testing.assert_array_equal(s.impropers.typeid,
                                     frame0.impropers.typeid)
    numpy.testing.assert_array_equal(s.impropers.group, frame0.impropers.group)

    assert s.constraints.N == frame0.constraints.N
    numpy.testing.assert_array_equal(s.constraints.value,
                                     frame0.constraints.value)
    numpy.testing.assert_array_equal(s.constraints.group,
                                     frame0.constraints.group)

    assert s.pairs.N == frame0.pairs.N
    assert s.pairs.types == frame0.pairs.types
    numpy.testing.assert_array_equal(s.pairs.typeid, frame0.pairs.typeid)
    numpy.testing.assert_array_equal(s.pairs.group, frame0.pairs.group)


def test_fallback(tmp_path, open_mode):
    """Test that properties fall back to defaults when the N changes."""
    frame0 = make_nondefault_frame()

    frame1 = gsd.hoomd.Frame()
    frame1.particles.N = 2
    frame1.particles.position = [[-2, -1, 0], [1, 3.0, 0.5]]
    frame1.bonds.N = None
    frame1.angles.N = None
    frame1.dihedrals.N = None
    frame1.impropers.N = None
    frame1.constraints.N = None
    frame1.pairs.N = None

    frame2 = gsd.hoomd.Frame()
    frame2.particles.N = 3
    frame2.particles.types = ['q', 's']
    frame2.particles.type_shapes = \
        [{}, {"type": "Ellipsoid", "a": 7.0, "b": 5.0, "c": 3.0}]
    frame2.bonds.N = 3
    frame2.angles.N = 4
    frame2.dihedrals.N = 5
    frame2.impropers.N = 6
    frame2.constraints.N = 4
    frame2.pairs.N = 7

    with gsd.hoomd.open(name=tmp_path / "test_fallback.gsd",
                        mode=open_mode.write) as hf:
        hf.extend([frame0, frame1, frame2])

    with gsd.hoomd.open(name=tmp_path / "test_fallback.gsd",
                        mode=open_mode.read) as hf:
        assert len(hf) == 3
        s = hf[0]

        assert_frames_equal(s, frame0)
        assert 'value' in s.log
        numpy.testing.assert_array_equal(s.log['value'], frame0.log['value'])

        # test that everything but position remained the same in frame 1
        s = hf[1]

        assert_frames_equal(s, frame0, check_position=False)
        assert 'value' in s.log
        numpy.testing.assert_array_equal(s.log['value'], frame0.log['value'])

        # check that the third frame goes back to defaults because it has a
        # different N
        s = hf[2]

        assert s.particles.N == 3
        assert s.particles.types == ['q', 's']
        assert s.particles.type_shapes == frame2.particles.type_shapes
        numpy.testing.assert_array_equal(
            s.particles.typeid, numpy.array([0, 0, 0], dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.particles.mass, numpy.array([1, 1, 1], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.diameter, numpy.array([1, 1, 1], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.body, numpy.array([-1, -1, -1], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.charge, numpy.array([0, 0, 0], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.moment_inertia,
            numpy.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.position,
            numpy.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.orientation,
            numpy.array([[1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]],
                        dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.velocity,
            numpy.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.angmom,
            numpy.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
                        dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.particles.image,
            numpy.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=numpy.int32))

        assert s.bonds.N == 3
        assert s.bonds.types == frame0.bonds.types
        numpy.testing.assert_array_equal(
            s.bonds.typeid, numpy.array([0, 0, 0], dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.bonds.group,
            numpy.array([[0, 0], [0, 0], [0, 0]], dtype=numpy.uint32))

        assert s.angles.N == 4
        assert s.angles.types == frame0.angles.types
        numpy.testing.assert_array_equal(
            s.angles.typeid, numpy.array([0, 0, 0, 0], dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.angles.group,
            numpy.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
                        dtype=numpy.uint32))

        assert s.dihedrals.N == 5
        assert s.dihedrals.types == frame0.dihedrals.types
        numpy.testing.assert_array_equal(
            s.dihedrals.typeid, numpy.array([0, 0, 0, 0, 0],
                                            dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.dihedrals.group,
            numpy.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                         [0, 0, 0, 0]],
                        dtype=numpy.uint32))

        assert s.impropers.N == 6
        assert s.impropers.types == frame0.impropers.types
        numpy.testing.assert_array_equal(
            s.impropers.typeid,
            numpy.array([0, 0, 0, 0, 0, 0], dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.impropers.group,
            numpy.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0],
                         [0, 0, 0, 0], [0, 0, 0, 0]],
                        dtype=numpy.uint32))

        assert s.constraints.N == 4
        numpy.testing.assert_array_equal(
            s.constraints.value, numpy.array([0, 0, 0, 0], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            s.constraints.group,
            numpy.array([[0, 0], [0, 0], [0, 0], [0, 0]], dtype=numpy.uint32))

        assert s.pairs.N == 7
        assert s.pairs.types == frame0.pairs.types
        numpy.testing.assert_array_equal(
            s.pairs.typeid, numpy.array([0] * 7, dtype=numpy.uint32))
        numpy.testing.assert_array_equal(
            s.pairs.group, numpy.array([[0, 0]] * 7, dtype=numpy.uint32))

        assert 'value' in s.log
        numpy.testing.assert_array_equal(s.log['value'], frame0.log['value'])


def test_fallback_to_frame0(tmp_path, open_mode):
    """Test that missing entries fall back to data in frame when N matches."""
    frame0 = make_nondefault_frame()

    frame1 = gsd.hoomd.Frame()
    frame1.configuration.step = 200000
    frame1.particles.N = None
    frame1.bonds.N = None
    frame1.angles.N = None
    frame1.dihedrals.N = None
    frame1.impropers.N = None
    frame1.constraints.N = None
    frame1.pairs.N = None

    with gsd.hoomd.open(name=tmp_path / "test_fallback2.gsd",
                        mode=open_mode.write) as hf:
        hf.extend([frame0, frame1])

    with gsd.hoomd.open(name=tmp_path / "test_fallback2.gsd",
                        mode=open_mode.read) as hf:
        assert len(hf) == 2

        s = hf[1]
        assert s.configuration.step == frame1.configuration.step
        assert_frames_equal(s, frame0, check_step=False)
        assert 'value' in s.log
        numpy.testing.assert_array_equal(s.log['value'], frame0.log['value'])


def test_no_fallback(tmp_path, open_mode):
    """Test that writes of default quantities do not fall back to frame 0."""
    frame0 = make_nondefault_frame()

    frame1 = gsd.hoomd.Frame()
    frame1.configuration.step = 200000
    frame1.configuration.dimensions = 3
    frame1.configuration.box = [1, 1, 1, 0, 0, 0]
    frame1.particles.N = frame0.particles.N
    frame1.particles.types = ['A']
    frame1.particles.typeid = [0] * frame0.particles.N
    frame1.particles.type_shapes = [{}]
    frame1.particles.mass = [1.0] * frame0.particles.N
    frame1.particles.charge = [0.0] * frame0.particles.N
    frame1.particles.diameter = [1.0] * frame0.particles.N
    frame1.particles.body = [-1] * frame0.particles.N
    frame1.particles.moment_inertia = [[0, 0, 0]] * frame0.particles.N
    frame1.particles.position = [[0, 0, 0]] * frame0.particles.N
    frame1.particles.orientation = [[1, 0, 0, 0]] * frame0.particles.N
    frame1.particles.velocity = [[0, 0, 0]] * frame0.particles.N
    frame1.particles.angmom = [[0, 0, 0, 0]] * frame0.particles.N
    frame1.particles.image = [[0, 0, 0]] * frame0.particles.N

    frame1.bonds.N = frame0.bonds.N
    frame1.bonds.types = ['A']
    frame1.bonds.typeid = [0] * frame0.bonds.N
    frame1.bonds.group = [[0, 0]] * frame0.bonds.N

    frame1.angles.N = frame0.angles.N
    frame1.angles.types = ['A']
    frame1.angles.typeid = [0] * frame0.angles.N
    frame1.angles.group = [[0, 0, 0]] * frame0.angles.N

    frame1.dihedrals.N = frame0.dihedrals.N
    frame1.dihedrals.types = ['A']
    frame1.dihedrals.typeid = [0] * frame0.dihedrals.N
    frame1.dihedrals.group = [[0, 0, 0, 0]] * frame0.dihedrals.N

    frame1.impropers.N = frame0.impropers.N
    frame1.impropers.types = ['A']
    frame1.impropers.typeid = [0] * frame0.impropers.N
    frame1.impropers.group = [[0, 0, 0, 0]] * frame0.impropers.N

    frame1.constraints.N = frame0.constraints.N
    frame1.constraints.value = [0] * frame0.constraints.N
    frame1.constraints.group = [[0, 0]] * frame0.constraints.N

    frame1.pairs.N = frame0.pairs.N
    frame1.pairs.types = ['A']
    frame1.pairs.typeid = [0] * frame0.pairs.N
    frame1.pairs.group = [[0, 0]] * frame0.pairs.N

    with gsd.hoomd.open(name=tmp_path / "test_no_fallback.gsd",
                        mode=open_mode.write) as hf:
        hf.extend([frame0, frame1])

    with gsd.hoomd.open(name=tmp_path / "test_no_fallback.gsd",
                        mode=open_mode.read) as hf:
        assert len(hf) == 2

        s = hf[1]
        assert s.configuration.step == frame1.configuration.step
        assert_frames_equal(s, frame1)


def test_iteration(tmp_path, open_mode):
    """Test the iteration protocols for hoomd trajectories."""
    with gsd.hoomd.open(name=tmp_path / "test_iteration.gsd",
                        mode=open_mode.write) as hf:
        hf.extend((create_frame(i) for i in range(20)))

    with gsd.hoomd.open(name=tmp_path / "test_iteration.gsd",
                        mode=open_mode.read) as hf:
        step = hf[-1].configuration.step
        assert step == 20

        step = hf[-2].configuration.step
        assert step == 19

        step = hf[-3].configuration.step
        assert step == 18

        step = hf[0].configuration.step
        assert step == 1

        step = hf[-20].configuration.step
        assert step == 1

        with pytest.raises(IndexError):
            step = hf[-21].configuration.step

        with pytest.raises(IndexError):
            step = hf[20]

        frames = hf[5:10]
        steps = [frame.configuration.step for frame in frames]
        assert steps == [6, 7, 8, 9, 10]

        frames = hf[15:50]
        steps = [frame.configuration.step for frame in frames]
        assert steps == [16, 17, 18, 19, 20]

        frames = hf[15:-3]
        steps = [frame.configuration.step for frame in frames]
        assert steps == [16, 17]


def test_slicing_and_iteration(tmp_path, open_mode):
    """Test that hoomd trajectories can be sliced."""
    with gsd.hoomd.open(name=tmp_path / "test_slicing.gsd",
                        mode=open_mode.write) as hf:
        hf.extend((create_frame(i) for i in range(20)))

    with gsd.hoomd.open(name=tmp_path / "test_slicing.gsd",
                        mode=open_mode.read) as hf:
        # Test len()-function on trajectory and sliced trajectory.
        assert len(hf) == 20
        assert len(hf[:10]) == 10

        # Test len()-function with explicit iterator.
        assert len(iter(hf)) == len(hf)
        assert len(iter(hf[:10])) == len(hf[:10])

        # Test iteration with implicit iterator.
        # All iterations are run twice to check for issues
        # with iterator exhaustion.
        assert len(list(hf)) == len(hf)
        assert len(list(hf)) == len(hf)
        assert len(list(hf[:10])) == len(hf[:10])
        assert len(list(hf[:10])) == len(hf[:10])

        # Test iteration with explicit iterator.
        hf_iter = iter(hf)
        assert len(hf_iter) == len(hf)  # sanity check
        assert len(list(hf_iter)) == len(hf)
        assert len(list(hf_iter)) == len(hf)

        # Test iteration with explicit sliced iterator.
        hf_iter = iter(hf[:10])
        assert len(hf_iter) == 10  # sanity check
        assert len(list(hf_iter)) == 10
        assert len(list(hf_iter)) == 10

        # Test frame selection
        with pytest.raises(IndexError):
            hf[len(hf)]
        assert hf[0].configuration.step == hf[0].configuration.step
        assert hf[len(hf) - 1].configuration.step == hf[-1].configuration.step


def test_view_slicing_and_iteration(tmp_path, open_mode):
    """Test that trajectories can be sliced."""
    with gsd.hoomd.open(name=tmp_path / "test_slicing.gsd",
                        mode=open_mode.write) as hf:
        hf.extend((create_frame(i) for i in range(40)))

    with gsd.hoomd.open(name=tmp_path / "test_slicing.gsd",
                        mode=open_mode.read) as hf:
        view = hf[::2]

        # Test len()-function on trajectory and sliced view.
        assert len(view) == 20
        assert len(view[:10]) == 10
        assert len(view[::2]) == 10

        # Test len()-function with explicit iterator.
        assert len(iter(view)) == len(view)
        assert len(iter(view[:10])) == len(view[:10])

        # Test iteration with implicit iterator.
        # All iterations are run twice to check for issues
        # with iterator exhaustion.
        assert len(list(view)) == len(view)
        assert len(list(view)) == len(view)
        assert len(list(view[:10])) == len(view[:10])
        assert len(list(view[:10])) == len(view[:10])
        assert len(list(view[::2])) == len(view[::2])
        assert len(list(view[::2])) == len(view[::2])

        # Test iteration with explicit iterator.
        view_iter = iter(view)
        assert len(view_iter) == len(view)  # sanity check
        assert len(list(view_iter)) == len(view)
        assert len(list(view_iter)) == len(view)

        # Test iteration with explicit sliced iterator.
        view_iter = iter(view[:10])
        assert len(view_iter) == 10  # sanity check
        assert len(list(view_iter)) == 10
        assert len(list(view_iter)) == 10

        # Test frame selection
        with pytest.raises(IndexError):
            view[len(view)]
        assert view[0].configuration.step == view[0].configuration.step
        assert view[len(view)
                    - 1].configuration.step == view[-1].configuration.step


def test_truncate(tmp_path):
    """Test the truncate API."""
    with gsd.hoomd.open(name=tmp_path / "test_iteration.gsd", mode='wb') as hf:
        hf.extend((create_frame(i) for i in range(20)))

        assert len(hf) == 20
        s = hf[10]  # noqa
        assert hf._initial_frame is not None

        hf.truncate()
        assert len(hf) == 0
        assert hf._initial_frame is None


def test_state(tmp_path, open_mode):
    """Test the state chunks."""
    frame0 = gsd.hoomd.Frame()

    frame0.state['hpmc/sphere/radius'] = [2.0]
    frame0.state['hpmc/sphere/orientable'] = [1]

    frame1 = gsd.hoomd.Frame()

    frame1.state['hpmc/convex_polyhedron/N'] = [3]
    frame1.state['hpmc/convex_polyhedron/vertices'] = [[-1, -1, -1], [0, 1, 1],
                                                       [1, 0, 0]]

    with gsd.hoomd.open(name=tmp_path / "test_state.gsd",
                        mode=open_mode.write) as hf:
        hf.extend([frame0, frame1])

    with gsd.hoomd.open(name=tmp_path / "test_state.gsd",
                        mode=open_mode.read) as hf:
        assert len(hf) == 2
        s = hf[0]

        numpy.testing.assert_array_equal(s.state['hpmc/sphere/radius'],
                                         frame0.state['hpmc/sphere/radius'])
        numpy.testing.assert_array_equal(s.state['hpmc/sphere/orientable'],
                                         frame0.state['hpmc/sphere/orientable'])

        s = hf[1]

        numpy.testing.assert_array_equal(
            s.state['hpmc/convex_polyhedron/N'],
            frame1.state['hpmc/convex_polyhedron/N'])
        numpy.testing.assert_array_equal(
            s.state['hpmc/convex_polyhedron/vertices'],
            frame1.state['hpmc/convex_polyhedron/vertices'])


def test_log(tmp_path, open_mode):
    """Test the log chunks."""
    frame0 = gsd.hoomd.Frame()

    frame0.log['particles/net_force'] = [[1, 2, 3], [4, 5, 6]]
    frame0.log['particles/pair_lj_energy'] = [0, -5, -8, -3]
    frame0.log['value/potential_energy'] = [10]
    frame0.log['value/pressure'] = [-3]

    frame1 = gsd.hoomd.Frame()

    frame1.log['particles/pair_lj_energy'] = [1, 2, -4, -10]
    frame1.log['value/pressure'] = [5]

    with gsd.hoomd.open(name=tmp_path / "test_log.gsd",
                        mode=open_mode.write) as hf:
        hf.extend([frame0, frame1])

    with gsd.hoomd.open(name=tmp_path / "test_log.gsd",
                        mode=open_mode.read) as hf:
        assert len(hf) == 2
        s = hf[0]

        numpy.testing.assert_array_equal(s.log['particles/net_force'],
                                         frame0.log['particles/net_force'])
        numpy.testing.assert_array_equal(s.log['particles/pair_lj_energy'],
                                         frame0.log['particles/pair_lj_energy'])
        numpy.testing.assert_array_equal(s.log['value/potential_energy'],
                                         frame0.log['value/potential_energy'])
        numpy.testing.assert_array_equal(s.log['value/pressure'],
                                         frame0.log['value/pressure'])

        s = hf[1]

        # unspecified entries pull from frame 0
        numpy.testing.assert_array_equal(s.log['particles/net_force'],
                                         frame0.log['particles/net_force'])
        numpy.testing.assert_array_equal(s.log['value/potential_energy'],
                                         frame0.log['value/potential_energy'])

        # specified entries are different in frame 1
        numpy.testing.assert_array_equal(s.log['particles/pair_lj_energy'],
                                         frame1.log['particles/pair_lj_energy'])
        numpy.testing.assert_array_equal(s.log['value/pressure'],
                                         frame1.log['value/pressure'])


def test_pickle(tmp_path, open_mode):
    """Test that hoomd trajectory objects can be pickled."""
    with gsd.hoomd.open(name=tmp_path / "test_pickling.gsd",
                        mode=open_mode.write) as traj:
        traj.extend((create_frame(i) for i in range(20)))
        with pytest.raises(pickle.PickleError):
            pkl = pickle.dumps(traj)
    with gsd.hoomd.open(name=tmp_path / "test_pickling.gsd",
                        mode=open_mode.read) as traj:
        pkl = pickle.dumps(traj)
        with pickle.loads(pkl) as hf:
            assert len(hf) == 20


@pytest.mark.parametrize(
    'container',
    ['particles', 'bonds', 'angles', 'dihedrals', 'impropers', 'pairs'])
def test_no_duplicate_types(tmp_path, container):
    """Test that duplicate types raise an error."""
    with gsd.hoomd.open(name=tmp_path / "test_create.gsd", mode='wb') as hf:

        frame = gsd.hoomd.Frame()

        getattr(frame, container).types = ['A', 'B', 'B', 'C']

        with pytest.raises(ValueError):
            hf.append(frame)


def test_read_log(tmp_path):
    """Test that data logged in gsd files are read correctly."""
    frame0 = gsd.hoomd.Frame()
    frame0.log['particles/pair_lj_energy'] = [0, -5, -8, -3]
    frame0.log['particles/pair_lj_force'] = [
        (0, 0, 0),
        (1, 1, 1),
        (2, 2, 2),
        (3, 3, 3),
    ]
    frame0.log['value/potential_energy'] = [10]
    frame0.log['value/pressure'] = [-3]

    frame1 = gsd.hoomd.Frame()
    frame1.configuration.step = 1
    frame1.log['particles/pair_lj_energy'] = [1, 2, -4, -10]
    frame1.log['particles/pair_lj_force'] = [
        (1, 1, 1),
        (2, 2, 2),
        (3, 3, 3),
        (4, 4, 4),
    ]
    frame1.log['value/pressure'] = [5]

    with gsd.hoomd.open(name=tmp_path / "test_log.gsd", mode='wb') as hf:
        hf.extend([frame0, frame1])

    # Test scalar_only = False
    logged_data_dict = gsd.hoomd.read_log(name=tmp_path / "test_log.gsd",
                                          scalar_only=False)

    assert len(logged_data_dict) == 5
    assert list(logged_data_dict.keys()) == [
        'configuration/step', 'log/particles/pair_lj_energy',
        'log/particles/pair_lj_force', 'log/value/potential_energy',
        'log/value/pressure'
    ]

    numpy.testing.assert_array_equal(logged_data_dict['configuration/step'],
                                     [0, 1])
    numpy.testing.assert_array_equal(
        logged_data_dict['log/particles/pair_lj_energy'], [
            frame0.log['particles/pair_lj_energy'],
            frame1.log['particles/pair_lj_energy']
        ])
    numpy.testing.assert_array_equal(
        logged_data_dict['log/particles/pair_lj_force'], [
            frame0.log['particles/pair_lj_force'],
            frame1.log['particles/pair_lj_force']
        ])
    numpy.testing.assert_array_equal(
        logged_data_dict['log/value/potential_energy'], [
            *frame0.log['value/potential_energy'],
            *frame0.log['value/potential_energy']
        ])
    numpy.testing.assert_array_equal(
        logged_data_dict['log/value/pressure'],
        [*frame0.log['value/pressure'], *frame1.log['value/pressure']])

    # Test scalar_only = True
    logged_data_dict = gsd.hoomd.read_log(name=tmp_path / "test_log.gsd",
                                          scalar_only=True)
    assert len(logged_data_dict) == 3
    assert list(logged_data_dict.keys()) == [
        'configuration/step', 'log/value/potential_energy', 'log/value/pressure'
    ]
    numpy.testing.assert_array_equal(logged_data_dict['configuration/step'],
                                     [0, 1])
    numpy.testing.assert_array_equal(
        logged_data_dict['log/value/potential_energy'], [
            *frame0.log['value/potential_energy'],
            *frame0.log['value/potential_energy']
        ])
    numpy.testing.assert_array_equal(
        logged_data_dict['log/value/pressure'],
        [*frame0.log['value/pressure'], *frame1.log['value/pressure']])


def test_read_log_warning(tmp_path):
    """Test that read_log issues a warning."""
    frame = gsd.hoomd.Frame()

    with gsd.hoomd.open(name=tmp_path / "test_log.gsd", mode='wb') as hf:
        hf.extend([frame])

    with pytest.warns(RuntimeWarning):
        log = gsd.hoomd.read_log(tmp_path / "test_log.gsd")

    assert list(log.keys()) == ['configuration/step']
