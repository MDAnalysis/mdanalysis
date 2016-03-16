# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from MDAnalysis.lib import autocorrelation

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_
from six.moves import range



class TestAutoCorrelation(object):
    @staticmethod
    def reference(vecs):
        # This is a numpy version of what the cython code should do
        n = vecs.shape[0]
        result = np.zeros(n, dtype=np.float32)
        for i in range(n):
            pos = 0
            for j in range(i, n):
                result[pos] += np.dot(vecs[i], vecs[j])
                pos += 1
        for i in range(n):
            result[i] /= n - i
        return result

    def test_autocorrelation(self):
        vecs = np.arange(30).reshape(10, 3).astype(np.float32)

        result = autocorrelation.vector_autocorrelation(vecs)

        assert_(result.shape == (10,))
        assert_(result.dtype == np.float32)
        assert_array_almost_equal(self.reference(vecs), result)


class TestWindowedAutoCorrelation(object):
    @staticmethod
    def reference(vecs, window):
        # numpy version of cython implementation
        result = np.zeros(window)

        n = vecs.shape[0] - window + 1

        for i in range(n):
            for j in range(window):
                result[j] += np.dot(vecs[i], vecs[i + j])

        for i in range(window):
            result[i] /= n

        return result

    def test_windowed_autocorrelation(self):
        vecs = np.arange(30).reshape(10, 3).astype(np.float32)
        window = 10

        result = autocorrelation.windowed_vector_autocorrelation(vecs, window)

        assert_(result.shape == (window,))
        assert_(result.dtype == np.float32)
        assert_array_almost_equal(self.reference(vecs, window), result)
