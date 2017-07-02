from unittest import TestCase

import pytest


class _IgnoreClass(TestCase):
    def test_number(self):
        assert 2 == 2

    def test_bool(self):
        assert True is True

    def dont_collect(self):
        assert 1 == 1


class IgnoreClass2(object):
    def test_number(self):
        assert 2 == 2

    def test_bool(self):
        assert True is False


class ParentClass(TestCase):
    __test__ = False

    def test_value(self):
        print('chala')
        assert self.a == 2


class TestChildClass(ParentClass):
    __test__ = True

    a = 2


def check_even(n, nn):
    print(n)
    assert n % 2 == 0 or nn % 2 == 0


def test_evens():
    for i in range(0, 5):
        yield check_even, i, i * 3


@pytest.fixture
def fx1():
    print('Worked')


class TestFixture:
    def test_fx(self, fx1):
        print(fx1)
        print(fx)
