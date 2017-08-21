import pytest


class A(object):

    def test_value(self, request):
        request.cls.filename == "foo"

@pytest.mark.parametrize('filename', (
    'foo',
    'bar'
))
class TestA(A):
    pass