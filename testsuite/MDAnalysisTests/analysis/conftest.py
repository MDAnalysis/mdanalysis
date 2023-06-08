import pytest

@pytest.fixture(params=[None, 'localdask'])
def scheduler(request):
    """
    Fixture for testing all possible schedulers.
    If scheduler raises an exception, you should explicitly check it
    using `with pytest.raises(NotImplementedError)`
    """
    return request.param
