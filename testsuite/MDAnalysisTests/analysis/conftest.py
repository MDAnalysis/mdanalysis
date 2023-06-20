import pytest

@pytest.fixture(params=[None, 
                       'dask:1', 'dask:2', 
                    #    'multiprocessing:1', 'multiprocessing:2'
                       ], 
                scope='module')
def scheduler(request):
    """
    Fixture for testing all possible schedulers.
    If scheduler raises an exception, you should explicitly check it
    using `with pytest.raises(NotImplementedError)`
    """
    if request.param is None:
        return {'scheduler':None}
    else:
        scheduler, n_workers = request.param.split(':')    
        n_workers = int(n_workers)
        return {'scheduler':scheduler, 'n_workers':n_workers}

@pytest.fixture(params=[None])
def localscheduler(request):
    """
    Fixture for testing only local, i.e. before-dask, codepath.
    """
    return request.param