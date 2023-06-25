import pytest

@pytest.fixture(params=[None, 
                       'dask:1', 'dask:2', 
                       'multiprocessing:1', 'multiprocessing:2'
                       ], 
                scope='module')
def schedulers_all(request):
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
def scheduler_only_current_process(request):
    """
    Fixture for testing only current process.
    Made to allow for easier replacement of fixture in the future
    upon implementation of parallelization for certain classes.

    If scheduler raises an exception, you should explicitly check it
    using `with pytest.raises(NotImplementedError)`
    """
    return {'scheduler':None}

@pytest.fixture(params=[None, 'dask:1', 'dask:2'])
def scheduler_current_or_dask(request):
    """
    Fixture for testing dask or current process (because multiprocessing is worse with serialization).

    If scheduler raises an exception, you should explicitly check it
    using `with pytest.raises(NotImplementedError)`
    """
    if request.param is None:
        return {'scheduler':None}
    else:
        scheduler, n_workers = request.param.split(':')    
        n_workers = int(n_workers)
        return {'scheduler':scheduler, 'n_workers':n_workers}
