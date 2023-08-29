from typing import Callable
class Multiprocessing:
     
    def __init__(self):
        pass
    
    @classmethod
    def apply(self, func: Callable, computations: list, n_workers: int) -> list:
        from multiprocessing import Pool

        with Pool(processes=n_workers) as pool:
            results = pool.map(func, computations)
        return results


class Dask:

    def __init__(self):
        pass

    @classmethod
    def apply(self, func: Callable, computations: list, n_workers: int) -> list:
        from dask.delayed import delayed
        import dask

        computations = [delayed(func)(task) for task in computations]
        results = dask.compute(computations, scheduler='processes', chunksize=1, n_workers=n_workers)[0]
        return results


class DaskDistributed:

    def __init__(self, client):
        self.client = client
    
    def apply(self, func: Callable, computations: list, n_workers: int) -> list:
        """Apply the function to the list of callables.
        n_workers is ignored - all workers from the client will be used.
        """
        from dask.delayed import delayed

        computations = [delayed(func)(task) for task in computations]
        return [obj.result() for obj in self.client.compute(computations)]
