import pickle
from collections import UserDict

import numpy as np
import pytest
from MDAnalysis.analysis import results as results_module
from numpy.testing import assert_equal


class Test_Results:
    @pytest.fixture
    def results(self):
        return results_module.Results(a=1, b=2)

    def test_get(self, results):
        assert results.a == results["a"] == 1

    def test_no_attr(self, results):
        msg = "'Results' object has no attribute 'c'"
        with pytest.raises(AttributeError, match=msg):
            results.c

    def test_set_attr(self, results):
        value = [1, 2, 3, 4]
        results.c = value
        assert results.c is results["c"] is value

    def test_set_key(self, results):
        value = [1, 2, 3, 4]
        results["c"] = value
        assert results.c is results["c"] is value

    @pytest.mark.parametrize("key", dir(UserDict) + ["data"])
    def test_existing_dict_attr(self, results, key):
        msg = f"'{key}' is a protected dictionary attribute"
        with pytest.raises(AttributeError, match=msg):
            results[key] = None

    @pytest.mark.parametrize("key", dir(UserDict) + ["data"])
    def test_wrong_init_type(self, key):
        msg = f"'{key}' is a protected dictionary attribute"
        with pytest.raises(AttributeError, match=msg):
            results_module.Results(**{key: None})

    @pytest.mark.parametrize("key", ("0123", "0j", "1.1", "{}", "a b"))
    def test_weird_key(self, results, key):
        msg = f"'{key}' is not a valid attribute"
        with pytest.raises(ValueError, match=msg):
            results[key] = None

    def test_setattr_modify_item(self, results):
        mylist = [1, 2]
        mylist2 = [3, 4]
        results.myattr = mylist
        assert results.myattr is mylist
        results["myattr"] = mylist2
        assert results.myattr is mylist2
        mylist2.pop(0)
        assert len(results.myattr) == 1
        assert results.myattr is mylist2

    def test_setitem_modify_item(self, results):
        mylist = [1, 2]
        mylist2 = [3, 4]
        results["myattr"] = mylist
        assert results.myattr is mylist
        results.myattr = mylist2
        assert results.myattr is mylist2
        mylist2.pop(0)
        assert len(results["myattr"]) == 1
        assert results["myattr"] is mylist2

    def test_delattr(self, results):
        assert hasattr(results, "a")
        delattr(results, "a")
        assert not hasattr(results, "a")

    def test_missing_delattr(self, results):
        assert not hasattr(results, "d")
        msg = "'Results' object has no attribute 'd'"
        with pytest.raises(AttributeError, match=msg):
            delattr(results, "d")

    def test_pop(self, results):
        assert hasattr(results, "a")
        results.pop("a")
        assert not hasattr(results, "a")

    def test_update(self, results):
        assert not hasattr(results, "spudda")
        results.update({"spudda": "fett"})
        assert results.spudda == "fett"

    def test_update_data_fail(self, results):
        msg = f"'data' is a protected dictionary attribute"
        with pytest.raises(AttributeError, match=msg):
            results.update({"data": 0})

    def test_pickle(self, results):
        results_p = pickle.dumps(results)
        results_new = pickle.loads(results_p)

    @pytest.mark.parametrize(
        "args, kwargs, length",
        [
            (({"darth": "tater"},), {}, 1),
            ([], {"darth": "tater"}, 1),
            (({"darth": "tater"},), {"yam": "solo"}, 2),
            (({"darth": "tater"},), {"darth": "vader"}, 1),
        ],
    )
    def test_initialize_arguments(self, args, kwargs, length):
        results = results_module.Results(*args, **kwargs)
        ref = dict(*args, **kwargs)
        assert ref == results
        assert len(results) == length

    def test_different_instances(self, results):
        new_results = results_module.Results(darth="tater")
        assert new_results.data is not results.data


class Test_ResultsGroup:
    @pytest.fixture
    def results_0(self):
        return results_module.Results(
            sequence=[0],
            ndarray_mean=np.array([0, 0, 0]),
            ndarray_sum=np.array([0, 0, 0, 0]),
            float=0.0,
            float_sum=0.0,
        )

    @pytest.fixture
    def results_1(self):
        return results_module.Results(
            sequence=[1],
            ndarray_mean=np.array([1, 1, 1]),
            ndarray_sum=np.array([1, 1, 1, 1]),
            float=1.0,
            float_sum=1.0,
        )

    @pytest.fixture
    def merger(self):
        RG = results_module.ResultsGroup
        lookup = {
            "sequence": RG.flatten_sequence,
            "ndarray_mean": RG.ndarray_mean,
            "ndarray_sum": RG.ndarray_sum,
            "float": RG.float_mean,
            "float_sum": lambda floats: sum(floats),
        }
        return RG(lookup=lookup)

    @pytest.mark.parametrize("n", [1, 2, 5, 14])
    def test_all_results(self, results_0, results_1, merger, n):
        from itertools import cycle

        objects = [obj for obj, _ in zip(cycle([results_0, results_1]), range(n))]

        arr = [i for _, i in zip(range(n), cycle([0, 1]))]
        answers = {
            "sequence": arr,
            "ndarray_mean": [np.mean(arr) for _ in range(3)],
            "ndarray_sum": [np.sum(arr) for _ in range(4)],
            "float": np.mean(arr),
            "float_sum": np.sum(arr),
        }

        results = merger.merge(objects)
        for attr, merged_value in results.items():
            assert_equal(merged_value, answers.get(attr), err_msg=f"{attr=}, {merged_value=}, {arr=}, {objects=}")
