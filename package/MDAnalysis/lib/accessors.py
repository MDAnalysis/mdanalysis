from functools import partial

from .. import _CONVERTERS
from ..core._get_readers import get_converter_for


class Accessor:
    def __init__(self, accessor):
        self._accessor = accessor

    def __get__(self, obj, cls):
        return self._accessor(obj)


class ConverterAccessor:
    def __init__(self, ag):
        self._ag = ag
        for lib, converter in _CONVERTERS.items():
            method_name = lib.lower()
            fconvert = partial(converter().convert, self._ag)
            setattr(self, method_name, fconvert)

    def __call__(self, package, **kwargs):
        converter = get_converter_for(package.upper())
        return converter().convert(self._ag, **kwargs)