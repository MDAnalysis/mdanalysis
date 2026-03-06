# Copyright (c) 2006-2021  Andrey Golovizin
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


"""Miscellaneous small utils."""
from __future__ import print_function, unicode_literals

import itertools
from collections import OrderedDict, deque
from functools import wraps
from types import GeneratorType
try:
    from collections.abc import MutableMapping, MutableSet, Sequence
except ImportError:
    from collections import MutableMapping, MutableSet, Sequence

from itertools import zip_longest


def deprecated(since, reason=None):
    def decorator(f):
        @wraps(f)
        def new_f(*args, **kwargs):
            import warnings
            message = u'{0}() is deprecated since {1}'.format(f.__name__, since)
            if reason:
                message += ': {0}'.format(reason)
            warnings.warn(message, DeprecationWarning, stacklevel=2)
            return f(*args, **kwargs)
        return new_f
    return decorator


def memoize(f, capacity=1024):
    memory = {}
    history = deque()

    @wraps(f)
    def new_f(*args):
        if args not in memory:
            if len(history) >= capacity:
                del memory[history.popleft()]
            memory[args] = f(*args)
            history.append(args)
        return memory[args]
    return new_f


def collect_iterable(f):
    @wraps(f)
    def new_f(*args, **kwargs):
        return list(f(*args, **kwargs))
    return new_f


def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip_longest(a, b)


class CaseInsensitiveDict(MutableMapping):
    """A dict with case-insensitive lookup.

    >>> d = CaseInsensitiveDict(TesT='passed')
    >>> d
    CaseInsensitiveDict({'TesT': 'passed'})
    >>> d.lower()
    CaseInsensitiveDict({'test': 'passed'})
    >>> print(d['TesT'])
    passed
    >>> print(d['test'])
    passed
    >>> print(d['Test'])
    passed
    >>> print(d.get('test'))
    passed
    >>> print(d.get('Test'))
    passed

    >>> d['Test'] = 'passed again'
    >>> print(d['test'])
    passed again
    >>> 'test' in d
    True
    >>> 'Test' in d
    True
    >>> list(d.keys())
    ['Test']
    >>> for key in d:
    ...     print(key)
    Test
    >>> for key, value in d.items():
    ...     print(key, value)
    Test passed again
    >>> bool(d)
    True
    >>> len(d)
    1

    >>> del d['test']
    >>> len(d)
    0
    >>> bool(d)
    False
    >>> 'test' in d
    False
    >>> 'Test' in d
    False
    >>> print(d.get('test'))
    None
    >>> print(d.get('Test'))
    None
    >>> print(d.get('Test', 'failed'))
    failed

    >>> CaseInsensitiveDict(
    ...     (key, value) for key, value in [('a', 'b')]
    ... )
    CaseInsensitiveDict({'a': 'b'})

    """

    def __init__(self, *args, **kwargs):
        initial = dict(*args, **kwargs)
        self._dict = dict((key.lower(), value) for key, value in initial.items())
        self._keys = dict((key.lower(), key) for key in initial)

    def __len__(self):
        return len(self._dict)

    def __iter__(self):
        return iter(self._keys.values())

    def __setitem__(self, key, value):
        """To implement lowercase keys."""
        key_lower = key.lower()
        self._dict[key_lower] = value
        self._keys[key_lower] = key

    def __getitem__(self, key):
        return self._dict[key.lower()]

    def __delitem__(self, key):
        key_lower = key.lower()
        del self._dict[key_lower]
        del self._keys[key_lower]

    def __contains__(self, key):
        return key.lower() in self._dict

    def __repr__(self):
        """A caselessDict version of __repr__ """
        dct = dict((key, self[key]) for key in self)
        return '{0}({1})'.format(
            type(self).__name__, repr(dct),
        )

    def items_lower(self):
        return ((key.lower(), value) for key, value in self.items())

    def lower(self):
        return type(self)(self.items_lower())


class CaseInsensitiveDefaultDict(CaseInsensitiveDict):
    """CaseInseisitiveDict with default factory, like collections.defaultdict

    >>> d = CaseInsensitiveDefaultDict(int)
    >>> d['a']
    0
    >>> d['a'] += 1
    >>> d['a']
    1
    >>> d['A']
    1
    >>> d['a'] = 3
    >>> d['a']
    3
    >>> d['B'] += 10
    >>> d['b']
    10

    """
    def __init__(self, default_factory):
        super(CaseInsensitiveDefaultDict, self).__init__()
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return super(CaseInsensitiveDefaultDict, self).__getitem__(key)
        except KeyError:
            return self.default_factory()


class OrderedCaseInsensitiveDict(CaseInsensitiveDict):
    """ An (incomplete) ordered case-insensitive dict.

    >>> d = OrderedCaseInsensitiveDict([
    ...     ('Uno', 1),
    ...     ('Dos', 2),
    ...     ('Tres', 3),
    ... ])
    >>> d
    OrderedCaseInsensitiveDict([('Uno', 1), ('Dos', 2), ('Tres', 3)])
    >>> d.lower()
    OrderedCaseInsensitiveDict([('uno', 1), ('dos', 2), ('tres', 3)])
    >>> list(d.keys())
    ['Uno', 'Dos', 'Tres']
    >>> list(d.items())
    [('Uno', 1), ('Dos', 2), ('Tres', 3)]
    >>> list(d.values())
    [1, 2, 3]
    >>> d['Cuatro'] = 4
    >>> list(d.keys())
    ['Uno', 'Dos', 'Tres', 'Cuatro']
    >>> list(d.items())
    [('Uno', 1), ('Dos', 2), ('Tres', 3), ('Cuatro', 4)]
    >>> list(d.values())
    [1, 2, 3, 4]
    >>> list(d)
    ['Uno', 'Dos', 'Tres', 'Cuatro']
    >>> 'Uno' in d
    True
    >>> 'uno' in d
    True
    >>> d['Uno']
    1
    >>> d['uno']
    1
    >>> d['UNO']
    1
    >>> 'Cuatro' in d
    True
    >>> 'CUATRO' in d
    True
    >>> d['Cuatro']
    4
    >>> d['cuatro']
    4
    >>> d['UNO'] = 'one'
    >>> d['uno']
    'one'
    >>> d['Uno']
    'one'
    >>> list(d.keys())
    ['UNO', 'Dos', 'Tres', 'Cuatro']
    >>> d['cuatro'] = 'four'
    >>> d['Cuatro']
    'four'
    >>> d['cuatro']
    'four'
    >>> list(d.keys())
    ['UNO', 'Dos', 'Tres', 'cuatro']
    >>> list(d.values())
    ['one', 2, 3, 'four']
    >>> del d['dos']
    >>> list(d.keys())
    ['UNO', 'Tres', 'cuatro']
    >>> list(d.values())
    ['one', 3, 'four']
    """

    def __init__(self, *args, **kwargs):
        initial = OrderedDict(*args, **kwargs)
        self._dict = dict((key.lower(), value) for key, value in initial.items())
        self._keys = OrderedDict((key.lower(), key) for key in initial)

    def __repr__(self):
        return '{0}({1})'.format(
            type(self).__name__, list(self.items()),
        )


class CaseInsensitiveSet(MutableSet):
    """A very basic case-insensitive set.

    >>> s = CaseInsensitiveSet()
    >>> len(s)
    0
    >>> 'a' in s
    False
    >>> list(CaseInsensitiveSet(['aaa', 'Aaa', 'AAA']))
    ['aaa']
    >>> s = CaseInsensitiveSet(['Aaa', 'Bbb'])
    >>> s
    CaseInsensitiveSet(['Aaa', 'Bbb'])
    >>> s.lower()
    CaseInsensitiveSet(['aaa', 'bbb'])
    >>> len(s)
    2
    >>> 'aaa' in s
    True
    >>> 'Aaa' in s
    True
    >>> 'AAA' in s
    True
    >>> 'bbb' in s
    True
    >>> 'Bbb' in s
    True
    >>> 'abc' in s
    False
    >>> s.add('ccc')
    >>> len(s)
    3
    >>> 'aaa' in s
    True
    >>> 'ccc' in s
    True
    >>> s.remove('AAA')
    >>> len(s)
    2
    >>> 'aaa' in s
    False

    >>> bool(CaseInsensitiveSet(['a']))
    True
    >>> bool(CaseInsensitiveSet([]))
    False
    >>> bool(CaseInsensitiveSet())
    False

    """

    def __init__(self, iterable=()):
        self._set = set()
        self._keys = dict()
        for item in iterable:
            self.add(item)

    def __contains__(self, key):
        return key.lower() in self._set

    def __iter__(self):
        return iter(self._set)

    def __len__(self):
        return len(self._set)

    def __repr__(self):
        """A caselessDict version of __repr__ """
        return '{0}({1})'.format(
            type(self).__name__, repr(sorted(self._keys.values()))
        )

    def add(self, key):
        key_lower = key.lower()
        self._set.add(key_lower)
        self._keys[key_lower] = key

    def discard(self, key):
        key_lower = key.lower()
        self._set.discard(key_lower)
        self._keys.pop(key_lower, None)

    def get_canonical_key(self, key):
        return self._keys[key.lower()]

    def lower(self):
        return type(self)(self._set)
