# Implementation of Scientific.Geometry.Vector in Pyrex
#
# Written by Konrad Hinsen
# last revision: 2005-6-1
#


cdef extern from "math.h":

    double sqrt(double x)
    double acos(double x)

import Numeric

#
# For efficiency reasons (calling __init__ makes the creation of a vector
# rather expensive), most of the operations happen in class "vector", which
# is not meant to be used directly in application code. Objects of this
# class are initialized with zero values and then initialized explicitly
# by calling the method "set".
# Application code should use the derived class "Vector", which adds
# only the __init__ method.
#
cdef class vector:

    cdef double xv, yv, zv

    property is_vector:
        def __get__(self):
            return 1

    property array:
        def __get__(self):
            return Numeric.array([self.xv, self.yv, self.zv])

    def __copy__(self, memo = None):
        return self

    def __deepcopy__(self, memo = None):
        return self

    def __getstate__(self):
        return [self.xv, self.yv, self.zv]

    def __setstate__(self, state):
        self.xv, self.yv, self.zv = state

    def __reduce__(self):
        return (Vector, (self.xv, self.yv, self.zv))

    cdef void set(self, double x, double y, double z):
        self.xv = x
        self.yv = y
        self.zv = z

    def x(self):
        "Returns the x coordinate."
        return self.xv

    def y(self):
        "Returns the y coordinate."
        return self.yv

    def z(self):
        "Returns the z coordinate."
        return self.zv

    def __repr__(self):
        return 'Vector(%f,%f,%f)' % (self.xv, self.yv, self.zv)

    def __str__(self):
        return str([self.xv, self.yv, self.zv])

    def __add__(vector self, vector other):
        result = vector()
        vector.set(result,
                   self.xv+other.xv, self.yv+other.yv, self.zv+other.zv)
        return result

    def __neg__(vector self):
        result = vector()
        vector.set(result, -self.xv, -self.yv, -self.zv)
        return result

    def __sub__(vector self, vector other):
        result = vector()
        vector.set(result,
                   self.xv-other.xv, self.yv-other.yv, self.zv-other.zv)
        return result

    def __mul__(x, y):
        cdef vector v1, v2
        from Scientific import Geometry
        if isinstance(y, vector):
            if isinstance(x, vector):
                v1 = x
                v2 = y
                return v1.xv*v2.xv+v1.yv*v2.yv+v1.zv*v2.zv
            else:
                x, y = y, x
        if Geometry.isTensor(y):
            product = Geometry.Tensor(x.array).dot(y)
            if product.rank == 1:
                result = vector()
                vector.set(result, product.array[0],
                           product.array[1], product.array[2])
                return result
            else:
                return product
        elif hasattr(y, "_product_with_vector"):
            return y._product_with_vector(self)
        else:
            v1 = x
            result = vector()
            vector.set(result, v1.xv*y, v1.yv*y, v1.zv*y)
            return result

    def __div__(vector self, double factor):
        result = vector()
        vector.set(result, self.xv/factor, self.yv/factor, self.zv/factor)
        return result

    def __richcmp__(vector self, vector other, int op):
        if op != 2 and op != 3:
            return NotImplemented
        eq = self.xv == other.xv and self.yv == other.yv \
             and self.zv == other.zv
        if op == 2:
            return eq
        else:
            return  not eq

    def __len__(self):
        return 3

    def __getitem__(self, int index):
        if index == 0:
            return self.xv
        elif index == 1:
            return self.yv
        elif index == 2:
            return self.zv
        raise IndexError

    def length(self):
        "Returns the length (norm)."
        return sqrt(self.xv*self.xv+self.yv*self.yv+self.zv*self.zv)

    def normal(self):
        "Returns a normalized copy."
        cdef double len
        len = sqrt(self.xv*self.xv+self.yv*self.yv+self.zv*self.zv)
        if len == 0:
            raise ZeroDivisionError, "Can't normalize a zero-length vector"
        result = vector()
        vector.set(result, self.xv/len, self.yv/len, self.zv/len)
        return result

    def cross(self, vector other):
        "Returns the cross product with vector |other|."
        result = vector()
        vector.set(result, self.yv*other.zv-self.zv*other.yv,
                           self.zv*other.xv-self.xv*other.zv,
                           self.xv*other.yv-self.yv*other.xv)
        return result

    def angle(self, vector other):
        "Returns the angle to vector |other|."
        cdef double cosa
        cosa = (self.xv*other.xv+self.yv*other.yv+self.zv*other.zv) / \
               sqrt((self.xv*self.xv+self.yv*self.yv+self.zv*self.zv)*
                    (other.xv*other.xv+other.yv*other.yv+other.zv*other.zv))
        if cosa > 1.:
            cosa = 1.
        if cosa < -1.:
            cosa = -1.
        return acos(cosa)

    def asTensor(self):
        "Returns an equivalent tensor object of rank 1."
        from Scientific import Geometry
        return Geometry.Tensor(self.array, 1)

    def dyadicProduct(self, other):
        "Returns the dyadic product with vector or tensor |other|."
        from Scientific import Geometry
        if isinstance(other, vector):
            return Geometry.Tensor(self.array, 1) * \
                   Geometry.Tensor(other.array, 1)
        elif Geometry.isTensor(other):
            return Geometry.Tensor(self.array, 1)*other
        else:
            raise TypeError, "Dyadic product with non-vector"
        

cdef class Vector(vector):

    property __safe_for_unpickling__:
        def __get__(self):
            return 1

    def __init__(self, x=None, y=None, z=None):
        if x is None:
            pass  # values are pre-initialized to zero
        elif y is None and z is None:
            self.xv, self.yv, self.zv = x
        else:
            self.xv = x
            self.yv = y
            self.zv = z

#
# For compatibility reasons, this routine works like its predecessor
# by testing the attribute is_vector. However, isinstance() works
# as well and is probably more efficient.
#
def isVector(x):
    try:
        return x.is_vector
    except AttributeError:
        return 0
