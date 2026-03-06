"""__init__.pyi - Type stubs for the netCDF4 Python package"""
# Notes:
#
# - The stubs in this file are manually-generated and must be updated if and when the API is changed.
# - The following **ruff** commands may be used to format this file according to
#   https://typing.readthedocs.io/en/latest/source/stubs.html
#
#   ruff format --line-length 130 src/netCDF4/__init__.pyi  # format code
#   ruff check --line-length 130 --select I --fix src/netCDF4/__init__.pyi  # sort imports
#
# - The Variable class is a generic and may thus be statically typed, but this has limited utility for the following reasons:
#   - The return type of `Variable.__getitem__()` (and `Variable.getValue()`) depends on a number of factors (e.g. variable
#     shape, key shape, whether masking is enabled) that cannot be easily determined statically.
#   - Similarly, the types and shapes of data that `Variable.__setitem__()` may accept varies widely depending on many factors
#     and is intractable to determine statically.
#   - Some facility for automatically typing a Variable on creation has been provided, however it is not exhaustive as a variable
#     may created with a string literal indicating its type and it would require an excessive number of overloads to enumerate
#     each of these cases.
#   - It is not possible to statically type a Variable of any user-defined type (CompoundType, EnumType, VLType) as these types
#     are created dynamically.
#   Thus it is most often best for the user to implement TypeGuards and/or perform other mixed static/runtime type-checking to
#   ensure the type and shape of data retrieved from this library.
# - The return type of some functions or properties (such as `Dataset.__getitem__()`) may one of a number of types. Where it is
#   not possible to narrow the type with overloads, the authors of these stubs have generally chosen to use `Any` as the return
#   type rather than a union of the possible types.
# - `MFDataset.dimensions` returns `dict[str, Dimension]` and `MFDataset.variables` returns `dict[str, Variable]` even though the
#   dict value types may actually be `_Dimension` and `_Variable`, respectively. The original authors of this stubfile have
#   elected to do this for simplicity's sake, but it may make sense to change this in the future, or just return `dict[str, Any]`.

import datetime as dt
import os
import sys
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    final,
    Final,
    Generic,
    Iterable,
    Iterator,
    Literal,
    Mapping,
    NoReturn,
    Sequence,
    TypedDict,
    TypeVar,
    Union,
    overload,
)

import cftime
import numpy as np
import numpy.typing as npt
from typing_extensions import Buffer, Self, TypeAlias, disjoint_base

__all__ = [
    "Dataset",
    "Variable",
    "Dimension",
    "Group",
    "MFDataset",
    "MFTime",
    "CompoundType",
    "VLType",
    "date2num",
    "num2date",
    "rc_get",
    "rc_set",
    "date2index",
    "stringtochar",
    "chartostring",
    "stringtoarr",
    "getlibversion",
    "EnumType",
    "get_chunk_cache",
    "set_chunk_cache",
    "set_alignment",
    "get_alignment",
]
__pdoc__ = {"utils": False}

# string type specifiers
# fmt: off
RealTypeLiteral: TypeAlias = Literal[
    "i1", "b", "B", "int8",        # NC_BYTE
    "u1", "uint8",                 # NC_UBYTE
    "i2", "h", "s", "int16",       # NC_SHORT
    "u2", "uint16",                # NC_USHORT
    "i4", "i", "l", "int32",       # NC_INT
    "u4", "uint32",                # NC_UINT
    "i8", "int64", "int",          # NC_INT64
    "u8", "uint64",                # NC_UINT64
    "f4", "f", "float32",          # NC_FLOAT
    "f8", "d", "float64", "float"  # NC_DOUBLE
]
# fmt: on
ComplexTypeLiteral: TypeAlias = Literal["c8", "c16", "complex64", "complex128"]
NumericTypeLiteral: TypeAlias = RealTypeLiteral | ComplexTypeLiteral
CharTypeLiteral: TypeAlias = Literal["S1", "c"]  # NC_CHAR
TypeLiteral: TypeAlias = NumericTypeLiteral | CharTypeLiteral

# Numpy types
NumPyRealType: TypeAlias = (
    np.int8 | np.uint8 | np.int16 | np.uint16 | np.int32 | np.uint32 | np.int64 | np.uint64 | np.float16 | np.float32 | np.float64
)
NumPyComplexType: TypeAlias = np.complex64 | np.complex128
NumPyNumericType: TypeAlias = NumPyRealType | NumPyComplexType
# Classes that can create instances of NetCDF user-defined types
NetCDFUDTClass: TypeAlias = CompoundType | VLType | EnumType
# Possible argument types for the datatype argument used in Variable creation. At this time, it is not possible to allow unknown
# strings arguments in the datatype field but exclude and string literals that are not one of `TypeLiteral`, so really
# `TypeLiteral` is made irrelevant, except for anyone who looks at this file.
DatatypeType: TypeAlias = (
    TypeLiteral
    | str
    | np.dtype[NumPyNumericType | np.str_]
    | type[int | float | NumPyNumericType | str | np.str_]
    | NetCDFUDTClass
)

VarT = TypeVar("VarT")
RealVarT = TypeVar("RealVarT", bound=NumPyRealType)
ComplexVarT = TypeVar("ComplexVarT", bound=NumPyComplexType)
NumericVarT = TypeVar("NumericVarT", bound=NumPyNumericType)

DimensionsType: TypeAlias = Union[str, bytes, Dimension, Iterable[Union[str, bytes, Dimension]]]
CompressionType: TypeAlias = Literal[
    "zlib", "szip", "zstd", "bzip2", "blosc_lz", "blosc_lz4", "blosc_lz4hc", "blosc_zlib", "blosc_zstd"
]
CompressionLevel: TypeAlias = Literal[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
AccessMode: TypeAlias = Literal["r", "w", "r+", "a", "x", "rs", "ws", "r+s", "as"]
Format: TypeAlias = Literal["NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_CLASSIC", "NETCDF3_64BIT_OFFSET", "NETCDF3_64BIT_DATA"]
DiskFormat: TypeAlias = Literal["NETCDF3", "HDF5", "HDF4", "PNETCDF", "DAP2", "DAP4", "UNDEFINED"]
QuantizeMode: TypeAlias = Literal["BitGroom", "BitRound", "GranularBitRound"]
EndianType: TypeAlias = Literal["native", "little", "big"]
CalendarType: TypeAlias = Literal[
    "standard", "gregorian", "proleptic_gregorian", "noleap", "365_day", "360_day", "julian", "all_leap", "366_day"
]
BoolInt: TypeAlias = Literal[0, 1]

DateTimeArray: TypeAlias = npt.NDArray[np.object_]
"""numpy array of datetime.datetime or cftime.datetime"""

class BloscInfo(TypedDict):
    compressor: Literal["blosc_lz", "blosc_lz4", "blosc_lz4hc", "blosc_zlib", "blosc_zstd"]
    shuffle: Literal[0, 1, 2]

class SzipInfo(TypedDict):
    coding: Literal["nn", "ec"]
    pixels_per_block: Literal[4, 8, 16, 32]

class FiltersDict(TypedDict):
    """Dict returned from netCDF4.Variable.filters()"""

    zlib: bool
    szip: Literal[False] | SzipInfo
    zstd: bool
    bzip2: bool
    blosc: Literal[False] | BloscInfo
    shuffle: bool
    complevel: int
    fletcher32: bool

__version__: str
__netcdf4libversion__: str
__hdf5libversion__: str
__has_rename_grp__: BoolInt
__has_nc_inq_path__: BoolInt
__has_nc_inq_format_extended__: BoolInt
__has_nc_open_mem__: BoolInt
__has_nc_create_mem__: BoolInt
__has_cdf5_format__: BoolInt
__has_parallel4_support__: BoolInt
__has_pnetcdf_support__: BoolInt
__has_parallel_support__: BoolInt
__has_quantization_support__: BoolInt
__has_zstandard_support__: BoolInt
__has_bzip2_support__: BoolInt
__has_blosc_support__: BoolInt
__has_szip_support__: BoolInt
__has_set_alignment__: BoolInt
__has_ncfilter__: BoolInt
__has_nc_rc_set__: BoolInt
is_native_little: bool
is_native_big: bool
default_encoding: Final = "utf-8"
unicode_error: Final = "replace"
default_fillvals: dict[str, int | float | str]

# date2index, date2num, and num2date are actually provided by cftime and if stubs for
# cftime are completed these should be removed.
def date2index(
    dates: dt.datetime | cftime.datetime | Sequence[dt.datetime | cftime.datetime] | DateTimeArray,
    nctime: Variable,
    calendar: CalendarType | str | None = None,
    select: Literal["exact", "before", "after", "nearest"] = "exact",
    has_year_zero: bool | None = None,
) -> int | npt.NDArray[np.int_]: ...
def date2num(
    dates: dt.datetime | cftime.datetime | Sequence[dt.datetime | cftime.datetime] | DateTimeArray,
    units: str,
    calendar: CalendarType | str | None = None,
    has_year_zero: bool | None = None,
    longdouble: bool = False,
) -> np.number | npt.NDArray[np.number]: ...
def num2date(
    times: Sequence[int | float | np.number] | npt.NDArray[np.number],
    units: str,
    calendar: CalendarType | str = "standard",
    only_use_cftime_datetimes: bool = True,
    only_use_python_datetimes: bool = False,
    has_year_zero: bool | None = None,
) -> dt.datetime | DateTimeArray: ...

class NetCDF4MissingFeatureException(Exception):
    def __init__(self, feature: str, version: str): ...

def dtype_is_complex(dtype: str) -> bool: ...

@disjoint_base
class Dataset:
    def __init__(
        self,
        filename: str | os.PathLike,
        mode: AccessMode = "r",
        clobber: bool = True,
        format: Format = "NETCDF4",
        diskless: bool = False,
        persist: bool = False,
        keepweakref: bool = False,
        memory: Buffer | int | None = None,
        encoding: str | None = None,
        parallel: bool = False,
        comm: Any = None,
        info: Any = None,
        auto_complex: bool = False,
        **kwargs: Any,
    ): ...
    @property
    def name(self) -> str: ...
    @property
    def groups(self) -> dict[str, Group]: ...
    @property
    def dimensions(self) -> dict[str, Dimension]: ...
    @property
    def variables(self) -> dict[str, Variable[Any]]: ...
    @property
    def cmptypes(self) -> dict[str, CompoundType]: ...
    @property
    def vltypes(self) -> dict[str, VLType]: ...
    @property
    def enumtypes(self) -> dict[str, EnumType]: ...
    @property
    def data_model(self) -> Format: ...
    @property
    def file_format(self) -> Format: ...
    @property
    def disk_format(self) -> DiskFormat: ...
    @property
    def parent(self) -> Dataset | None: ...
    @property
    def path(self) -> str: ...
    @property
    def keepweakref(self) -> bool: ...
    @property
    def auto_complex(self) -> bool: ...
    @property
    def _ncstring_attrs__(self) -> bool: ...
    @property
    def __orthogonal_indexing__(self) -> bool: ...
    def filepath(self, encoding: str | None = None) -> str: ...
    def isopen(self) -> bool: ...
    def close(self) -> memoryview: ...  # only if writing and memory != None, but otherwise people ignore the return None anyway
    def sync(self) -> None: ...
    def set_fill_on(self) -> None: ...
    def set_fill_off(self) -> None: ...
    def createDimension(self, dimname: str, size: int | None = None) -> Dimension: ...
    def renameDimension(self, oldname: str, newname: str) -> None: ...
    @overload
    def createVariable(
        self,
        varname: str,
        datatype: np.dtype[NumericVarT] | type[NumericVarT],
        dimensions: DimensionsType = (),
        compression: CompressionType | None = None,
        zlib: bool = False,
        complevel: CompressionLevel | None = 4,
        shuffle: bool = True,
        szip_coding: Literal["nn", "ec"] = "nn",
        szip_pixels_per_block: Literal[4, 8, 16, 32] = 8,
        blosc_shuffle: Literal[0, 1, 2] = 1,
        fletcher32: bool = False,
        contiguous: bool = False,
        chunksizes: Sequence[int] | None = None,
        endian: EndianType = "native",
        least_significant_digit: int | None = None,
        significant_digits: int | None = None,
        quantize_mode: QuantizeMode = "BitGroom",
        fill_value: int | float | np.generic | str | bytes | Literal[False] | np.ndarray | None = None,
        chunk_cache: int | None = None,
    ) -> Variable[NumericVarT]: ...
    @overload
    def createVariable(
        self,
        varname: str,
        datatype: np.dtype[np.str_] | type[str | np.str_],
        dimensions: DimensionsType = (),
        compression: CompressionType | None = None,
        zlib: bool = False,
        complevel: CompressionLevel | None = 4,
        shuffle: bool = True,
        szip_coding: Literal["nn", "ec"] = "nn",
        szip_pixels_per_block: Literal[4, 8, 16, 32] = 8,
        blosc_shuffle: Literal[0, 1, 2] = 1,
        fletcher32: bool = False,
        contiguous: bool = False,
        chunksizes: Sequence[int] | None = None,
        endian: EndianType = "native",
        least_significant_digit: int | None = None,
        significant_digits: int | None = None,
        quantize_mode: QuantizeMode = "BitGroom",
        fill_value: int | float | np.generic | str | bytes | Literal[False] | np.ndarray | None = None,
        chunk_cache: int | None = None,
    ) -> Variable[str]: ...
    @overload
    def createVariable(
        self,
        varname: str,
        datatype: DatatypeType,
        dimensions: DimensionsType = (),
        compression: CompressionType | None = None,
        zlib: bool = False,
        complevel: CompressionLevel | None = 4,
        shuffle: bool = True,
        szip_coding: Literal["nn", "ec"] = "nn",
        szip_pixels_per_block: Literal[4, 8, 16, 32] = 8,
        blosc_shuffle: Literal[0, 1, 2] = 1,
        fletcher32: bool = False,
        contiguous: bool = False,
        chunksizes: Sequence[int] | None = None,
        endian: EndianType = "native",
        least_significant_digit: int | None = None,
        significant_digits: int | None = None,
        quantize_mode: QuantizeMode = "BitGroom",
        fill_value: int | float | np.generic | str | bytes | Literal[False] | np.ndarray | None = None,
        chunk_cache: int | None = None,
    ) -> Variable: ...
    def renameVariable(self, oldname: str, newname: str) -> None: ...
    def createGroup(self, groupname: str) -> Group: ...
    def renameGroup(self, oldname: str, newname: str) -> None: ...
    def renameAttribute(self, oldname: str, newname: str) -> None: ...
    def createCompoundType(
        self, datatype: npt.DTypeLike | Sequence[tuple[str, npt.DTypeLike]], datatype_name: str
    ) -> CompoundType: ...
    def createVLType(self, datatype: npt.DTypeLike, datatype_name: str) -> VLType: ...
    def createEnumType(
        self,
        datatype: np.dtype[np.integer] | type[np.integer] | type[int] | str,
        datatype_name: str,
        enum_dict: Mapping[str, int | np.integer],
    ) -> EnumType: ...
    def ncattrs(self) -> list[str]: ...
    def setncattr_string(self, name: str, value: Any) -> None: ...
    def setncattr(self, name: str, value: Any) -> None: ...
    def setncatts(self, attdict: Mapping[str, Any]) -> None: ...
    def getncattr(self, name: str, encoding: str = "utf-8") -> Any: ...
    def delncattr(self, name: str) -> None: ...
    def set_auto_chartostring(self, value: bool) -> None: ...
    def set_auto_maskandscale(self, value: bool) -> None: ...
    def set_auto_mask(self, value: bool) -> None: ...
    def set_auto_scale(self, value: bool) -> None: ...
    def set_always_mask(self, value: bool) -> None: ...
    def set_ncstring_attrs(self, value: bool) -> None: ...
    def get_variables_by_attributes(self, **kwargs: Callable[[Any], bool] | Any) -> list[Variable]: ...
    @staticmethod
    def fromcdl(
        cdlfilename: str | os.PathLike, ncfilename: str | os.PathLike | None = None, mode: AccessMode = "a", format: Format = "NETCDF4"
    ) -> Dataset: ...
    @overload
    def tocdl(self, coordvars: bool = False, data: bool = False, outfile: None = None) -> str: ...
    @overload
    def tocdl(self, coordvars: bool = False, data: bool = False, *, outfile: str | os.PathLike) -> None: ...
    def has_blosc_filter(self) -> bool: ...
    def has_zstd_filter(self) -> bool: ...
    def has_bzip2_filter(self) -> bool: ...
    def has_szip_filter(self) -> bool: ...
    def __getitem__(self, elem: str) -> Any: ...  # should be Group | Variable, but this causes too many problems
    # __iter__ and __contains__ always error because iteration and membership ops are not allowed
    def __iter__(self) -> NoReturn: ...
    def __contains__(self, key) -> NoReturn: ...
    def __setattr__(self, name: str, value: Any) -> None: ...
    def __getattr__(self, name: str) -> Any: ...
    def __delattr__(self, name: str): ...
    def __reduce__(self) -> NoReturn: ...
    def __enter__(self) -> Self: ...
    def __exit__(self, atype, value, traceback) -> None: ...

class Group(Dataset):
    def __init__(self, parent: Dataset, name: str, **kwargs: Any) -> None: ...
    def close(self) -> NoReturn: ...

@final
class Dimension:
    def __init__(self, grp: Dataset, name: str, size: int | None = None, **kwargs: Any) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def size(self) -> int: ...
    def group(self) -> Dataset: ...
    def isunlimited(self) -> bool: ...
    def __len__(self) -> int: ...

class _VarDatatypeProperty:
    # A faux descriptor definition of the property to allow overloads
    @overload
    def __get__(self, instance: Variable[RealVarT], owner: Any) -> np.dtype[RealVarT]: ...
    @overload
    def __get__(self, instance: Variable[ComplexVarT], owner: Any) -> CompoundType: ...
    @overload
    def __get__(self, instance: Variable[str], owner: Any) -> VLType: ...
    @overload
    def __get__(
        self, instance: Variable, owner: Any
    ) -> Any: ...  # actual return type np.dtype | CompoundType | VLType | EnumType

class _VarDtypeProperty:
    # A faux descriptor definition of the property to allow overloads
    @overload
    def __get__(self, instance: Variable[NumericVarT], owner: Any) -> np.dtype[NumericVarT]: ...
    @overload
    def __get__(self, instance: Variable[str], owner: Any) -> type[str]: ...
    @overload
    def __get__(self, instance: Variable, owner: Any) -> Any: ...  # actual return type np.dtype | Type[str]

@final
class Variable(Generic[VarT]):
    # Overloads of __new__ are provided for some cases where the Variable's type may be statically inferred from the datatype arg
    @overload
    def __new__(
        cls,
        grp: Dataset,
        name: str,
        datatype: np.dtype[NumericVarT] | type[NumericVarT],
        dimensions: DimensionsType = (),
        compression: CompressionType | None = None,
        zlib: bool = False,
        complevel: CompressionLevel | None = 4,
        shuffle: bool = True,
        szip_coding: Literal["nn", "ec"] = "nn",
        szip_pixels_per_block: Literal[4, 8, 16, 32] = 8,
        blosc_shuffle: Literal[0, 1, 2] = 1,
        fletcher32: bool = False,
        contiguous: bool = False,
        chunksizes: Sequence[int] | None = None,
        endian: EndianType = "native",
        least_significant_digit: int | None = None,
        significant_digits: int | None = None,
        quantize_mode: QuantizeMode = "BitGroom",
        fill_value: int | float | np.generic | str | bytes | Literal[False] | np.ndarray | None = None,
        chunk_cache: int | None = None,
        **kwargs: Any,
    ) -> Variable[NumericVarT]: ...
    @overload
    def __new__(
        cls,
        grp: Dataset,
        name: str,
        datatype: np.dtype[np.str_] | type[str | np.str_],
        dimensions: DimensionsType = (),
        compression: CompressionType | None = None,
        zlib: bool = False,
        complevel: CompressionLevel | None = 4,
        shuffle: bool = True,
        szip_coding: Literal["nn", "ec"] = "nn",
        szip_pixels_per_block: Literal[4, 8, 16, 32] = 8,
        blosc_shuffle: Literal[0, 1, 2] = 1,
        fletcher32: bool = False,
        contiguous: bool = False,
        chunksizes: Sequence[int] | None = None,
        endian: EndianType = "native",
        least_significant_digit: int | None = None,
        significant_digits: int | None = None,
        quantize_mode: QuantizeMode = "BitGroom",
        fill_value: int | float | np.generic | str | bytes | Literal[False] | np.ndarray | None = None,
        chunk_cache: int | None = None,
        **kwargs: Any,
    ) -> Variable[str]: ...
    @overload
    def __new__(
        cls,
        grp: Dataset,
        name: str,
        datatype: DatatypeType,
        dimensions: DimensionsType = (),
        compression: CompressionType | None = None,
        zlib: bool = False,
        complevel: CompressionLevel | None = 4,
        shuffle: bool = True,
        szip_coding: Literal["nn", "ec"] = "nn",
        szip_pixels_per_block: Literal[4, 8, 16, 32] = 8,
        blosc_shuffle: Literal[0, 1, 2] = 1,
        fletcher32: bool = False,
        contiguous: bool = False,
        chunksizes: Sequence[int] | None = None,
        endian: EndianType = "native",
        least_significant_digit: int | None = None,
        significant_digits: int | None = None,
        quantize_mode: QuantizeMode = "BitGroom",
        fill_value: int | float | np.generic | str | bytes | Literal[False] | np.ndarray | None = None,
        chunk_cache: int | None = None,
        **kwargs: Any,
    ) -> Variable: ...
    def __init__(
        self,
        grp: Dataset,
        name: str,
        datatype: DatatypeType,
        dimensions: DimensionsType = (),
        compression: CompressionType | None = None,
        zlib: bool = False,
        complevel: CompressionLevel | None = 4,
        shuffle: bool = True,
        szip_coding: Literal["nn", "ec"] = "nn",
        szip_pixels_per_block: Literal[4, 8, 16, 32] = 8,
        blosc_shuffle: Literal[0, 1, 2] = 1,
        fletcher32: bool = False,
        contiguous: bool = False,
        chunksizes: Sequence[int] | None = None,
        endian: EndianType = "native",
        least_significant_digit: int | None = None,
        significant_digits: int | None = None,
        quantize_mode: QuantizeMode = "BitGroom",
        fill_value: int | float | np.generic | str | bytes | Literal[False] | np.ndarray | None = None,
        chunk_cache: int | None = None,
        **kwargs: Any,
    ) -> None: ...
    datatype: _VarDatatypeProperty
    dtype: _VarDtypeProperty
    @property
    def name(self) -> str: ...
    @property
    def shape(self) -> tuple[int, ...]: ...
    @property
    def size(self) -> int: ...
    @property
    def dimensions(self) -> tuple[str, ...]: ...
    @property
    def ndim(self) -> int: ...
    @property
    def scale(self) -> bool: ...
    @property
    def mask(self) -> bool: ...
    @property
    def chartostring(self) -> bool: ...
    @property
    def always_mask(self) -> bool: ...
    @property
    def __orthogonal_indexing__(self) -> bool: ...
    def group(self) -> Dataset: ...
    def ncattrs(self) -> list[str]: ...
    def setncattr(self, name: str, value: Any) -> None: ...
    def setncattr_string(self, name: str, value: Any) -> None: ...
    def setncatts(self, attdict: Mapping[str, Any]) -> None: ...
    def getncattr(self, name: str, encoding="utf-8"): ...
    def delncattr(self, name: str) -> None: ...
    def filters(self) -> FiltersDict: ...
    def quantization(self) -> tuple[int, QuantizeMode] | None: ...
    def endian(self) -> EndianType: ...
    def chunking(self) -> Literal["contiguous"] | list[int]: ...
    def get_var_chunk_cache(self) -> tuple[int, int, float]: ...
    def set_var_chunk_cache(
        self, size: int | None = None, nelems: int | None = None, preemption: float | None = None
    ) -> None: ...
    def renameAttribute(self, oldname: str, newname: str) -> None: ...
    def assignValue(self, val: Any) -> None: ...
    def getValue(self) -> Any: ...
    def get_fill_value(self) -> Any: ...
    def set_auto_chartostring(self, chartostring: bool) -> None: ...
    def use_nc_get_vars(self, use_nc_get_vars: bool) -> None: ...
    def set_auto_maskandscale(self, maskandscale: bool) -> None: ...
    def set_auto_scale(self, scale: bool) -> None: ...
    def set_auto_mask(self, mask: bool) -> None: ...
    def set_always_mask(self, always_mask: bool) -> None: ...
    def set_ncstring_attrs(self, ncstring_attrs: bool) -> None: ...
    def set_collective(self, value: bool) -> None: ...
    def get_dims(self) -> tuple[Dimension, ...]: ...
    def __delattr__(self, name: str) -> None: ...
    def __setattr__(self, name: str, value: Any) -> None: ...
    def __getattr__(self, name: str) -> Any: ...
    def __getitem__(self, elem: Any) -> Any: ...
    def __setitem__(self, elem: Any, data: npt.ArrayLike) -> None: ...
    def __array__(self) -> np.ndarray: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[Any]: ...  # faux method so mypy believes Variable is iterable

@final
class CompoundType:
    dtype: np.dtype
    dtype_view: np.dtype
    name: str

    def __init__(
        self, grp: Dataset, dt: npt.DTypeLike | Sequence[tuple[str, npt.DTypeLike]], dtype_name: str, **kwargs: Any
    ) -> None: ...
    def __reduce__(self) -> NoReturn: ...

@final
class VLType:
    dtype: np.dtype
    name: str | None

    def __init__(self, grp: Dataset, dt: npt.DTypeLike, dtype_name: str, **kwargs: Any) -> None: ...
    def __reduce__(self) -> NoReturn: ...

@final
class EnumType:
    dtype: np.dtype[np.integer]
    name: str
    enum_dict: Mapping[str, int]

    def __init__(
        self,
        grp: Dataset,
        dt: np.dtype[np.integer] | type[np.integer] | type[int] | str,
        dtype_name: str,
        enum_dict: Mapping[str, int | np.integer],
        **kwargs: Any,
    ) -> None: ...
    def __reduce__(self) -> NoReturn: ...

class MFDataset(Dataset):
    def __init__(
        self,
        files: str | Sequence[str | os.PathLike],
        check: bool = False,
        aggdim: str | None = None,
        exclude: Sequence[str] = [],
        master_file: str | os.PathLike | None = None,
    ) -> None: ...
    @property
    def dimensions(self) -> dict[str, Dimension]: ...  # this should be: dict[str, Dimension | _Dimension]
    @property
    def variables(self) -> dict[str, Variable[Any]]: ...  # this should be: dict[str, _Variable[Any] | _Variable]

class _Dimension:
    dimlens: list[int]
    dimtolen: int

    def __init__(self, dimname: str, dim: Dimension, dimlens: list[int], dimtotlen: int) -> None: ...
    def __len__(self) -> int: ...
    def isunlimited(self) -> Literal[True]: ...

class _Variable:
    dimensions: tuple[str, ...]
    dtype: np.dtype | type[str]

    def __init__(self, dset: Dataset, varname: str, var: Variable[Any], recdimname: str) -> None: ...

    # shape, ndim, and name actually come from __getattr__
    @property
    def shape(self) -> tuple[int, ...]: ...
    @property
    def ndim(self) -> int: ...
    @property
    def name(self) -> str: ...
    def typecode(self) -> np.dtype | type[str]: ...
    def ncattrs(self) -> list[str]: ...
    def _shape(self) -> tuple[int, ...]: ...
    def set_auto_chartostring(self, val: bool) -> None: ...
    def set_auto_maskandscale(self, val: bool) -> None: ...
    def set_auto_mask(self, val: bool) -> None: ...
    def set_auto_scale(self, val: bool) -> None: ...
    def set_always_mask(self, val: bool) -> None: ...
    def __getattr__(self, name: str) -> Any: ...
    def __getitem__(self, elem: Any) -> Any: ...
    def __len__(self) -> int: ...

class MFTime(_Variable):
    calendar: CalendarType | None
    units: str | None

    def __init__(self, time: Variable, units: str | None = None, calendar: CalendarType | str | None = None): ...
    def __getitem__(self, elem: Any) -> np.ndarray: ...

@overload
def stringtoarr(
    string: str,
    NUMCHARS: int,
    dtype: Literal["S"] | np.dtype[np.bytes_] = "S",
) -> npt.NDArray[np.bytes_]: ...
@overload
def stringtoarr(
    string: str,
    NUMCHARS: int,
    dtype: Literal["U"] | np.dtype[np.str_],
) -> npt.NDArray[np.str_]: ...
@overload
def stringtochar(
    a: npt.NDArray[np.character],
    encoding: Literal["none", "None", "bytes"],
    n_strlen: int | None = None,
) -> npt.NDArray[np.bytes_]: ...
@overload
def stringtochar(
    a: npt.NDArray[np.character],
    encoding: str = ...,
) -> npt.NDArray[np.str_] | npt.NDArray[np.bytes_]: ...
@overload
def chartostring(
    b: npt.NDArray[np.character],
    encoding: Literal["none", "None", "bytes"] = ...,
) -> npt.NDArray[np.bytes_]: ...
@overload
def chartostring(
    b: npt.NDArray[np.character],
    encoding: str = ...,
) -> npt.NDArray[np.str_] | npt.NDArray[np.bytes_]: ...
def getlibversion() -> str: ...
def rc_get(key: str) -> str | None: ...
def rc_set(key: str, value: str) -> None: ...
def set_alignment(threshold: int, alignment: int): ...
def get_alignment() -> tuple[int, int]: ...
def set_chunk_cache(size: int | None = None, nelems: int | None = None, preemption: float | None = None) -> None: ...
def get_chunk_cache() -> tuple[int, int, float]: ...
