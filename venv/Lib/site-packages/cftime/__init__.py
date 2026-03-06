from ._cftime import (datetime, real_datetime,
                      _parse_date, _dateparse, _datesplit, is_leap_year)
from ._cftime import num2date, date2num, date2index, time2index, num2pydate, to_tuple
from ._cftime import (microsec_units, millisec_units,
                     sec_units, hr_units, day_units, min_units,
                     UNIT_CONVERSION_FACTORS)
from ._cftime import CFWarning
# these will be removed in a future release
from ._cftime import (DatetimeNoLeap, DatetimeAllLeap, Datetime360Day,
                     Datetime360Day, DatetimeJulian, 
                     DatetimeGregorian, DatetimeProlepticGregorian)


def __getattr__(item: str):
    if item == "__version__":
        from importlib.metadata import version
        return version("cftime")
    
    raise AttributeError(f"module 'cftime' has no attribute {item!r}")
