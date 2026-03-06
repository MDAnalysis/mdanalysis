"""
Performs conversions of netCDF time coordinate data to/from datetime objects.
"""

from cpython.object cimport (PyObject_RichCompare, Py_LT, Py_LE, Py_EQ,
                             Py_NE, Py_GT, Py_GE)
from numpy cimport int64_t, int32_t
import cython
import numpy as np
cimport numpy as np
import re
import time
from datetime import datetime as datetime_python
from datetime import timedelta, MINYEAR, MAXYEAR
import warnings
from ._strptime import _strptime

np.import_array()

microsec_units = ['microseconds','microsecond', 'microsec', 'microsecs']
millisec_units = ['milliseconds', 'millisecond', 'millisec', 'millisecs', 'msec', 'msecs', 'ms']
sec_units =      ['second', 'seconds', 'sec', 'secs', 's']
min_units =      ['minute', 'minutes', 'min', 'mins']
hr_units =       ['hour', 'hours', 'hr', 'hrs', 'h']
day_units =      ['day', 'days', 'd']
month_units =    ['month', 'months'] # only allowed for 360_day calendar
year_units =     ['common_year', 'common_years'] # only allowed for 365_day and noleap calendars
_units = microsec_units+millisec_units+sec_units+min_units+hr_units+day_units
# supported calendars. Includes synonyms ('standard'=='gregorian',
# '366_day'=='all_leap','365_day'=='noleap')
# see http://cfconventions.org/cf-conventions/cf-conventions#calendar
# for definitions.
_calendars = ['standard', 'gregorian', 'proleptic_gregorian',
              'noleap', 'julian', 'all_leap', '365_day', '366_day', '360_day']
_idealized_calendars= ['all_leap','noleap','366_day','365_day','360_day']
# Following are number of days per month
cdef int[12] _dayspermonth      = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
cdef int[12] _dayspermonth_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
# same as above, but including accumulated days of previous months.
cdef int[13] _cumdayspermonth = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
cdef int[13] _cumdayspermonth_leap = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]

# Adapted from http://delete.me.uk/2005/03/iso8601.html
# Note: This regex ensures that all ISO8601 timezone formats are accepted - but, due to legacy support for other timestrings, not all incorrect formats can be rejected.
#       For example, the TZ spec "+01:0" will still work even though the minutes value is only one character long.
ISO8601_REGEX = re.compile(r"(?P<year>[+-]?[0-9]+)(-(?P<month>[0-9]{1,2})(-(?P<day>[0-9]{1,2})"
                           r"(((?P<separator1>.)(?P<hour>[0-9]{1,2}):(?P<minute>[0-9]{1,2})(:(?P<second>[0-9]{1,2})(\.(?P<fraction>[0-9]+))?)?)?"
                           r"((?P<separator2>.?)(?P<timezone>Z|(([-+])([0-9]{2})((:([0-9]{2}))|([0-9]{2}))?)))?)?)?)?"
                           )
# Note: The re module apparently does not support branch reset groups that allow redefinition of the same group name in alternative branches as PCRE does.
#       Using two different group names is also somewhat ugly, but other solutions might hugely inflate the expression. feel free to contribute a better solution.
TIMEZONE_REGEX = re.compile(
    "(?P<prefix>[+-])(?P<hours>[0-9]{2})(?:(?::(?P<minutes1>[0-9]{2}))|(?P<minutes2>[0-9]{2}))?")

class real_datetime(datetime_python):
    """add dayofwk, dayofyr, daysinmonth attributes to python datetime instance"""
    @property
    def dayofwk(self):
        # 0=Monday, 6=Sunday
        return self.weekday()
    @property
    def dayofyr(self):
        return self.timetuple().tm_yday
    @property
    def daysinmonth(self):
        if _is_leap(self.year,'proleptic_gregorian'):
            return _dayspermonth_leap[self.month-1]
        else:
            return _dayspermonth[self.month-1]
    nanosecond = 0 # workaround for pandas bug (cftime issue #77)

def _datesplit(timestr):
    """split a time string into two components, units and the remainder
    after 'since'
    """
    try:
        (units, sincestring, remainder) = timestr.split(None,2)
    except ValueError as e:
        raise ValueError('Incorrectly formatted CF date-time unit_string')

    if sincestring.lower() != 'since':
        raise ValueError("no 'since' in unit_string")

    return units.lower(), remainder

def _dateparse(timestr,calendar,has_year_zero=None):
    """parse a string of the form time-units since yyyy-mm-dd hh:mm:ss,
    return a datetime instance"""
    # same as version in cftime, but returns a timezone naive
    # python datetime instance with the utc_offset included.

    calendar = calendar.lower()
    # set calendar-specific defaults for has_year_zero
    if has_year_zero is None:
        has_year_zero = _year_zero_defaults(calendar)
    (units, isostring) = _datesplit(timestr)
    if not ((units in month_units and calendar=='360_day') or (units in year_units and calendar in {'365_day', 'noleap'}) or units in _units):
        if units in month_units and calendar != '360_day':
            raise ValueError("'months since' units only allowed for '360_day' calendar")
        if units in year_units and calendar not in {'365_day', 'noleap'}:    
            raise ValueError("'%s' units only allowed for '365_day' and 'noleap' calendars" % units) 
        else:
            raise ValueError(
            "In general, units must be one of 'microseconds', 'milliseconds', "
            "'seconds', 'minutes', 'hours', or 'days' (or select abbreviated "
            "versions of these).  For the '360_day' calendar, "
            "'months' can also be used, or for the 'noleap' calendar 'common_years' "
            "can also be used. Got '%s' instead, which are not recognized." % units)
    # parse the date string.
    year, month, day, hour, minute, second, microsecond, utc_offset =\
        _parse_date( isostring.strip() )
    if year == 0 and not has_year_zero and calendar in ['julian', 'standard', 'gregorian', 'proleptic_gregorian']:
        msg='zero not allowed as a reference year when has_year_zero=False'
        raise ValueError(msg)
    if calendar in ['noleap', '365_day'] and month == 2 and day == 29:
        raise ValueError(
            'cannot specify a leap day as the reference time with the noleap calendar')
    if calendar == '360_day' and day > 30:
        raise ValueError(
            'there are only 30 days in every month with the 360_day calendar')
    basedate = datetime(year, month, day, hour, minute, second,
                        microsecond,calendar=calendar,has_year_zero=has_year_zero)
    # subtract utc_offset from basedate time instance (which is timezone naive)
    if utc_offset:
        basedate -= timedelta(days=utc_offset/1440.)
    return basedate

def _can_use_python_datetime(date,calendar):
    #gregorian = datetime(1582,10,15,calendar=calendar,has_year_zero=date.has_year_zero)
    #return ((calendar == 'proleptic_gregorian' and date.year >= MINYEAR and date.year <= MAXYEAR) or \
    #       (calendar in ['gregorian','standard'] and date > gregorian and date.year <= MAXYEAR))
    return  (calendar == 'proleptic_gregorian' and date.year >= MINYEAR and date.year <= MAXYEAR) or \
            ((calendar in ['gregorian','standard'] and date.year <= MAXYEAR) and (date.year > 1582 or \
            (date.year == 1582 and date.month >= 10 and date.day > 15)))

@cython.embedsignature(True)
def date2num(dates, units, calendar=None, has_year_zero=None, longdouble=False):
    """
    Return numeric time values given datetime objects. The units
    of the numeric time values are described by the **units** argument
    and the **calendar** keyword. The datetime objects must
    be in UTC with no time-zone offset.  If there is a
    time-zone offset in **units**, it will be applied to the
    returned numeric values.

    **dates**: A datetime object or a sequence of datetime objects.
    The datetime objects should not include a time-zone offset. They
    can be either native python datetime instances (which use
    the proleptic gregorian calendar) or cftime.datetime instances.

    **units**: a string of the form **<time units> since <reference time>**
    describing the time units. **<time units>** can be days, hours, minutes,
    seconds, milliseconds or microseconds. **<reference time>** is the time
    origin. **months since** is allowed *only* for the **360_day** calendar
    and **common_years since** is allowed *only* for the **365_day** calendar.

    **calendar**: describes the calendar to be used in the time calculations.
    All the values currently defined in the
    `CF metadata convention <http://cfconventions.org/cf-conventions/cf-conventions#calendar>`__ are supported.
    Valid calendars **'standard', 'gregorian', 'proleptic_gregorian'
    'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day'**.
    Default is `None` which means the calendar associated with the first
    input datetime instance will be used.

    **has_year_zero**: If set to True, astronomical year numbering
    is used and the year zero exists.  If set to False for real-world
    calendars, then historical year numbering is used and the year 1 is
    preceded by year -1 and no year zero exists.
    The defaults are set to conform with
    CF version 1.9 conventions (False for 'julian', 'gregorian'/'standard', True
    for 'proleptic_gregorian' (ISO 8601) and True for the idealized
    calendars 'noleap'/'365_day', '360_day', 366_day'/'all_leap')
    Note that CF v1.9 does not specifically mention whether year zero
    is allowed in the proleptic_gregorian calendar, but ISO-8601 has
    a year zero so we have adopted this as the default.
    The defaults can only be over-ridden for the real-world calendars,
    for the the idealized calendars the year zero 
    always exists and the has_year_zero kwarg is ignored.
    This kwarg is not needed to define calendar systems allowed by CF
    (the calendar-specific defaults do this).

    **longdouble**: If set True, output is in the long double float type
    (numpy.float128) instead of float (numpy.float64), allowing microsecond
    accuracy when converting a time value to a date and back again. Otherwise
    this is only possible if the discretization of the time variable is an
    integer multiple of the units.

    returns a numeric time value, or an array of numeric time values
    with approximately 1 microsecond accuracy.
    """

    # input a scale or array-like?
    isscalar = False
    try:
        dates[0]
    except:
        if not dates:
            # if empty list or array input, return empty array (issue #315)
            return np.array([],dtype=float)
        else:
            isscalar = True

    # masked array input?
    ismasked = False
    if np.ma.isMA(dates) and np.ma.is_masked(dates):
        mask = dates.mask
        ismasked = True

    # are all input dates 'real' python datetime objects?
    dates = np.asanyarray(dates) # convert to numpy array
    shape = dates.shape # save shape of input
    all_python_datetimes = True
    for date in dates.flat:
        if not isinstance(date,datetime_python):
           all_python_datetimes = False
           break

    # if has_year_zero is None and calendar is None, use first input cftime.datetime instance
    # to determine if year zero is to be included.
    # If calendar is specified, use calendar specific defaults
    if has_year_zero is None:
        if calendar is None:
            d0 = dates.item(0)
            if isinstance(d0,datetime_python):
                has_year_zero = False
            else:
                try:
                    has_year_zero = d0.has_year_zero
                except AttributeError:
                    raise ValueError('has_year_zero not specified',type(d0))
        else:
            has_year_zero = _year_zero_defaults(calendar)

    # if calendar is None or '', use calendar of first input cftime.datetime instances.
    # if inputs are 'real' python datetime instances, use proleptic gregorian.
    if not calendar:
        if all_python_datetimes:
            calendar = 'proleptic_gregorian'
        else:
            d0 = dates.item(0)
            if isinstance(d0,datetime_python):
                calendar = 'proleptic_gregorian' 
            else:
                try:
                    calendar = d0.calendar
                except AttributeError:
                    raise ValueError('no calendar specified',type(d0))

    calendar = calendar.lower()
    basedate = _dateparse(units,calendar=calendar,has_year_zero=has_year_zero)
    (unit, isostring) = _datesplit(units)
    # real-world calendars cannot have zero as a reference year.
    if calendar in ['julian', 'standard', 'gregorian', 'proleptic_gregorian'] and not has_year_zero:
        if basedate.year == 0:
            msg='zero not allowed as a reference year, does not exist in Julian or Gregorian calendars'
            raise ValueError(msg)
    if unit not in UNIT_CONVERSION_FACTORS:
        raise ValueError("Unsupported time units provided, {!r}.".format(unit))
    if unit in ["months", "month"] and calendar != "360_day":
        raise ValueError("Units of months only valid for 360_day calendar.")
    unit_timedelta = timedelta(microseconds=UNIT_CONVERSION_FACTORS[unit])
    can_use_python_basedatetime = calendar=='proleptic_gregorian' and _can_use_python_datetime(basedate,calendar)

    if can_use_python_basedatetime and all_python_datetimes:
        use_python_datetime = True
        if not isinstance(basedate, datetime_python):
            basedate = real_datetime(basedate.year, basedate.month, basedate.day, 
                       basedate.hour, basedate.minute, basedate.second,
                       basedate.microsecond)
    else:
        use_python_datetime = False
        # convert basedate to specified calendar
        basedate =  to_calendar_specific_datetime(basedate, calendar, False, has_year_zero=has_year_zero)
    times = []
    for n, date in enumerate(dates.flat):
        if ismasked and mask.flat[n]:
            times.append(None)
            continue

        # use python datetime if possible.
        if use_python_datetime:
            # remove time zone offset
            if getattr(date, 'tzinfo',None) is not None:
                date = date.replace(tzinfo=None) - date.utcoffset()
        else: # convert date to same calendar specific cftime.datetime instance
            date = to_calendar_specific_datetime(date, calendar, False, has_year_zero=has_year_zero)

        td = date - basedate
        if td % unit_timedelta == timedelta(0):
            # Explicitly cast result to np.int64 for Windows compatibility
            quotient = np.int64(td // unit_timedelta)
            times.append(quotient)
        else:
            if longdouble:
                # Division of timedelta's is in float64 precision,
                # i.e. losing microsecond precision.
                # Conversion to float128 helps but can still lead to imprecision
                # of +-1 microsecond in division:
                #   quotient = (np.longdouble(td.total_seconds()) /
                #               np.longdouble(unit_timedelta.total_seconds()))
                # -> Convert to (64-bit) integers of microseconds
                mtd = (td.days * 86400000000 +
                       td.seconds * 1000000 +
                       td.microseconds)
                munit = (unit_timedelta.days * 86400000000 +
                         unit_timedelta.seconds * 1000000 +
                         unit_timedelta.microseconds)
                quotient = np.longdouble(mtd) / np.longdouble(munit)
            else:
                quotient = td / unit_timedelta
            times.append(quotient)

    if ismasked: # convert to masked array if input was masked array
        times = np.array(times, dtype=float)  # None -> nan
        times = np.ma.masked_invalid(times)
        if isscalar:
            return times[0]
        else:
            return np.reshape(times, shape)
    if isscalar:
        return times[0]
    else:
        return np.reshape(np.array(times), shape)


@cython.embedsignature(True)
def num2pydate(times,units,calendar='standard'):
    """
    Always returns python datetime.datetime
    objects and raise an error if this is not possible.

    Same as
    num2date(times,units,calendar,only_use_cftime_datetimes=False,only_use_python_datetimes=True)
    """
    return num2date(times,units,calendar,only_use_cftime_datetimes=False,only_use_python_datetimes=True)


UNIT_CONVERSION_FACTORS = {
    "microseconds": 1,
    "microsecond": 1,
    "microsec": 1,
    "microsecs": 1,
    "milliseconds": 1000,
    "millisecond": 1000,
    "millisec": 1000,
    "millisecs": 1000,
    "msec": 1000,
    "msecs": 1000,
    "ms": 1000,
    "seconds": 1000000,
    "second": 1000000,
    "sec": 1000000,
    "secs": 1000000,
    "s": 1000000,
    "minutes": 60 * 1000000,
    "minute": 60 * 1000000,
    "min": 60 * 1000000,
    "mins": 60 * 1000000,
    "hours": 3600 * 1000000,
    "hour": 3600 * 1000000,
    "hr": 3600 * 1000000,
    "hrs": 3600 * 1000000,
    "h": 3600 * 1000000,
    "day": 86400 * 1000000,
    "days": 86400 * 1000000,
    "d": 86400 * 1000000,
    "month": 30 * 86400 * 1000000,  # Only allowed for 360_day calendar
    "months": 30 * 86400 * 1000000,
    "common_year": 365 * 86400 * 1000000, # Only allowed for 365_day and no_leap calendars
    "common_years": 365 * 86400 * 1000000 # Only allowed for 365_day and no_leap calendars
    
}

DATE_TYPES = {
     "proleptic_gregorian": DatetimeProlepticGregorian,
     "standard": DatetimeGregorian,
     "noleap": DatetimeNoLeap,
     "365_day": DatetimeNoLeap,
     "all_leap": DatetimeAllLeap,
     "366_day": DatetimeAllLeap,
     "julian": DatetimeJulian,
     "360_day": Datetime360Day,
     "gregorian": DatetimeGregorian
 }

def to_calendar_specific_datetime(dt, calendar, use_python_datetime,has_year_zero=None):
    if use_python_datetime:
        return real_datetime(
               dt.year,
               dt.month,
               dt.day,
               dt.hour,
               dt.minute,
               dt.second,
               dt.microsecond)
    else:
        date_type = DATE_TYPES[calendar]
        return date_type(
               dt.year,
               dt.month,
               dt.day,
               dt.hour,
               dt.minute,
               dt.second,
               dt.microsecond,
               calendar=calendar,has_year_zero=has_year_zero)

_MAX_INT64 = np.iinfo("int64").max
_MIN_INT64 = np.iinfo("int64").min


def cast_to_int(num, units=None):
    if num.dtype.kind in "iu":
        return num
    else:
        #if np.any(num < _MIN_INT64) or np.any(num > _MAX_INT64):
        # use this instead to avoid test failures on windows with numpy 1.23.0 (issue #279)
        if np.any(np.less(num,_MIN_INT64,casting='same_kind')) or np.any(np.greater(num,_MAX_INT64,casting='same_kind')):
            raise OverflowError('time values outside range of 64 bit signed integers')
        if isinstance(num, np.ma.core.MaskedArray):
            int_num = np.ma.masked_array(np.rint(num), dtype=np.int64)
            # use ceil instead of rint if 1 microsec less than a second
            # or floor if 1 microsec greater than a second (issue #187)
            if units not in microsec_units and units not in millisec_units:
                int_num = np.ma.where(int_num%1000000 == 1, \
                          np.ma.masked_array(np.floor(num),dtype=np.int64), int_num)
                int_num = np.ma.where(int_num%1000000 == 999999, \
                          np.ma.masked_array(np.ceil(num),dtype=np.int64), int_num)
        else:
            int_num = np.array(np.rint(num), dtype=np.int64)
            if units not in microsec_units and units not in millisec_units:
                int_num = np.where(int_num%1000000 == 1, \
                          np.array(np.floor(num),dtype=np.int64), int_num)
                int_num = np.where(int_num%1000000 == 999999, \
                          np.array(np.ceil(num),dtype=np.int64), int_num)
        return int_num


def upcast_times(num):
    """Cast times array to larger float or integer dtype before scaling
    to units of microseconds.
    """
    if num.dtype.kind == "f":
        return num.astype(np.longdouble)
    else:
        return num.astype(np.int64)


def scale_times(num, factor):
    """Scale times by a factor, raise an error if values will not fit in np.int64"""
    if num.dtype.kind == "f":
        return factor * num
    else:
        if num.size == 0: # empty array (issue #287)
            return num
        # Python integers have arbitrary precision, so convert min and max
        # returned by NumPy functions through item, prior to multiplying by
        # factor.
        minimum = np.min(num).item() * factor
        maximum = np.max(num).item() * factor
        if minimum < _MIN_INT64 or maximum > _MAX_INT64:
            # allowable time range = 292,471 years
            raise OverflowError('time values outside range of 64 bit signed integers')
        else:
            return num * factor


def decode_date_from_scalar(time_in_microseconds, basedate):
    """Decode a date from a scalar input."""
    delta = time_in_microseconds.astype("timedelta64[us]").astype(timedelta)
    try:
        return basedate + delta
    except OverflowError:
        raise ValueError("OverflowError in datetime, possibly because year < datetime.MINYEAR")


def decode_dates_from_array(times_in_microseconds, basedate):
    """Decode values encoded by an integer array in units of microseconds to dates.
    
    This is an optimized algorithm that operates by flattening and sorting the input
    array of integers, decoding the first date using the original base date, and then
    incrementally adding timedeltas to decode the rest of the dates in the array.  This
    is an optimal approach, because it minimizes the length of the timedeltas used in
    each addition operation required to decode the times (timedelta addition is the rate
    limiting step in the process).  The original order of the elements and shape of the
    array are restored at the end.  The sorting and unsorting steps add only a small
    overhead.  See discussion and timing results in GitHub issue 269.
    """
    original_shape = times_in_microseconds.shape
    times_in_microseconds = times_in_microseconds.ravel()

    sort_indices = np.argsort(times_in_microseconds)
    unsort_indices = np.argsort(sort_indices)
    times_in_microseconds = times_in_microseconds[sort_indices]

    # We first cast to the np.timedelta64[us] dtype out of convenience, but ultimately
    # cast to datetime.timedelta objects for operations with cftime objects (we cannot
    # cast from integers to datetime.timedelta objects directly).
    deltas = times_in_microseconds.astype("timedelta64[us]")
    differential_deltas = np.diff(deltas).astype(timedelta)

    dates = np.empty(times_in_microseconds.shape, dtype="O")
    try:
        dates[0] = basedate + deltas[0].astype(timedelta)
        for i in range(len(differential_deltas)):
            dates[i + 1] = dates[i] + differential_deltas[i]
    except OverflowError:
        raise ValueError("OverflowError in datetime, possibly because year < datetime.MINYEAR")

    return dates[unsort_indices].reshape(original_shape)


@cython.embedsignature(True)
def num2date(
    times,
    units,
    calendar='standard',
    only_use_cftime_datetimes=True,
    only_use_python_datetimes=False,
    has_year_zero=None
):
    """
    Return datetime objects given numeric time values. The units
    of the numeric time values are described by the **units** argument
    and the **calendar** keyword. The returned datetime objects represent
    UTC with no time-zone offset, even if the specified
    **units** contain a time-zone offset.

    **times**: numeric time values.

    **units**: a string of the form **<time units> since <reference time>**
    describing the time units. **<time units>** can be days, hours, minutes,
    seconds, milliseconds or microseconds. **<reference time>** is the time
    origin. **months since** is allowed *only* for the **360_day** calendar
    and **common_years since** is allowed *only* for the **365_day** calendar.

    **calendar**: describes the calendar used in the time calculations.
    All the values currently defined in the
    `CF metadata convention <http://cfconventions.org/cf-conventions/cf-conventions#calendar>`__ are supported.
    Valid calendars **'standard', 'gregorian', 'proleptic_gregorian'
    'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day'**.
    Default is **'standard'**, which is a mixed Julian/Gregorian calendar.

    **only_use_cftime_datetimes**: if False, python datetime.datetime
    objects are returned from num2date where possible; if True dates which
    subclass cftime.datetime are returned for all calendars. Default **True**.

    **only_use_python_datetimes**: always return python datetime.datetime
    objects and raise an error if this is not possible. Ignored unless
    **only_use_cftime_datetimes=False**. Default **False**.

    **has_year_zero**: if set to True, astronomical year numbering
    is used and the year zero exists.  If set to False for real-world
    calendars, then historical year numbering is used and the year 1 is
    preceded by year -1 and no year zero exists.
    The defaults are set to conform with
    CF version 1.9 conventions (False for 'julian', 'gregorian'/'standard', True
    for 'proleptic_gregorian' (ISO 8601) and True for the idealized
    calendars 'noleap'/'365_day', '360_day', 366_day'/'all_leap')
    The defaults can only be over-ridden for the real-world calendars,
    for the the idealized calendars the year zero 
    always exists and the has_year_zero kwarg is ignored.
    This kwarg is not needed to define calendar systems allowed by CF
    (the calendar-specific defaults do this).

    returns a datetime instance, or an array of datetime instances with
    microsecond accuracy, if possible.

    ***Note***: If only_use_cftime_datetimes=False and
    use_only_python_datetimes=False, the datetime instances
    returned are 'real' python datetime
    objects if **calendar='proleptic_gregorian'**, or
    **calendar='standard'** or **'gregorian'**
    and the date is after the breakpoint between the Julian and
    Gregorian calendars (1582-10-15). Otherwise, they are ctime.datetime
    objects which support some but not all the methods of native python
    datetime objects. The datetime instances
    do not contain a time-zone offset, even if the specified **units**
    contains one.
    """
    calendar = calendar.lower()
    # set calendar-specific defaults for has_year_zero
    if has_year_zero is None:
        has_year_zero = _year_zero_defaults(calendar)
    basedate = _dateparse(units,calendar=calendar,has_year_zero=has_year_zero)

    can_use_python_datetime=_can_use_python_datetime(basedate,calendar)
    if not only_use_cftime_datetimes and only_use_python_datetimes:
        if not can_use_python_datetime:
            msg='illegal calendar or reference date for python datetime'
            raise ValueError(msg)

    unit, ignore = _datesplit(units)

    # real-world calendars limited to positive reference years.
    if calendar in ['julian', 'standard', 'gregorian', 'proleptic_gregorian'] and not has_year_zero:
        if basedate.year == 0:
            msg='zero not allowed as a reference year unless has_year_zero=True'
            raise ValueError(msg)

    use_python_datetime = False
    if only_use_python_datetimes and not only_use_cftime_datetimes:
        # only_use_cftime_datetimes takes precedence
        use_python_datetime = True
    if not only_use_python_datetimes and not only_use_cftime_datetimes and can_use_python_datetime:
        # if only_use_cftimes_datetimes and only_use_python_datetimes are False
        # return python datetime if possible.
        use_python_datetime = True

    basedate = to_calendar_specific_datetime(basedate, calendar, use_python_datetime, has_year_zero=has_year_zero)

    if unit not in UNIT_CONVERSION_FACTORS:
        raise ValueError("Unsupported time units provided, {!r}.".format(unit))

    if unit in ["months", "month"] and calendar != "360_day":
        raise ValueError("Units of months only valid for 360_day calendar.")

    factor = UNIT_CONVERSION_FACTORS[unit]
    times = np.asanyarray(times)  # Allow list as input
    times = upcast_times(times)
    # convert to masked array if any nan or inf values present
    if not np.isfinite(times).all():
        times = np.ma.masked_invalid(times)
    scaled_times = scale_times(times, factor)
    scaled_times = cast_to_int(scaled_times,units=unit)

    if scaled_times.ndim == 0 or scaled_times.size == 0:
        return decode_date_from_scalar(scaled_times, basedate)
    else:
        if isinstance(scaled_times, np.ma.MaskedArray):
            # The algorithm requires data be present for all values. To handle this, we fill 
            # masked values with 0 temporarily and then restore the mask at the end.
            original_mask = np.ma.getmask(scaled_times)
            scaled_times = scaled_times.filled(0)
            dates = decode_dates_from_array(scaled_times, basedate)
            return np.ma.MaskedArray(dates, mask=original_mask)
        else:
            return decode_dates_from_array(scaled_times, basedate)


@cython.embedsignature(True)
def date2index(dates, nctime, calendar=None, select='exact', has_year_zero=None):
    """
    Return indices of a netCDF time variable corresponding to the given dates.

    **dates**: A datetime object or a sequence of datetime objects.
    The datetime objects should not include a time-zone offset.

    **nctime**: A netCDF time variable object. The nctime object must have a
    **units** attribute.

    **calendar**: describes the calendar to be used in the time calculations.
    All the values currently defined in the
    `CF metadata convention <http://cfconventions.org/cf-conventions/cf-conventions#calendar>`__ are supported.
    Valid calendars **'standard', 'gregorian', 'proleptic_gregorian'
    'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day'**.
    Default is `None` which means the calendar associated with the first
    input datetime instance will be used.

    **select**: **'exact', 'before', 'after', 'nearest'**
    The index selection method. **exact** will return the indices perfectly
    matching the dates given. **before** and **after** will return the indices
    corresponding to the dates just before or just after the given dates if
    an exact match cannot be found. **nearest** will return the indices that
    correspond to the closest dates.

    **has_year_zero**: if set to True, astronomical year numbering
    is used and the year zero exists.  If set to False for real-world
    calendars, then historical year numbering is used and the year 1 is
    preceded by year -1 and no year zero exists.
    The defaults are set to conform with
    CF version 1.9 conventions (False for 'julian', 'gregorian'/'standard', True
    for 'proleptic_gregorian' (ISO 8601) and True for the idealized
    calendars 'noleap'/'365_day', '360_day', 366_day'/'all_leap')
    The defaults can only be over-ridden for the real-world calendars,
    for the the idealized calendars the year zero 
    always exists and the has_year_zero kwarg is ignored.
    This kwarg is not needed to define calendar systems allowed by CF
    (the calendar-specific defaults do this).

    returns an index (indices) of the netCDF time variable corresponding
    to the given datetime object(s).
    """
    try:
        nctime.units
    except AttributeError:
        raise AttributeError("netcdf time variable is missing a 'units' attribute")

    dates_test = np.asanyarray(dates) # convert to numpy array

    # if has_year_zero is None and calendar is None, use first input cftime.datetime instance
    # to determine if year zero is to be included.
    # If calendar is specified, use calendar specific defaults
    if has_year_zero is None:
        if calendar is None:
            d0 = dates_test.item(0)
            if isinstance(d0,datetime_python):
                has_year_zero = False
            else:
                try:
                    has_year_zero = d0.has_year_zero
                except AttributeError:
                    raise ValueError('has_year_zero not specified',type(d0))
        else:
            has_year_zero = _year_zero_defaults(calendar)

    # if calendar is None or '', use calendar of first input cftime.datetime instances.
    # if inputs are 'real' python datetime instances, use proleptic gregorian.
    if not calendar:
        d0 = dates_test.item(0)
        if isinstance(d0,datetime_python):
            calendar = 'proleptic_gregorian' 
        else:
            try:
                calendar = d0.calendar
            except AttributeError:
                raise ValueError('no calendar specified',type(d0))

    calendar = calendar.lower()
    basedate = _dateparse(nctime.units,calendar=calendar,has_year_zero=has_year_zero)
    # real-world calendars limited to positive reference years.
    if calendar in ['julian', 'standard', 'gregorian', 'proleptic_gregorian'] and not has_year_zero:
        if basedate.year == 0:
            msg='zero not allowed as a reference year, does not exist in Julian or Gregorian calendars'
            raise ValueError(msg)

    if _can_use_python_datetime(basedate,calendar):
        # use python datetime
        times = date2num(dates,nctime.units,calendar=calendar,has_year_zero=has_year_zero)
        return time2index(times, nctime, calendar, select)
    else: # use cftime module for other cases
        return _date2index(dates, nctime, calendar, select, has_year_zero=has_year_zero)


cdef _parse_timezone(tzstring):
    """Parses ISO 8601 time zone specs into tzinfo offsets

    Adapted from pyiso8601 (http://code.google.com/p/pyiso8601/)
    """
    if tzstring == "Z":
        return 0
    # This isn't strictly correct, but it's common to encounter dates without
    # time zones so I'll assume the default (which defaults to UTC).
    if tzstring is None:
        return 0
    m = TIMEZONE_REGEX.match(tzstring)
    prefix, hours, minutes1, minutes2 = m.groups()
    hours = int(hours)
    # Note: Minutes don't have to be specified in tzstring, so if the group is not found it means minutes is 0.
    #       Also, due to the timezone regex definition, there are two mutually exclusive groups that might hold the minutes value, so check both.
    minutes = int(minutes1) if minutes1 is not None else int(minutes2) if minutes2 is not None else 0
    if prefix == "-":
        hours = -hours
        minutes = -minutes
    return minutes + hours * 60.


cpdef _parse_date(datestring):
    """Parses ISO 8601 dates into datetime objects

    The timezone is parsed from the date string, assuming UTC
    by default.

    Note that a seconds element with a fractional component
    (e.g. 12.5) is converted into integer seconds and integer
    microseconds.

    Adapted from pyiso8601 (http://code.google.com/p/pyiso8601/)

    """
    if not isinstance(datestring, str) and not isinstance(datestring, unicode):
        raise ValueError("Expecting a string %r" % datestring)
    m = ISO8601_REGEX.match(datestring.strip())
    if not m:
        raise ValueError("Unable to parse date string %r" % datestring)
    groups = m.groupdict()
    tzoffset_mins = _parse_timezone(groups["timezone"])
    if groups["hour"] is None:
        groups["hour"] = 0
    if groups["minute"] is None:
        groups["minute"] = 0
    if groups["second"] is None:
        groups["second"] = 0
    if groups["fraction"] is None:
        groups["fraction"] = 0
    else:
        groups["fraction"] = int(float("0.%s" % groups["fraction"]) * 1e6)
    iyear = int(groups["year"])
    return iyear, int(groups["month"]), int(groups["day"]),\
        int(groups["hour"]), int(groups["minute"]), int(groups["second"]),\
        int(groups["fraction"]),\
        tzoffset_mins

cdef _check_index(indices, times, nctime, select):
    """Return True if the time indices given correspond to the given times,
    False otherwise.

    Parameters:

    indices : sequence of integers
    Positive integers indexing the time variable.

    times : sequence of times.
    Reference times.

    nctime : netCDF Variable object
    NetCDF time object.

    select : string
    Index selection method.
    """
    N = nctime.shape[0]
    if (indices < 0).any():
        return False

    if (indices >= N).any():
        return False

    try:
        t = nctime[indices]
        nctime = nctime
    # WORKAROUND TO CHANGES IN SLICING BEHAVIOUR in 1.1.2
    # this may be unacceptably slow...
    # if indices are unsorted, or there are duplicate
    # values in indices, read entire time variable into numpy
    # array so numpy slicing rules can be used.
    except IndexError:
        nctime = nctime[:]
        t = nctime[indices]
# if fancy indexing not available, fall back on this.
#   t=[]
#   for ind in indices:
#       t.append(nctime[ind])

    if select == 'exact':
        return np.all(t == times)

    elif select == 'before':
        ta = nctime[np.clip(indices + 1, 0, N - 1)]
        return np.all(t <= times) and np.all(ta > times)

    elif select == 'after':
        tb = nctime[np.clip(indices - 1, 0, N - 1)]
        return np.all(t >= times) and np.all(tb < times)

    elif select == 'nearest':
        ta = nctime[np.clip(indices + 1, 0, N - 1)]
        tb = nctime[np.clip(indices - 1, 0, N - 1)]
        delta_after = ta - t
        delta_before = t - tb
        delta_check = np.abs(times - t)
        return np.all(delta_check <= delta_after) and np.all(delta_check <= delta_before)


def _date2index(dates, nctime, calendar=None, select='exact', has_year_zero=None):
    """
    Return indices of a netCDF time variable corresponding to the given dates.

    **dates**: A datetime object or a sequence of datetime objects.
    The datetime objects should not include a time-zone offset.

    **nctime**: A netCDF time variable object. The nctime object must have a
    `units` attribute. The entries are assumed to be stored in increasing
    order.

    **calendar**: Describes the calendar used in the time calculation.
    Valid calendars 'standard', 'gregorian', 'proleptic_gregorian'
    'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day'.
    Default is 'standard', which is a mixed Julian/Gregorian calendar
    If `calendar` is None, its value is given by `nctime.calendar` or
    `standard` if no such attribute exists.

    **select**: 'exact', 'before', 'after', 'nearest'
    The index selection method. `exact` will return the indices perfectly
    matching the dates given. `before` and `after` will return the indices
    corresponding to the dates just before or just after the given dates if
    an exact match cannot be found. `nearest` will return the indices that
    correspond to the closest dates.

    **has_year_zero**: if set to True, astronomical year numbering
    is used and the year zero exists.  If set to False for real-world
    calendars, then historical year numbering is used and the year 1 is
    preceded by year -1 and no year zero exists.
    The defaults are set to conform with
    CF version 1.9 conventions (False for 'julian', 'gregorian'/'standard', True
    for 'proleptic_gregorian' (ISO 8601) and True for the idealized
    calendars 'noleap'/'365_day', '360_day', 366_day'/'all_leap')
    The defaults can only be over-ridden for the real-world calendars,
    for the the idealized calendars the year zero 
    always exists and the has_year_zero kwarg is ignored.
    This kwarg is not needed to define calendar systems allowed by CF
    (the calendar-specific defaults do this).
    """
    try:
        nctime.units
    except AttributeError:
        raise AttributeError("netcdf time variable is missing a 'units' attribute")
    times = date2num(dates,nctime.units,calendar=calendar, has_year_zero=has_year_zero)
    return time2index(times, nctime, calendar=calendar, select=select)


@cython.embedsignature(True)
def time2index(times, nctime, calendar=None, select='exact'):
    """
    Return indices of a netCDF time variable corresponding to the given times.

    **times**: A numeric time or a sequence of numeric times.

    **nctime**: A netCDF time variable object. The nctime object must have a
    `units` attribute. The entries are assumed to be stored in increasing
    order.

    **calendar**: Describes the calendar used in the time calculation.
    Valid calendars 'standard', 'gregorian', 'proleptic_gregorian'
    'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day'.
    Default is `standard`, which is a mixed Julian/Gregorian calendar
    If `calendar` is None, its value is given by `nctime.calendar` or
    `standard` if no such attribute exists.

    **select**: **'exact', 'before', 'after', 'nearest'**
    The index selection method. `exact` will return the indices perfectly
    matching the times given. `before` and `after` will return the indices
    corresponding to the times just before or just after the given times if
    an exact match cannot be found. `nearest` will return the indices that
    correspond to the closest times.
    """
    try:
        nctime.units
    except AttributeError:
        raise AttributeError("netcdf time variable is missing a 'units' attribute")
    # Setting the calendar.
    if calendar == None:
        calendar = getattr(nctime, 'calendar', 'standard')

    if select != 'exact':
        # if select works, then 'nearest' == 'exact', 'before' == 'exact'-1 and
        # 'after' == 'exact'+1
        try:
            index = time2index(times, nctime, calendar=calendar, select='exact')
            if select == 'nearest':
                return index
            elif select == 'before':
                return index-1
            else:
                return index+1
        except ValueError:
            pass

    num = np.atleast_1d(times)
    N = len(nctime)

    # Trying to infer the correct index from the starting time and the stride.
    # This assumes that the times are increasing uniformly.
    if len(nctime) >= 2:
        t0, t1 = nctime[:2]
        dt = t1 - t0
    else:
        t0 = nctime[0]
        dt = 1.
    if select in ['exact', 'before']:
        index = np.array((num - t0) / dt, int)
    elif select == 'after':
        index = np.array(np.ceil((num - t0) / dt), int)
    else:
        index = np.array(np.around((num - t0) / dt), int)

    # Checking that the index really corresponds to the given time.
    # If the times do not correspond, then it means that the times
    # are not increasing uniformly and we try the bisection method.
    if not _check_index(index, times, nctime, select):

        # Use the bisection method. Assumes nctime is ordered.
        import bisect
        index = np.array([bisect.bisect_right(nctime, n) for n in num], int)
        before = index == 0

        index = np.array([bisect.bisect_left(nctime, n) for n in num], int)
        after = index == N

        if select in ['before', 'exact'] and np.any(before):
            raise ValueError(
                'Some of the times given are before the first time in **nctime**.')

        if select in ['after', 'exact'] and np.any(after):
            raise ValueError(
                'Some of the times given are after the last time in **nctime**.')

        # Find the times for which the match is not perfect.
        # Use list comprehension instead of the simpler **nctime[index]** since
        # not all time objects support numpy integer indexing (eg dap).
        index[after] = N - 1
        ncnum = np.squeeze([nctime[i] for i in index])
        mismatch = np.nonzero(ncnum != num)[0]

        if select == 'exact':
            if len(mismatch) > 0:
                raise ValueError(
                    'Some of the times specified were not found in the **nctime** variable.')

        elif select == 'before':
            index[after] = N
            index[mismatch] -= 1

        elif select == 'after':
            pass

        elif select == 'nearest':
            nearest_to_left = num[mismatch] < np.array(
                [float(nctime[i - 1]) + float(nctime[i]) for i in index[mismatch]]) / 2.
            index[mismatch] = index[mismatch] - 1 * nearest_to_left

        else:
            raise ValueError(
                "%s is not an option for the **select** argument." % select)

        # Correct for indices equal to -1
        index[before] = 0

    # convert numpy scalars or single element arrays to python ints.
    return _toscalar(index)


cdef _toscalar(a):
    if a.shape in [(), (1,)]:
        return a.item()
    else:
        return a

def to_tuple(dt):
    """Turn a datetime.datetime instance into a tuple of integers. Elements go
    in the order of decreasing significance, making it easy to compare
    datetime instances. Parts of the state that don't affect ordering
    are omitted. Compare to datetime.timetuple()."""
    return (dt.year, dt.month, dt.day, dt.hour, dt.minute,
            dt.second, dt.microsecond)

cdef _year_zero_defaults(calendar):
    if calendar: calendar = calendar.lower()
    if calendar in ['standard','gregorian','julian']:
       return False
    elif calendar in ['proleptic_gregorian']:
       return True # ISO 8601 year zero=1 BC
    elif calendar in _idealized_calendars:
       return True
    else:
       return False

# factory function without optional kwargs that can be used in datetime.__reduce__
def _create_datetime(date_type, args, kwargs): return date_type(*args, **kwargs)
# custorm warning for invalid CF dates.
cfwarnmsg="this date/calendar/year zero convention is not supported by CF"
class CFWarning(UserWarning):
    pass

@cython.embedsignature(True)
cdef class datetime(object):
    """
This class mimics datetime.datetime but support calendars other than the proleptic
Gregorian calendar.

Supports timedelta operations by overloading +/-, and
comparisons with other instances (even if they use different calendars).
When comparing instances with different calendars, the second instance in the comparison (RHS) is
converted to the calendar of the LHS instance.  When comparing a list or array of instances (all using
the same calendar) to a single scalar instance,
it is much faster to convert the single instance to the calendar of the array before doing
the comparison.

Comparison with native python datetime instances is possible
for cftime.datetime instances using
'gregorian' and 'proleptic_gregorian' calendars.

All the calendars currently defined in the
`CF metadata convention <http://cfconventions.org/cf-conventions/cf-conventions#calendar>`__ are supported.
Valid calendars are 'standard', 'gregorian', 'proleptic_gregorian'
'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day'.
Default is 'standard', which is a mixed Julian/Gregorian calendar.
'standard' and 'gregorian' are synonyms, as are 'all_leap'/'366_day'
and 'noleap'/'365_day'.

If the calendar kwarg is set to a blank string ('') or None (the default is 'standard') the 
instance will not be calendar-aware and some methods will not work.

If the has_year_zero kwarg is set to True, astronomical year numbering
is used and the year zero exists.  If set to False for real-world
calendars, then historical year numbering is used and the year 1 is
preceded by year -1 and no year zero exists.
The defaults are set to conform with
CF version 1.9 conventions (False for 'julian', 'gregorian'/'standard', True
for 'proleptic_gregorian' (ISO 8601) and True for the idealized
calendars 'noleap'/'365_day', '360_day', 366_day'/'all_leap')
The defaults can only be over-ridden for the real-world calendars,
for the the idealized calendars the year zero 
always exists and the has_year_zero kwarg is ignored.
This kwarg is not needed to define calendar systems allowed by CF
(the calendar-specific defaults do this).

Has isoformat, strftime, timetuple, replace, dayofwk, dayofyr, daysinmonth,
__repr__,__format__, __add__, __sub__, __str__ and comparison methods. 

dayofwk, dayofyr, daysinmonth, __add__ and __sub__ only work for calendar-aware
instances.

The default format of the string produced by strftime is controlled by self.format
(default %Y-%m-%d %H:%M:%S).
    """
    cdef readonly int year, month, day, hour, minute
    cdef readonly int second, microsecond
    cdef readonly str calendar
    cdef readonly int _dayofwk, _dayofyr
    cdef readonly bint has_year_zero
    cdef readonly object tzinfo

    # Python's datetime.datetime uses the proleptic Gregorian
    # calendar. This boolean is used to decide whether a
    # cftime.datetime instance can be converted to
    # datetime.datetime.
    cdef readonly bint datetime_compatible

    def __init__(self, int year, int month, int day, int hour=0, int minute=0,
                       int second=0, int microsecond=0, int dayofwk=-1, 
                       int dayofyr = -1, calendar='standard', has_year_zero=None):

        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second
        self.microsecond = microsecond
        self._dayofwk = dayofwk
        self._dayofyr = dayofyr
        self.tzinfo = None
        if calendar:
            calendar = calendar.lower()
        # set calendar-specific defaults for has_year_zero
        if has_year_zero is None:
            if year == 0:
                # assume if user sets year to zero, the calendar should
                # include the year zero (issue #248)
                # warn if calendar is being set to non-CF calendar
                msg="year=0 was specified - this date/calendar/year zero convention is not supported by CF"
                if calendar is not None and not _year_zero_defaults(calendar):
                    warnings.warn(msg,category=CFWarning)
                has_year_zero=True
            else:
                has_year_zero = _year_zero_defaults(calendar)
        # warn if requested date not allowed by CF
        # (no years < 1 in mixed Julian/Gregorian calendar).
        # CF version 1.9 does not specify whether a year zero should exist
        # for the proleptic_gregorian calendar. IS0 8601 uses proleptic_gregorian
        # and has a year zero, so for now this is the default in cftime.
        if calendar in ['julian','gregorian','standard'] and year <= 0:
            warnings.warn(cfwarnmsg,category=CFWarning)
        # raise exception if year zero requested but has_year_zero set
        # to False (issue #248).
        if year == 0 and has_year_zero==False:
            msg='year zero requested, but has_year_zero=False'
            raise ValueError(msg)
        if not has_year_zero and calendar in _idealized_calendars:
            warnings.warn('has_year_zero kwarg ignored for idealized calendars (always True)')
        self.has_year_zero = has_year_zero
        if calendar == 'gregorian' or calendar == 'standard':
            # dates after 1582-10-15 can be converted to and compared to
            # proleptic Gregorian dates
            self.calendar = 'standard'
            if self.to_tuple() >= (1582, 10, 15, 0, 0, 0, 0):
                self.datetime_compatible = True
            else:
                self.datetime_compatible = False
            assert_valid_date(self, is_leap_gregorian, True, has_year_zero=has_year_zero)
        elif calendar == 'noleap' or calendar == '365_day':
            self.calendar = 'noleap'
            self.datetime_compatible = False
            assert_valid_date(self, no_leap, False, has_year_zero=True)
        elif calendar == 'all_leap' or calendar == '366_day':
            self.calendar = 'all_leap'
            self.datetime_compatible = False
            assert_valid_date(self, all_leap, False, has_year_zero=True)
        elif calendar == '360_day':
            self.calendar = calendar
            self.datetime_compatible = False
            assert_valid_date(self, no_leap, False, has_year_zero=True, is_360_day=True)
        elif calendar == 'julian':
            self.calendar = calendar
            self.datetime_compatible = False
            assert_valid_date(self, is_leap_julian, False, has_year_zero=has_year_zero)
        elif calendar == 'proleptic_gregorian':
            self.calendar = calendar
            self.datetime_compatible = True
            assert_valid_date(self, is_leap_proleptic_gregorian, False, has_year_zero=has_year_zero)
        elif calendar == '' or calendar is None:
            # instance not calendar-aware, some method will not work
            self.calendar = ''
            self.datetime_compatible = False
        else:
            raise ValueError(
                "calendar must be one of %s, got '%s'" % (str(_calendars), calendar))

    @property
    def format(self):
        return '%Y-%m-%d %H:%M:%S'

    @property
    def dayofwk(self):
        if self._dayofwk < 0 and self.calendar:
            jd = self.toordinal()
            dayofwk = (jd + 1) % 7
            # convert to ISO 8601 (0 = Monday, 6 = Sunday), like python datetime
            dayofwk -= 1
            if dayofwk == -1: dayofwk = 6
            # cache results for dayofwk
            self._dayofwk = dayofwk
            return dayofwk
        else:
            return self._dayofwk

    @property
    def dayofyr(self):
        if self._dayofyr < 0 and self.calendar:
            if self.calendar == '360_day':
                dayofyr = (self.month-1)*30+self.day
            else:
                if _is_leap(self.year,self.calendar,has_year_zero=self.has_year_zero):
                    dayofyr = _cumdayspermonth_leap[self.month-1]+self.day
                else:
                    dayofyr = _cumdayspermonth[self.month-1]+self.day
            # cache results for dayofyr
            self._dayofyr = dayofyr
            return dayofyr
        else:
            return self._dayofyr

    @property
    def daysinmonth(self):
        if self.calendar == '360_day':
            return 30
        else:
            if _is_leap(self.year,self.calendar,
                    has_year_zero=self.has_year_zero):
                return _dayspermonth_leap[self.month-1]
            else:
                return _dayspermonth[self.month-1]

    def strftime(self, format=None):
        """
        Return a string representing the date, controlled by an explicit format
        string. For a complete list of formatting directives, see section
        'strftime() and strptime() Behavior' in the base Python documentation.
        """
        if format is None:
            format = self.format
        return _strftime(self, format)

    @staticmethod
    def strptime(datestring, format, calendar='standard', has_year_zero=None):
        """
        Return a datetime corresponding to date_string, parsed according to format,
        with a specified calendar and year zero convention.
        The format directives 'y','Y','m','B','b','d','H','M','S' and 'f'
        are supported for all calendars and dates.  If the date is valid
        in the python 'proleptic_gregorian' calendar, then python's
        datetime.strptime is used. For a complete list of formatting directives
        supported in python's datetime.strptime, see section
        'strftime() and strptime() Behavior' in the base Python documentation.
        """
        # if possible use python's datetime.strptime to get a python datetime instance
        # (works for dates in proleptic_gregorian calendar)
        fd = [d[0] for d in format.split('%') if d] # extract format descriptors
        # calendar specific format descriptors that won't work will all calendars
        special_fd = ['a', 'A', 'w', 'j', 'U', 'W', 'G', 'u', 'V']
        try:
            pydatetime = datetime_python.strptime(datestring, format)
            # remove time zone offset
            if getattr(pydatetime, 'tzinfo',None) is not None:
                 pydatetime = pydatetime.replace(tzinfo=None) - pydatetime.utcoffset()
            compatible_date =\
            calendar == 'proleptic_gregorian' or \
            (calendar in ['gregorian','standard'] and (pydatetime.year > 1582 or \
             (pydatetime.year == 1582 and pydatetime.month > 10) or \
             (pydatetime.year == 1582 and pydatetime.month == 10 and pydatetime.day > 15)))
            if not compatible_date and any(x in special_fd for x in fd):
                msg='one of the supplied format directives may not be consistent with the chosen calendar'
                raise KeyError(msg)
            # convert the cftime datetime instance
            return datetime(pydatetime.year, pydatetime.month, pydatetime.day,
                            pydatetime.hour, pydatetime.minute, pydatetime.second,
                            pydatetime.microsecond, calendar=calendar, has_year_zero=has_year_zero)
        # otherwise use a stripped-down version of C-python's _strptime.py
        # (doesn't understand all possible formats, just
        # 'y','Y','m','B','b','d','H','M','S' and 'f')
        except ValueError:
            year,month,day,hour,minute,second,microsecond = _strptime(datestring,format)
            return datetime(year,month,day,hour,minute,second,microsecond,
                            calendar=calendar,has_year_zero=has_year_zero)

    def __format__(self, format):
        # the string format "{t_obj}".format(t_obj=t_obj)
        # without an explicit format gives an empty string (format='')
        # so set this to None to get the default strftime behaviour
        if format == '':
            format = None
        return self.strftime(format)

    def replace(self, **kwargs):
        """Return datetime with new specified fields."""
        args = {"year": self.year,
                "month": self.month,
                "day": self.day,
                "hour": self.hour,
                "minute": self.minute,
                "second": self.second,
                "microsecond": self.microsecond,
                "has_year_zero": self.has_year_zero,
                "calendar": self.calendar}

        if 'dayofyr' in kwargs or 'dayofwk' in kwargs:
            raise ValueError('Replacing the dayofyr or dayofwk of a datetime is '
                             'not supported.')

        if 'calendar' in kwargs:
            raise ValueError('Replacing the calendar of a datetime is '
                             'not supported.')

        # if attempting to set year to zero, also set has_year_zero=True
        # (issue #248)
        if 'year' in kwargs:
            if kwargs['year']==0 and 'has_year_zero' not in kwargs:
                kwargs['has_year_zero']=True

        for name, value in kwargs.items():
            args[name] = value

        return self.__class__(**args)

    def timetuple(self):
        """
        Return a time.struct_time such as returned by time.localtime().
        The DST flag is -1. d.timetuple() is equivalent to
        time.struct_time((d.year, d.month, d.day, d.hour, d.minute,
        d.second, d.weekday(), yday, dst)), where yday is the
        day number within the current year starting with 1 for January 1st.
        """
        return time.struct_time((self.year, self.month, self.day, self.hour,
                self.minute, self.second, self.dayofwk, self.dayofyr, -1))

    cpdef _to_real_datetime(self):
        return real_datetime(self.year, self.month, self.day,
                             self.hour, self.minute, self.second,
                             self.microsecond)

    def __repr__(self):
        if self.__class__.__name__ != 'datetime': # a calendar-specific sub-class
            return "{0}.{1}({2}, {3}, {4}, {5}, {6}, {7}, {8}, has_year_zero={9})".format('cftime',
            self.__class__.__name__,
            self.year,self.month,self.day,self.hour,self.minute,self.second,self.microsecond,self.has_year_zero)
        if self.calendar == None:
            return "{0}.{1}({2}, {3}, {4}, {5}, {6}, {7}, {8}, calendar={9}, has_year_zero={10})".format('cftime',
            self.__class__.__name__,
            self.year,self.month,self.day,self.hour,self.minute,self.second,self.microsecond,self.calendar,self.has_year_zero)
        else:
            return "{0}.{1}({2}, {3}, {4}, {5}, {6}, {7}, {8}, calendar='{9}', has_year_zero={10})".format('cftime',
            self.__class__.__name__,
            self.year,self.month,self.day,self.hour,self.minute,self.second,self.microsecond,self.calendar,self.has_year_zero)

    def __str__(self):
        return self.isoformat(' ')

    def isoformat(self, sep='T', timespec='auto'):
        """
        ISO date representation

        """
        if self.year < 0:
            form0 = '{:05d}-{:02d}-{:02d}'
        else:
            form0 = '{:04d}-{:02d}-{:02d}'
        if timespec == 'days':
            form = form0
            return form.format(self.year, self.month, self.day)
        elif timespec == 'hours':
            form = form0 + '{:s}{:02d}'
            return form.format(self.year, self.month, self.day, sep,
                               self.hour)
        elif timespec == 'minutes':
            form = form0 + '{:s}{:02d}:{:02d}'
            return form.format(self.year, self.month, self.day, sep,
                               self.hour, self.minute)
        elif timespec == 'seconds':
            form = form0 + '{:s}{:02d}:{:02d}:{:02d}'
            return form.format(self.year, self.month, self.day, sep,
                               self.hour, self.minute, self.second)
        elif timespec in ['auto', 'microseconds', 'milliseconds']:
            second = '{:02d}'.format(self.second)
            if timespec == 'milliseconds':
                millisecs = int(round(self.microsecond / 1000, 0))
                second += '.{:03d}'.format(millisecs)
            elif timespec == 'microseconds':
                second += '.{:06d}'.format(self.microsecond)
            else:
                if self.microsecond > 0:
                    second += '.{:06d}'.format(self.microsecond)
            form = form0 + '{:s}{:02d}:{:02d}:{:s}'
            return form.format(self.year, self.month, self.day, sep,
                               self.hour, self.minute, second)
        else:
            raise ValueError('illegal timespec')

    def __hash__(self):
        try:
            d = self._to_real_datetime()
        except ValueError:
            return hash(self.timetuple())
        return hash(d)

    def to_tuple(self):
        return (self.year, self.month, self.day, self.hour, self.minute,
                self.second, self.microsecond)

    def __richcmp__(self, other, int op):
        cdef datetime dt, dt_other, d1, d2
        dt = self
        if isinstance(other, datetime):
            dt_other = other
            # comparing two datetime instances
            if dt.calendar == dt_other.calendar and dt.has_year_zero == dt_other.has_year_zero:
                return PyObject_RichCompare(dt.to_tuple(), dt_other.to_tuple(), op)
            else:
                # raise an error if either is an idealized calendar (comparison not valid)
                if dt.calendar in _idealized_calendars or other.calendar in _idealized_calendars:
                    raise TypeError("cannot compare {0!r} and {1!r}".format(dt, other))
                # convert one to match calendar of the other.
                other2 = other.change_calendar(dt.calendar,has_year_zero=dt.has_year_zero)
                return PyObject_RichCompare(dt.to_tuple(), other2.to_tuple(), op)
        elif isinstance(other, datetime_python):
            # comparing datetime and real_datetime
            if not dt.datetime_compatible:
                raise TypeError("cannot compare {0!r} and {1!r} (different calendars)".format(self, other))
            return PyObject_RichCompare(dt.to_tuple(), to_tuple(other), op)
        else:
            # Delegate responsibility of comparison to the other object.
            return NotImplemented

    cdef _getstate(self):
        """return args and kwargs needed to create class instance"""
        args = (self.year, self.month, self.day)
        kwargs = {'hour': self.hour,
                  'minute': self.minute,
                  'second': self.second,
                  'microsecond': self.microsecond,
                  'dayofwk': self._dayofwk,
                  'dayofyr': self._dayofyr,
                  'calendar': self.calendar,
                  'has_year_zero': self.has_year_zero}
        return args, kwargs

    def __reduce__(self):
        """special method that allows instance to be pickled"""
        args, kwargs = self._getstate()
        date_type = type(self)
        return (_create_datetime, (date_type, args, kwargs))

    cdef _add_timedelta(self, other):
        return NotImplemented

    @staticmethod
    def fromordinal(jday,calendar='standard',has_year_zero=None):
        """Create a datetime instance from a julian day ordinal, calendar
        and (optionally) year zero convention (inverse of toordinal). The
        Julian day number is the number of days since noon UTC January 1, -4713
        in the proleptic julian calendar with no year zero  (November 24, -4713 
        in the proleptic gregorian calendar that includes the year zero). For
        idealized calendars, the origin is noon UTC of the year zero."""
        calendar = calendar.lower()
        # set calendar-specific defaults for has_year_zero
        if has_year_zero is None:
            has_year_zero = _year_zero_defaults(calendar)
        if calendar in ['standard','julian','gregorian']:
            if has_year_zero:
               units = 'days since -4712-1-1-12'
            else:
               units = 'days since -4713-1-1-12'
        elif calendar == 'proleptic_gregorian':
            if has_year_zero:
                units = 'days since -4713-11-24-12'
            else:
                units = 'days since -4714-11-24-12'
        else:
            units = 'days since 0-1-1-12'
        # suppress warning about invalid CF date (year <= 0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=CFWarning)
            jd = num2date(jday,units=units,calendar=calendar,has_year_zero=has_year_zero)
        return jd

    def toordinal(self,fractional=False):
        """Return (integer) julian day ordinal.

        Day 0 starts at noon January 1 of the year -4713 for the
        julian, gregorian and standard calendars (year -4712 if year
        zero allowed).

        Day 0 starts at noon on November 24 of the year -4714 for the
        proleptic gregorian calendar (year -4713 if year zero allowed).

        Day 0 starts at noon on January 1 of the year zero is for the
        360_day, 365_day, 366_day, all_leap and noleap calendars.
        
        If fractional=True, fractional part of day is included (default
        False)."""
        ijd = _IntJulianDayFromDate(self.year, self.month, self.day, self.calendar,
               skip_transition=False,has_year_zero=self.has_year_zero)
        if fractional:
            fracday = self.hour / np.array(24.,np.longdouble) + \
                      self.minute / np.array(1440.0,np.longdouble) + \
                      (self.second + self.microsecond/(np.array(1.e6,np.longdouble)))/\
                      np.array(86400.0,np.longdouble)
            # at this point jd is an integer representing noon UTC on the given
            # year,month,day.
            # compute fractional day from hour,minute,second,microsecond
            return ijd - 0.5 + fracday
        else:
            return ijd

    def change_calendar(self,calendar,has_year_zero=None):
        cdef datetime dt
        """return a new cftime.datetime instance with a different 'real-world' calendar."""
        if calendar in _idealized_calendars or self.calendar in _idealized_calendars:
            raise ValueError('change_calendar only works for real-world calendars')
        # fixed frame of reference is days since -4713-1-1 in the Julian calendar with no year zero
        return self.fromordinal(self.toordinal(fractional=True),
                                calendar=calendar,has_year_zero=has_year_zero)

    def __add__(self, other):
        cdef datetime dt
        cdef bint has_year_zero
        if isinstance(self, datetime) and isinstance(other, timedelta):
            dt = self
            calendar = self.calendar
            has_year_zero = self.has_year_zero
            delta = other
        elif isinstance(self, timedelta) and isinstance(other, datetime):
            dt = other
            calendar = other.calendar
            has_year_zero = other.has_year_zero
            delta = self
        else:
            return NotImplemented
        # return calendar-specific subclasses for backward compatibility,
        # even though after 1.3.0 this is no longer necessary.
        if calendar == '360_day':
            #return dt.__class__(*add_timedelta_360_day(dt, delta),calendar=calendar)
            return Datetime360Day(*add_timedelta_360_day(dt, delta))
        elif calendar == 'noleap':
            #return dt.__class__(*add_timedelta(dt, delta, no_leap, False, True),calendar=calendar)
            return DatetimeNoLeap(*add_timedelta(dt, delta, no_leap, False, True))
        elif calendar == 'all_leap':
            #return dt.__class__(*add_timedelta(dt, delta, all_leap, False, True),calendar=calendar)
            return DatetimeAllLeap(*add_timedelta(dt, delta, all_leap, False, True))
        elif calendar == 'julian':
            #return dt.__class__(*add_timedelta(dt, delta, is_leap_julian, False, has_year_zero),calendar=calendar,has_year_zero=has_year_zero)
            return DatetimeJulian(*add_timedelta(dt, delta, is_leap_julian, False, has_year_zero),has_year_zero=has_year_zero)
        elif calendar == 'standard':
            #return dt.__class__(*add_timedelta(dt, delta, is_leap_gregorian, True, has_year_zero),calendar=calendar,has_year_zero=has_year_zero)
            return DatetimeGregorian(*add_timedelta(dt, delta, is_leap_gregorian, True, has_year_zero),has_year_zero=has_year_zero)
        elif calendar == 'proleptic_gregorian':
            #return dt.__class__(*add_timedelta(dt, delta, is_leap_proleptic_gregorian, False, has_year_zero),calendar=calendar,has_year_zero=has_year_zero)
            return DatetimeProlepticGregorian(*add_timedelta(dt, delta, is_leap_proleptic_gregorian, False, has_year_zero),has_year_zero=has_year_zero)
        else:
            return NotImplemented

    def __sub__(self, other):
        cdef datetime dt
        cdef bint has_year_zero
        if isinstance(self, datetime): # left arg is a datetime instance
            dt = self
            if isinstance(other, datetime):
                # datetime - datetime
                if dt.calendar != other.calendar:
                    raise TypeError("cannot compute the time difference between dates with different calendars")
                if dt.calendar == "":
                    raise TypeError("cannot compute the time difference between dates that are not calendar-aware")
                if dt.has_year_zero != other.has_year_zero:
                    raise TypeError("cannot compute the time difference between dates with year zero conventions")
                ordinal_self = self.toordinal() # julian day
                ordinal_other = other.toordinal()
                days = ordinal_self - ordinal_other
                seconds_self = dt.second + 60 * dt.minute + 3600 * dt.hour
                seconds_other = other.second + 60 * other.minute + 3600 * other.hour
                seconds = seconds_self - seconds_other
                microseconds = dt.microsecond - other.microsecond
                return timedelta(days, seconds, microseconds)
            elif isinstance(other, datetime_python):
                # datetime - real_datetime
                if not dt.datetime_compatible:
                    msg="""
Cannot compute the time difference between dates with different calendars.
One of the datetime objects may have been converted to a native python
datetime instance.  Try using only_use_cftime_datetimes=True when creating the
datetime object."""
                    raise TypeError(msg)
                return dt._to_real_datetime() - other
            elif isinstance(other, timedelta):
                # datetime - timedelta
                # return calendar-specific subclasses for backward compatibility,
                # even though after 1.3.0 this is no longer necessary.
                has_year_zero=self.has_year_zero
                if self.calendar == '360_day':
                    #return self.__class__(*add_timedelta_360_day(self, -other),calendar=self.calendar)
                    return Datetime360Day(*add_timedelta_360_day(self, -other))
                elif self.calendar == 'noleap':
                    #return self.__class__(*add_timedelta(self, -other, no_leap, False, True),calendar=self.calendar)
                    return DatetimeNoLeap(*add_timedelta(self, -other, no_leap, False, True))
                elif self.calendar == 'all_leap':
                    #return self.__class__(*add_timedelta(self, -other, all_leap, False, True),calendar=self.calendar)
                    return DatetimeAllLeap(*add_timedelta(self, -other, all_leap, False, True))
                elif self.calendar == 'julian':
                    #return self.__class__(*add_timedelta(self, -other, 
                    #     is_leap_julian, False, has_year_zero),calendar=self.calendar,has_year_zero=self.has_year_zero)
                    return DatetimeJulian(*add_timedelta(self, -other, is_leap_julian, False, has_year_zero),has_year_zero=self.has_year_zero)
                elif self.calendar == 'standard':
                    #return self.__class__(*add_timedelta(self, -other, 
                    #     is_leap_gregorian, True, has_year_zero),calendar=self.calendar,has_year_zero=self.has_year_zero)
                    return DatetimeGregorian(*add_timedelta(self, -other, is_leap_gregorian, True, has_year_zero),has_year_zero=self.has_year_zero)
                elif self.calendar == 'proleptic_gregorian':
                    #return self.__class__(*add_timedelta(self, -other,
                    #    is_leap_proleptic_gregorian, False, has_year_zero),calendar=self.calendar,has_year_zero=self.has_year_zero)
                    return DatetimeProlepticGregorian(*add_timedelta(self, -other, is_leap_proleptic_gregorian, False, has_year_zero),has_year_zero=self.has_year_zero)
                else:
                    return NotImplemented
            else:
                return NotImplemented
        else:
            if isinstance(self, datetime_python):
                # real_datetime - datetime
                if not other.datetime_compatible:
                    msg="""
Cannot compute the time difference between dates with different calendars.
One of the datetime objects may have been converted to a native python
datetime instance.  Try using only_use_cftime_datetimes=True when creating the
datetime object."""
                    raise TypeError(msg)
                return self - other._to_real_datetime()
            else:
                return NotImplemented

_illegal_s = re.compile(r"((^|[^%])(%%)*%s)")


cdef _findall(text, substr):
    # Also finds overlaps
    sites = []
    i = 0
    while 1:
        j = text.find(substr, i)
        if j == -1:
            break
        sites.append(j)
        i = j + 1
    return sites

# Every 28 years the calendar repeats, except through century leap
# years where it's 6 years.  But only if you're using the Gregorian
# calendar.  ;)
# Make also 4-digit negative years
# Allow .%f for microseconds
cdef _strftime(datetime dt, fmt):
    if _illegal_s.search(fmt):
        raise TypeError("This strftime implementation does not handle %s")
    # don't use strftime method at all.
    # if dt.year > 1900:
    #    return dt.strftime(fmt)
    if '%f' in fmt:
        if not fmt.endswith('.%f'):
            raise TypeError('If %f is used for microseconds it must be the'
                            ' at the end as .%f')
        else:
            ihavems = True
            fmt1 = fmt[:-3]
    else:
        ihavems = False
        fmt1 = fmt

    year = dt.year
    # For every non-leap year century, advance by
    # 6 years to get into the 28-year repeat cycle
    delta = 2000 - year
    off = 6 * (delta // 100 + delta // 400)
    year = year + off

    # Move to around the year 2000
    year = year + ((2000 - year) // 28) * 28
    timetuple = dt.timetuple()
    s1 = time.strftime(fmt1, (year,) + timetuple[1:])
    twodigityear = 'y' in fmt1
    if twodigityear:
        sites1 = _findall(s1, str(year)[-2:])
    else:
        sites1 = _findall(s1, str(year))

    s2 = time.strftime(fmt1, (year + 28,) + timetuple[1:])
    if twodigityear:
       sites2 = _findall(s2, str(year + 28)[-2:])
    else:
       sites2 = _findall(s2, str(year + 28))

    sites = []
    for site in sites1:
        if site in sites2:
            sites.append(site)

    s = s1
    if dt.year < 0:
        syear = "%05d" % (dt.year,)
    else:
        syear = "%04d" % (dt.year,)
    n=4
    if twodigityear:
        syear = syear[-2:]
        if dt.year < 0:
            syear = '-'+syear
        n=2
    for site in sites:
        s = s[:site] + syear + s[site + n:]
    if ihavems:
        s = s + '.{:06d}'.format(dt.microsecond)
    return s

cdef bint is_leap_julian(int year, bint has_year_zero):
    "Return 1 if year is a leap year in the Julian calendar, 0 otherwise."
    return _is_leap(year, calendar='julian',has_year_zero=has_year_zero)

cdef bint is_leap_proleptic_gregorian(int year, bint has_year_zero):
    "Return 1 if year is a leap year in the Proleptic Gregorian calendar, 0 otherwise."
    return _is_leap(year, calendar='proleptic_gregorian',has_year_zero=has_year_zero)

cdef bint is_leap_gregorian(int year, bint has_year_zero):
    "Return 1 if year is a leap year in the Gregorian calendar, 0 otherwise."
    return _is_leap(year, calendar='standard',has_year_zero=has_year_zero)

cdef bint all_leap(int year, bint has_year_zero):
    "Return True for all years."
    return True

cdef bint no_leap(int year, bint has_year_zero):
    "Return False for all years."
    return False

cdef int * month_lengths(bint (*is_leap)(int,bint), int year, bint has_year_zero):
    if is_leap(year,has_year_zero):
        return _dayspermonth_leap
    else:
        return _dayspermonth

cdef void assert_valid_date(datetime dt, bint (*is_leap)(int, bint),
                            bint julian_gregorian_mixed,
                            bint has_year_zero=False,
                            bint is_360_day=False) except *:
    cdef int[12] month_length

    if not has_year_zero:
        if dt.year == 0:
            raise ValueError("invalid year provided in {0!r}".format(dt))
    if is_360_day:
        month_length = 12*[30]
    else:
        month_length = month_lengths(is_leap, dt.year, has_year_zero)


    if dt.month < 1 or dt.month > 12:
        raise ValueError("invalid month provided in {0!r}".format(dt))

    if dt.day < 1 or dt.day > month_length[dt.month-1]:
        raise ValueError("invalid day number provided in {0!r}".format(dt))

    if julian_gregorian_mixed and dt.year == 1582 and dt.month == 10 and dt.day > 4 and dt.day < 15:
        raise ValueError("{0!r} is not present in the mixed Julian/Gregorian calendar".format(dt))

    if dt.hour < 0 or dt.hour > 23:
        raise ValueError("invalid hour provided in {0!r}".format(dt))

    if dt.minute < 0 or dt.minute > 59:
        raise ValueError("invalid minute provided in {0!r}".format(dt))

    if dt.second < 0 or dt.second > 59:
        raise ValueError("invalid second provided in {0!r}".format(dt))

    if dt.microsecond < 0 or dt.microsecond > 999999:
        raise ValueError("invalid microsecond provided in {0!r}".format(dt))

# Add a datetime.timedelta to a cftime.datetime instance. Uses
# integer arithmetic to avoid rounding errors and preserve
# microsecond accuracy.
#
# The argument is_leap is the pointer to a function returning 1 for leap years and 0 otherwise.
#
# This implementation supports 365_day (no_leap), 366_day (all_leap),
# julian, proleptic_gregorian, and the mixed Julian/Gregorian
# (standard, gregorian) calendars by using different is_leap and
# julian_gregorian_mixed arguments.
#
# The date of the transition from the Julian to Gregorian calendar and
# the number of invalid dates are hard-wired (1582-10-4 is the last day
# of the Julian calendar, after which follows 1582-10-15).
cdef tuple add_timedelta(datetime dt, delta, bint (*is_leap)(int,bint), bint julian_gregorian_mixed, bint has_year_zero):
    cdef int microsecond, second, minute, hour, day, month, year
    cdef int delta_microseconds, delta_seconds, delta_days
    cdef int* month_length
    cdef int extra_days, n_invalid_dates

    # extract these inputs here to avoid type conversion in the code below
    delta_microseconds = delta.microseconds
    delta_seconds = delta.seconds
    delta_days = delta.days

    # shift microseconds, seconds, days
    microsecond = dt.microsecond + delta_microseconds
    second = dt.second + delta_seconds
    minute = dt.minute
    hour = dt.hour
    day = dt.day
    month = dt.month
    year = dt.year

    month_length = month_lengths(is_leap, year, has_year_zero)

    n_invalid_dates = 10 if julian_gregorian_mixed else 0

    # Normalize microseconds, seconds, minutes, hours.
    second += microsecond // 1000000
    microsecond = microsecond % 1000000
    minute += second // 60
    second = second % 60
    hour += minute // 60
    minute = minute % 60
    extra_days = hour // 24
    hour = hour % 24

    delta_days += extra_days

    while delta_days < 0:
        if year == 1582 and month == 10 and day > 14 and day + delta_days < 15:
            delta_days -= n_invalid_dates    # skip over invalid dates
        if day + delta_days < 1:
            delta_days += day
            # decrement month
            month -= 1
            if month < 1:
                month = 12
                year -= 1
                if year == 0 and not has_year_zero:
                    year = -1
                month_length = month_lengths(is_leap, year, has_year_zero)
            day = month_length[month-1]
        else:
            day += delta_days
            delta_days = 0

    while delta_days > 0:
        if year == 1582 and month == 10 and day < 5 and day + delta_days > 4:
            delta_days += n_invalid_dates    # skip over invalid dates
        if day + delta_days > month_length[month-1]:
            delta_days -= month_length[month-1] - (day - 1)
            # increment month
            month += 1
            if month > 12:
                month = 1
                year += 1
                if year == 0 and not has_year_zero:
                    year = 1
                month_length = month_lengths(is_leap, year, has_year_zero)
            day = 1
        else:
            day += delta_days
            delta_days = 0

    return (year, month, day, hour, minute, second, microsecond, -1, -1)

# Add a datetime.timedelta to a cftime.datetime instance with the 360_day calendar.
#
# Assumes that the 360_day,365_day and 366_day calendars (unlike the rest of supported
# calendars) have the year 0. Also, there are no leap years and all
# months are 30 days long, so we can compute month and year by using
# "//" and "%".
cdef tuple add_timedelta_360_day(datetime dt, delta):
    cdef int microsecond, second, minute, hour, day, month, year
    cdef int delta_microseconds, delta_seconds, delta_days

    # extract these inputs here to avoid type conversion in the code below
    delta_microseconds = delta.microseconds
    delta_seconds = delta.seconds
    delta_days = delta.days

    # shift microseconds, seconds, days
    microsecond = dt.microsecond + delta_microseconds
    second = dt.second + delta_seconds
    minute = dt.minute
    hour = dt.hour
    day = dt.day + delta_days
    month = dt.month
    year = dt.year

    # Normalize microseconds, seconds, minutes, hours, days, and months.
    second += microsecond // 1000000
    microsecond = microsecond % 1000000
    minute += second // 60
    second = second % 60
    hour += minute // 60
    minute = minute % 60
    day += hour // 24
    hour = hour % 24
    # day and month are counted from 1; all months have 30 days
    month += (day - 1) // 30
    day = (day - 1) % 30 + 1
    # all years have 12 months
    year += (month - 1) // 12
    month = (month - 1) % 12 + 1

    return (year, month, day, hour, minute, second, microsecond, -1, -1)

@cython.embedsignature(True)
def is_leap_year(year, calendar, has_year_zero=None):
    """returns `True` if specified year in specified calendar is
    a leap year. Optional kwarg `has_year_zero`
    controls whether astronomical year numbering
    is used and the year zero exists.  If not specified,
    calendar-specific default is assumed."""
    # public version of _is_leap (issue #259)
    return _is_leap(year, calendar, has_year_zero=has_year_zero)

cdef _is_leap(int year, calendar, has_year_zero=None):
    cdef int tyear
    cdef bint leap
    calendar = _check_calendar(calendar)
    # set calendar-specific defaults for has_year_zero
    if has_year_zero is None:
        has_year_zero = _year_zero_defaults(calendar)
    if year == 0 and not has_year_zero:
        raise ValueError('year zero does not exist in the %s calendar' %\
                calendar)
    # Because there is no year 0 in the Julian calendar, years -1, -5, -9, etc
    # are leap years. year zero is a leap year if it exists.
    if year < 0 and not has_year_zero:
        tyear = year + 1
    else:
        tyear = year
    if calendar == 'proleptic_gregorian' or (calendar == 'standard' and year > 1581):
        if tyear % 4: # not divisible by 4
            leap = False
        elif tyear % 100: # not divisible by 100
            leap = True
        elif tyear % 400: # not divisible by 400
            leap = False
        else:
            leap = True
    elif calendar == 'julian' or (calendar == 'standard' and year < 1582):
        leap = tyear % 4 == 0
    elif calendar == '366_day':
        leap = True
    else:
        leap = False
    return leap

cdef _check_calendar(calendar):
    """validate calendars, convert to subset of names to get rid of synonyms"""
    calendar = calendar.lower()
    if calendar not in _calendars:
        raise ValueError('unsupported calendar')
    calout = calendar
    # remove 'gregorian','noleap','all_leap'
    if calendar == 'gregorian':
        calout = 'standard'
    if calendar == 'noleap':
        calout = '365_day'
    if calendar == 'all_leap':
        calout = '366_day'
    return calout

# The following function (_IntJulianDayFromDate) is based on
# algorithms described in the book
# "Calendrical Calculations" by Dershowitz and Rheingold, 3rd edition, Cambridge University Press, 2007
# and the C implementation provided at https://reingold.co/calendar.C
# with modifications to handle non-real-world calendars and negative years.


cdef _IntJulianDayFromDate(int year,int month,int day,calendar,skip_transition=False,has_year_zero=None):
    """Compute integer Julian Day from year,month,day and calendar.

    Allowed calendars are 'standard', 'gregorian', 'julian',
    'proleptic_gregorian','360_day', '365_day', '366_day', 'noleap',
    'all_leap'.

    'noleap' is a synonym for '365_day'
    'all_leap' is a synonym for '366_day'
    'gregorian' is a synonym for 'standard'

    Integer julian day is number of days since noon UTC -4713-1-1
    in the julian or mixed julian/gregorian calendar, or noon UTC
    -4714-11-24 in the proleptic_gregorian calendar (without year zero).
    Reference date is noon UTC 0-1-1 for other calendars.

    If the has_year_zero kwarg is set to True, astronomical year numbering
    is used and the year zero exists.  If set to False for real-world
    calendars, then historical year numbering is used and the year 1 is
    preceded by year -1 and no year zero exists.
    The defaults are set to conform with
    CF version 1.9 conventions (False for 'julian', 'gregorian'/'standard', True
    for 'proleptic_gregorian' (ISO 8601) and True for the idealized
    calendars 'noleap'/'365_day', '360_day', 366_day'/'all_leap')
    The defaults can only be over-ridden for the real-world calendars,
    for the the idealized calendars the year zero 
    always exists and the has_year_zero kwarg is ignored.
    This kwarg is not needed to define calendar systems allowed by CF
    (the calendar-specific defaults do this).

    Subtract 0.5 to get 00 UTC on that day.

    optional kwarg 'skip_transition':  When True, leave a 10-day
    gap in Julian day numbers between Oct 4 and Oct 15 1582 (the transition
    from Julian to Gregorian calendars).  Default False, ignored
    unless calendar = 'standard'."""
    cdef int jday, jday_jul, jday_greg
    cdef bint leap

    # validate inputs.
    calendar = _check_calendar(calendar)
    if month < 1 or month > 12 or day < 1 or day > 31:
        msg = "date %04d-%02d-%02d does not exist in the %s calendar" %\
        (year,month,day,calendar)
        raise ValueError(msg)
    # set calendar-specific defaults for has_year_zero
    if has_year_zero is None:
        has_year_zero = _year_zero_defaults(calendar)
    if calendar == '360_day':
        return year*360 + (month-1)*30 + day - 1
    elif calendar == '365_day':
        if month == 2 and day == 29:
            raise ValueError('no leap days in 365_day calendar')
        return year*365 + _cumdayspermonth[month-1] + day - 1
    elif calendar == '366_day':
        return year*366 + _cumdayspermonth_leap[month-1] + day - 1

    # handle standard, julian, proleptic_gregorian calendars.
    if year == 0 and not has_year_zero:
        raise ValueError('year zero does not exist in the %s calendar' %\
                calendar)
    leap = _is_leap(year,calendar,has_year_zero=has_year_zero)
    if not leap and month == 2 and day == 29:
        raise ValueError('%s is not a leap year' % year)

    if leap:
        jday = day + _cumdayspermonth_leap[month-1]
    else:
        jday = day + _cumdayspermonth[month-1]
    if year < 0 and not has_year_zero: year += 1
    year += 4800 # add offset so -4800 is year 0.
    # 1st term is the number of days in the last year
    # 2nd term is the number of days in each preceding non-leap year
    # last terms are the number of preceding leap years since -4800
    jday_jul = jday + 365*(year-1) + (year-1)//4
    # remove offset for 87 years before -4713 (including leap days)
    jday_jul -= 31777
    jday_greg = jday + 365*(year-1) + (year-1)//4 - (year-1)//100 + (year-1)//400
    # remove offset, and account for the fact that -4713/1/1 is jday=38 in
    # gregorian calendar.
    jday_greg -= 31739
    if calendar == 'julian':
        return jday_jul
    elif calendar == 'proleptic_gregorian':
        return jday_greg
    elif calendar in ['standard','gregorian']:
        # check for invalid days in mixed calendar (there are 10 missing)
        if jday_jul >= 2299161 and jday_jul < 2299171:
            raise ValueError('invalid date in mixed calendar')
        if jday_jul < 2299161: # 1582 October 15
            return jday_jul
        else:
            if skip_transition:
                return jday_greg+10
            else:
                return jday_greg

# legacy calendar specific sub-classes (will be removed in a future release).

@cython.embedsignature(True)
cdef class DatetimeNoLeap(datetime):
    """
Phony datetime object which mimics the python datetime object,
but uses the "noleap" ("365_day") calendar.
    """
    def __init__(self, *args, **kwargs):
        kwargs['calendar']='noleap'
        super().__init__(*args, **kwargs)

@cython.embedsignature(True)
cdef class DatetimeAllLeap(datetime):
    """
Phony datetime object which mimics the python datetime object,
but uses the "all_leap" ("366_day") calendar.
    """
    def __init__(self, *args, **kwargs):
        kwargs['calendar']='all_leap'
        super().__init__(*args, **kwargs)

@cython.embedsignature(True)
cdef class Datetime360Day(datetime):
    """
Phony datetime object which mimics the python datetime object,
but uses the "360_day" calendar.
    """
    def __init__(self, *args, **kwargs):
        kwargs['calendar']='360_day'
        super().__init__(*args, **kwargs)

@cython.embedsignature(True)
cdef class DatetimeJulian(datetime):
    """
Phony datetime object which mimics the python datetime object,
but uses the "julian" calendar.
    """
    def __init__(self, *args, **kwargs):
        kwargs['calendar']='julian'
        super().__init__(*args, **kwargs)

@cython.embedsignature(True)
cdef class DatetimeGregorian(datetime):
    """
Phony datetime object which mimics the python datetime object,
but uses the mixed Julian-Gregorian ("standard") calendar.
    """
    def __init__(self, *args, **kwargs):
        kwargs['calendar']='standard'
        super().__init__(*args, **kwargs)

@cython.embedsignature(True)
cdef class DatetimeProlepticGregorian(datetime):
    """
Phony datetime object which mimics the python datetime object,
but allows for dates that don't exist in the proleptic gregorian calendar.
    """
    def __init__(self, *args, **kwargs):
        kwargs['calendar']='proleptic_gregorian'
        super().__init__( *args, **kwargs)
