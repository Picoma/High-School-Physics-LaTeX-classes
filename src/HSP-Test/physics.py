# -*- coding: utf-8 -*-
"""
IPython 5+ extension for physical quantity input.
See README.txt for usage examples.

Author: Georg Brandl <georg@python.org>.
This file has been placed in the public domain.
"""

import re
from functools import reduce, total_ordering

import numpy as np
from numpy import pi

def valuetype(val_unit):
    return val_unit[0]

class UnitError(ValueError):
    pass

# Adapted from ScientificPython:
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr>
# with contributions from Greg Ward
# last revision: 2007-5-25

class NumberDict(dict):
    """Dictionary storing numerical values.

    An instance of this class acts like an array of number with generalized
    (non-integer) indices. A value of zero is assumed for undefined
    entries. NumberDict instances support addition, and subtraction with other
    NumberDict instances, and multiplication and division by scalars.
    """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            return 0

    def __add__(self, other):
        sum_dict = NumberDict()
        for key in self:
            sum_dict[key] = self[key]
        for key in other:
            sum_dict[key] = sum_dict[key] + other[key]
        return sum_dict

    def __radd__(self, other):
        return NumberDict(other) + self

    def __sub__(self, other):
        sum_dict = NumberDict()
        for key in self:
            sum_dict[key] = self[key]
        for key in other:
            sum_dict[key] = sum_dict[key] - other[key]
        return sum_dict

    def __rsub__(self, other):
        return NumberDict(other) - self

    def __mul__(self, other):
        new = NumberDict()
        for key in self:
            new[key] = other*self[key]
        return new
    
    __rmul__ = __mul__

    def __truediv__(self, other):
        new = NumberDict()
        for key in self:
            new[key] = self[key]/other
        return new


# Type checks

def isPhysicalUnit(x):
    return hasattr(x, 'factor') and hasattr(x, 'powers')


def isPhysicalQuantity(x):
    return hasattr(x, 'value') and hasattr(x, 'unit')


@total_ordering
class PhysicalUnit(object):
    """Physical unit.

    A physical unit is defined by a name (possibly composite), a scaling factor,
    and the exponentials of each of the SI base units that enter into it. Units
    can be multiplied, divided, and raised to integer powers.
    """

    def __init__(self, names, factor, powers, offset=0):
        """
        @param names: a dictionary mapping each name component to its
                      associated integer power (e.g. C{{'m': 1, 's': -1}})
                      for M{m/s}). As a shorthand, a string may be passed
                      which is assigned an implicit power 1.
        @param factor: a scaling factor
        @param powers: the integer powers for each of the nine base units
        @param offset: an additive offset to the base unit (used only for
                       temperatures)
        """
        if isinstance(names, str):
            self.names = NumberDict()
            self.names[names] = 1
        else:
            self.names = names
        self.factor = factor
        self.offset = offset
        self.powers = powers

    def set_name(self, name):
        self.names = NumberDict()
        self.names[name] = 1

    def name(self):
        num = ''
        denom = ''
        for unit in self.names:
            power = self.names[unit]
            if power < 0:
                denom = denom + '/' + unit
                if power < -1:
                    denom = denom + '**' + str(-power)
            elif power > 0:
                num = num + '*' + unit
                if power > 1:
                    num = num + '**' + str(power)
        if len(num) == 0:
            num = '1'
            if denom == '':
                return ''
        else:
            num = num[1:]
        return num + denom

    def latex(self):
        num = ''
        denom = ''
        for unit in self.names:
            power = self.names[unit]
            if power < 0:
                denom = denom + _latex_table[unit]
#                if power != -1:
                denom = denom + r'\tothe{' + str(power) + r"}"
            elif power > 0:
                num = num + _latex_table[unit]
                if power != 1:
                    num = num + r'\tothe{' + str(power) + r"}"
        return num + denom

    @property
    def is_dimensionless(self):
        return not reduce(lambda a, b: a or b, self.powers)

    @property
    def is_angle(self):
        return self.powers[7] == 1 and \
            reduce(lambda a, b: a + b, self.powers) == 1

    def __str__(self):
        return self.name()

    def __repr__(self):
        return '<PhysicalUnit ' + self.name() + '>'

    def __eq__(self, other):
        if self.powers != other.powers:
            raise UnitError('Incompatible units')
        return self.factor == other.factor

    def __lt__(self, other):
        if self.powers != other.powers:
            raise UnitError('Incompatible units')
        return self.factor < other.factor

    def __mul__(self, other):
        if self.offset != 0 or (isPhysicalUnit(other) and other.offset != 0):
            raise UnitError('Cannot multiply units with non-zero offset')
        if isPhysicalUnit(other):
            return PhysicalUnit(self.names + other.names,
                                self.factor * other.factor,
                                [a+b for (a,b) in zip(self.powers, other.powers)])
        else:
            return PhysicalUnit(self.names + {str(other): 1},
                                self.factor * other, self.powers,
                                self.offset * other)
    __rmul__ = __mul__

    def __truediv__(self, other):
        if self.offset != 0 or (isPhysicalUnit(other) and other.offset != 0):
            raise UnitError('Cannot divide units with non-zero offset')
        if isPhysicalUnit(other):
            return PhysicalUnit(self.names - other.names,
                                self.factor / other.factor,
                                [a-b for (a,b) in zip(self.powers, other.powers)])
        else:
            return PhysicalUnit(self.names+{str(other): -1},
                                self.factor/other, self.powers)

    def __rtruediv__(self, other):
        if self.offset != 0 or (isPhysicalUnit(other) and other.offset != 0):
            raise UnitError('Cannot divide units with non-zero offset')
        if isPhysicalUnit(other):
            return PhysicalUnit(other.names - self.names,
                                other.factor/self.factor,
                                [a-b for (a,b) in zip(other.powers, self.powers)])
        else:
            return PhysicalUnit({str(other): 1} - self.names,
                                other / self.factor,
                                [-x for x in self.powers])

    def __pow__(self, other):
        if self.offset != 0:
            raise UnitError('Cannot exponentiate units with non-zero offset')
        if isinstance(other, int):
            return PhysicalUnit(other*self.names, pow(self.factor, other),
                                [x*other for x in self.powers])
        if isinstance(other, float):
            inv_exp = 1./other
            rounded = int(np.floor(inv_exp + 0.5))
            if abs(inv_exp-rounded) < 1.e-10:
                if reduce(lambda a, b: a and b,
                          map(lambda x, e=rounded: x%e == 0, self.powers)):
                    f = pow(self.factor, other)
                    p = [x/rounded for x in self.powers]
                    if reduce(lambda a, b: a and b,
                              map(lambda x, e=rounded: x%e == 0,
                                  self.names.values())):
                        names = self.names/rounded
                    else:
                        names = NumberDict()
                        if f != 1.:
                            names[str(f)] = 1
                        for i in range(len(p)):
                            names[_base_names[i]] = p[i]
                    return PhysicalUnit(names, f, p)
                else:
                    raise UnitError('Illegal exponent %f' % other)
        raise UnitError('Only integer and inverse integer exponents allowed')

    def conversion_factor_to(self, other):
        """Return conversion factor to another unit."""
        if self.powers != other.powers:
            raise UnitError('Incompatible units')
        if self.offset != other.offset and self.factor != other.factor:
            raise UnitError(('Unit conversion (%s to %s) cannot be expressed ' +
                             'as a simple multiplicative factor') %
                             (self.name(), other.name()))
        return self.factor/other.factor

    def conversion_tuple_to(self, other):
        """Return conversion factor and offset to another unit."""
        if self.powers != other.powers:
            raise UnitError('Incompatible units')

        # let (s1,d1) be the conversion tuple from 'self' to base units
        #   (ie. (x+d1)*s1 converts a value x from 'self' to base units,
        #   and (x/s1)-d1 converts x from base to 'self' units)
        # and (s2,d2) be the conversion tuple from 'other' to base units
        # then we want to compute the conversion tuple (S,D) from
        #   'self' to 'other' such that (x+D)*S converts x from 'self'
        #   units to 'other' units
        # the formula to convert x from 'self' to 'other' units via the
        #   base units is (by definition of the conversion tuples):
        #     ( ((x+d1)*s1) / s2 ) - d2
        #   = ( (x+d1) * s1/s2) - d2
        #   = ( (x+d1) * s1/s2 ) - (d2*s2/s1) * s1/s2
        #   = ( (x+d1) - (d1*s2/s1) ) * s1/s2
        #   = (x + d1 - d2*s2/s1) * s1/s2
        # thus, D = d1 - d2*s2/s1 and S = s1/s2
        factor = self.factor / other.factor
        offset = self.offset - (other.offset * other.factor / self.factor)
        return (factor, offset)


# Helper functions

def _findUnit(unit):
    if isinstance(unit, str):
        name = unit.strip().replace('^', '**').replace('µ', 'mu').replace('°', 'deg')
        try:
            unit = eval(name, _unit_table)
        except NameError:
            raise UnitError('Invalid or unknown unit in %r' % unit)
        for cruft in ['__builtins__', '__args__']:
            try: del _unit_table[cruft]
            except: pass
    if not isPhysicalUnit(unit):
        raise UnitError(str(unit) + ' is not a unit')
    return unit


def _convertValue(value, src_unit, target_unit):
    (factor, offset) = src_unit.conversion_tuple_to(target_unit)
    return (value + offset) * factor


@total_ordering
class PhysicalQuantity(object):
    """Physical quantity with units.

    PhysicalQuantity instances allow addition, subtraction, multiplication, and
    division with each other as well as multiplication, division, and
    exponentiation with numbers.  Addition and subtraction check that the units
    of the two operands are compatible and return the result in the units of the
    first operand. A limited set of mathematical functions (from numpy) is
    applicable as well.
    """

    _number = re.compile(r'([+-]?[0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)'
                         r'(?:\s+\+\/-\s+([+-]?[0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?))?')

    def __init__(self, value, unit=None, stdev=None, prec=3, force_scientific=False):
        """There are two constructor calling patterns:

        1. PhysicalQuantity(value, unit), where value is any number and unit is
           a string defining the unit

        2. PhysicalQuantity(value_with_unit), where value_with_unit is a string
           that contains both the value and the unit, i.e. '1.5 m/s'. This form
           is provided for more convenient interactive use.
        REMOVE THE FORMATTER OPTION ITS USELESS NOW
        """
        self.force_scientific = force_scientific
        self.prec = prec
        self.stdev = stdev
        if unit is not None:
            self.value = valuetype((value, stdev or 0))
            self.unit = _findUnit(unit)
        else:
            s = value.strip()
            match = self._number.match(s)
            if match is None:
                raise UnitError('No number found in %r' % value)
            self.value = valuetype((float(match.group(1)),
                                    float(match.group(2) or 0)))
            self.unit = _findUnit(s[match.end(0):])
    
    def _get_formated_value_string(self):
        formatter = '{val:.{prec}{c}}'
        formatted_num = ''
        if (self.value < 1e2 and self.value > 1e-2) and not self.force_scientific:
            formatted_num = formatter.format(val=self.value, prec=self.prec, c='g')
        else :
            formatted_num = formatter.format(val=self.value, prec=self.prec-1, c='e')
        return(formatted_num)
    
    def __str__(self):
        prec = self.prec
        unit = self.unit.name().replace('**', '^')
        return '{val:.{prec}g} {unit}'.format(val=self.value, prec=prec, unit=unit)
    
    def latex(self, gray=False, addition=''):
        formatted_num = self._get_formated_value_string()
        # On ajoute l'incertitude s'il y a :
        if self.stdev != None and ('e' in formatted_num ):
            # Si il y a incertitude ET un exposant : on intercale "(stdev)" entre la mantisse et l'exposant
            left, right = formatted_num.split("e")
            formatted_num = left + "({})e".format(self.stdev) + right
        elif self.stdev != None:
            # Pas d'exposant qui fait chier : on ajoute direct l'incertitude
            formatted_num += "({})".format(self.stdev)
        
        prec = self.prec
        unit = self.unit.latex()
        if gray :
            return(r'\SI[unit-color=gray]{' \
                   + formatted_num \
                   + addition \
                   + r'}{' \
                   + '{unit}'.format(unit=unit) \
                   + r'}') # addition serves if i want to manually add significant zeros, that are removed by the %g format
        else :
            return(r'\SI{' \
                   + formatted_num \
                   + addition \
                   + r'}{' \
                   + '{unit}'.format(unit=unit) \
                   + r'}')

    def __repr__(self):
        return self.__str__()

    def _sum(self, other, sign1, sign2):
        if not isPhysicalQuantity(other):
            raise UnitError('Incompatible types')
        new_value = sign1 * self.value + \
            sign2 * other.value * other.unit.conversion_factor_to(self.unit)
        new_sig_fig = min(self.prec, other.prec)
        return self.__class__(new_value, self.unit, prec=new_sig_fig)

    def __add__(self, other):
        return self._sum(other, 1, 1)

    __radd__ = __add__

    def __sub__(self, other):
        return self._sum(other, 1, -1)

    def __rsub__(self, other):
        return self._sum(other, -1, 1)

    def __eq__(self, other):
        diff = self._sum(other, 1, -1)
        return diff.value == 0

    def __lt__(self, other):
        diff = self._sum(other, 1, -1)
        return diff.value < 0

    def __mul__(self, other):
        if not isPhysicalQuantity(other):
            return self.__class__(self.value * other, self.unit, \
                prec = self.prec)
        value = self.value * other.value
        unit = self.unit * other.unit
        sig_fig = min(self.prec, other.prec)
        if unit.is_dimensionless:
            return value * unit.factor
        else:
            return self.__class__(value, unit, prec=sig_fig)

    __rmul__ = __mul__

    def __div__(self, other):
        if not isPhysicalQuantity(other):
            return self.__class__(self.value / other, self.unit, \
                prec = self.prec)
        value = self.value / other.value
        unit = self.unit / other.unit
        sig_fig = min(self.prec, other.prec)
        if unit.is_dimensionless:
            return value * unit.factor
        else:
            return self.__class__(value, unit, prec = sig_fig)

    def __rdiv__(self, other):
        if not isPhysicalQuantity(other):
            return self.__class__(other / self.value, pow(self.unit, -1), \
                prec = self.prec)
        value = other.value / self.value
        unit = other.unit / self.unit
        sig_fig = min(self.prec, other.prec)
        if unit.is_dimensionless:
            return value * unit.factor
        else:
            return self.__class__(value, unit, prec = sig_fig)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __pow__(self, other):
        if isPhysicalQuantity(other):
            raise UnitError('Exponents must be dimensionless')
        return self.__class__(pow(self.value, other), pow(self.unit, other), \
            prec = self.prec)

    def __rpow__(self, other):
        raise UnitError('Exponents must be dimensionless')

    def __abs__(self):
        return self.__class__(abs(self.value), self.unit, prec = self.prec)

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.value, self.unit, prec = self.prec)

    def __nonzero__(self):
        return self.value != 0

    def __format__(self, *args, **kw):
        return "{1:{0}} {2}".format(args[0],self.value, self.unit)

    def convert(self, unit):
        """Change the unit and adjust the value such that the combination is
        equivalent to the original one. The new unit must be compatible with the
        previous unit of the object.
        """
        unit = _findUnit(unit)
        self.value = _convertValue(self.value, self.unit, unit)
        self.unit = unit

    def _round(self, x):
        if np.greater(x, 0.):
            return np.floor(x)
        else:
            return np.ceil(x)

    def to(self, *units):
        """Express the quantity in different units. If one unit is specified, a
        new PhysicalQuantity object is returned that expresses the quantity in
        that unit. If several units are specified, the return value is a tuple
        of PhysicalObject instances with with one element per unit such that the
        sum of all quantities in the tuple equals the the original quantity and
        all the values except for the last one are integers. This is used to
        convert to irregular unit systems like hour/minute/second.
        """
        units = [_findUnit(x) for x in units]
        if len(units) == 1:
            unit = units[0]
            value = _convertValue(self.value, self.unit, unit)
            return self.__class__(value, unit)
        else:
            units.sort()
            result = []
            value = self.value
            unit = self.unit
            for i in range(len(units)-1,-1,-1):
                value = value*unit.conversion_factor_to(units[i])
                if i == 0:
                    rounded = value
                else:
                    rounded = self._round(value)
                result.append(self.__class__(rounded, units[i]))
                value = value - rounded
                unit = units[i]
            return tuple(result)

    @staticmethod
    def any_to(qty, unit):
        if not isPhysicalQuantity(qty):
            qty = PhysicalQuantity(qty, 'rad')
        return qty.to(unit)

    @property
    def base(self):
        """Returns the same quantity converted to base units."""
        new_value = self.value * self.unit.factor
        num = ''
        denom = ''
        for i in range(9):
            unit = _base_names[i]
            power = self.unit.powers[i]
            if power < 0:
                denom += '/' + unit
                if power < -1:
                    denom += '**' + str(-power)
            elif power > 0:
                num += '*' + unit
                if power > 1:
                    num += '**' + str(power)
        if len(num) == 0:
            num = '1'
            if denom == '':
                return new_value
        else:
            num = num[1:]
        return self.__class__(new_value, num + denom)

    @property
    def cgs(self):
        """Returns the same quantity converted to cgs units."""
        new_value = self.value * self.unit.factor
        num = ''
        denom = ''
        for i in range(9):

            unit_name = _base_names[i]
            cgs_name = _cgs_names[i]
            power = self.unit.powers[i]

            conversion_factor = Q('1 '+unit_name).to(cgs_name).value
            new_value *= conversion_factor**power

            if power < 0:
                denom += '/' + cgs_name
                if power < -1:
                    denom += '**' + str(-power)
            elif power > 0:
                num += '*' + cgs_name
                if power > 1:
                    num += '**' + str(power)
        if len(num) == 0:
            num = '1'
            if denom == '':
                return new_value
        else:
            num = num[1:]

        return self.__class__(new_value, num + denom)

    # implementations of special functions, used by numpy ufuncs

    def sqrt(self):
        return pow(self, 0.5)

    def sin(self):
        if self.unit.is_angle:
            return np.sin(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of sin must be an angle')

    def cos(self):
        if self.unit.is_angle:
            return np.cos(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of cos must be an angle')

    def tan(self):
        if self.unit.is_angle:
            return np.tan(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of tan must be an angle')


Q = PhysicalQuantity


# SI unit definitions

_base_names = ['m', 'kg', 's', 'A', 'K', 'mol', 'cd', 'rad', 'sr']

_base_units = [
    ('m',   PhysicalUnit('m',   1.,    [1,0,0,0,0,0,0,0,0]), r'\meter'    ),
    ('g',   PhysicalUnit('g',   0.001, [0,1,0,0,0,0,0,0,0]), r'\gram'     ),
    ('s',   PhysicalUnit('s',   1.,    [0,0,1,0,0,0,0,0,0]), r'\second'   ),
    ('A',   PhysicalUnit('A',   1.,    [0,0,0,1,0,0,0,0,0]), r'\ampere'   ),
    ('K',   PhysicalUnit('K',   1.,    [0,0,0,0,1,0,0,0,0]), r'\kelvin'   ),
    ('mol', PhysicalUnit('mol', 1.,    [0,0,0,0,0,1,0,0,0]), r'\mole'     ),
    ('cd',  PhysicalUnit('cd',  1.,    [0,0,0,0,0,0,1,0,0]), r'\candela'  ),
    ('rad', PhysicalUnit('rad', 1.,    [0,0,0,0,0,0,0,1,0]), r'\radian'   ),
    ('sr',  PhysicalUnit('sr',  1.,    [0,0,0,0,0,0,0,0,1]), r'\steradian'),
]

_cgs_names = ['cm', 'g', 's', 'abA', 'K', 'mol', 'cd', 'rad', 'sr']

_prefixes = [
    ('Y',  1.e24,  r'\yotta'),
    ('Z',  1.e21,  r'\zetta'),
    ('E',  1.e18,  r'\exa'  ),
    ('P',  1.e15,  r'\peta' ),
    ('T',  1.e12,  r'\tera' ),
    ('G',  1.e9,   r'\giga' ),
    ('M',  1.e6,   r'\mega' ),
    ('k',  1.e3,   r'\kilo' ),
    ('h',  1.e2,   r'\hecto'),
    ('da', 1.e1,   r'\deca' ),
    ('d',  1.e-1,  r'\deci' ),
    ('c',  1.e-2,  r'\centi'),
    ('m',  1.e-3,  r'\milli'),
    ('mu', 1.e-6,  r'\micro'),
    ('n',  1.e-9,  r'\nano' ),
    ('p',  1.e-12, r'\pico' ),
    ('f',  1.e-15, r'\femto'),
    ('a',  1.e-18, r'\atto' ),
    ('z',  1.e-21, r'\zepto'),
    ('y',  1.e-24, r'\yocto'),
]

_unit_table = {}

_latex_table = {}

for unit in _base_units:
    _unit_table[unit[0]] = unit[1]
    _latex_table[unit[0]] = unit[2]

#TODO : modifier les fonctions _addUnit, _addPrefixed pour prendre en compte les commandes LaTeX des unités de siunitx

def _addUnit(name, unit, latex_code, comment=''):
    if name in _unit_table:
        raise KeyError('Unit ' + name + ' already defined')
    if type(unit) == type(''):
        unit = eval(unit,{}, _unit_table)
    unit.set_name(name)
    _unit_table[name] = unit
    _latex_table[name] = latex_code

def _addPrefixed(unit : str):
    _prefixed_names = []
    for prefix in _prefixes:
        name = prefix[0] + unit
        _addUnit(name, prefix[1]*_unit_table[unit], prefix[2] + _latex_table[unit])
        _prefixed_names.append(name)


# SI derived units; these automatically get prefixes
_unit_table['kg'] = PhysicalUnit('kg',   1., [0,1,0,0,0,0,0,0,0])

_addUnit('Hz' , '1/s'      , r'\hertz'    , 'Hertz'    )
_addUnit('N'  , 'm*kg/s**2', r'\newton'   , 'Newton'   )
_addUnit('Pa' , 'N/m**2'   , r'\pascal'   , 'Pascal'   )
_addUnit('J'  , 'N*m'      , r'\joule'    , 'Joule'    )
_addUnit('W'  , 'J/s'      , r'\watt'     , 'Watt'     )
_addUnit('C'  , 's*A'      , r'\coulomb'  , 'Coulomb'  )
_addUnit('V'  , 'W/A'      , r'\volt'     , 'Volt'     )
_addUnit('F'  , 'C/V'      , r'\farad'    , 'Farad'    )
_addUnit('ohm', 'V/A'      , r'\ohm'      , 'Ohm'      )
_addUnit('S'  , 'A/V'      , r'\siemens'  , 'Siemens'  )
_addUnit('Wb' , 'V*s'      , r'\weber'    , 'Weber'    )
_addUnit('T'  , 'Wb/m**2'  , r'\tesla'    , 'Tesla'    )
_addUnit('H'  , 'Wb/A'     , r'\henry'    , 'Henry'    )
_addUnit('lm' , 'cd*sr'    , r'\lumen'    , 'Lumen'    )
_addUnit('lx' , 'lm/m**2'  , r'\lux'      , 'Lux'      )
_addUnit('Bq' , '1/s'      , r'\becquerel', 'Becquerel')
_addUnit('Gy' , 'J/kg'     , r'\gray'     , 'Gray'     )
_addUnit('Sv' , 'J/kg'     , r'\sievert'  , 'Sievert'  )
_addUnit('kat', 'mol/s'    , r'\katal'    , 'Katal'    )

del _unit_table['kg']

for unit in list(_unit_table):
    _addPrefixed(unit)

# Fundamental constants, as far as needed to define other units
_unit_table['pi'] = np.pi
_addUnit('c0'     , '299792458.*m/s'    , '', 'speed of light'        )
_addUnit('mu0'    , '4.e-7*pi*N/A**2'   , '', 'permeability of vacuum')
_addUnit('eps0'   , '1/mu0/c0**2'       , '', 'permittivity of vacuum')
_addUnit('hplanck', '6.62606957e-34*J*s', '', 'Planck constant'       )
_addUnit('hbar'   , 'hplanck/(2*pi)'    , '', 'Planck constant / 2pi' )
_addUnit('e0'     , '1.602176565e-19*C' , '', 'elementary charge'     )
_addUnit('me'     , '9.10938291e-31*kg' , '', 'electron mass'         )
_addUnit('kb'     , '1.3806488e-23*J/K' , '', 'Boltzmann constant'    )

# Time units
_addUnit('min'      , '60*s'     , r'\minute', 'minute')
_addUnit('h'        , '60*min'   , r'\hour'  , 'hour'  )
_addUnit('d'        , '24*h'     , r'\day'   , 'day'   )
#_addUnit('wk'       , '7*d'      , r'WEEK', 'week')
#_addUnit('yr'       , '365.25*d' , r'YEAR', 'year')
#_addPrefixed('yr')
#_addUnit('fortnight', '1209600*s', r'', '14 days')

# Length units
_addUnit('inch'   , '2.54*cm'                   , r'in', 'inch')
_addUnit('ft'     , '12*inch'                   , r'ft', 'foot')
_addUnit('yd'     , '3*ft'                      , r'yd', 'yard')
_addUnit('mi'     , '5280.*ft'                  , r'mi', '(British) mile')
_addUnit('nmi'    , '1852.*m'                   , r'\nauticalmile', 'Nautical mile')
_addUnit('Ang'    , '1.e-10*m'                  , r'\angstrom', 'Angstrom')
_addUnit('AA'     , '1.e-10*m'                  , r'\angstrom', 'Angstrom')
_addUnit('Bohr'   , '4*pi*eps0*hbar**2/me/e0**2', r'\bohr', 'Bohr radius')
_addUnit('au'     , '149597870691*m'            , r'\astronomicalunit', 'astronomical unit')

# Area units
_addUnit('ha'   , '10000*m**2', r'\hectare', 'hectare')
_addUnit('acres', 'mi**2/640' , r'acres'   , 'acre')
_addUnit('b'    , '1.e-28*m'  , r'\barn'   , 'barn')

# Volume units
_addUnit('L'    , 'dm**3'           , r'\liter', 'liter')
_addPrefixed('L')

# Mass units
_addUnit('t'  , '1000*kg'           , r'\tonne'         , 'Metric ton')
_addUnit('amu', '1.660538921e-27*kg', r'\atomicmassunit', 'atomic mass units')
_addUnit('Da' , '1*amu'             , r'\dalton'        , 'Dalton')

# Energy units
_addUnit('eV'     , 'e0*V'    , r'\electronvolt', 'electron volt')
_addUnit('cal'    , '4.184*J' , r'cal'          , 'thermochemical calorie')
_addUnit('kcal'   , '1000*cal', r'kcal'         , 'thermochemical kilocalorie')

_addPrefixed('eV')

# Pressure units
_addUnit('bar', '1.e5*Pa'         , r'\bar'      , 'bar (cgs unit)')
_addUnit('mbar', '1.e2*Pa'        , r'\milli\bar', 'millibar')
_addUnit('kbar', '1.e8*Pa'        , r'\kilo\bar' , 'kilobar')
_addUnit('atm', '101325.*Pa'      , r'atm'       , 'standard atmosphere')
_addUnit('torr', 'atm/760'        , r'\mmHg'     , 'torr = mm of mercury')
_addUnit('psi', '6894.75729317*Pa', r'psi'       , 'pounds per square inch')

# Angle units
_addUnit('deg'   , 'pi*rad/180'     , r'\degree'   , 'degrees')
_addUnit('arcmin', 'pi*rad/180/60'  , r'\arcminute', 'minutes of arc')
_addUnit('arcsec', 'pi*rad/180/3600', r'\arcsecond', 'seconds of arc')
_unit_table['cycles'] = 2*np.pi

# Temperature units -- can't use the 'eval' trick that _addUnit provides
# for degC and degF because you can't add units
kelvin = _findUnit('K')
_addUnit('degR', '(5./9.)*K', 'degrees Rankine')
_addUnit('degC', PhysicalUnit(None, 1.0, kelvin.powers, 273.15), r'\degreeCelsius',
         'degrees Celcius')
_addUnit('degF', PhysicalUnit(None, 5./9., kelvin.powers, 459.67), r'\degree F',
         'degree Fahrenheit')
del kelvin

# Important physical constants
_constants = [
    ('pi'    , np.pi),
    ('e'     , np.e),
    ('c0'    , Q('299792458. m/s')),
    ('mu0'   , Q('4.e-7 pi*N/A**2').base),
    ('eps0'  , Q('1 1/mu0/c0**2').base),
    ('Grav'  , Q('6.67384e-11 m**3/kg/s**2')),
    ('hpl'   , Q('6.62606957e-34 J*s')),
    ('hbar'  , Q('6.62606957e-34 J*s')/(2*pi)),
    ('e0'    , Q('1.602176565e-19 C')),
    ('me'    , Q('9.10938291e-31 kg')),
    ('mp'    , Q('1.672621777e-27 kg')),
    ('mn'    , Q('1.674927351e-27 kg')),
    ('NA'    , Q('6.02214129e23 1/mol')),
    ('kb'    , Q('1.3806488e-23 J/K')),
    ('g0'    , Q('9.80665 m/s**2')),
    ('R'     , Q('8.3144621 J/mol/K')),
    ('alpha' , 7.2973525698e-3),
    ('Ry'    , Q('10973731.568539 1/m')),
    ('mu_n'  , Q('-0.96623647e-26 J/T')),
    ('gamma' , Q('183.247179 MHz/T')),
    ('h0'    , 0.704),  # WMAP-7 + BAO constraint
    ('sigmaT', Q('6.652453e-29 m**2')),
]
