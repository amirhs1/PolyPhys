"""\
==========================================================
:mod:`polyphys.tests.strategies`
==========================================================

Hypothesis strategies shared by the package's tests and the study tests in
``projects/``.

:func:`lineage_names` draws an artifact name of one lineage of a parser
class from the parser's own declaration: one value per keyword of the
lineage's pattern, written in the text a parser rebuilds a short name with,
so that parsing the name and rebuilding it is the identity.

Examples
========
>>> from hypothesis import find
>>> from polyphys.tests.strategies import lineage_names
>>> from polyphys.tests.synthetic import ChainCyl
>>> name, values = find(lineage_names(ChainCyl, 'space'), lambda _: True)
>>> sorted(values)
['D', 'N', 'ac']
"""
from decimal import Decimal

from hypothesis import assume
from hypothesis import strategies as st

from polyphys.manage.parser import KeywordSpec, ParserBase
from polyphys.manage.types import LineageT

__all__ = ['keyword_text', 'lineage_names', 'stored_values']


def keyword_text(value: int | Decimal) -> str:
    """
    Return the text of a keyword value as a parser rebuilds it in a short
    name: an `int` unchanged, and a float that is a whole number as an
    integer, any other as its :func:`repr`.

    Examples
    --------
    >>> from decimal import Decimal
    >>> keyword_text(200), keyword_text(Decimal('2.000')), \\
    ...     keyword_text(Decimal('1.250'))
    ('200', '2', '1.25')
    """
    if isinstance(value, int):
        return str(value)
    number = float(value)
    return str(int(number)) if number.is_integer() else repr(number)


def _values(spec: KeywordSpec) -> st.SearchStrategy[int | Decimal]:
    """The values of one keyword: integers from 0 to 100000, or decimals
    with up to three places from 0.001 to 1000."""
    if spec.kind is int:
        return st.integers(min_value=0, max_value=100_000)
    return st.decimals(min_value=Decimal('0.001'), max_value=Decimal('1000'),
                       places=3, allow_nan=False, allow_infinity=False)


@st.composite
def lineage_names(
    draw: st.DrawFn, parser_class: type[ParserBase], lineage: LineageT
) -> tuple[str, dict[str, int | Decimal]]:
    """
    Draw a name of `lineage` for `parser_class` and the value of each of
    its keywords.

    The keywords follow the lineage's pattern, each followed by the text of
    its value (:func:`keyword_text`); in a ``segment`` name the segment
    keyword ``j`` follows a dot, as the probes write it. In a cylinder,
    a name is kept only when the cylinder size exceeds every species' size,
    so that each species' centers have room and the parser's bulk values
    exist (:func:`~polyphys.manage.utils.number_density_cylinder`).

    Returns
    -------
    tuple of (str, dict)
        The name, and each keyword of the pattern mapped to its value.
    """
    keywords = parser_class._keywords
    patterns = parser_class._patterns
    assert keywords is not None and patterns is not None
    values = {keyword: draw(_values(keywords[keyword]))
              for keyword in patterns[lineage]}
    if parser_class._geometry == 'cylindrical':
        stored = {**parser_class._defaults,
                  **stored_values(parser_class, values)}
        if 'dcyl' in stored:
            assume(all(stored.get(species.size, 0.0) < stored['dcyl']
                       for species in parser_class._species))
    parts = [('.' if keyword == 'j' else '') + keyword + keyword_text(value)
             for keyword, value in values.items()]
    return ''.join(parts), values


def stored_values(
    parser_class: type[ParserBase], values: dict[str, int | Decimal]
) -> dict[str, int | float]:
    """
    Map each keyword's attribute to the value a parser stores for it,
    :math:`y = s x + c` (:class:`~polyphys.manage.parser.KeywordSpec`), in
    float arithmetic.

    Examples
    --------
    >>> from decimal import Decimal
    >>> from polyphys.tests.synthetic import ChainCyl
    >>> stored_values(ChainCyl, {'N': 200, 'r': Decimal('10.5')})
    {'nmon': 200, 'dcyl': 20.0}
    """
    keywords = parser_class._keywords
    assert keywords is not None
    stored: dict[str, int | float] = {}
    for keyword, value in values.items():
        spec = keywords[keyword]
        stored[spec.attr] = (value if isinstance(value, int)
                             else spec.scale * float(value) + spec.offset)
    return stored
