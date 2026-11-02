"""\
=====================================================================
Shared LAMMPS pieces --- :mod:`polyphys.simulate.core`
=====================================================================

The :mod:`~polyphys.simulate.core` module holds what the LAMMPS log, dump,
and data-file modules of :mod:`polyphys.simulate` share: the simulation box,
:class:`Box`, and :class:`LammpsFormatError`, the base class of their
errors.

.. versionadded:: 0.6.0

The box
=======

A LAMMPS box has an origin :math:`\\mathbf{o} = (x_{lo}, y_{lo}, z_{lo})`,
lengths :math:`l_x = x_{hi} - x_{lo}`, :math:`l_y`, and :math:`l_z`, and, when
it is a restricted triclinic box, the tilt factors :math:`xy`, :math:`xz`, and
:math:`yz`. Its edge vectors are

.. math::

    \\mathbf{A} = (l_x, 0, 0), \\quad \\mathbf{B} = (xy, l_y, 0), \\quad
    \\mathbf{C} = (xz, yz, l_z).

Scaled coordinates :math:`(s_a, s_b, s_c)`, the fractions of the edge
vectors that dumps write as ``xs ys zs``, map to Cartesian coordinates
:math:`\\mathbf{r} = \\mathbf{o} + s_a \\mathbf{A} + s_b \\mathbf{B} + s_c
\\mathbf{C}`, and the image flags :math:`(i_x, i_y, i_z)` unwrap a wrapped
position to :math:`\\mathbf{r} + i_x \\mathbf{A} + i_y \\mathbf{B} + i_z
\\mathbf{C}`. In an orthogonal box the tilt factors are zero. Lengths are in
the length unit of the run; nothing is converted.

References
==========

The LAMMPS documentation of `Howto_triclinic
<https://docs.lammps.org/Howto_triclinic.html>`_, and these lines of the
LAMMPS source, release ``stable_30Sep2026_update1``:

- The edge vectors: `domain.cpp, lines 256 to 270
  <https://github.com/lammps/lammps/blob/098ffdefa525cabc82272925de02e9884f58c657/src/domain.cpp#L256-L270>`_.
- Scaled to Cartesian coordinates: `domain.cpp, lines 2336 to 2345
  <https://github.com/lammps/lammps/blob/098ffdefa525cabc82272925de02e9884f58c657/src/domain.cpp#L2336-L2345>`_
  (``Domain::lamda2x``).
- Unwrapping by image flags: `domain.cpp, lines 1919 to 1934
  <https://github.com/lammps/lammps/blob/098ffdefa525cabc82272925de02e9884f58c657/src/domain.cpp#L1919-L1934>`_
  (``Domain::unmap``).

Dependencies
============
- Python standard library: `dataclasses`, `gzip`, `pathlib`
- Third-party: `numpy`
"""

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Any

import numpy as np
import numpy.typing as npt

__all__ = ["Box", "LammpsFormatError"]

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]


class LammpsFormatError(ValueError):
    """A LAMMPS file is cut off or does not follow its format.

    The base class of :class:`~polyphys.simulate.log.BrokenLogError`,
    :class:`~polyphys.simulate.dump.BrokenDumpError`, and
    :class:`~polyphys.simulate.data.BrokenDataError`, so one ``except``
    clause catches a broken file of any of the three kinds. It is a
    :class:`ValueError`.

    .. versionadded:: 0.6.0
    """


@dataclass(frozen=True, eq=False)
class Box:
    r"""A LAMMPS simulation box, orthogonal or restricted triclinic.

    .. versionadded:: 0.6.0

    Attributes
    ----------
    lo, hi : numpy.ndarray
        Shape ``(3,)``, float64: :math:`(x_{lo}, y_{lo}, z_{lo})` and
        :math:`(x_{hi}, y_{hi}, z_{hi})` of the box itself, with `lo` below
        `hi`, in the length unit of the run. Read-only.
    tilt : numpy.ndarray
        Shape ``(3,)``, float64: the tilt factors :math:`(xy, xz, yz)`;
        zeros by default. Read-only.
    boundary : tuple of str or None
        The three boundary flags, such as ``('pp', 'pp', 'ff')``, as a dump
        gives them; None when the file does not say, as in a data file.
    triclinic : bool
        True for a restricted triclinic box, whose tilt factors LAMMPS
        stores even when they are zero.

    Raises
    ------
    ValueError
        If `lo`, `hi`, or `tilt` does not have shape ``(3,)``, `lo` is not
        below `hi`, or an orthogonal box has a non-zero tilt factor.

    Examples
    --------
    >>> box = Box(lo=np.zeros(3), hi=np.full(3, 6.0),
    ...           tilt=np.array([1.0, 0.5, -0.5]), triclinic=True)
    >>> box.lengths
    array([6., 6., 6.])
    >>> box.edges()
    array([[ 6. ,  0. ,  0. ],
           [ 1. ,  6. ,  0. ],
           [ 0.5, -0.5,  6. ]])
    """

    lo: FloatArray
    hi: FloatArray
    tilt: FloatArray = field(default_factory=lambda: np.zeros(3))
    boundary: tuple[str, str, str] | None = None
    triclinic: bool = False

    def __post_init__(self) -> None:
        for name in ("lo", "hi", "tilt"):
            value = _frozen(getattr(self, name), np.float64)
            if value.shape != (3,):
                raise ValueError(
                    f"'{name}' must have shape (3,), got {value.shape}"
                )
            object.__setattr__(self, name, value)
        if not np.all(self.lo < self.hi):
            raise ValueError(f"'lo' {self.lo} must be below 'hi' {self.hi}")
        if not self.triclinic and np.any(self.tilt != 0):
            raise ValueError(
                f"an orthogonal box has zero tilt factors, got {self.tilt}"
            )

    @property
    def lengths(self) -> FloatArray:
        """The box lengths :math:`(l_x, l_y, l_z)`, `hi` minus `lo`."""
        return self.hi - self.lo

    def edges(self) -> FloatArray:
        r"""Return the box edge vectors as the rows of a matrix.

        Returns
        -------
        numpy.ndarray
            Shape ``(3, 3)``, float64: the rows are
            :math:`\mathbf{A} = (l_x, 0, 0)`,
            :math:`\mathbf{B} = (xy, l_y, 0)`, and
            :math:`\mathbf{C} = (xz, yz, l_z)`.
        """
        lx, ly, lz = self.lengths
        xy, xz, yz = self.tilt
        return np.array(
            [[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]], dtype=np.float64
        )

    def to_cartesian(self, scaled: npt.ArrayLike) -> FloatArray:
        r"""Return Cartesian coordinates of scaled coordinates.

        Parameters
        ----------
        scaled : array_like
            Shape ``(n_atoms, 3)``: the fractions :math:`(s_a, s_b, s_c)` of
            the edge vectors, as in the ``xs ys zs`` or ``xsu ysu zsu``
            columns of a dump.

        Returns
        -------
        numpy.ndarray
            Shape ``(n_atoms, 3)``, float64: :math:`\mathbf{o} + s_a
            \mathbf{A} + s_b \mathbf{B} + s_c \mathbf{C}`. Time and memory
            are :math:`O(n_a)` for :math:`n_a` atoms.

        Examples
        --------
        >>> box = Box(lo=np.array([1.0, 0.0, 0.0]),
        ...           hi=np.array([5.0, 4.0, 4.0]))
        >>> box.to_cartesian([[0.5, 0.25, 0.0]])
        array([[3., 1., 0.]])
        """
        cartesian: FloatArray = self.lo + np.asarray(
            scaled, dtype=np.float64) @ self.edges()
        return cartesian

    def unwrap(
        self, positions: npt.ArrayLike, image: npt.ArrayLike
    ) -> FloatArray:
        r"""Return wrapped positions unwrapped by their image flags.

        Parameters
        ----------
        positions : array_like
            Shape ``(n_atoms, 3)``: Cartesian coordinates wrapped into the
            box.
        image : array_like
            Shape ``(n_atoms, 3)``: the integer image flags
            :math:`(i_x, i_y, i_z)`.

        Returns
        -------
        numpy.ndarray
            Shape ``(n_atoms, 3)``, float64: :math:`\mathbf{r} + i_x
            \mathbf{A} + i_y \mathbf{B} + i_z \mathbf{C}`. Time and memory
            are :math:`O(n_a)`.

        Examples
        --------
        >>> box = Box(lo=np.zeros(3), hi=np.full(3, 10.0))
        >>> box.unwrap([[9.5, 5.0, 5.0]], [[-1, 0, 0]])
        array([[-0.5,  5. ,  5. ]])
        """
        unwrapped: FloatArray = np.asarray(positions, dtype=np.float64) + (
            np.asarray(image, dtype=np.float64) @ self.edges())
        return unwrapped


def _frozen(value: Any, dtype: type[np.generic]) -> npt.NDArray[Any]:
    """Return `value` as a read-only array copy of type `dtype`."""
    array = np.array(value, dtype=dtype)
    array.flags.writeable = False
    return array


def _open_text(path: Path, errors: str = "strict") -> IO[str]:
    """Open a file as UTF-8 text, through gzip when its name ends in .gz."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors=errors)
    return path.open("r", encoding="utf-8", errors=errors)
