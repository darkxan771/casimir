from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction
from typing import Protocol

import numpy as np

from ..lie import LieGroup, SymmetricSpace
from ..systems import Weight


class TooLongError(Exception):
    def __init__(self):
        super().__init__("The length of the list is too large")


def complete_list(n: int, L: Sequence[int]) -> list[int]:
    if len(L) > n:
        raise TooLongError
    return list(L) + [0] * (n - len(L))


class Label(Protocol):
    """
    A generic class for the labelling of irreducible representations
    of compact Lie groups. Signatures are used for irreducibles of
    unitary groups and even special orthogonal groups, while integer
    partitions are used for special unitary groups, compact symplectic
    groups and odd special orthogonal groups.
    """

    terms: list[int]

    def __eq__(self, other: object) -> bool:
        from .partitions import Partition
        from .signatures import Signature

        if isinstance(other, Partition) or isinstance(other, Signature):
            return self.terms == other.terms
        else:
            return False

    def __len__(self) -> int:
        return len(self.terms)

    def __abs__(self) -> int:
        return sum(self.terms)

    def __repr__(self) -> str:
        if self.terms[-1] != 0:
            return repr(self.terms)
        else:
            i = self.terms.index(0)
            if i < len(self.terms) - 1:
                res = repr(self.terms[:i])[:-1]
                res += ", 0, ..]"
            else:
                res = repr(self.terms)
            return res

    def __iter__(self):
        return iter(self.terms)

    def __hash__(self):
        return hash(tuple(self.terms))

    @property
    def length(self) -> int:
        """
        The number of terms of the label.
        """
        return len(self)

    @property
    def size(self) -> int:
        """
        The size of the label (sum of terms).
        """
        return abs(self)

    def highest_weight(self, G: LieGroup) -> Weight: ...

    def highest_K_type(self, X: SymmetricSpace) -> list[Label]: ...

    def dimension(self, G: LieGroup) -> int: ...

    def casimir(self, G: LieGroup, coeff: None | int) -> Fraction: ...

    def hypocasimir(self, X: SymmetricSpace) -> Fraction: ...

    def up_term(self, S: LieGroup | SymmetricSpace) -> float:
        """
        Returns the term of the upper bound on dtv(mu_t, Haar)^2 at cutoff time.
        """
        if isinstance(S, LieGroup):
            return float(
                (self.dimension(S) ** 2)
                * np.exp(self.casimir(S) * S.mixing_time)
            )
        else:
            return float(
                (self.dimension(S.isometry_group) ** 2)
                * np.exp(self.hypocasimir(S) * S.mixing_time)
            )


def general_casimir(L: list[Label], G: list[LieGroup], coeff: int) -> Fraction:
    return Fraction(sum([L[i].casimir(G[i], coeff) for i in range(len(L))]))
