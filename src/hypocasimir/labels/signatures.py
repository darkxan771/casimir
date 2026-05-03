from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction
from math import prod
from typing import Callable

import matplotlib.pyplot as plt

from ..lie import LieGroup, SymmetricSpace
from ..systems import Weight
from .labels import Label, general_casimir


class NotLabelError(Exception):
    def __init__(self):
        super().__init__(
            "This is not the label of an irreducible representation of G"
        )


def cht_SO2U(n: int, L: Sequence[int]) -> list[Label]:
    S = Signature(L)
    if len(S) != n or not S.is_signed_partition:
        raise NotLabelError
    return [S]


def cht_GrR(p: int, q: int, L: Sequence[int]) -> list[Label]:
    n = p + q
    S = Signature(L)
    if len(S) != int(n / 2) or not S.is_signed_partition:
        raise NotLabelError
    if p % 2 == 0:
        return [Signature(L[: int(p / 2)]), Signature(L[int(p / 2) :])]
    else:
        from .partitions import Partition

        return [
            Partition(L[: int((p - 1) / 2)]),
            Partition(L[int((p + 1) / 2) :]),
        ]


compute_highest_type_structure: dict[
    str, Callable[[int, Sequence[int]], list[Label]]
] = {}
compute_highest_type_structure["SO2/U"] = cht_SO2U
compute_highest_type_grassmann: dict[
    str, Callable[[int, int, Sequence[int]], list[Label]]
] = {}
compute_highest_type_grassmann["GrR"] = cht_GrR


class Signature(Label):
    """
    A class for the manipulation of signatures. They label
    the irreducible representations of unitary groups and even special
    orthogonal groups. For special unitary groups, odd special orthogonal
    groups and compact symplectic groups, use partitions.
    """

    def __init__(self, L: Sequence[int]):
        if any(L[i] < L[i + 1] for i in range(len(L) - 1)):
            raise ValueError("The terms should be in nonincreasing order.")
        self.terms = list(L)

    @property
    def length(self) -> int:
        """
        The number of terms of the signature.
        """
        return len(self)

    @property
    def is_signed_partition(self) -> bool:
        """
        Checks whether the signature is a signed partition (thus, corresponding
        to an irreducible representation of SO(2n)).
        """
        return all(x >= 0 for x in self.terms[:-1]) and self.terms[-2] >= abs(
            self.terms[-1]
        )

    def highest_weight(self, G: LieGroup) -> Weight:
        """
        If the signature is a signed partition with length n and if G=SO(2n),
        converts the signed partition to the corresponding highest weight
        of the irreducible representation of G.
        """
        if not (
            self.is_signed_partition and G == LieGroup("SO", 2 * len(self))
        ):
            raise NotLabelError
        return Weight("D", self.terms)

    def highest_K_type(self, X: SymmetricSpace) -> list[Label]:
        """
        Given a symmetric space X=G/K, computes the highest K-type of the
        restriction of the irreducible representation of G labelled
        by the signature.
        """
        if X.stype == "structure":
            return compute_highest_type_structure[X.type](X.n, self.terms)
        else:
            return compute_highest_type_grassmann[X.type](X.p, X.q, self.terms)

    def dimension(self, G: LieGroup) -> int:
        """
        Computes the dimension of the irreducible representation of G
        labelled by the signature.
        """
        M = self.terms
        n = len(self)
        if G == LieGroup("U", n):
            num = prod(
                M[i] - M[j] + j - i for i in range(n) for j in range(i + 1, n)
            )
            denom = prod(j - i for i in range(n) for j in range(i + 1, n))
            return int(num / denom)
        elif G == LieGroup("SO", 2 * n):
            return self.highest_weight(G).weyl_dimension()
        else:
            raise NotLabelError

    def casimir(self, G: LieGroup, coeff: None | int = None) -> Fraction:
        """
        Computes the Casimir coefficient of the action of the elliptic
        Laplace-Beltrami operator on the irreducible representation of
        G labelled by the signature.
        """
        n = len(self)
        M = self.terms
        if G == LieGroup("U", n):
            if coeff is None:
                coe = 2 * n
            else:
                coe = coeff
            num = sum(M[i] * (M[i] + n - 1 - 2 * i) for i in range(n))
            return Fraction(-num, coe)
        elif G == LieGroup("SO", 2 * n):
            return self.highest_weight(G).casimir(coeff)
        else:
            raise NotLabelError

    def hypocasimir(self, X: SymmetricSpace) -> Fraction:
        """
        Computes the maximal Casimir coefficient of the action of the
        hypoelliptic Laplace-Beltrami operator on the irreducible
        representation of G labelled by the signature.
        """
        n = len(self)
        G = X.isometry_group
        if not (G == LieGroup("SO", 2 * n) and self.is_signed_partition):
            raise NotLabelError
        coeff = G.weight_space.killing_coefficient
        cG = self.casimir(G)
        K = X.isotropy_group
        L = self.highest_K_type(X)
        if X.type == "SO2/U":
            coeff *= 2
        cK = general_casimir(L, K, coeff)
        return cG - cK

    def draw(self) -> None:
        """
        Draws the Young diagram of the signature.
        """
        _, ax0 = plt.subplots()
        from .plot import draw_signature_on_ax

        draw_signature_on_ax(self.terms, ax0)
        plt.show()
