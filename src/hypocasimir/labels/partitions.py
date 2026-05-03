from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction
from typing import Callable

import matplotlib.pyplot as plt

from ..lie import LieGroup, SymmetricSpace
from ..systems import Weight
from .labels import Label, TooLongError, complete_list, general_casimir
from .signatures import NotLabelError, Signature

convert_highest_weight: dict[str, Callable[[int, list[int]], Weight]] = {}
convert_highest_weight["A"] = lambda n, L: Weight(
    "A", [x - Fraction(sum(L), n + 1) for x in complete_list(n, L) + [0]]
)
convert_highest_weight["B"] = lambda n, L: Weight("B", complete_list(n, L))
convert_highest_weight["C"] = lambda n, L: Weight("C", complete_list(n, L))


def cht_SUSO(n: int, L: Sequence[int]) -> list[Label]:
    if len(L) >= n:
        raise TooLongError
    M = complete_list(n, L)
    if n % 2 == 1:
        return [
            Partition([M[i] - M[n - 1 - i] for i in range(int((n - 1) / 2))])
        ]
    else:
        return [Signature([M[i] - M[n - 1 - i] for i in range(int(n / 2))])]


def cht_SU2SP(n: int, L: Sequence[int]) -> list[Label]:
    if len(L) >= 2 * n:
        raise TooLongError
    M = complete_list(2 * n, L)
    return [Partition([M[i] - M[2 * n - 1 - i] for i in range(n)])]


def cht_SPU(n: int, L: Sequence[int]) -> list[Label]:
    return [Signature(complete_list(n, L))]


def cht_GrR(p: int, q: int, L: Sequence[int]) -> list[Label]:
    M = complete_list(int((p + q - 1) / 2), L)
    if p % 2 == 0:
        return [Signature(M[: int(p / 2)]), Partition(M[int(p / 2) :])]
    else:
        return [
            Partition(M[: int((p - 1) / 2)]),
            Signature(M[int((p - 1) / 2) :]),
        ]


def cht_GrC(p: int, q: int, L: Sequence[int]) -> list[Label]:
    M = complete_list(p + q - 1, L)
    return [Signature(M[:p]), Partition(M[p:])]


def cht_GrH(p: int, q: int, L: Sequence[int]) -> list[Label]:
    M = complete_list(p + q, L)
    return [Partition(M[:p]), Partition(M[p:])]


compute_highest_type_structure: dict[
    str, Callable[[int, Sequence[int]], list[Label]]
] = {}
compute_highest_type_structure["SU/SO"] = cht_SUSO
compute_highest_type_structure["SU2/SP"] = cht_SU2SP
compute_highest_type_structure["SP/U"] = cht_SPU

compute_highest_type_grassmann: dict[
    str, Callable[[int, int, Sequence[int]], list[Label]]
] = {}
compute_highest_type_grassmann["GrR"] = cht_GrR
compute_highest_type_grassmann["GrC"] = cht_GrC
compute_highest_type_grassmann["GrH"] = cht_GrH


class Partition(Label):
    """
    A class for the manipulation of integer partitions. They label
    the irreducible representations of special unitary groups, odd
    special orthogonal groups and compact symplectic groups. For
    unitary groups and even special orthogonal groups, use signatures.
    """

    def __init__(self, L: Sequence[int]):
        if any(x < 0 for x in L):
            raise ValueError("All parts should be positive integers.")
        if any(L[i] < L[i + 1] for i in range(len(L) - 1)):
            raise ValueError("The parts should be in nonincreasing order.")
        self.terms: list[int] = list(filter((0).__ne__, L))

    def highest_weight(self, G: LieGroup) -> Weight:
        """
        Converts the integer partition to a highest weight of an irreducible
        representation of a compact Lie group.
        """
        if G.type == "U" or (G.type == "SO" and G.size % 2 == 0):
            raise NotLabelError
        RW = G.weight_space
        return convert_highest_weight[RW.system](RW.rank, self.terms)

    def highest_K_type(self, X: SymmetricSpace) -> list[Label]:
        """
        Given a symmetric space X=G/K, computes the highest K-type of the
        restriction of the irreducible representation of G labelled
        by the integer partition.
        """
        if X.stype == "structure":
            return compute_highest_type_structure[X.type](X.n, self.terms)
        else:
            return compute_highest_type_grassmann[X.type](X.p, X.q, self.terms)

    def dimension(self, G: LieGroup) -> int:
        """
        Computes the dimension of the irreducible representation of G
        labelled by the integer partition.
        """
        return self.highest_weight(G).weyl_dimension()

    def casimir(self, G: LieGroup, coeff: None | int = None) -> Fraction:
        """
        Computes the Casimir coefficient of the action of the elliptic
        Laplace-Beltrami operator on the irreducible representation of
        G labelled by the integer partition.
        """
        return self.highest_weight(G).casimir(coeff)

    def hypocasimir(self, X: SymmetricSpace) -> Fraction:
        """
        Computes the maximal Casimir coefficient of the action of the
        hypoelliptic Laplace-Beltrami operator on the irreducible
        representation of G labelled by the integer partition.
        """
        G = X.isometry_group
        L = self.highest_K_type(X)
        coeff = G.weight_space.killing_coefficient
        cG = self.casimir(G)
        if X.type == "GrC":
            n, p, q = X.n, X.p, X.q
            K1, K2 = LieGroup("SU", p), LieGroup("SU", q)
            L1 = Partition([a - L[0].terms[-1] for a in L[0].terms[:-1]])
            L2 = L[1]
            cK = general_casimir([L1, L2], [K1, K2], coeff)
            cK -= Fraction(
                (q * L[0].size - p * L2.size) ** 2, 2 * n * n * p * q
            )
        else:
            K = X.isotropy_group
            if X.type == "SP/U":
                coeff *= 2
            cK = general_casimir(L, K, coeff)
        return cG - cK

    def draw(self, style: str = "french") -> None:
        """
        Draws the Young diagram of the integer partition.
        Available options (styles): "french", "english", "russian".
        """
        _, ax0 = plt.subplots()
        from .plot import draw_partition_on_ax

        draw_partition_on_ax(self.terms, ax0, style)
        plt.show()
