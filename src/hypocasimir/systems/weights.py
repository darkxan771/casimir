from __future__ import annotations

from math import prod
from typing import Callable

from sympy import Rational

# DOMINANT WEIGHTS
is_dominant: dict[str, Callable[[int, list], bool]] = {}
is_dominant["A"] = lambda n, L: all(L[i] >= L[i + 1] for i in range(n))
is_dominant["B"] = lambda n, L: (
    all(L[i] >= L[i + 1] for i in range(n - 1)) and L[n - 1] >= 0
)
is_dominant["C"] = is_dominant["B"]
is_dominant["D"] = lambda n, L: (
    all(L[i] >= L[i + 1] for i in range(n - 2)) and L[n - 2] >= abs(L[n - 1])
)


# INTEGRAL WEIGHTS
is_integral: dict[str, Callable[[int, list], bool]] = {}
is_integral["A"] = lambda n, L: all((L[i] - L[-1]).is_integer for i in range(n))
is_integral["B"] = lambda n, L: all(L[i].is_integer for i in range(n))
is_integral["C"] = is_integral["B"]
is_integral["D"] = is_integral["B"]


# LABELS
compute_label: dict[str, Callable[[list], list]] = {}
compute_label["A"] = lambda L: [L[i] - L[-1] for i in range(len(L) - 1)]
compute_label["B"] = lambda L: L
compute_label["C"] = compute_label["B"]
compute_label["D"] = compute_label["B"]


class Weight:
    """
    A class for the manipulation of weights in classical weight spaces
    (type A/B/C/D).
    """

    def __init__(self, system: str = "A", L: list = []):
        self.system = system
        self.coordinates = list(Rational(x) for x in L)
        self.rank = len(self.coordinates)
        if self.system == "A":
            self.rank -= 1

    def __repr__(self):
        return f"{self.system}{self.coordinates}"

    def __add__(self, other) -> Weight:
        if not (self.system == other.system and self.rank == other.rank):
            raise ValueError("The two weights cannot be added.")
        r = len(self.coordinates)
        L = [self.coordinates[i] + other.coordinates[i] for i in range(r)]
        return Weight(self.system, L)

    def __neg__(self) -> Weight:
        return Weight(self.system, [-x for x in self.coordinates])

    @property
    def space(self):
        """
        Returns the weight space containing the weight.
        """
        from .spaces import WeightSpace

        return WeightSpace(self.system, self.rank)

    @property
    def is_dominant(self) -> bool:
        """
        Checks whether the weight belongs to the Weyl chamber determined
        by the choice of positive roots in the weight space.
        """
        return is_dominant[self.system](self.rank, self.coordinates)

    @property
    def is_integral(self) -> bool:
        """
        Checks whether the weight belongs to the lattice spanned by
        the weights of the representations of the compact Lie group
        associated to the weight space. Beware that in the case of
        special orthogonal groups (type B and D), this is strictly
        stronger than being an integral linear combination of the
        fundamental weights.
        """
        return is_integral[self.system](self.rank, self.coordinates)

    def _check_hw(self) -> None:
        nhw = "This is not the highest weight of an irreducible representation"
        if not self.is_dominant and self.is_integral:
            raise ValueError(nhw)
        return None

    def label(self):
        """
        Checks whether the weight is dominant and integral, and returns
        the corresponding integer partition or signed integer partition
        (in the case of even special orthogonal groups).
        """
        _ = self._check_hw()
        from ..labels import Partition, Signature

        if self.system == "D":
            return Signature(compute_label["D"](self.coordinates))
        else:
            return Partition(compute_label[self.system](self.coordinates))

    def weyl_dimension(self) -> int:
        """
        Checks whether the weight is dominant and integral, and then
        computes the dimension of the corresponding irreducible
        representation, by using the Weyl formula.
        """
        _ = self._check_hw()
        RW = self.space
        rho = RW.rho
        aug_self = self + rho
        Phi = RW.positive_roots()
        num = prod([RW.scalar_product(aug_self, alpha) for alpha in Phi])
        denom = prod([RW.scalar_product(rho, alpha) for alpha in Phi])
        return int(num / denom)

    def casimir(self, coeff: None | int = None):
        """
        Computes the Casimir coefficient of the action of the elliptic
        Laplace-Beltrami operator on the irreducible representation
        corresponding to the weight.
        """
        _ = self._check_hw()
        RW = self.space
        rho = RW.rho
        rho2 = rho + rho
        if coeff is None:
            coe = RW.killing_coefficient
        else:
            coe = coeff
        return RW.scalar_product(-self, self + rho2, coe)
