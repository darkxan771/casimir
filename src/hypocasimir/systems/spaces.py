from typing import Callable, Container

from sympy import Rational

from ..lie import LieGroup
from .named import (
    long_root_B,
    long_root_C,
    root_A,
    root_D,
    short_root_B,
    short_root_C,
    weight_A,
    weight_B,
    weight_C,
    weight_D,
)
from .weights import Weight

# SIMPLE ROOTS
simple_roots: dict[str, Callable[[int], list]] = {}
simple_roots["A"] = lambda n: [root_A(n, i, i + 1) for i in range(n)]
simple_roots["B"] = lambda n: (
    [long_root_B(n, i, i + 1, -1) for i in range(n - 1)]
    + [short_root_B(n, n - 1)]
)
simple_roots["C"] = lambda n: (
    [short_root_C(n, i, i + 1, -1) for i in range(n - 1)]
    + [long_root_C(n, n - 1)]
)
simple_roots["D"] = lambda n: (
    [root_D(n, i, i + 1, -1) for i in range(n - 1)]
    + [root_D(n, n - 2, n - 1, +1)]
)


# POSITIVE ROOTS
positive_roots: dict[str, Callable[[int], list]] = {}
positive_roots["A"] = lambda n: [
    root_A(n, i, j) for i in range(n) for j in range(i + 1, n + 1)
]
positive_roots["B"] = lambda n: (
    [
        long_root_B(n, i, j, sign)
        for i in range(n - 1)
        for j in range(i + 1, n)
        for sign in [-1, +1]
    ]
    + [short_root_B(n, i) for i in range(n)]
)
positive_roots["C"] = lambda n: (
    [
        short_root_C(n, i, j, sign)
        for i in range(n - 1)
        for j in range(i + 1, n)
        for sign in [-1, +1]
    ]
    + [long_root_C(n, i) for i in range(n)]
)
positive_roots["D"] = lambda n: [
    root_D(n, i, j, sign)
    for i in range(n)
    for j in range(i + 1, n)
    for sign in [-1, +1]
]


# ROOT NUMBERS
root_number: dict[str, Callable[[int], int]] = {}
root_number["A"] = lambda n: n * (n + 1)
root_number["B"] = lambda n: 2 * n * n
root_number["C"] = root_number["B"]
root_number["D"] = lambda n: 2 * n * (n - 1)


# FUNDAMENTAL WEIGHTS
fundamental_weights: dict[str, Callable[[int, int], Weight]] = {}
fundamental_weights["A"] = weight_A
fundamental_weights["B"] = weight_B
fundamental_weights["C"] = weight_C
fundamental_weights["D"] = weight_D


# LIE GROUPS
lie_groups: dict[str, Callable[[int], LieGroup]] = {}
lie_groups["A"] = lambda n: LieGroup("SU", n + 1)
lie_groups["B"] = lambda n: LieGroup("SO", 2 * n + 1)
lie_groups["C"] = lambda n: LieGroup("SP", n)
lie_groups["D"] = lambda n: LieGroup("SO", 2 * n)


# RHO
rho_vector: dict[str, Callable[[int], Weight]] = {}
rho_vector["A"] = lambda n: Weight(
    "A", [Rational(n - 2 * i, 2) for i in range(n + 1)]
)
rho_vector["B"] = lambda n: Weight(
    "B", [Rational(2 * n - 1 - 2 * i, 2) for i in range(n)]
)
rho_vector["C"] = lambda n: Weight("C", [n - i for i in range(n)])
rho_vector["D"] = lambda n: Weight("D", [n - 1 - i for i in range(n)])


# KILLING COEFFICIENT
killing_coefficient: dict[str, Callable[[int], int]] = {}
killing_coefficient["A"] = lambda n: 2 * n + 2
killing_coefficient["B"] = lambda n: 2 * n - 1
killing_coefficient["C"] = lambda n: 2 * n + 2
killing_coefficient["D"] = lambda n: 2 * n - 2


class WeightSpace(Container):
    """
    A class for the manipulation of the weight space associated
    to a compact Lie group and its root system.
    """

    def __init__(self, system: str = "A", rank: int = 1):
        self.system = system
        self.rank = rank

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, WeightSpace)
            and self.system == other.system
            and self.rank == other.rank
        )

    def __repr__(self):
        return f"Weight space {self.system}({self.rank})"

    def __contains__(self, obj) -> bool:
        return (
            isinstance(obj, Weight)
            and self.system == obj.system
            and self.rank == obj.rank
        )

    def simple_roots(self) -> list[Weight]:
        """
        Returns the list of simple roots of the weight space.
        """
        return simple_roots[self.system](self.rank)

    def positive_roots(self) -> list[Weight]:
        """
        Returns the list of positive roots of the weight space.
        """
        return positive_roots[self.system](self.rank)

    def number_of_roots(self) -> int:
        """
        The number of roots of the weight space. It is also equal to
        dim G - rank G, where G is the compact Lie group associated
        to the weight space.
        """
        return root_number[self.system](self.rank)

    def fundamental_weights(self) -> list[Weight]:
        """
        Returns the list of fundamental weights of the weight space.
        """
        return [
            fundamental_weights[self.system](self.rank, i)
            for i in range(self.rank)
        ]

    def scalar_product(
        self, w1: Weight, w2: Weight, coeff: int = 1
    ) -> Rational:
        """
        Computes the scalar product of the weights w1 and w2, with respect
        to the bilinear form on the weight space which is dual to the form
        - coeff*tr(XY) on the Lie algebra g of the compact Lie group G
        associated to the weight space.
        """
        if not (w1 in self and w2 in self):
            raise ValueError("The weights do not belong to this weight space")
        res = Rational(1, coeff) * sum(
            [
                w1.coordinates[i] * w2.coordinates[i]
                for i in range(len(w1.coordinates))
            ]
        )
        if self.system == "A":
            return res
        else:
            return Rational(1, 2) * res

    def lie_group(self) -> LieGroup:
        """
        Returns the compact Lie group canonically attached to the
        weight space and its root system.
        """
        return lie_groups[self.system](self.rank)

    @property
    def rho(self) -> Weight:
        """
        Returns the half-sum of positive roots, also equal to the sum of all the
        fundamental weights.
        """
        return rho_vector[self.system](self.rank)

    @property
    def killing_coefficient(self) -> int:
        """
        Returns the coefficient c such that the Killing form of the associated
        compact Lie group G is - c*tr(XY).
        """
        return killing_coefficient[self.system](self.rank)
