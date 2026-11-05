from collections.abc import Callable

import numpy as np
from scipy.linalg import block_diag

from .bases_lie_algebra import bases_lie_algebra

names: dict[str, str] = {}
names["U"] = "Unitary Group"
names["SU"] = "Special Unitary Group"
names["SO"] = "Special Orthogonal Group"
names["SP"] = "Compact Symplectic Group"

group_dimension: dict[str, Callable[[int], int]] = {}
group_dimension["U"] = lambda n: n * n
group_dimension["SU"] = lambda n: n * n - 1
group_dimension["SO"] = lambda n: int(n * (n - 1) / 2)
group_dimension["SP"] = lambda n: n * (2 * n + 1)

group_rank: dict[str, Callable[[int], int]] = {}
group_rank["U"] = lambda n: n
group_rank["SU"] = lambda n: n - 1
group_rank["SO"] = lambda n: int(n / 2)
group_rank["SP"] = lambda n: n

compute_system: dict[str, Callable[[int], str]] = {}
compute_system["SU"] = lambda _: "A"
compute_system["SO"] = lambda n: "B" if n % 2 == 1 else "D"
compute_system["SP"] = lambda _: "C"


def Jmatrix(n: int) -> np.array:
    J = np.array([[0, 1], [-1, 0]])
    return block_diag(*([J] * n))


def involution_sp(M: np.array) -> np.array:
    n = int(M.shape[0] / 2)
    Jn = Jmatrix(n)
    res = -(Jn @ M) @ Jn
    return res.transpose()


defining_involution: dict[str, Callable[[np.array], np.array]] = {}
defining_involution["SU"] = lambda M: (M + 0j).conjugate().transpose()
defining_involution["SO"] = lambda M: M.transpose()
defining_involution["SP"] = involution_sp


class LieGroup:
    """
    A class for the manipulation of classical compact Lie groups
    and their representations.
    """

    def __init__(self, type: str = "SU", size: int = 1):
        if type not in {"U", "SU", "SO", "SP"}:
            raise NotImplementedError
        self.type: str = type
        if size <= 0:
            raise ValueError("The size should be a positive integer")
        self.size: int = size

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, LieGroup)
            and self.type == other.type
            and self.size == other.size
        )

    def __repr__(self):
        return f"{self.type}({self.size})"

    @property
    def description(self) -> str:
        """
        An explicit description of the compact Lie group.
        """
        return f"{names[self.type]} {repr(self)}"

    @property
    def dimension(self) -> int:
        """
        Returns the real dimension of the Lie group.
        """
        return group_dimension[self.type](self.size)

    @property
    def rank(self) -> int:
        """
        Returns the dimension of a maximal torus inside the Lie group.
        """
        return group_rank[self.type](self.size)

    @property
    def mixing_time(self) -> float:
        """
        Returns the mixing time of the standard Brownian motion on the Lie
        group.
        """
        return float(4 * np.log(self.size))

    @property
    def weight_space(self):
        """
        Returns the weight space associated to the compact Lie group.
        """
        if self.type == "U":
            raise NotImplementedError("The unitary group is not semisimple")
        from ..systems import WeightSpace

        return WeightSpace(compute_system[self.type](self.size), self.rank)

    def basis_lie_algebra(self, normalised: bool = True) -> list[np.array]:
        """
        Returns an orthogonal basis of the Lie algebra of the compact Lie
        group.
        """
        return bases_lie_algebra[self.type](self.size, normalised)

    def defining_involution(self, g: np.array) -> np.array:
        """
        The defining involution of the compact Lie Group.
        """
        return defining_involution[self.type](g)

    @property
    def scalar_product(self) -> Callable[[np.array, np.array], float]:
        """
        The scalar product given by the opposite of the Killing form (or by
        -2n * tr(XY) in the case of unitary groups).
        """
        if self.type == "U":
            c = 2 * self.size
        else:
            c = self.weight_space.killing_coefficient
        return lambda X, Y: -c * (X @ Y).trace().real

    def casimir(self) -> np.array:
        """
        Computes the Casimir element of the group G, which is the sum

        C = sum X otimes X

        in the tensor product g otimes g of the Lie algebra g of G by itself.
        The sum runs over an orthonormal basis of elements X of g.
        """
        b = self.basis_lie_algebra()
        return sum([np.kron(X, X) for X in b])
