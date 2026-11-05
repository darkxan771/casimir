from collections.abc import Callable

import numpy as np

from .bases_lie_algebra import (
    bases_k_grassmann,
    bases_k_structure,
    bases_m_grassmann,
    bases_m_structure,
)
from .groups import LieGroup

space_dimension_structure: dict[str, Callable[[int], int]] = {}
space_dimension_structure["SU/SO"] = lambda n: int((n - 1) * (n + 2) / 2)
space_dimension_structure["SU2/SP"] = lambda n: int(2 * n * n - 1 - n)
space_dimension_structure["SP/U"] = lambda n: int(n * (n + 1))
space_dimension_structure["SO2/U"] = lambda n: int(n * (n - 1))
space_dimension_grassmann: dict[str, Callable[[int, int], int]] = {}
space_dimension_grassmann["GrR"] = lambda p, q: int(p * q)
space_dimension_grassmann["GrC"] = lambda p, q: int(2 * p * q)
space_dimension_grassmann["GrH"] = lambda p, q: int(4 * p * q)

compute_repr_structure: dict[str, Callable[[int], str]] = {}
compute_repr_structure["SU/SO"] = lambda n: f"SU({n})/SO({n})"
compute_repr_structure["SU2/SP"] = lambda n: f"SU({2 * n})/SP({n})"
compute_repr_structure["SP/U"] = lambda n: f"SP({n})/U({n})"
compute_repr_structure["SO2/U"] = lambda n: f"SO({2 * n})/U({n})"
compute_repr_grassmann: dict[str, Callable[[int, int], str]] = {}
compute_repr_grassmann["GrR"] = lambda p, q: f"SO({p + q})/(SO({p}) x SO({q}))"
compute_repr_grassmann["GrC"] = lambda p, q: f"SU({p + q})/S(U({p}) x U({q}))"
compute_repr_grassmann["GrH"] = lambda p, q: f"SP({p + q})/(SP({p}) x SP({q}))"

compute_field: dict[str, str] = {}
compute_field["SU/SO"] = "real"
compute_field["SU2/SP"] = "quaternionic"
compute_field["SP/U"] = "complex"
compute_field["SO2/U"] = "complex"
compute_field["GrR"] = "real"
compute_field["GrC"] = "complex"
compute_field["GrH"] = "quaternionic"

compute_G: dict[str, Callable[[int], LieGroup]] = {}
compute_G["SU/SO"] = lambda n: LieGroup("SU", n)
compute_G["SU2/SP"] = lambda n: LieGroup("SU", 2 * n)
compute_G["SP/U"] = lambda n: LieGroup("SP", n)
compute_G["SO2/U"] = lambda n: LieGroup("SO", 2 * n)
compute_G["GrR"] = lambda n: LieGroup("SO", n)
compute_G["GrC"] = lambda n: LieGroup("SU", n)
compute_G["GrH"] = lambda n: LieGroup("SP", n)

compute_K_structure: dict[str, Callable[[int], list[LieGroup]]] = {}
compute_K_structure["SU/SO"] = lambda n: [LieGroup("SO", n)]
compute_K_structure["SU2/SP"] = lambda n: [LieGroup("SP", n)]
compute_K_structure["SP/U"] = lambda n: [LieGroup("U", n)]
compute_K_structure["SO2/U"] = lambda n: [LieGroup("U", n)]
compute_K_grassmann: dict[str, Callable[[int, int], list[LieGroup]]] = {}
compute_K_grassmann["GrR"] = lambda p, q: [LieGroup("SO", p), LieGroup("SO", q)]
compute_K_grassmann["GrH"] = lambda p, q: [LieGroup("SP", p), LieGroup("SP", q)]


class SymmetricSpace:
    """
    A class for the manipulation of classical compact symmetric spaces.
    """

    def __init__(self, type: str = "SU/SO", size: int | tuple[int, int] = 2):
        self.type: str = type
        self.stype: str = ""
        self.size: int = 0
        if type in {"SU/SO", "SU2/SP", "SP/U", "SO2/U"} and isinstance(
            size, int
        ):
            self.stype = "structure"
            self.size = size
        elif type in {"GrR", "GrC", "GrH"} and isinstance(size, tuple):
            self.stype = "grassmann"
            self.p: int = max(size)
            self.q: int = min(size)
            self.size = self.p + self.q
        else:
            raise NotImplementedError

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, SymmetricSpace)
            and self.type == other.type
            and self.size == other.size
            and ((self.stype == "structure") or (self.p == other.p))
        )

    def __repr__(self):
        if self.stype == "structure":
            return compute_repr_structure[self.type](self.size)
        else:
            return compute_repr_grassmann[self.type](self.p, self.q)

    @property
    def description(self) -> str:
        """
        An explicit description of the compact symmetric space.
        """
        if self.stype == "structure":
            return (
                f"Space of {compute_field[self.type]} structures {repr(self)}"
            )
        else:
            res = compute_field[self.type].capitalize()
            res += " Grassmannian manifold "
            res += repr(self)
            return res

    @property
    def dimension(self) -> int:
        """
        Returns the real dimension of the symmetric space.
        """
        if self.stype == "structure":
            return space_dimension_structure[self.type](self.size)
        else:
            return space_dimension_grassmann[self.type](self.p, self.q)

    @property
    def mixing_time(self) -> float:
        """
        Returns the mixing time of the hypoelliptic Brownian motion on the
        isometry group of the symmetric space.
        """
        if self.stype == "structure":
            return float(8 * np.log(self.isometry_group.size))
        else:
            return float(
                4 * self.size / self.q * np.log(self.isometry_group.size)
            )

    @property
    def isometry_group(self) -> LieGroup:
        """
        Returns the simple simply connected compact Lie group G such
        that X=G/K.
        """
        return compute_G[self.type](self.size)

    @property
    def isotropy_group(self) -> list[LieGroup]:
        """
        Returns the stabilizer K of a point for the action of the
        isometry group G on X=G/K. Not implemented for complex Grassmannian
        manifolds.
        """
        if self.stype == "structure":
            return compute_K_structure[self.type](self.size)
        else:
            if self.type == "GrC":
                raise NotImplementedError
            return compute_K_grassmann[self.type](self.p, self.q)

    def basis_k(self) -> list[np.array]:
        """
        Returns an orthonormal basis of the Lie subalgebra k of g.
        """
        if self.stype == "structure":
            return bases_k_structure[self.type](self.size)
        else:
            return bases_k_grassmann[self.type](self.p, self.q)

    def basis_m(self) -> list[np.array]:
        """
        Returns an orthonormal basis of the orthogonal complement m of the
        Lie subalgebra k of g.
        """
        if self.stype == "structure":
            return bases_m_structure[self.type](self.size)
        else:
            return bases_m_grassmann[self.type](self.p, self.q)

    def casimir(self) -> np.array:
        """
        Computes the hypoelliptic Casimir element of the symmetric space
        X = G/K, which is the sum

        C = sum V otimes V

        in the tensor product g otimes g of the Lie algebra g of G by itself.
        The sum runs over an orthonormal basis of elements V of m, the
        orthogonal complement of k in g.
        """
        b = self.basis_m()
        return sum([np.kron(X, X) for X in b])

    def completed_casimir(self) -> np.array:
        """
        Computes the completed hypoelliptic Casimir element of the symmetric
        space X = G/K, which is the sum

        C = sum (I/2 otimes V^2 + V otimes V + V^2 otimes I/2)

        in the tensor product g otimes g of the Lie algebra g of G by itself.
        The sum runs over an orthonormal basis of elements V of m, the
        orthogonal complement of k in g.
        """
        b = self.basis_m()
        if self.isometry_group.type == "SP":
            d = 2 * self.isometry_group.size
        else:
            d = self.isometry_group.size
        A = (1 / 2) * sum([np.kron(np.eye(d), X @ X) for X in b])
        B = sum([np.kron(X, X) for X in b])
        C = (1 / 2) * sum([np.kron(X @ X, np.eye(d)) for X in b])
        return A + B + C

    def conj_completed_casimir(self) -> np.array:
        b = self.basis_m()
        if self.isometry_group.type == "SP":
            d = 2 * self.isometry_group.size
        else:
            d = self.isometry_group.size
        A = (1 / 2) * sum(
            [np.kron(np.eye(d), X.conjugate() @ X.conjugate()) for X in b]
        )
        B = sum([np.kron(X, X.conjugate()) for X in b])
        C = (1 / 2) * sum([np.kron(X @ X, np.eye(d)) for X in b])
        return A + B + C
