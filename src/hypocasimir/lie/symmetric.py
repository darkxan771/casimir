from typing import Callable

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
        if type in {"SU/SO", "SU2/SP", "SP/U", "SO2/U"} and isinstance(
            size, int
        ):
            self.type = type
            self.stype = "structure"
            self.n = size
        elif type in {"GrR", "GrC", "GrH"} and isinstance(size, tuple):
            self.type = type
            self.stype = "grassmann"
            self.p = max(size)
            self.q = min(size)
            self.n = self.p + self.q
        else:
            raise NotImplementedError

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, SymmetricSpace)
            and self.type == other.type
            and self.n == other.n
            and ((self.stype == "structure") or (self.p == other.p))
        )

    def __repr__(self):
        if self.stype == "structure":
            return compute_repr_structure[self.type](self.n)
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
        if self.stype == "structure":
            return space_dimension_structure[self.type](self.n)
        else:
            return space_dimension_grassmann[self.type](self.p, self.q)

    @property
    def isometry_group(self) -> LieGroup:
        """
        Returns the simple simply connected compact Lie group G such
        that X=G/K.
        """
        return compute_G[self.type](self.n)

    @property
    def isotropy_group(self) -> list[LieGroup]:
        """
        Returns the stabilizer K of a point for the action of the
        isometry group G on X=G/K. Not implemented for complex Grassmannian
        manifolds.
        """
        if self.stype == "structure":
            return compute_K_structure[self.type](self.n)
        else:
            if self.type == "GrC":
                raise NotImplementedError
            return compute_K_grassmann[self.type](self.p, self.q)
