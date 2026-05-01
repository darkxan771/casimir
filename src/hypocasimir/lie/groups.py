from typing import Callable

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


class LieGroup:
    """
    A class for the manipulation of classical compact Lie groups
    and their representations.
    """

    def __init__(self, type: str = "SU", size: int = 1):
        if type not in {"U", "SU", "SO", "SP"}:
            raise NotImplementedError
        self.type = type
        if size <= 0:
            raise ValueError("The size should be a positive integer")
        self.size = size

    def __eq__(self, other) -> bool:
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
    def weight_space(self):
        """
        Returns the weight space associated to the compact Lie group.
        """
        if self.type == "U":
            raise NotImplementedError("The unitary group is not semisimple")
        from ..systems import WeightSpace

        return WeightSpace(compute_system[self.type](self.size), self.rank)
