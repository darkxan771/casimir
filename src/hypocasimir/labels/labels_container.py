from collections.abc import Callable, Container, Iterator
from itertools import chain, count

from .partitions import Partition


class InfiniteSetError(Exception):
    pass


def generating_series_P(N: int) -> list[int]:
    """
    Computes the N first terms of the generating series of
    integer partitions.
    """
    divs = [[] for _ in range(N + 1)]
    for i in range(1, N + 1):
        for j in range(i, N + 1, i):
            divs[j].append(i)
    res = [0] * (N + 1)
    sigma = [sum(x, 0) for x in divs]
    res[0] = 1
    for n in range(1, N + 1):
        res[n] = sum([sigma[n - k] * res[k] for k in range(n)]) // n
    return res


class _PartitionsIterator_n(Iterator):
    def __init__(self, n: int):
        self.a = [0] * (n + 1)
        self.k = 1
        self.x = 1
        self.y = n - 1
        self.end_of_cycle = True

    def __iter__(self):
        return self

    def __next__(self):
        if self.k == 0 and self.end_of_cycle:
            raise StopIteration
        else:
            if self.end_of_cycle:
                self.end_of_cycle = False
                self.x = self.a[self.k - 1] + 1
                self.k -= 1
            while 2 * self.x <= self.y:
                self.a[self.k] = self.x
                self.y -= self.x
                self.k += 1
            if self.x <= self.y:
                self.a[self.k] = self.x
                self.a[self.k + 1] = self.y
                res = Partition(self.a[self.k + 1 :: -1])
                self.x += 1
                self.y -= 1
                return res
            else:
                self.a[self.k] = self.x + self.y
                self.y = self.x + self.y - 1
                res = Partition(self.a[self.k :: -1])
                self.end_of_cycle = True
                return res


class Partitions(Container):
    """
    A container for integer partitions.

    If a size n is given, the container is restricted to integer
    partitions with size n. In any case, one can iterate upon the
    container.
    """

    def __init__(self, n: int | None = None):
        self.category = "partition"
        self.order = n

    def __repr__(self) -> str:
        res = "Partitions"
        if self.order is not None:
            res += f" with size {self.order}"
        return res

    def __len__(self) -> int:
        if self.order is None:
            raise InfiniteSetError("Infinite set")
        else:
            return int(type(self).generating_series(self.order)[-1])

    @property
    def cardinality(self) -> int:
        return len(self)

    @classmethod
    def generating_series(cls, N: int) -> list[int]:
        return generating_series_P(N)

    @classmethod
    def iter_n(cls) -> Callable[[int], Iterator]:
        return lambda n: _PartitionsIterator_n(n)

    def __iter__(self) -> Iterator:
        if self.order is not None:
            return type(self).iter_n()(self.order)
        else:
            return chain.from_iterable(
                type(self).iter_n()(n) for n in count(1)
            ).__iter__()

    def __contains__(self, obj: object) -> bool:
        A = isinstance(obj, Partition)
        B = True
        if self.order is not None:
            B = self.order == obj.size

        return A and B
