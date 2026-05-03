from fractions import Fraction

from .weights import Weight


# type A
def root_A(n: int, i: int, j: int) -> Weight:
    L = [0] * (n + 1)
    L[i] = 1
    L[j] = -1
    return Weight("A", L)


def weight_A(n, i: int) -> Weight:
    L = [Fraction(n - i, n + 1)] * (i + 1) + [Fraction(-i - 1, n + 1)] * (n - i)
    return Weight("A", L)


# type B
def long_root_B(n: int, i: int, j: int, sign: int) -> Weight:
    L = [0] * n
    L[i] = 1
    L[j] = sign
    return Weight("B", L)


def short_root_B(n: int, i: int) -> Weight:
    L = [0] * n
    L[i] = 1
    return Weight("B", L)


def weight_B(n, i: int) -> Weight:
    if i < n - 1:
        L = [1] * (i + 1) + [0] * (n - i - 1)
        return Weight("B", L)
    else:
        return Weight("B", [Fraction(1, 2) for _ in range(n)])


# type C
def short_root_C(n: int, i: int, j: int, sign: int) -> Weight:
    L = [0] * n
    L[i] = 1
    L[j] = sign
    return Weight("C", L)


def long_root_C(n: int, i: int) -> Weight:
    L = [0] * n
    L[i] = 2
    return Weight("C", L)


def weight_C(n: int, i: int) -> Weight:
    L = [1] * (i + 1) + [0] * (n - i - 1)
    return Weight("C", L)


# type D
def root_D(n: int, i: int, j: int, sign: int) -> Weight:
    L = [0] * n
    L[i] = 1
    L[j] = sign
    return Weight("D", L)


def weight_D(n, i: int) -> Weight:
    if i < n - 2:
        L = [1] * (i + 1) + [0] * (n - i - 1)
        return Weight("D", L)
    elif i == n - 1:
        return Weight(
            "D", [Fraction(1, 2) for _ in range(n - 1)] + [Fraction(-1, 2)]
        )
    else:
        return Weight("D", [Fraction(1, 2) for _ in range(n)])
