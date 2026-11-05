from collections.abc import Callable

import numpy as np


def basis_e(n: int, i: int, j: int) -> np.array:
    """
    The e(i,j) matrix in size n*n.
    """
    res = np.zeros((n, n), dtype=np.int64)
    res[i, j] = 1
    return res


def sp_matrix_1(n: int, i: int, j: int) -> np.array:
    """
    The quaternionic matrix e(i,j) represented in a 2n*2n complex matrix.
    """
    return basis_e(2 * n, 2 * i, 2 * j) + basis_e(2 * n, 2 * i + 1, 2 * j + 1)


def sp_matrix_i(n: int, i: int, j: int) -> np.array:
    """
    The quaternionic matrix i*e(i,j) represented in a 2n*2n complex matrix.
    """
    return (1j) * (
        basis_e(2 * n, 2 * i, 2 * j) - basis_e(2 * n, 2 * i + 1, 2 * j + 1)
    )


def sp_matrix_j(n: int, i: int, j: int) -> np.array:
    """
    The quaternionic matrix j*e(i,j) represented in a 2n*2n complex matrix.
    """
    return basis_e(2 * n, 2 * i, 2 * j + 1) - basis_e(2 * n, 2 * i + 1, 2 * j)


def sp_matrix_k(n: int, i: int, j: int) -> np.array:
    """
    The quaternionic matrix k*e(i,j) represented in a 2n*2n complex matrix.
    """
    return (1j) * (
        basis_e(2 * n, 2 * i, 2 * j + 1) + basis_e(2 * n, 2 * i + 1, 2 * j)
    )


def basis_su(n: int, normalised: bool = True) -> list[np.array]:
    """
    Returns an orthogonal basis of the Lie algebra su(n).
    """
    res = [
        (basis_e(n, i, j) - basis_e(n, j, i))
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        (1j) * (basis_e(n, i, j) + basis_e(n, j, i))
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    if normalised:
        res = [M / np.sqrt(4 * n) for M in res]
        res += [
            (1j)
            * (sum([basis_e(n, i, i) for i in range(k)]) - k * basis_e(n, k, k))
            / np.sqrt(2 * n * k * (k + 1))
            for k in range(1, n)
        ]
    else:
        res += [
            (1j)
            * (sum([basis_e(n, i, i) for i in range(k)]) - k * basis_e(n, k, k))
            for k in range(1, n)
        ]
    return res


def basis_u(n: int, normalised: bool = True) -> list[np.array]:
    """
    Returns an orthogonal basis of the Lie algebra u(n).
    """
    res = basis_su(n, normalised)[: n * (n - 1)]
    if normalised:
        res += [(1j) / (np.sqrt(2 * n)) * basis_e(n, i, i) for i in range(n)]
    else:
        res += [(1j) * basis_e(n, i, i) for i in range(n)]
    return res


def basis_so(n: int, normalised: bool = True) -> list[np.array]:
    """
    Returns an orthogonal basis of the Lie algebra so(n).
    """
    res = [
        (basis_e(n, i, j) - basis_e(n, j, i))
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    if normalised:
        res = [M / np.sqrt(2 * n - 4) for M in res]
    return res


def basis_sp(n: int, normalised: bool = True) -> list[np.array]:
    """
    Returns an orthogonal basis of the Lie algebra sp(n).
    """
    res = [
        sp_matrix_1(n, i, j) - sp_matrix_1(n, j, i)
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        sp_matrix_i(n, i, j) + sp_matrix_i(n, j, i)
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        sp_matrix_j(n, i, j) + sp_matrix_j(n, j, i)
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        sp_matrix_k(n, i, j) + sp_matrix_k(n, j, i)
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [sp_matrix_i(n, i, i) for i in range(n)]
    res += [sp_matrix_j(n, i, i) for i in range(n)]
    res += [sp_matrix_k(n, i, i) for i in range(n)]
    if normalised:
        res = [
            res[i] / np.sqrt(8 * (n + 1)) for i in range(2 * n * (n - 1))
        ] + [
            res[i] / np.sqrt(4 * (n + 1))
            for i in range(2 * n * (n - 1), n * (2 * n + 1))
        ]
    return res


bases_lie_algebra: dict[str, Callable[[int, bool], list[np.array]]] = {}
bases_lie_algebra["U"] = basis_u
bases_lie_algebra["SU"] = basis_su
bases_lie_algebra["SO"] = basis_so
bases_lie_algebra["SP"] = basis_sp


def double_matrix_spu(m: np.array) -> np.array:
    """
    Takes a complex n*n matrix and returns the complex 2n*2n matrix
    corresponding to the embedding of H into C^(2*2), restricted to C.
    """
    n = m.shape[0]
    res = np.zeros((2 * n, 2 * n)) + 0j
    for i in range(n):
        for j in range(n):
            res[2 * i, 2 * j] = m[i, j]
            res[2 * i + 1, 2 * j + 1] = (m[i, j] + 0j).conjugate()
    return res


def double_matrix_so2u(m: np.array) -> np.array:
    """
    Takes a complex n*n matrix and returns the real 2n*2n matrix
    corresponding to the embedding of C into R^(2*2).
    """
    n = m.shape[0]
    res = np.zeros((2 * n, 2 * n))
    for i in range(n):
        for j in range(n):
            res[2 * i, 2 * j] = (m[i, j] + 0).real
            res[2 * i + 1, 2 * j + 1] = (m[i, j] + 0).real
            res[2 * i, 2 * j + 1] = (m[i, j] + 0).imag
            res[2 * i + 1, 2 * j] = -(m[i, j] + 0).imag
    return res


def basis_m_su2sp(n: int) -> list[np.array]:
    """
    Returns an orthonormal basis of the complement of sp(n) in su(2n).
    """
    res = []
    res += [
        (
            basis_e(2 * n, 2 * i, 2 * j)
            - basis_e(2 * n, 2 * i + 1, 2 * j + 1)
            - basis_e(2 * n, 2 * j, 2 * i)
            + basis_e(2 * n, 2 * j + 1, 2 * i + 1)
        )
        / np.sqrt(16 * n)
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        (1j)
        * (
            basis_e(2 * n, 2 * i, 2 * j)
            + basis_e(2 * n, 2 * i + 1, 2 * j + 1)
            + basis_e(2 * n, 2 * j, 2 * i)
            + basis_e(2 * n, 2 * j + 1, 2 * i + 1)
        )
        / np.sqrt(16 * n)
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        (
            basis_e(2 * n, 2 * i, 2 * j + 1)
            + basis_e(2 * n, 2 * i + 1, 2 * j)
            - basis_e(2 * n, 2 * j + 1, 2 * i)
            - basis_e(2 * n, 2 * j, 2 * i + 1)
        )
        / np.sqrt(16 * n)
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        (1j)
        * (
            basis_e(2 * n, 2 * i, 2 * j + 1)
            - basis_e(2 * n, 2 * i + 1, 2 * j)
            + basis_e(2 * n, 2 * j + 1, 2 * i)
            - basis_e(2 * n, 2 * j, 2 * i + 1)
        )
        / np.sqrt(16 * n)
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        (1j)
        * np.diag([1] * (2 * k) + [-k, -k] + [0] * (2 * n - 2 * k - 2))
        / np.sqrt(8 * n * k * (k + 1))
        for k in range(1, n)
    ]
    return res


def basis_m_spu(n: int) -> list[np.array]:
    """
    Returns an orthonormal basis of the complement of u(n) in sp(n).
    """
    res = []
    res += [
        (sp_matrix_j(n, i, j) + sp_matrix_j(n, j, i)) / np.sqrt(8 * (n + 1))
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        (sp_matrix_k(n, i, j) + sp_matrix_k(n, j, i)) / np.sqrt(8 * (n + 1))
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [sp_matrix_j(n, i, i) / np.sqrt(4 * (n + 1)) for i in range(n)]
    res += [sp_matrix_k(n, i, i) / np.sqrt(4 * (n + 1)) for i in range(n)]
    return res


def basis_m_so2u(n: int) -> list[np.array]:
    """
    Returns an orthonormal basis of the complement of u(n) in so(2n).
    """
    res = []
    res += [
        (
            basis_e(2 * n, 2 * i, 2 * j)
            - basis_e(2 * n, 2 * i + 1, 2 * j + 1)
            - basis_e(2 * n, 2 * j, 2 * i)
            + basis_e(2 * n, 2 * j + 1, 2 * i + 1)
        )
        / np.sqrt(8 * (n - 1))
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        (
            basis_e(2 * n, 2 * i, 2 * j + 1)
            + basis_e(2 * n, 2 * i + 1, 2 * j)
            - basis_e(2 * n, 2 * j + 1, 2 * i)
            - basis_e(2 * n, 2 * j, 2 * i + 1)
        )
        / np.sqrt(8 * (n - 1))
        for i in range(n - 1)
        for j in range(i + 1, n)
    ]
    return res


def basis_k_grr(p: int, q: int) -> list[np.array]:
    """
    Returns an orthonormal basis of so(p) * so(q) in so(n), n = p + q.
    """
    n = p + q
    res = [
        (basis_e(n, i, j) - basis_e(n, j, i))
        for i in range(p - 1)
        for j in range(i + 1, p)
    ]
    res += [
        (basis_e(n, i, j) - basis_e(n, j, i))
        for i in range(p, n)
        for j in range(i + 1, n)
    ]
    res = [M / np.sqrt(2 * n - 4) for M in res]
    return res


def basis_m_grr(p: int, q: int) -> list[np.array]:
    """
    Returns an orthonormal basis of the complement of so(p) * so(q)
    in so(n), n = p + q.
    """
    n = p + q
    res = [
        (basis_e(n, i, j) - basis_e(n, j, i))
        for i in range(p)
        for j in range(p, n)
    ]
    res = [M / np.sqrt(2 * n - 4) for M in res]
    return res


def basis_k_grc(p: int, q: int) -> list[np.array]:
    """
    Returns an orthonormal basis of s(u(p) * u(q)) in su(n), n = p + q.
    """
    n = p + q
    res = [
        (basis_e(n, i, j) - basis_e(n, j, i))
        for i in range(p - 1)
        for j in range(i + 1, p)
    ]
    res += [
        (1j) * (basis_e(n, i, j) + basis_e(n, j, i))
        for i in range(p - 1)
        for j in range(i + 1, p)
    ]
    res += [
        (basis_e(n, i, j) - basis_e(n, j, i))
        for i in range(p, n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        (1j) * (basis_e(n, i, j) + basis_e(n, j, i))
        for i in range(p, n - 1)
        for j in range(i + 1, n)
    ]
    res = [M / np.sqrt(4 * n) for M in res]
    res += [
        (1j)
        * (sum([basis_e(n, i, i) for i in range(k)]) - k * basis_e(n, k, k))
        / np.sqrt(2 * n * k * (k + 1))
        for k in range(1, n)
    ]
    return res


def basis_m_grc(p: int, q: int) -> list[np.array]:
    """
    Returns an orthonormal basis of the complement of s(u(p) * u(q))
    in su(n), n = p + q.
    """
    n = p + q
    res = [
        (basis_e(n, i, j) - basis_e(n, j, i))
        for i in range(p)
        for j in range(p, n)
    ]
    res += [
        (1j) * (basis_e(n, i, j) + basis_e(n, j, i))
        for i in range(p)
        for j in range(p, n)
    ]
    res = [M / np.sqrt(4 * n) for M in res]
    return res


def basis_k_grh(p: int, q: int) -> list[np.array]:
    """
    Returns an orthonormal basis of sp(p) * sp(q) in sp(n), n = p + q.
    """
    n = p + q
    res = [
        sp_matrix_1(n, i, j) - sp_matrix_1(n, j, i)
        for i in range(p - 1)
        for j in range(i + 1, p)
    ]
    res += [
        sp_matrix_i(n, i, j) + sp_matrix_i(n, j, i)
        for i in range(p - 1)
        for j in range(i + 1, p)
    ]
    res += [
        sp_matrix_j(n, i, j) + sp_matrix_j(n, j, i)
        for i in range(p - 1)
        for j in range(i + 1, p)
    ]
    res += [
        sp_matrix_k(n, i, j) + sp_matrix_k(n, j, i)
        for i in range(p - 1)
        for j in range(i + 1, p)
    ]
    res += [
        sp_matrix_1(n, i, j) - sp_matrix_1(n, j, i)
        for i in range(p, n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        sp_matrix_i(n, i, j) + sp_matrix_i(n, j, i)
        for i in range(p, n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        sp_matrix_j(n, i, j) + sp_matrix_j(n, j, i)
        for i in range(p, n - 1)
        for j in range(i + 1, n)
    ]
    res += [
        sp_matrix_k(n, i, j) + sp_matrix_k(n, j, i)
        for i in range(p, n - 1)
        for j in range(i + 1, n)
    ]
    res = [M / np.sqrt(8 * (n + 1)) for M in res]
    res += [sp_matrix_i(n, i, i) / np.sqrt(4 * (n + 1)) for i in range(n)]
    res += [sp_matrix_j(n, i, i) / np.sqrt(4 * (n + 1)) for i in range(n)]
    res += [sp_matrix_k(n, i, i) / np.sqrt(4 * (n + 1)) for i in range(n)]
    return res


def basis_m_grh(p: int, q: int) -> list[np.array]:
    """
    Returns an orthonormal basis of the complement of sp(p) * sp(q)
    in sp(n), n = p + q.
    """
    n = p + q
    res = [
        sp_matrix_1(n, i, j) - sp_matrix_1(n, j, i)
        for i in range(p)
        for j in range(p, n)
    ]
    res += [
        sp_matrix_i(n, i, j) + sp_matrix_i(n, j, i)
        for i in range(p)
        for j in range(p, n)
    ]
    res += [
        sp_matrix_j(n, i, j) + sp_matrix_j(n, j, i)
        for i in range(p)
        for j in range(p, n)
    ]
    res += [
        sp_matrix_k(n, i, j) + sp_matrix_k(n, j, i)
        for i in range(p)
        for j in range(p, n)
    ]
    res = [M / np.sqrt(8 * (n + 1)) for M in res]
    return res


bases_k_structure: dict[str, Callable[[int], list[np.array]]] = {}
bases_k_structure["SU/SO"] = lambda n: basis_su(n)[0 : int(n * (n - 1) / 2)]
bases_k_structure["SU2/SP"] = lambda n: [
    np.sqrt((n + 1) / (2 * n)) * M for M in basis_sp(n)
]

bases_k_structure["SP/U"] = lambda n: [
    np.sqrt(n / (2 * n + 2)) * double_matrix_spu(M) for M in basis_u(n)
]
bases_k_structure["SO2/U"] = lambda n: [
    np.sqrt(n / (2 * n - 2)) * double_matrix_so2u(M) for M in basis_u(n)
]

bases_m_structure: dict[str, Callable[[int], list[np.array]]] = {}
bases_m_structure["SU/SO"] = lambda n: basis_su(n)[int(n * (n - 1) / 2) :]
bases_m_structure["SU2/SP"] = basis_m_su2sp
bases_m_structure["SP/U"] = basis_m_spu
bases_m_structure["SO2/U"] = basis_m_so2u
bases_k_grassmann: dict[str, Callable[[int, int], list[np.array]]] = {}
bases_k_grassmann["GrR"] = basis_k_grr
bases_k_grassmann["GrC"] = basis_k_grc
bases_k_grassmann["GrH"] = basis_k_grh
bases_m_grassmann: dict[str, Callable[[int, int], list[np.array]]] = {}
bases_m_grassmann["GrR"] = basis_m_grr
bases_m_grassmann["GrC"] = basis_m_grc
bases_m_grassmann["GrH"] = basis_m_grh
