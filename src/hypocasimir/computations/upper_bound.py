from collections.abc import Callable, Iterable
from itertools import chain

from ..labels import Label, Partitions, Signature
from ..labels.labels import complete_list
from ..lie import LieGroup, SymmetricSpace


def ci_SO(N: int, n: int) -> Iterable:
    if n % 2 == 1:
        return (P for P in Partitions(N) if len(P) <= int((n - 1) / 2))
    else:
        m = int(n / 2)
        gen = (P for P in Partitions(N) if len(P) <= m)
        gen1 = [Signature(complete_list(m, P.terms)) for P in gen]
        gen2 = [S.switch_signed_partition() for S in gen1 if S.terms[-1] != 0]
        return chain(gen1, gen2)


compute_iterable: dict[tuple[str], Callable[[int, int], Iterable]] = {}
compute_iterable["SU"] = lambda N, n: (P for P in Partitions(N) if len(P) < n)
compute_iterable["SP"] = lambda N, n: (P for P in Partitions(N) if len(P) <= n)
compute_iterable["SO"] = ci_SO


def list_terms(N: int, S: LieGroup | SymmetricSpace) -> dict[Label, float]:
    """
    Returns a dictionary whose keys are the labels with size smaller than N,
    and whose values are the corresponding terms in the upper bound on
    dtv(mu_t, Haar)^2 at cutoff time.
    """
    if isinstance(S, LieGroup):
        CI = compute_iterable[S.type]
        return {L: L.up_term(S) for k in range(1, N + 1) for L in CI(k, S.size)}
    else:
        G = S.isometry_group
        CI = compute_iterable[G.type]
        return {L: L.up_term(S) for k in range(1, N + 1) for L in CI(k, G.size)}
