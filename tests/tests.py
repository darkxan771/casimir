# Tests

import sys

sys.path.insert(0, "../src/")

from sympy import Rational

from hypocasimir.labels import Partition, Signature
from hypocasimir.lie import LieGroup, SymmetricSpace
from hypocasimir.systems import WeightSpace


def test_SUSO():
    RW = WeightSpace("A", 4)
    G = LieGroup("SU", 5)
    K = LieGroup("SO", 5)
    X = SymmetricSpace("SU/SO", 5)
    P = Partition([5, 3, 2])
    W = P.highest_weight(G)
    assert G.description == "Special Unitary Group SU(5)"
    assert X.description == "Space of real structures SU(5)/SO(5)"
    assert G.dimension == 24
    assert X.dimension == 14
    assert X.isometry_group == G
    assert X.isotropy_group == [K]
    assert G.rank == 4
    assert K.rank == 2
    assert G.weight_space == RW
    assert repr(RW) == "Weight space A(4)"
    assert RW.killing_coefficient == 10
    assert W.coordinates == [3, 1, 0, -2, -2]
    assert P.size == 10
    assert P.dimension(G) == 2700
    assert P.casimir(G) == Rational(-22, 5)
    assert P.highest_K_type(X) == [Partition([5, 3])]
    assert P.hypocasimir(X) == Rational(-9, 5)
    G2 = LieGroup("SU", 6)
    K2 = LieGroup("SO", 6)
    X2 = SymmetricSpace("SU/SO", 6)
    assert X2.isometry_group == G2
    assert X2.isotropy_group == [K2]
    assert K2.rank == 3
    assert P.dimension(G2) == 15750
    assert P.casimir(G2) == Rational(-43, 9)
    assert P.highest_K_type(X2) == [Signature([5, 3, 2])]
    assert P.hypocasimir(X2) == Rational(-19, 9)


def test_SU2SP():
    G = LieGroup("SU", 10)
    K = LieGroup("SP", 5)
    X = SymmetricSpace("SU2/SP", 5)
    P = Partition([5, 3, 2, 1, 1])
    assert X.description == "Space of quaternionic structures SU(10)/SP(5)"
    assert G.dimension == 99
    assert X.dimension == 44
    assert X.isometry_group == G
    assert X.isotropy_group == [K]
    assert P.dimension(G) == 11561550
    assert P.casimir(G) == Rational(-132, 25)
    assert P.highest_K_type(X) == [P]
    assert P.hypocasimir(X) == Rational(-99, 50)


def test_SPU():
    RW = WeightSpace("C", 5)
    G = LieGroup("SP", 5)
    K = LieGroup("U", 5)
    X = SymmetricSpace("SP/U", 5)
    P = Partition([5, 3, 2])
    assert G.description == "Compact Symplectic Group SP(5)"
    assert X.description == "Space of complex structures SP(5)/U(5)"
    assert G.dimension == 55
    assert X.dimension == 30
    assert X.isometry_group == G
    assert X.isotropy_group == [K]
    assert G.rank == 5
    assert K.rank == 5
    assert G.weight_space == RW
    assert repr(RW) == "Weight space C(5)"
    assert RW.killing_coefficient == 12
    assert P.dimension(G) == 1514700
    assert P.casimir(G) == Rational(-31, 6)
    assert P.highest_K_type(X) == [Signature([5, 3, 2, 0, 0])]
    assert P.hypocasimir(X) == Rational(-5, 2)


def test_SO2U():
    RW = WeightSpace("D", 5)
    G = LieGroup("SO", 10)
    K = LieGroup("U", 5)
    X = SymmetricSpace("SO2/U", 5)
    P = Signature([5, 3, 2, 0, 0])
    assert G.description == "Special Orthogonal Group SO(10)"
    assert X.description == "Space of complex structures SO(10)/U(5)"
    assert G.dimension == 45
    assert X.dimension == 20
    assert X.isometry_group == G
    assert X.isotropy_group == [K]
    assert G.rank == 5
    assert K.rank == 5
    assert G.weight_space == RW
    assert repr(RW) == "Weight space D(5)"
    assert RW.killing_coefficient == 8
    assert P.dimension(G) == 1316250
    assert P.casimir(G) == Rational(-13, 2)
    assert P.highest_K_type(X) == [Signature([5, 3, 2, 0, 0])]
    assert P.hypocasimir(X) == Rational(-5, 2)


def test_GrR():
    G1 = LieGroup("SO", 10)
    X1 = SymmetricSpace("GrR", (6, 4))
    (K11, K12) = (LieGroup("SO", 6), LieGroup("SO", 4))
    P = Signature([5, 3, 2, 1, 1])
    assert G1.description == "Special Orthogonal Group SO(10)"
    assert X1.description == "Real Grassmannian manifold SO(10)/(SO(6) x SO(4))"
    assert G1.dimension == 45
    assert X1.dimension == 24
    assert X1.isometry_group == G1
    assert X1.isotropy_group == [K11, K12]
    assert G1.rank == 5
    assert P.dimension(G1) == 2502500
    assert P.casimir(G1) == Rational(-27, 4)
    assert P.highest_K_type(X1) == [Signature([5, 3, 2]), Signature([1, 1])]
    assert P.hypocasimir(X1) == Rational(-5, 2)
    X2 = SymmetricSpace("GrR", (7, 3))
    assert P.highest_K_type(X2) == [Partition([5, 3, 2]), Partition([1])]
    assert P.hypocasimir(X2) == -2
    X3 = SymmetricSpace("GrR", (7, 4))
    P = Partition([5, 3, 2, 1, 1])
    assert P.highest_K_type(X3) == [Partition([5, 3, 2]), Signature([1, 1])]
    assert P.hypocasimir(X3) == Rational(-7, 3)


def test_GrC():
    G = LieGroup("SU", 7)
    X = SymmetricSpace("GrC", (4, 3))
    P = Partition([5, 3, 2, 1, 1])
    assert G.description == "Special Unitary Group SU(7)"
    assert X.description == "Complex Grassmannian manifold SU(7)/S(U(4) x U(3))"
    assert G.dimension == 48
    assert X.dimension == 24
    assert X.isometry_group == G
    assert G.rank == 6
    assert P.dimension(G) == 107800
    assert P.casimir(G) == Rational(-222, 49)
    assert P.highest_K_type(X) == [Signature([5, 3, 2, 1]), Partition([1])]
    assert P.hypocasimir(X) == Rational(-29, 14)


def test_GrH():
    G = LieGroup("SP", 7)
    K = [LieGroup("SP", 4), LieGroup("SP", 3)]
    X = SymmetricSpace("GrH", (4, 3))
    P = Partition([5, 3, 2, 1, 1])
    assert G.description == "Compact Symplectic Group SP(7)"
    assert (
        X.description
        == "Quaternionic Grassmannian manifold SP(7)/(SP(4) x SP(3))"
    )
    assert G.dimension == 105
    assert X.dimension == 48
    assert X.isometry_group == G
    assert X.isotropy_group == K
    assert G.rank == 7
    assert P.dimension(G) == 506970464
    assert P.casimir(G) == Rational(-45, 8)
    assert P.highest_K_type(X) == [Partition([5, 3, 2, 1]), Partition([1])]
    assert P.hypocasimir(X) == Rational(-33, 16)
