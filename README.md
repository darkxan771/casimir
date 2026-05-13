# casimir
Python code for the manipulation of integer partitions, Casimir coefficients 
of elliptic and hypoelliptic Brownian motions, and dimensions of representations
of compact Lie groups.

## Installation

```bash
git clone https://github.com/darkxan771/casimir.git
```

In the jupyter/ folder, one can launch `jupyter lab`, 
then create a Python notebook and execute in the first cell: 
```python
import sys
sys.path.insert(0, "../src/")
from hypocasimir import *
```


## Features and usage
The classical simple and simply connected compact Lie groups and symmetric
spaces are accessed by:

```python
G = LieGroup("SO", 10)
X = SymmetricSpace("GrC", (7, 4))
G.description
>>> 'Special Orthogonal Group SO(10)'
X.description
>>> 'Complex Grassmannian manifold SU(11)/S(U(7) x U(4))'
```

For the moment they are not defined as containers, although this might be
possible in future versions (in particular for the Lie groups). Various
properties and methods are defined: `dimension`, `rank`, `isometry_group`, 
etc.

On the other hand, labels of irreducible representations are encoded as 
Partition (for SU, SOodd, SP) or Signature (for U, SOeven) instances. 
Then, various quantities can be computed: `dimension` (of an irreducible
representation of G), `casimir` (for the action of the Laplace-Beltrami 
operator on the irreducible representation space), `highest_K_type` (restriction
from G to K, with X=G/K symmetric space), and `hypocasimir` (for the action
of the hypoelliptic operator associated to a symmetric space).

```python
P = Partition([5, 3, 2])
G = LieGroup("SU", 5)
X = SymmetricSpace("SU/SO", 5)
P.dimension(G)
>>> 2700
P.casimir(G)
>>> -22/5
P.highest_K_type(X)
>>> [5, 3]
P.hypocasimir(X)
>>> -9/5
```

A dictionary with all terms with label smaller than N in the upper bound
can be obtained with `list_terms(N, S)` where S is a Lie group or a symmetric
space.

```python
for _, (k, v) in enumerate(
    list_terms(5, SymmetricSpace("SU/SO", 30)).items()
    ):
    print(k, v)
>>> 
[1] 0.8092640218929349
[1, 1] 0.24817237016913793
[2] 0.11449376029790077
[1, 1, 1] 0.051148586919033956
[2, 1] 0.06433576277557028
[3] 0.005023980233364324
[1, 1, 1, 1] 0.008943907752393036
[2, 1, 1] 0.017296876219644392
[3, 1] 0.003682591691227675
[2, 2] 0.003562973603797588
[4] 8.636642458499071e-05
[1, 1, 1, 1, 1] 0.0015055652820075143
[2, 1, 1, 1] 0.0035468985168519663
[3, 1, 1] 0.0011610641242476506
[2, 2, 1] 0.0017552349905100195
[4, 1] 7.424000659718629e-05
[3, 2] 0.00037369833200651484
[5] 6.60590546404583e-07
```
