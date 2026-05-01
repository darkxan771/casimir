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


