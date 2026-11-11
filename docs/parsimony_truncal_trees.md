# Parsimony Truncal Trees and Constraint Matrices

This note summarizes the binary tree representation of the truncal part of the
parsimony routes extracted from the hand-drawn notes.

CN states use total:minor notation. Internally, a CN state has

```text
major = total - minor
minor = minor
```

The current SPUR tree model represents the major and minor allelic
gain-history trees separately:

```text
truncal tree = major tree + minor tree
```

For a binary tree, each row of the constraint matrix corresponds to one extant
copy. Each row contains `1` on the timing nodes along the path from the tree root
to that copy.

## Tree Shapes

### One Copy

```text
m0
```

Shape:

```text
L
```

Constraint block:

```text
[1]
```

### Two Copies

```text
      m0
     /  \
   m1    m2
```

Shape:

```text
(L,L)
```

Constraint block, columns `[m0, m1, m2]`:

```text
[[1, 1, 0],
 [1, 0, 1]]
```

### Three Copies

The canonical three-copy tree is:

```text
        m0
       /  \
     m1    m2
          /  \
        m3    m4
```

Shape:

```text
(L,(L,L))
```

Constraint block, columns `[m0, m1, m2, m3, m4]`:

```text
[[1, 1, 0, 0, 0],
 [1, 0, 1, 1, 0],
 [1, 0, 1, 0, 1]]
```

## Unique Truncal CN States

### Truncal 2:1

Major tree:

```text
M0
```

Minor tree:

```text
B0
```

Canonical key:

```text
major=L;minor=L
```

Columns:

```text
[M0, B0]
```

Constraint matrix:

```text
[[1, 0],
 [0, 1]]
```

### Truncal 3:1

Major tree:

```text
      M0
     /  \
   M1    M2
```

Minor tree:

```text
B0
```

Canonical key:

```text
major=(L,L);minor=L
```

Columns:

```text
[M0, M1, M2, B0]
```

Constraint matrix:

```text
[[1, 1, 0, 0],
 [1, 0, 1, 0],
 [0, 0, 0, 1]]
```

### Truncal 4:2

Major tree:

```text
      M0
     /  \
   M1    M2
```

Minor tree:

```text
      B0
     /  \
   B1    B2
```

Canonical key:

```text
major=(L,L);minor=(L,L)
```

Columns:

```text
[M0, M1, M2, B0, B1, B2]
```

Constraint matrix:

```text
[[1, 1, 0, 0, 0, 0],
 [1, 0, 1, 0, 0, 0],
 [0, 0, 0, 1, 1, 0],
 [0, 0, 0, 1, 0, 1]]
```

If the truncal 4:2 state is explicitly WGD-derived, add a WGD equality
constraint tying the major and minor root timing variables:

```text
M0 = B0
```

or

```text
[1, 0, 0, -1, 0, 0]
```

So the WGD-aware matrix is:

```text
[[ 1, 1, 0,  0, 0, 0],
 [ 1, 0, 1,  0, 0, 0],
 [ 0, 0, 0,  1, 1, 0],
 [ 0, 0, 0,  1, 0, 1],
 [ 1, 0, 0, -1, 0, 0]]
```

### Truncal 4:1

Major tree:

```text
        M0
       /  \
     M1    M2
          /  \
        M3    M4
```

Minor tree:

```text
B0
```

Canonical key:

```text
major=(L,(L,L));minor=L
```

Columns:

```text
[M0, M1, M2, M3, M4, B0]
```

Constraint matrix:

```text
[[1, 1, 0, 0, 0, 0],
 [1, 0, 1, 1, 0, 0],
 [1, 0, 1, 0, 1, 0],
 [0, 0, 0, 0, 0, 1]]
```

For an explicitly WGD-derived 4:1 route, the hand-written notes indicate that
the ordinary 4:1 matrix can be reused:

```text
A(4:1) = A_W(4:1)
```

This means the WGD event is encoded by the allowed route construction, not by an
extra visible equality row in the reduced surviving-copy matrix.

### Truncal 5:2

Major tree:

```text
        M0
       /  \
     M1    M2
          /  \
        M3    M4
```

Minor tree:

```text
      B0
     /  \
   B1    B2
```

Canonical key:

```text
major=(L,(L,L));minor=(L,L)
```

Columns:

```text
[M0, M1, M2, M3, M4, B0, B1, B2]
```

Constraint matrix:

```text
[[1, 1, 0, 0, 0, 0, 0, 0],
 [1, 0, 1, 1, 0, 0, 0, 0],
 [1, 0, 1, 0, 1, 0, 0, 0],
 [0, 0, 0, 0, 0, 1, 1, 0],
 [0, 0, 0, 0, 0, 1, 0, 1]]
```

If explicitly WGD-derived and both homolog roots are retained in the reduced
tree, the analogous WGD equality is:

```text
M0 = B0
```

or

```text
[1, 0, 0, 0, 0, -1, 0, 0]
```

## Route-to-Truncal Mapping

### WGD Parsimony Routes

| Route | Truncal tree used for shared timing | Notes |
|---|---|---|
| `3:1 vs 2:1` | `3:1`, `major=(L,L);minor=L` | WGD followed by loss; exact reduced WGD constraint is ambiguous in the notes. |
| `4:2 vs 3:1` | `4:2`, `major=(L,L);minor=(L,L)` | Use WGD-aware 4:2 matrix with `M0=B0`. |
| `5:2 vs 4:2` | `4:2` before branch-specific gain to `5:2` | Shared WGD state is 4:2. |
| `5:2 vs 3:1` | `4:2` before branch-specific gain/loss | Shared WGD state is 4:2. |
| `4:1 vs 3:1` | `3:1` before branch-specific gain to `4:1` | Notes say the constraint matrix is the same as 3:1. |
| `3:1 vs 3:1` | `3:1` | Both samples share 3:1. |
| `4:2 vs 4:2` | `4:2` | Both samples share WGD 4:2. |
| `4:1 vs 4:1` | `4:1` | Notes say `A(4:1)=A_W(4:1)`. |
| `5:2 vs 5:2` | `5:2` | WGD 4:2 followed by one shared gain. |

### Non-WGD Parsimony Routes

| Route | Truncal tree used for shared timing | Notes |
|---|---|---|
| `3:1 vs 2:1` | `2:1` or `3:1`, depending on whether the gain is truncal or branch-specific | Hand notes show `[2,1]` vs `[1,1]`. |
| `4:2 vs 3:1` | `3:1` before branch-specific gain to `4:2` | Notes show `A(3:1)=A_W(3:1)`. |
| `4:1 vs 3:1` | `3:1` before branch-specific gain to `4:1` | Fate vector `[2,1,1]`. |
| `4:1 vs 2:1` | `2:1` before branch-specific gains to `4:1` | Fate vector `[3,1]`. |
| `3:1 vs 3:1` | `3:1` | Both samples share 3:1. |
| `4:2 vs 4:2` | `4:2` | Both samples share 4:2; notes say WGD constraint can be removed for non-WGD. |

## Important Modeling Choice

For implementation, it is better to keep two layers separate:

1. The ordinary surviving-copy constraint matrix from the binary tree.
2. A WGD annotation or WGD equality constraint, when the route explicitly
   requires a WGD event.

This avoids forcing WGD into the generic `TreeNode` data structure, which
currently only knows tree shape and copy paths.
