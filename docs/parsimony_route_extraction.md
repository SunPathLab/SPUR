# Parsimony Route Extraction

Source files:

- `/home/yunong-xia/Downloads/Parsimony.pdf`
- `/home/yunong-xia/Downloads/T2352.pdf`

This document transcribes the hand-drawn parsimony route sketches. Some pages in
`Parsimony.pdf` are faint, so entries from those pages are marked with a
confidence label.

## Notation

- CN states use total:minor notation.
- `|` separates the two allelic homologs.
- `WGD` duplicates both homologs:

  ```text
  F1 | F2  --WGD-->  F3 F4 | F5 F6
  ```

- A route is written as:

  ```text
  truncal history -> sample 1 state
                  -> sample 2 state
  ```

- Fate vectors are written over the post-WGD copies, usually in the order:

  ```text
  [F3, F4, F5, F6]
  ```

  For non-WGD routes with three copies, the order is usually:

  ```text
  [F3, F4, F5]
  ```

## WGD Parsimony Routes

### 3:1 vs 2:1

Confidence: inferred from `Parsimony.pdf` page 1.

Sketch:

```text
F1 | F2
  --WGD-->
F3 F4 | F5 F6
  --delete one copy-->
3:1 state
  -> sample 1: 3:1
  -> sample 2: 2:1
```

The page shows a WGD node followed by a 3:1 MRCA-like state and a 2:1 branch.
The exact deletion identity is not fully clear in the scan.

### 4:2 vs 3:1

Confidence: clear from `Parsimony.pdf` page 1.

Sketch:

```text
F1 | F2
  --WGD-->
F3 F4 | F5 F6
  -> sample 1: 4:2
  -> sample 2: 3:1
```

Sample 1 keeps the WGD state:

```text
[1, 1, 1, 1]
```

Sample 2 deletes one of the four WGD-derived copies:

```text
[1, 1, 1, 0]
[1, 1, 0, 1]
[1, 0, 1, 1]
[0, 1, 1, 1]
```

### 5:2 vs 4:2

Confidence: clear from `T2352.pdf` page 1.

Sketch:

```text
F1 | F2
  --WGD-->
F3 F4 | F5 F6
  -> sample 1: 5:2
  -> sample 2: 4:2
```

Sample 2 keeps the WGD state:

```text
[1, 1, 1, 1]
```

Sample 1 gains one additional copy from the WGD-derived copies:

```text
[2, 1, 1, 1]
[1, 2, 1, 1]
[1, 1, 2, 1]
[1, 1, 1, 2]
```

### 5:2 vs 3:1

Confidence: clear from `T2352.pdf` page 2.

Sketch:

```text
F1 | F2
  --WGD-->
F3 F4 | F5 F6
  -> sample 1: 5:2
  -> sample 2: 3:1
```

Sample 1 gains one additional copy from the WGD-derived copies:

```text
[2, 1, 1, 1]
[1, 2, 1, 1]
[1, 1, 2, 1]
[1, 1, 1, 2]
```

Sample 2 deletes one WGD-derived copy:

```text
[0, 1, 1, 1]
[1, 0, 1, 1]
[1, 1, 0, 1]
[1, 1, 1, 0]
```

The page also groups possible 3:1 descendants into three hand-written cases:

```text
I:   F3 F4 | F5
II:  F3 | F5 F6
III: F4 | F5 F6
```

### 4:1 vs 3:1

Confidence: clear from `Parsimony.pdf` page 4.

Sketch:

```text
F1 | F2
  --WGD-->
F3 F4 | F5 F6
  --delete one copy-->
3:1
  -> sample 1: 4:1
  -> sample 2: 3:1
```

The drawn example is consistent with a 3:1 state followed by one extra gain on
one branch. The page notes that the constraint matrix is the same as for 3:1.

### 3:1 vs 3:1

Confidence: clear from `Parsimony.pdf` page 5.

Sketch:

```text
F1 | F2
  --WGD-->
F3 F4 | F5 F6
  --delete one copy-->
3:1
  -> sample 1: 3:1
  -> sample 2: 3:1
```

Both samples share the same 3:1 state.

Fate vectors shown:

```text
[1, 1, 1] vs [1, 1, 1]
```

### 4:2 vs 4:2

Confidence: clear from `Parsimony.pdf` page 5.

Sketch:

```text
F1 | F2
  --WGD-->
F3 F4 | F5 F6
  -> sample 1: 4:2
  -> sample 2: 4:2
```

Both samples keep the WGD state:

```text
[1, 1, 1, 1] vs [1, 1, 1, 1]
```

### 4:1 vs 4:1

Confidence: inferred from `Parsimony.pdf` page 6 and page 10.

Sketch:

```text
F1 | F2
  --WGD-->
F3 F4 | F5 F6
  --delete one copy-->
3:1
  --gain one copy-->
4:1
  -> sample 1: 4:1
  -> sample 2: 4:1
```

The pages note:

```text
A(4:1) = A_W(4:1)
```

### 5:2 vs 5:2

Confidence: inferred from `Parsimony.pdf` page 6.

Sketch:

```text
F1 | F2
  --WGD-->
F3 F4 | F5 F6
  --gain one copy-->
5:2
  -> sample 1: 5:2
  -> sample 2: 5:2
```

The route is represented as a WGD state followed by one additional gain.

## Non-WGD Parsimony Routes

### 3:1 vs 2:1

Confidence: clear from `T2352.pdf` page 3.

Sketch:

```text
F1 | F2
  --gain one copy-->
3:1
  -> sample 1: 3:1
  -> sample 2: 2:1
```

Fate vectors shown:

```text
3:1: [2, 1]
2:1: [1, 1]
```

Constraint matrix shown:

```text
A = [[1, 0],
     [0, 1]]
```

### 4:2 vs 3:1

Confidence: clear from `T2352.pdf` page 3.

Sketch:

```text
F1 | F2
  --gain one copy-->
3:1
  --gain one copy on the other homolog-->
4:2
  -> sample 1: 4:2
  -> sample 2: 3:1
```

Fate vector shown for the 4:2 branch:

```text
[1, 1, 2]
```

The page notes:

```text
A(3:1) = A_W(3:1)
```

### 4:1 vs 3:1

Confidence: clear from `Parsimony.pdf` page 8.

Sketch:

```text
F1 | F2
  --gain one copy-->
3:1
  --gain one copy on the amplified homolog-->
4:1
  -> sample 1: 4:1
  -> sample 2: 3:1
```

Fate vector shown:

```text
[2, 1, 1]
```

The page notes:

```text
A(3:1) = A_W(3:1)
```

### 4:1 vs 2:1

Confidence: clear from `Parsimony.pdf` page 9.

Sketch:

```text
F1 | F2
  --gain one copy-->
3:1
  --gain one copy on the amplified homolog-->
4:1
  -> sample 1: 4:1
  -> sample 2: 2:1
```

Fate vector shown for the 4:1 branch:

```text
[3, 1]
```

The page notes:

```text
A(2:1) = [[1, 0],
          [0, 1]]
```

### 3:1 vs 3:1

Confidence: clear from `Parsimony.pdf` page 9.

Sketch:

```text
F1 | F2
  --gain one copy-->
3:1
  -> sample 1: 3:1
  -> sample 2: 3:1
```

The page notes:

```text
A(3:1) = A_W(3:1)
```

### 4:2 vs 4:2

Confidence: inferred from `Parsimony.pdf` page 9.

Sketch:

```text
F1 | F2
  --gain one copy-->
3:1
  --gain one copy on the other homolog-->
4:2
  -> sample 1: 4:2
  -> sample 2: 4:2
```

The page notes:

```text
A(4:2) = A_W(4:2) without W constraint
```

## Implementation Notes

The WGD routes can be represented as ordinary allelic gain-history trees plus a
node-level `WGD` annotation. In these sketches, WGD creates four post-WGD copies
from the ancestral diploid state. Subsequent parsimony routes are represented by
copy fates over those WGD-derived copies.

The non-WGD routes are standard gain-history routes without a WGD node. Several
pages explicitly state that the same constraint matrix can be reused as the
corresponding WGD-derived matrix after removing the WGD constraint.

## Unclear Entries

The following items are present but not fully legible in the scans:

- Some details of the WGD 3:1 vs 2:1 route on `Parsimony.pdf` page 1.
- Some intermediate diagrams on `Parsimony.pdf` pages 2, 3, and 7.
- Exact row ordering of several hand-written matrices. The route-level copy
  fates above are more reliable than the matrix row order.
