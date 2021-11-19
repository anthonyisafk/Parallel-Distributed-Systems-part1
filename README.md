# pds-part1
## Parallel and Distributed Systems - Assignment 1: Sparse matrices
### Aristotle University Thessaloniki - Electrical and Computer Engineering
#### Authored by:
###### Antonios Antoniou - aantonii@ece.auth.gr
###### Efthymios Grigorakis - eegrigor@ece.auth.gr

## The goal of this assignment
The first assignment of the Parallel and Distributed Systems course requires us to read an `.mtx` file, depiciting a square, symmetric, undirected, non-weighted graph, in the form of a triangular matrix, whose complete set of values has to be determined by the programmer. For this table, we have to calculate how many triangles (i.e subgraphs with 3 vertices and 3 edges forming a closed structure) each vertex is a member of.
\
\
This is realized by calculating the
\
<img src="https://render.githubusercontent.com/render/math?math=H = A \bigodot A^{2}"> matrix
\
(where <img src="https://render.githubusercontent.com/render/math?math=\bigodot"> denotes the Hadamard, element-wise multiplication of two matrices). Afterwards, we calculate the
\
<img src="https://render.githubusercontent.com/render/math?math=C = (H \cdot e) / 2"> vector,
\
(where is a vector filled with ones, the same size as A). Each element on the `C` matrix now represents how many triangles this vertex is in.
\
