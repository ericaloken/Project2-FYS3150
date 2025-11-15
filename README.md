# Project 2 – Solving Eigenvalue problem with Jacobi Rotation Algorithm
_This repository contains the work done by Erica Løken, Mads Mestl and Heine Husdal for Project 3 in the course 
FYS3150 – Computational Physics at the University of Oslo._

## Project introduction
In this project, we will solve an eigenvalue that arises as a special case of a one-dimensional buckling beam. We will study the shapes the beam can take when a compressive force is applied, using Jacobi's rotation algorithm. The buckling beam is modelled by the differential equation 

$$
\gamma \frac{d^2 u(x)}{dx^2} = -F u(x), \quad x \in [0, L],
$$
which, after discretization, corresponds to the tridiagonal eigenvalue problem 
$$
A \vec{v} = \lambda \vec{v}.
$$  

The Jacobi algortithm solves this problem by iteratively diagonalizing \(A\). The eigenvectors and eigenvalues from the Jacobi algorithm are compared to the analytical eigenvectors and eigenvecotors and the results from Armadillos eigensolver.

## Folder structure

- `src/` - C++ files
- `include/` - Header files
- `data/` - Text files
- `python/` - Python files

## Running the project
To compile and run the project, simply run
```bash
make
```
in the `/project2` folder and the Makefile will handle compiling and linking.


