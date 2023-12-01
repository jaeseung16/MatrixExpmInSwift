# MatrixExp
![](https://img.shields.io/badge/version-0.3.1-blue)
![](https://img.shields.io/badge/license-MIT-green)
![](https://img.shields.io/badge/last%20updated-October%202020-orange)

Copyright (c) 2020 Jae-Seung Lee

License: MIT License

MatrixExp is a Swift package for numerically evaluating matrix exponentiation.

### References

  1. N. J. Higham, The scaling and squaring method for the matrix exponential revisited. SIAM J. Matrix Anal. Appl., 26(4), (2005), pp. 1179-1193. [doi:10.1137/04061101X](https://doi.org/10.1137/04061101X)
  2. A. H. Al-Mohy and N. J. Higham, A new scaling and squaring algorithm for the matrix exponential, SIAM J. Matrix Anal. Appl., 31(3), (2009), pp. 970-989. [doi:10.1137/09074721X](https://doi.org/10.1137/09074721X)

### Dependencies

  1. [Swift Numerics](https://github.com/apple/swift-numerics)
  2. [LANumerics](https://github.com/phlegmaticprogrammer/LANumerics)

### Usage

Initialize a matrix using LANumerics' `Matrix<Element>` and instantiate an instance of `MatrixExponentiator<Element>`. Then, call compute() to get the result, which is either `Matrix<Element>` or `nil`.

```swift
let M1 = Matrix<Double>(rows: [[0, 1], [1, 0]])
let result = MatrixExponentiator<Double>(M1).compute()

let M2 = Matrix<Complex<Double>>(rows: [[Complex<Double>(0,0), Complex<Double>(0, -1)], [Complex<Double>(0, 1), Complex<Double>(0, 0)]])
let result = MatrixExponentiator<Complex>(M2).compute()
```

### Update History

#### Version 0.3.0 (11/29/2023)

- Update LANumerics to 0.1.12 and swift-numerics to 0.1.0
- Update the interface and usage

#### Version 0.2.0 (12/10/2020)

Implement a function corresponding to `normest1` in Matlab or [onenormest](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.onenormest.html) in Scipy to estimate the 1-norm of the powers of a given matrix [1]


[1] Nicholas J. Higham and Fran\c{c}oise Tisseur, A Block Algorithm for Matrix 1-Norm Estimation with an Application to 1-Norm Pseudospectra, SIAM J. Matrix Anal. App. 21, 1185-1201, 2000. [doi: 10.1137/S0895479899356080](https://doi.org/10.1137/S0895479899356080)

