//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 11/27/23.
//

import Foundation
import Numerics
import LANumerics

public class MatrixExponentiator<T> where T: Exponentiable & ElementaryFunctions, T.Magnitude: Real {
    public let matrix: Matrix<T>
    public var result: Matrix<T>?
    public var scaling: Int = 0
    public var orderPadeApproximant: PadeApproximantOrder?
    public var calculationType = MatrixExpCalculationType.notApplicable
    
    public init(_ matrix: Matrix<T>) {
        self.matrix = matrix
    }
    
    public func compute() -> Matrix<T>? {
        if matrix.isSquare {
            if matrix.isDiag {
                calculationType = .diag
                result = expDiag()
            } else if matrix.isHermitian {
                calculationType = .hermitian
                result = expHermitian()
            } else {
                calculationType = .pade
                (result, scaling, orderPadeApproximant) = exp()
            }
        } else {
            result = nil
        }
        return result
    }
    
    private func expDiag() -> Matrix<T> {
        let range = 0..<matrix.columns
        let diag = range.map { T.exp(matrix[$0,$0]) }
        return Matrix<T>(rows: matrix.rows, columns: matrix.columns, diagonal: diag)
    }
    
    private func expHermitian() -> Matrix<T> {
        let (_, schurForm, schurVectors) = matrix.schur()!
        let range = 0..<schurForm.columns
        let diag = range.map { T.exp(schurForm[$0,$0]) }
        let expSchurForm = Matrix<T>(rows: schurForm.rows, columns: schurForm.columns, diagonal: diag)
        return schurVectors * expSchurForm * schurVectors.adjoint
    }
    
    private func exp() -> (Matrix<T>, Int, PadeApproximantOrder) {
        var result = Matrix<T>.zeros(matrix.rows, matrix.columns)
        
        let (scaling, order, matrixPowers) = MatrixExponentiator<T>.getExpmParams(from: matrix)
        
        let factor = scaling > 0 ? T.raiseTwo(to: scaling) : T.one
        var scaledMatrix = scaling > 0 ? matrix.map { $0 / factor } : matrix
        let scaledMatrixPowers = (0..<matrixPowers.count).map { power in
            if scaling > 0 {
                let factor = T.raiseTwo(to: (power+1) * scaling)
                return matrixPowers[power].map { $0 / factor }
            } else {
                return matrixPowers[power]
            }
        }
        
        result = MatrixExponentiator<T>.padeApprox(for: scaledMatrix, evenPowersOfM: scaledMatrixPowers, order: order)!
        
        let blockFormat = matrix.isSchur ? matrix.quasiTrianglularStructure() : nil
        
        if let blockFormat = blockFormat {
            MatrixExponentiator<T>.recomputeBlockDiag(scaledMatrix, exponentiated: &result, structure: blockFormat)
        }
        
        if (scaling > 0) {
            for _ in 0..<scaling {
                result = result * result
                if let blockFormat = blockFormat {
                    scaledMatrix = 2.0 * scaledMatrix
                    MatrixExponentiator<T>.recomputeBlockDiag(scaledMatrix, exponentiated: &result, structure: blockFormat)
                }
            }
        }
        
        return (result, scaling, order)
    }

    static func getExpmParams(from M: Matrix<T>) -> (Int, PadeApproximantOrder, [Matrix<T>]) {
        let isSmall = M.rows < 150
        
        // Need only even powers
        // 2 -> 1 -> 0
        // 4 -> 3 -> 1
        // 6 -> 5 -> 2
        var evenPowers = [Matrix<T>](repeating: M, count: 3)
        evenPowers[0] = M * M
        evenPowers[1] = evenPowers[0] * evenPowers[0]
        evenPowers[2] = evenPowers[0] * evenPowers[1]
        
        var scaling = 0
        
        let d4 = T.raise(evenPowers[1].manhattanNorm, to: 1.0/4.0)
        let d6 = T.raise(evenPowers[2].manhattanNorm, to: 1.0/6.0)
        let η1 = T.greater(of: d4, and: d6)
        
        if (η1 <= MatrixExpConst<T>.theta(for: .three)! && ell(M, coeff: MatrixExpConst<T>.coefficientsOfBackwardsErrorFunction[0], order: 3) == 0.0) {
            return (scaling, .three, evenPowers)
        }
        if (η1 <= MatrixExpConst<T>.theta(for: .five)! && ell(M, coeff: MatrixExpConst<T>.coefficientsOfBackwardsErrorFunction[1], order: 5) == 0.0) {
            return (scaling, .five, evenPowers)
        }
        
        let d8 = isSmall ? T.raise((evenPowers[1] * evenPowers[1]).manhattanNorm, to: 1.0/8.0) : T.raise(MatrixPowerOneNormEstimator(A: evenPowers[3], order: 2).compute(), to: 1.0/8.0)
        
        let η3 = T.greater(of: d6, and: d8)
        if (η3 <= MatrixExpConst<T>.theta(for: .seven)! && ell(M, coeff: MatrixExpConst<T>.coefficientsOfBackwardsErrorFunction[2], order: 7) == 0.0) {
            return (scaling, .seven, evenPowers)
        }
        if (η3 <= MatrixExpConst<T>.theta(for: .nine)! && ell(M, coeff: MatrixExpConst<T>.coefficientsOfBackwardsErrorFunction[3], order: 9) == 0.0) {
            return (scaling, .nine, evenPowers)
        }
        
        let d10 = isSmall ? T.raise((evenPowers[1] * evenPowers[1]).manhattanNorm, to: 1.0/10.0) : T.raise(MatrixPowerOneNormEstimator(A: evenPowers[0], order: 5).compute(), to: 1.0/10.0)
        
        let η4 = T.greater(of: d8, and: d10)
        let η5 = T.greater(of: η3, and: η4)
        scaling = getScaling(η4: η4, η5: η5)
        
        let factor = T.raiseTwo(to: scaling)
        let scaledM = M.map { $0 / factor }
        let sFromEll = ell(scaledM, coeff: MatrixExpConst<T>.coefficientsOfBackwardsErrorFunction[4], order: 13)
        
        scaling = sFromEll.isNaN ? revertToOldEstimate(M) : scaling + sFromEll.toInt
        
        return (scaling, .thirteen, evenPowers)
    }
    
    private static func getScaling(η4: T.Magnitude, η5: T.Magnitude) -> Int {
        let temp = ceil( T.Magnitude.log2( η5 / MatrixExpConst<T>.theta(for: .thirteen)! ) ).toInt
        return Int(max(temp, 0))
    }
    
    private static func revertToOldEstimate(_ matrix: Matrix<T>) -> Int {
        let twoToThePowerOfScaling = matrix.manhattanNorm / MatrixExpConst<T>.theta(for: .thirteen)!
        let t = T.Magnitude.log2(twoToThePowerOfScaling.significand)
        return Int(twoToThePowerOfScaling.exponent) - (t == 0.5 ? 1 : 0)
    }
    
    static func padeApprox(for M: Matrix<T>, evenPowersOfM: [Matrix<T>], order: PadeApproximantOrder) -> Matrix<T>? {
        let coeffs = MatrixExpConst<T>.padeCoefficients(order)
        let I = Matrix<T>.eye(M.rows)
        
        var U = I
        var V = I
        
        if (order == .thirteen) {
            let matrix2 = evenPowersOfM[0]
            let matrix4 = evenPowersOfM[1]
            let matrix6 = evenPowersOfM[2]
            
            let U1 = coeffs[13] * matrix6 + coeffs[11] * matrix4 + coeffs[9] * matrix2
            let U2 = coeffs[7] * matrix6 + coeffs[5] * matrix4 + coeffs[3] * matrix2
            U = M * (matrix6 * U1 + U2 + coeffs[1] * I)
            
            let V1 = coeffs[12] * matrix6 + coeffs[10] * matrix4 + coeffs[8] * matrix2
            let V2 = coeffs[6] * matrix6 + coeffs[4] * matrix4 + coeffs[2] * matrix2
            V = matrix6 * V1 + V2 + coeffs[0] * I
        } else if (order == .three || order == .five || order == .seven || order == .nine) {
            var evenPowers = evenPowersOfM
            if (order == .nine && evenPowersOfM.count == 3) {
                evenPowers.append(evenPowersOfM[0] * evenPowersOfM[2])
            }
            
            U = coeffs[1] * I
            V = coeffs[0] * I
            
            // k = 9, 7, 5, 3 -> index for matrix powers = 3, 2, 1, 0
            for k in stride(from: order.rawValue, through: 3, by: -2) {
                let index = (k-3)/2
                U = U + coeffs[k] * evenPowers[index]
                V = V + coeffs[k-1] * evenPowers[index]
            }
            
            U = M * U
        } else {
            // Should not come here
            return nil
        }
        
        guard let solve = (V-U).solve(2.0 * U) else {
            // Should not come here
            return nil
        }
        return solve + I
    }
    
    static func ell(_ matrix: Matrix<T>, coeff: T.Magnitude, order: Int) -> T.Magnitude {
        let factor = T.Magnitude.pow(coeff, 1.0 / T.Magnitude(magnitude: T.Magnitude(2 * order + 1)))
        let scaledMatrix = matrix.map { T(magnitude: $0.length * factor)  }
        let oneNormEstimator = MatrixPowerOneNormEstimator(A: scaledMatrix, order: (2 * order + 1))
        let alpha = oneNormEstimator.compute() / matrix.manhattanNorm
        let u = T.raise(2.0, to: -52.0)
        return T.Magnitude.maximum(ceil( T.Magnitude.log2(2.0 * alpha / u)/T.Magnitude(2 * order) ), 0.0)
    }
    
    static func recomputeBlockDiag(_ matrix: Matrix<T>, exponentiated: inout Matrix<T>, structure: [Int]) {
        let n = matrix.rows - 1
        
        for k in 0..<n {
            switch structure[k] {
            case 1:
                let t11 = matrix[k,k]
                let t22 = matrix[k+1,k+1]
                
                let avg = (t11+t22)/T.two
                let df = (t11-t22)/T.two
                
                var x12: T
                // Compare lengths because T can be Complex
                if max(avg.length as! Double, df.length as! Double) < log(Double.greatestFiniteMagnitude) {
                    x12 = matrix[k, k+1] * T.exp(avg) * T.sinch((t22 - t11) / T.two)
                } else {
                    x12 = matrix[k, k+1] * (T.exp(t22) - T.exp(t11)) / (t22 - t11)
                }
                exponentiated[k,k] = T.exp(t11)
                exponentiated[k,k+1] = x12
                exponentiated[k+1, k+1] = T.exp(t22)
            case 2:
                let a = matrix[k, k]
                let b = matrix[k, k+1]
                let c = matrix[k+1, k]
                let d = matrix[k+1, k+1]
                
                let avg = (a+d)/T.two
                let df = (a-d)/T.two
                
                let μ = T.two * T.two * (df * df + b * c)
                let delta = T.sqrt(μ)/T.two
                
                let expad2 = T.exp(avg)
                let coshdelta = T.cosh(delta)
                let sinchdelta = T.sinch(delta)
                
                exponentiated[k,k] = expad2 * ( coshdelta + df * sinchdelta)
                exponentiated[k,k+1] = expad2 * b * sinchdelta
                exponentiated[k+1,k] = expad2 * c * sinchdelta
                exponentiated[k+1, k+1] = expad2 * ( coshdelta - df * sinchdelta)
            default:
                continue
            }
        }
    }
    
}

