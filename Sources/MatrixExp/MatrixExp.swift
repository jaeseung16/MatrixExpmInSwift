import Foundation
import Numerics
import LANumerics

public class MatrixExp<Type> where Type: Exponentiable, Type.Magnitude: Real {
    public let matrix: Matrix<Type>
    public let result: Matrix<Type>?
    public let scaling: Int
    public let orderPadeApproximant: Int
    public let calculationType: MatrixExpCalculationType
    
    public init(_ matrix: Matrix<Type>) {
        self.matrix = matrix
        
        if matrix.isSquare {
            if matrix.isDiag {
                self.result = MatrixExp<Type>.exp(diagMatrix: matrix)
                calculationType = .diag
                self.scaling = 0
                self.orderPadeApproximant = 0
            } else if matrix.isHermitian {
                self.result = MatrixExp<Type>.exp(hermitianMatrix: matrix)
                calculationType = .hermitian
                self.scaling = 0
                self.orderPadeApproximant = 0
            } else {
                (self.result, self.scaling, self.orderPadeApproximant) = MatrixExp<Type>.exp(matrix: matrix)
                calculationType = .pade
            }
        } else {
            self.result = nil
            self.scaling = 0
            self.orderPadeApproximant = 0
            calculationType = .notApplicable
        }
    }
    
    static func exp(diagMatrix: Matrix<Type>) -> Matrix<Type> {
        var result = Matrix<Type>.zeros(diagMatrix.rows, diagMatrix.columns)
        for k in 0..<diagMatrix.columns {
            result[k, k] = diagMatrix[k, k].exponentiation()
        }
        return result
    }
    
    static func exp(hermitianMatrix: Matrix<Type>) -> Matrix<Type> {
        let (eigenValues, schurForm, schurVectors) = hermitianMatrix.schur()!
        
        var expSchurForm = Matrix<Type>.eye(hermitianMatrix.rows)
        for k in 0..<schurForm.columns {
            expSchurForm[k, k] = schurForm[k, k].exponentiation()
        }
        
        return schurVectors * expSchurForm * schurVectors.adjoint
    }
    
    static func exp(matrix: Matrix<Type>) -> (Matrix<Type>, Int, Int) {
        let recomputeDiags = matrix.isSchur
        let copiedMatrix = matrix
        var result = Matrix<Type>.zeros(matrix.rows, matrix.columns)
        
        let blockFormat = recomputeDiags ? quasiTrianglularStructure(copiedMatrix) : nil
        
        let (scaling, order, Mpowers) = expmParams(for: copiedMatrix)
        
        let factor = scaling > 0 ? convertToType(floatLiteral: pow(2.0, Double(scaling))) : 1.0
        var scaledM = copiedMatrix.map { $0 / factor }
        
        var scaledMpowers = Mpowers
        if (scaling > 0) {
            for k in 0..<Mpowers.count {
                let factor = convertToType(floatLiteral: pow(2.0, Double((k+1) * scaling)))
                print("k = \(k), factor = \(factor)")
                scaledMpowers[k] = Mpowers[k].map { $0 / factor }
            }
        }
        
        result = padeApprox(for: scaledM, Mpowers: scaledMpowers, order: order)!
        
        if (recomputeDiags) {
            recomputeBlockDiag(scaledM, exponentiated: &result, structure: blockFormat!)
        }
        
        if (scaling > 0) {
            for _ in 0..<scaling {
                result = result * result
                if (recomputeDiags) {
                    scaledM = 2.0 * scaledM
                    recomputeBlockDiag(scaledM, exponentiated: &result, structure: blockFormat!)
                }
            }
        }
        
        return (result, scaling, order)
    }
    
    static func sinch(_ x: Type) -> Type {
        var value: Type
        if x == convertToType(floatLiteral: 0.0) {
            value = convertToType(floatLiteral: 1.0)
        } else {
            value = (x.exponentiation() - (-x).exponentiation()) / x / convertToType(floatLiteral: 2.0)
        }
        return value
    }
    
    static func expmParams(for M: Matrix<Type>) -> (Int, Int, [Matrix<Type>]) {
        let isSmall = M.rows < 150
        
        var Mpowers = [Matrix<Type>](repeating: M, count: 6)
        Mpowers[1] = M * M
        Mpowers[3] = Mpowers[1] * Mpowers[1]
        Mpowers[5] = Mpowers[1] * Mpowers[3]
        
        var scaling = 0
        var orderPadeApproximant = 0
        
        let d4 = pow(Mpowers[3].manhattanNorm as! Double, 1.0/4.0)
        let d6 = pow(Mpowers[5].manhattanNorm as! Double, 1.0/6.0)
        let eta1 = max(d4, d6)
        
        if (eta1 <= MatrixExpConst<Double>.theta(for: 3)! && ell(M, coeff: MatrixExpConst<Double>.coefficientsOfBackwardsErrorFunction[0], order: 3) == 0.0) {
            orderPadeApproximant = 3
            return (scaling, orderPadeApproximant, Mpowers)
        }
        if (eta1 <= MatrixExpConst<Double>.theta(for: 5)! && ell(M, coeff: MatrixExpConst<Double>.coefficientsOfBackwardsErrorFunction[1], order: 5) == 0.0) {
            orderPadeApproximant = 5
            return (scaling, orderPadeApproximant, Mpowers)
        }
        
        var d8: Double
        if (isSmall) {
            d8 = pow((Mpowers[3] * Mpowers[3]).manhattanNorm as! Double, 1.0/8.0)
        } else {
            let normest = NormEst1<Type>(A: Mpowers[3], order: 2)
            d8 = pow(normest.estimate, 1.0/8.0)
        }
        
        let eta3 = max(d6, d8)
        if (eta3 <= MatrixExpConst<Double>.theta(for: 7)! && ell(M, coeff: MatrixExpConst<Double>.coefficientsOfBackwardsErrorFunction[2], order: 7) == 0.0) {
            orderPadeApproximant = 7
            return (scaling, orderPadeApproximant, Mpowers)
        }
        if (eta3 <= MatrixExpConst<Double>.theta(for: 9)! && ell(M, coeff: MatrixExpConst<Double>.coefficientsOfBackwardsErrorFunction[3], order: 9) == 0.0) {
            orderPadeApproximant = 9
            return (scaling, orderPadeApproximant, Mpowers)
        }
        
        var d10: Double
        if (isSmall) {
            d10 = pow((Mpowers[3] * Mpowers[5]).manhattanNorm as! Double, 1.0/10.0)
        } else {
            let normest = NormEst1<Type>(A: Mpowers[1], order: 5)
            d10 = pow(normest.estimate, 1.0/10.0)
        }

        let eta4 = max(d8, d10)
        let eta5 = min(eta3, eta4)
        scaling = Int(max( ceil( log2( eta5 / MatrixExpConst<Double>.theta(for: 13)! ) ), 0))
        
        let factor = convertToType(floatLiteral: pow(2.0, Double(scaling)))
        let scaledM = M.map { $0 / factor }
        let sFromEll = ell(scaledM, coeff: MatrixExpConst<Double>.coefficientsOfBackwardsErrorFunction[4], order: 5)
        
        if (sFromEll.isNaN) {
            let norm1 = M.manhattanNorm as! Double
            let theta = MatrixExpConst<Type>.theta(for: 13)!
            
            let needAName = norm1/theta
            
            let t = log2(needAName.significand)
            scaling = needAName.exponent - (t == 0.5 ? 1 : 0)
        } else {
            scaling = scaling + Int(sFromEll)
        }
        orderPadeApproximant = 13

        return (scaling, orderPadeApproximant, Mpowers)
    }
    
    static func padeApprox(for M: Matrix<Type>, Mpowers: [Matrix<Type>], order: Int) -> Matrix<Type>? {
        let coeffs = MatrixExpConst<Type>.padeCoefficients(for: order)!
        let I = Matrix<Type>.eye(M.rows)
        var F = Matrix<Type>.eye(M.rows)
        var Mpowers2 = Mpowers
        
        var U = Matrix<Type>.eye(M.rows)
        var V = Matrix<Type>.eye(M.rows)
        
        if (order == 13) {
            let U1 = coeffs[13] * Mpowers[5] + coeffs[11] * Mpowers[3] + coeffs[9] * Mpowers[1]
            let U2 = coeffs[7] * Mpowers[5] + coeffs[5] * Mpowers[3] + coeffs[3] * Mpowers[1]
            U = M * (Mpowers[5] * U1 + U2 + coeffs[1] * I)
            
            let V1 = coeffs[12] * Mpowers[5] + coeffs[10] * Mpowers[3] + coeffs[8] * Mpowers[1]
            let V2 = coeffs[6] * Mpowers[5] + coeffs[4] * Mpowers[3] + coeffs[2] * Mpowers[1]
            V = Mpowers[5] * V1 + V2 + coeffs[0] * I
        } else if (order == 3 || order == 5 || order == 7 || order == 9) {
            if (order == 9 && Mpowers.count == 6) {
                Mpowers2.append(Matrix<Type>.eye(M.rows))
                Mpowers2.append(Mpowers[1] * Mpowers[5])
            }
            
            U = coeffs[1] * Matrix<Type>.eye(M.rows)
            V = coeffs[0] * Matrix<Type>.eye(M.rows)
            
            for k in stride(from: order, through: 3, by: -2) {
                U = U + coeffs[k] * Mpowers2[k-2]
                V = V + coeffs[k-1] * Mpowers2[k-2]
            }
            
            U = M * U
        } else {
            return nil
        }
        
        guard let solve = (V-U).solve(2.0 * U) else {
            return nil
        }
        
        F = solve + I
        
        return F
    }
    
    static func ell(_ matrix: Matrix<Type>, coeff: Double, order: Int) -> Double {
        let factor = pow(coeff, 1.0 / Double(2 * order + 1))
        let scaledMatrix = matrix.map { ($0.length as! Double) * factor}
        let alpha = NormOfMatrixPower(scaledMatrix, power: (2 * order + 1)).estimatedNorm / (matrix.manhattanNorm as! Double)
        let u = pow(2.0, -52.0)
        return max(ceil(log2(2.0 * alpha / u) / Double(2 * order)), 0.0)
    }
    
    static func recomputeBlockDiag(_ matrix: Matrix<Type>, exponentiated: inout Matrix<Type>, structure: [Int]) {
        let n = matrix.rows - 1
        let two = convertToType(floatLiteral: 2.0)
        
        for k in 0..<n {
            switch structure[k] {
            case 1:
                let t11 = matrix[k, k]
                let t22 = matrix[k+1, k+1]
                
                let ave = (t11 + t22) / two
                let df = (t11 - t22) / two
                
                var x12: Type
                // Compare lengths because Type can be Complex
                if (max(ave.length as! Double, df.length as! Double) < log(Double.greatestFiniteMagnitude)) {
                    let factor = ave.exponentiation() * sinch((t22 - t11) / two)
                    x12 = matrix[k, k+1] * factor
                } else {
                    let factor = (t22.exponentiation() - t11.exponentiation()) / (t22 - t11)
                    x12 = matrix[k, k+1] * factor
                }
                exponentiated[k,k] = t11.exponentiation()
                exponentiated[k,k+1] = x12
                exponentiated[k+1, k+1] = t22.exponentiation()
            case 2:
                let a = matrix[k, k]
                let b = matrix[k, k+1]
                let c = matrix[k+1, k]
                let d = matrix[k+1, k+1]
                
                let ave = (a + d) / two
                let df = (a - d) / two
                
                let μ = two * two * (df * df + b * c)
                let delta = μ.squareRoot() / two
                let expad2 = ave.exponentiation()
                let coshdelta = ( delta.exponentiation() + (-delta).exponentiation() ) / two
                let sinchdelta = sinch(delta)
                
                exponentiated[k,k] = expad2 * ( coshdelta + df * sinchdelta)
                exponentiated[k,k+1] = expad2 * b * sinchdelta
                exponentiated[k+1,k] = expad2 * c * sinchdelta
                exponentiated[k+1, k+1] = expad2 * ( coshdelta - df * sinchdelta)
                
            default:
                continue
            }
        }
    }
    
    static func quasiTrianglularStructure(_ matrix: Matrix<Type>) -> [Int] {
        let zero = convertToType(floatLiteral: 0.0)
        var structure = [Int]()
        
        print("matrix.rows = \(matrix.rows)")
        if (matrix.rows == 1) {
            structure.append(0)
            return structure
        } else if (matrix.rows == 2) {
            if (matrix[1, 0] == zero) {
                structure.append(1)
                return structure
            } else {
                structure.append(2)
                return structure
            }
        }
        
        var k = 0
        
        while (k < (matrix.rows - 2)) {
            if (matrix[k+1,k] != zero) {
                structure.append(2)
                structure.append(0)
                k = k + 2
            } else if (matrix[k+1, k] == zero && matrix[k+2, k+1] == zero) {
                structure.append(1)
                k = k + 1
            } else {
                structure.append(0)
                k = k + 1
            }
        }
        
        if (matrix[matrix.rows-1, matrix.rows-2] != zero) {
            structure.append(2)
        } else if (structure.last == 0 || structure.last == 1) {
            structure.append(1)
        }
        
        return structure
    }
    
    static func convertToType(floatLiteral: Double) -> Type {
        return Type(floatLiteral: floatLiteral as! Type.FloatLiteralType)
    }
}

