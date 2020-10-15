import Foundation
import Numerics
import RealModule
import LANumerics

class MatrixExp<Type> where Type: Exponentiable {
    var text = "Hello, World!"
    
    static func evaluate(for matrix: Matrix<Type>) -> Matrix<Type>? {
        guard isSquare(matrix) else {
            NSLog("Need a square matrix as an input: \(matrix) ")
            return nil
        }
        
        print("input = \(matrix)")
        var result = Matrix<Type>.zeros(matrix.rows, matrix.columns)
        print("\(result)")
        
        if (isDiag(matrix)) {
            print("diag")
            for k in 0..<matrix.columns {
                result[k, k] = matrix[k, k].exponentiation()
            }
        } else {
            let (scaling, order, Mpowers) = expmParams(for: matrix)
            print("order = \(order)")
            print("scaling = \(scaling)")
            print("Mpowers = \(Mpowers)")
            
            let factor = scaling > 0 ? Type(floatLiteral: pow(2.0, Double(scaling)) as! Type.FloatLiteralType) : 1.0
            let scaledM = matrix.map { $0 / factor }
            
            var scaledMpowers = Mpowers
            
            if (scaling > 0) {
                for k in 0..<Mpowers.count {
                    let factor = Type(floatLiteral: pow(2.0, Double((k+1) * scaling)) as! Type.FloatLiteralType)
                    print("k = \(k), factor = \(factor)")
                    scaledMpowers[k] = Mpowers[k].map { $0 / factor }
                }
            }
            
            print("scaledM = \(scaledM)")
            print("scaledMpowers = \(scaledMpowers)")
            
            result = padeApprox(for: scaledM, Mpowers: scaledMpowers, order: order)!
            print("result = \(result)")
            
            if (scaling > 0) {
                for _ in 0..<scaling {
                    result = result * result
                    print("result = \(result)")
                }
            }
        }
        
        print("output = \(result)")
        return result
    }
    
    static func isDiag(_ matrix: Matrix<Type>) -> Bool {
        var isDiag = false
        
        var diag = Vector<Type>()
        
        if (matrix.rows == matrix.columns) {
            for k in 0..<matrix.columns {
                diag.append(matrix[k,k])
            }
            isDiag = Matrix<Type>(diagonal: diag) == matrix
        }
        
        return isDiag
        
    }
    
    static func isSquare(_ matrix: Matrix<Type>) -> Bool {
        return matrix.rows == matrix.columns
    }
    
    static func sinch(_ x: Type) -> Type {
        var value: Type
        if (x == 0.0) {
            value = Type(floatLiteral: 1.0 as! Type.FloatLiteralType)
        } else {
            value = (x.exponentiation() - (-x).exponentiation()) / x / Type(floatLiteral: 2.0 as! Type.FloatLiteralType)
        }
        return x
    }
    
    static func expmParams(for M: Matrix<Type>) -> (Int, Int, [Matrix<Type>]) {
        // Use old estimate first
        // TODO: Implement the logic for the smaller order: 3, 5, 7, and 9
        
        var Mpowers = [Matrix<Type>](repeating: M, count: 6)
        Mpowers[1] = M * M
        Mpowers[3] = Mpowers[1] * Mpowers[1]
        Mpowers[5] = Mpowers[1] * Mpowers[3]
        
        var s = 0
        var order = 0
        
        let d4 = pow(Mpowers[3].manhattanNorm as! Double, 1.0/4.0)
        let d6 = pow(Mpowers[5].manhattanNorm as! Double, 1.0/6.0)
        let eta1 = max(d4, d6)
        
        if (eta1 <= MatrixExpConst<Double>.theta(for: 3)! && ell(M, coeff: MatrixExpConst<Double>.coefficientsOfBackwardsErrorFunction[0], order: 3)! == 0.0) {
            order = 3
            return (s, order, Mpowers)
        }
        if (eta1 <= MatrixExpConst<Double>.theta(for: 5)! && ell(M, coeff: MatrixExpConst<Double>.coefficientsOfBackwardsErrorFunction[1], order: 5)! == 0.0) {
            order = 5
            return (s, order, Mpowers)
        }
        
        let d8 = pow((Mpowers[3] * Mpowers[3]).manhattanNorm as! Double, 1.0/8.0)
        let eta3 = max(d6, d8)
        if (eta3 <= MatrixExpConst<Double>.theta(for: 7)! && ell(M, coeff: MatrixExpConst<Double>.coefficientsOfBackwardsErrorFunction[2], order: 7)! == 0.0) {
            order = 7
            return (s, order, Mpowers)
        }
        if (eta3 <= MatrixExpConst<Double>.theta(for: 9)! && ell(M, coeff: MatrixExpConst<Double>.coefficientsOfBackwardsErrorFunction[3], order: 9)! == 0.0) {
            order = 9
            return (s, order, Mpowers)
        }
        
        let d10 = pow((Mpowers[3] * Mpowers[5]).manhattanNorm as! Double, 1.0/10.0)
        let eta4 = max(d8, d10)
        let eta5 = min(eta3, eta4)
        s = Int(max( ceil( log2( eta5 / MatrixExpConst<Double>.theta(for: 13)! ) ), 0))
        
        let factor = Type(floatLiteral: pow(2.0, Double(s)) as! Type.FloatLiteralType)
        let scaledM = M.map { $0 / factor }
        let sFromEll = ell(scaledM, coeff: MatrixExpConst<Double>.coefficientsOfBackwardsErrorFunction[4], order: 5)!
        
        if (sFromEll.isNaN) {
            let norm1 = M.manhattanNorm as! Double
            let theta = MatrixExpConst<Type>.theta(for: 13)!
            
            let needAName = norm1/theta
            
            let t = log2(needAName.significand)
            s = needAName.exponent - (t == 0.5 ? 1 : 0)
        } else {
            s = s + Int(sFromEll)
        }
        order = 13

        return (s, order, Mpowers)
    }
    
    static func padeApprox(for M: Matrix<Type>, Mpowers: [Matrix<Type>], order: Int) -> Matrix<Type>? {
        let coeffs = MatrixExpConst<Type>.padeCoefficients(for: order)!
        let I = Matrix<Type>.eye(M.rows)
        var F = Matrix<Type>.eye(M.rows)
        var Mpowers2 = Mpowers
        
        var U = Matrix<Type>.eye(M.rows)
        var V = Matrix<Type>.eye(M.rows)
        
        print("F = \(F)")
        
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
        
        print("U = \(U)")
        print("V = \(V)")
        
        guard let solve = (V-U).solve(2.0 * U) else {
            print("returning nil")
            return nil
        }
        
        print("solve = \(solve)")
        F = solve + I
        
        return F
    }
    
    static func isSchur(_ matrix: Matrix<Type>) -> Bool {
        var result: Bool
        
        if (isScalar(matrix)) {
            result = true
        } else if (!isSquare(matrix)) {
            result = false
        } else if (matrix is Matrix<Double> || matrix is Matrix<Float>) {
            result = matrix.isQuasiUpperTriangle
        } else {
            result = matrix.isUpperTriangle
        }
        
        return result
    }
    
    static func isScalar(_ matrix: Matrix<Type>) -> Bool {
        return matrix.rows == 1 && matrix.columns == 1
    }
    
    static func isHermitian(_ matrix: Matrix<Type>) -> Bool {
        return matrix == matrix.adjoint
    }
    
    static func ell(_ matrix: Matrix<Type>, coeff: Double, order: Int) -> Double? {
        guard let realValueOfCoeff = coeff.length as? Double else {
            print("coeff should be real: coeff = \(coeff)")
            return nil
        }
        
        let factor = pow(realValueOfCoeff, 1.0 / Double(2 * order + 1))
        let scaledMatrix = matrix.map { ($0.length as! Double) * factor}
        let alpha = norm(of: scaledMatrix, power: (2 * order + 1))! / (matrix.manhattanNorm as! Double)
        
        let u = pow(2.0, -52.0)

        return max(ceil(log2(2.0 * alpha / u) / Double(2 * order)), 0.0)
    }
    
    static func norm(of M: Matrix<Double>, power: Int) -> Double? {
        var estimatedNorm: Double?
        var Mpower = Matrix<Double>.eye(M.rows)
        
        if (M.rows < 50 || !isNonNegative(M)) {
            for _ in 0..<power {
                Mpower = M * Mpower
            }
            estimatedNorm = Mpower.manhattanNorm
        } else if (isNonNegative(M)) {
            var e = Matrix<Double>(Vector<Double>(repeating: 1.0, count: M.rows))
            for k in 0..<power {
                e = M.transpose * e
            }
            estimatedNorm = e.infNorm
        } else {
            // MATLAB normAM has an implementation for this case
            // However, we know that M is real and positive for matrix exponentiation
            // For now, return nil
            estimatedNorm = nil
        }
        
        return estimatedNorm
    }
    
    static func isNonNegative(_ M: Matrix<Double>) -> Bool {
        return M.forall {$0 >= 0}
    }
    
    static func recomputeBlockDiag(_ matrix: Matrix<Type>, exponentiated: inout Matrix<Type>, structure: [Int]) {
        let n = matrix.rows - 1
        let two = Type(floatLiteral: 2.0 as! Type.FloatLiteralType)
        
        for k in 0..<n {
            switch structure[k] {
            case 1:
                let t11 = matrix[k, k]
                let t22 = matrix[k+1, k+1]
                
                let ave = (t11 + t22) / two
                let df = (t11 - t22) / two
                
                var x12: Type
                if (max(ave.length as! Double, df.length as! Double) < Double.greatestFiniteMagnitude) {
                    let factor = ave.exponentiation() * sinch((t22 - t11) / two)
                    x12 = matrix[k, k+1] * factor
                } else {
                    let factor = (t22.exponentiation() - t11.exponentiation()) / (t22 - t11)
                    x12 = matrix[k, k+1] * factor
                }
                exponentiated[k,k] = t11.exponentiation()
                exponentiated[k,k+1] = x12.exponentiation()
                exponentiated[k+1, k+1] = t22.exponentiation()
            case 2:
                let a = matrix[k, k]
                let b = matrix[k, k+1]
                let c = matrix[k+1, k]
                let d = matrix[k+1, k+1]
                
                let ave = (a + d) / two
                let df = (a - d) / two
                
                let μ = df * df + two * two * b * c
                let delta = μ.squareRoot() / two
                let expad2 = ( ave / two ).exponentiation()
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
        let zero = Type(floatLiteral: 0.0 as! Type.FloatLiteralType)
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
        
        print("structure.count = \(structure.count)")
        
        if (matrix[matrix.rows-1, matrix.rows-2] != zero) {
            structure.append(2)
        } else if (structure.last == 0 || structure.last == 1) {
            structure.append(1)
        }
        
        return structure
    }
}

