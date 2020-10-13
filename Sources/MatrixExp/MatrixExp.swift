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
            //print("scaling = \(scaling)")
            //print("order = \(order)")
            //print("Mpowers = \(Mpowers)")
            result = padeApprox(for: matrix, Mpowers: Mpowers, order: 13)!
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
    
    static func sinch(_ x: Double) -> Double {
        var value: Double
        if (x == 0.0) {
            value = 1
        } else {
            value = sinh(x) / x
        }
        return x
    }
    
    static func ell(T: Matrix<Double>, coeff: Double, order: Int) -> Double {
        return 0.0
    }
    
    static func norm(of M: Matrix<Double>, power: Int) -> (Double, Int) {
        var norm = 0.0
        var mv = 0
        var Mpower = M
        
        if (M.columns < 50) {
            for _ in 1..<power {
                Mpower = Mpower * M
            }
            norm = Mpower.manhattanNorm
        } else if (M.forall { $0 >= 0 }) {
            mv = power
        } else {
            
        }
        return (norm, mv)
    }
    
    static func expmParams(for M: Matrix<Type>) -> (Int, Int, [Matrix<Type>]) {
        // Use old estimate first
        // TODO: Implement the logic for the smaller order: 3, 5, 7, and 9
        
        var Mpowers = [Matrix<Type>](repeating: M, count: 6)
        Mpowers[1] = M * M
        Mpowers[3] = Mpowers[1] * Mpowers[1]
        Mpowers[5] = Mpowers[1] * Mpowers[3]
        
        let order = 13
        let norm1 = M.manhattanNorm as! Double
        let theta = MatrixExpConst<Type>.theta(for: order)!
        
        let needAName = norm1/theta
        
        let t = log2(needAName.significand)
        let s = needAName.exponent - (t == 0.5 ? 1 : 0)
        
        return (s, order, Mpowers)
    }
    
    static func padeApprox(for M: Matrix<Type>, Mpowers: [Matrix<Type>], order: Int) -> Matrix<Type>? {
        let coeffs = MatrixExpConst<Type>.padeCoefficients(for: order)!
        let I = Matrix<Type>.eye(M.rows)
        var F = Matrix<Type>.eye(M.rows)
        
        print("F = \(F)")
        
        if (order == 13) {
            let U1 = coeffs[13] * Mpowers[5] + coeffs[11] * Mpowers[3] + coeffs[9] * Mpowers[1]
            let U2 = coeffs[7] * Mpowers[5] + coeffs[5] * Mpowers[3] + coeffs[3] * Mpowers[1]
            let U = M * (Mpowers[5] * U1 + U2 + coeffs[1] * I)
            
            let V1 = coeffs[12] * Mpowers[5] + coeffs[10] * Mpowers[3] + coeffs[8] * Mpowers[1]
            let V2 = coeffs[6] * Mpowers[5] + coeffs[4] * Mpowers[3] + coeffs[2] * Mpowers[1]
            let V = Mpowers[5] * V1 + V2 + coeffs[0] * I
            
            print("U = \(U)")
            print("V = \(V)")
            
            guard let solve = (V-U).solve(2.0 * U) else {
                print("returning nil")
                return nil
            }
            
            print("solve = \(solve)")
            F = solve + I
        }
        
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
    
}

