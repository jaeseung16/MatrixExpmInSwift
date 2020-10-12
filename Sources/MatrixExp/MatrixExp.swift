import Foundation
import Numerics
import RealModule
import LANumerics

class MatrixExp {
    var text = "Hello, World!"
    
    static func evaluate(for matrix: Matrix<Double>) -> Matrix<Double>? {
        guard isSquare(matrix) else {
            NSLog("Need a square matrix as an input: \(matrix) ")
            return nil
        }
        
        var result = Matrix<Double>.zeros(matrix.rows, matrix.columns)
        
        if (isDiag(matrix)) {
            if (matrix.rows == matrix.columns) {
                for k in 0..<matrix.columns {
                    result[k, k] = exp(matrix[k, k])
                }
            }
        }
        
        return result
    }
    
    static func isDiag(_ matrix: Matrix<Double>) -> Bool {
        var isDiag = false
        
        var diag = Vector<Double>()
        
        if (matrix.rows == matrix.columns) {
            for k in 0..<matrix.columns {
                diag.append(matrix[k,k])
            }
            isDiag = Matrix<Double>(diagonal: diag) == matrix
        }
        
        return isDiag
        
    }
    
    static func isSquare(_ matrix: Matrix<Double>) -> Bool {
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
}

