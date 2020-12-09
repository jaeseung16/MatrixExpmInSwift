//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 12/9/20.
//

import Foundation
import Numerics
import LANumerics

class NormOfMatrixPower {
    let M: Matrix<Double>
    let power: Int
    let estimatedNorm: Double
    let numberOfMatrixVectorProduct: Int
    
    init(_ M: Matrix<Double>, power: Int) {
        self.M = M
        self.power = power
        
        (self.estimatedNorm, self.numberOfMatrixVectorProduct) = NormOfMatrixPower.norm(of: M, power: power)
    }
    
    static func norm(of M: Matrix<Double>, power: Int) -> (Double, Int) {
        var estimatedNorm: Double
        var numberOfMatrixVectorProduct: Int
        var Mpower = Matrix<Double>.eye(M.rows)
        
        if M.rows < 50 || !isNonNegative(M) {
            for _ in 0..<power {
                Mpower = M * Mpower
            }
            estimatedNorm = Mpower.manhattanNorm
            numberOfMatrixVectorProduct = 0
        } else if isNonNegative(M) {
            var e = Matrix<Double>(Vector<Double>(repeating: 1.0, count: M.rows))
            for _ in 0..<power {
                e = M.transpose * e
            }
            estimatedNorm = e.infNorm
            numberOfMatrixVectorProduct = power
        } else {
            let normEst = NormEst1<Double>(A: M, order: power)
            estimatedNorm = normEst.estimate
            numberOfMatrixVectorProduct = 2 * power * normEst.numberOfProducts
        }
        return (estimatedNorm, numberOfMatrixVectorProduct)
    }
    
    static func isNonNegative(_ M: Matrix<Double>) -> Bool {
        return M.forall {$0 >= 0}
    }
}
