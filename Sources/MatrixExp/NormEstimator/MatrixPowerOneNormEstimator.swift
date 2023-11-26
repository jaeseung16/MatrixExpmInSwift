//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 11/23/23.
//

import Foundation
import Numerics
import LANumerics

class MatrixPowerOneNormEstimator<T>: OneNormEstimator<T> where T: Exponentiable, T.Magnitude: Real {
    let order: Int
    
    init(A: Matrix<T>, order: Int = 1, t: Int = 2, X0: Matrix<T> = Matrix<T>(), toPrint: Bool = false) {
        self.order = order
        super.init(A: A, t: t, X0: X0, toPrint: toPrint)
    }
    
    override func compute() {
        if dimension < 50 {
            estimate = multiply(by: Matrix<T>.eye(dimension)).manhattanNorm
        } else if isNonNegative(A) {
            let e = Matrix<T>(Vector<T>(repeating: 1.0, count: dimension))
            estimate = conjugateTransposeAndMultiply(by: e).infNorm
        } else {
            super.compute()
        }
    }
    
    func isNonNegative(_ M: Matrix<T>) -> Bool {
        return M.forall {
            ($0 is Double && $0 as! Double >= 0) || ($0 is Float && $0 as! Float >= 0)
        }
    }
    
    override func computeExactly() {
        var M = A
        for _ in 0..<(order-1) {
            M = A * M
        }
        super.computeExactly(M)
    }
    
    override var dimension: Int {
        return n
    }
    
    override var isReal: Bool {
        return T.self is (any Real)
    }
    
    override func multiply(by S: Matrix<T>) -> Matrix<T> {
        var X = S
        for _ in 0..<order {
            X = A * X
        }
        return X
    }
    
    override func conjugateTransposeAndMultiply(by S: Matrix<T>) -> Matrix<T> {
        var X = S
        for _ in 0..<order {
            X = A.adjoint * X
        }
        return X
    }
}
