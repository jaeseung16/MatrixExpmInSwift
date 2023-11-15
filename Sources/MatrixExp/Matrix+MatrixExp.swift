//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 12/7/20.
//

import LANumerics

extension Matrix {
    public var isSquare: Bool {
        return self.rows == self.columns
    }
    
    public var isDiag: Bool {
        var isDiag = false
        if self.isSquare {
            var diag = Vector<Element>()
            for k in 0..<self.columns {
                diag.append(self[k,k])
            }
            isDiag = Matrix(diagonal: diag) == self
        }
        return isDiag
    }
    
    public var isScalar: Bool {
        return self.rows == 1 && self.columns == 1
    }
    
    public var isHermitian: Bool {
        return self == self.adjoint
    }
    
    public var isSchur: Bool {
        var result: Bool
        if (self.isScalar) {
            result = true
        } else if (!self.isSquare) {
            result = false
        } else if (self is Matrix<Double> || self is Matrix<Float>) {
            result = self.isQuasiUpperTriangle
        } else {
            result = self.isUpperTriangle
        }
        return result
    }
}
