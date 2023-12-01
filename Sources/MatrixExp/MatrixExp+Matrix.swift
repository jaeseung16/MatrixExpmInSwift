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
            let range = 0..<self.columns
            isDiag = Matrix(diagonal: range.map { self[$0,$0] }) == self
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

extension Matrix {
    func quasiTrianglularStructure() -> [Int] {
        let zero = Element.zero
        
        guard self.rows > 1 else {
            return [0]
        }
        
        guard self.rows > 2 else {
            return self[1,0] == zero ? [1] : [2]
        }
        
        var structure = [Int]()
        var k = 0
        while k < (self.rows - 2) {
            if self[k+1,k] != zero {
                structure.append(2)
                structure.append(0)
                k = k + 2
            } else if self[k+1, k] == zero && self[k+2, k+1] == zero {
                structure.append(1)
                k = k + 1
            } else {
                structure.append(0)
                k = k + 1
            }
        }
        
        if self[self.rows-1, self.rows-2] != zero {
            structure.append(2)
        } else if structure.last == 0 || structure.last == 1 {
            structure.append(1)
        }
        
        return structure
    }
}
