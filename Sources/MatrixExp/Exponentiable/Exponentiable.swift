//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 10/13/20.
//

import Foundation
import Numerics
import LANumerics

protocol Exponentiable: MatrixElement, LANumeric {
    func exponentiation() -> Self
}

extension Complex: Exponentiable where RealType: LANumeric {
    func exponentiation() -> Complex<RealType> {
        let exponentiationOfRealPart = Complex<RealType>(RealType.exp(self.real))
        let exponentiationOfImaginaryPart = Complex<RealType>(RealType.cos(self.imaginary), RealType.sin(self.imaginary))
        return exponentiationOfRealPart * exponentiationOfImaginaryPart
    }
}

extension Double: Exponentiable {
    func exponentiation() -> Double {
        return Double.exp(self)
    }
}

extension Float: Exponentiable {
    func exponentiation() -> Float {
        return Float.exp(self)
    }
}
