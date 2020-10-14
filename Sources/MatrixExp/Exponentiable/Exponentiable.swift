//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 10/13/20.
//

import Foundation
import Numerics
import LANumerics

protocol Exponentiable: LANumeric {
    func exponentiation() -> Self
    func squareRoot() -> Self
}

extension Complex: Exponentiable where RealType: LANumeric {
    func exponentiation() -> Complex<RealType> {
        let exponentiationOfRealPart = Complex<RealType>(RealType.exp(self.real))
        let exponentiationOfImaginaryPart = Complex<RealType>(RealType.cos(self.imaginary), RealType.sin(self.imaginary))
        return exponentiationOfRealPart * exponentiationOfImaginaryPart
    }
    
    func squareRoot() -> Complex<RealType> {
        return Complex<RealType>(length: self.length.squareRoot(), phase: self.phase/2.0)
    }
}

extension Double: Exponentiable {
    func exponentiation() -> Double {
        return Double.exp(self)
    }
    
    func squareRoot() -> Double {
        return Double.sqrt(self)
    }
}

extension Float: Exponentiable {
    func exponentiation() -> Float {
        return Float.exp(self)
    }
    
    func squareRoot() -> Float {
        return Float.sqrt(self)
    }
}
