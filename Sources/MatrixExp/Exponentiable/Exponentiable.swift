//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 10/13/20.
//

import Numerics
import LANumerics

public protocol Exponentiable: LANumeric {
    func exponentiation() -> Self
    func squareRoot() -> Self
}

extension Complex: Exponentiable where RealType: LANumeric {
    public func exponentiation() -> Complex<RealType> {
        let exponentiationOfRealPart = Complex<RealType>(RealType.exp(real))
        let exponentiationOfImaginaryPart = Complex<RealType>(RealType.cos(imaginary), RealType.sin(imaginary))
        return exponentiationOfRealPart * exponentiationOfImaginaryPart
    }
    
    public func squareRoot() -> Complex<RealType> {
        return Complex<RealType>(length: length.squareRoot(), phase: phase/2.0)
    }
}

extension Double: Exponentiable {
    public func exponentiation() -> Double {
        return Double.exp(self)
    }
    
    public func squareRoot() -> Double {
        return Double.sqrt(self)
    }
}

extension Float: Exponentiable {
    public func exponentiation() -> Float {
        return Float.exp(self)
    }
    
    public func squareRoot() -> Float {
        return Float.sqrt(self)
    }
}
