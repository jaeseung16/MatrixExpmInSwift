//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 10/13/20.
//

import Numerics
import LANumerics

public protocol Exponentiable: LANumeric, ElementaryFunctions {
    static var one: Self { get }
    static var two: Self { get }

    static func sinch(_ x: Self) -> Self
    static func raiseTwo(to: Int) -> Self
    
    static func raise(_ a: Self.Magnitude, to b: Self.Magnitude) -> Self.Magnitude
    static func greater(of a: Self.Magnitude, and b: Self.Magnitude) -> Self.Magnitude
}

extension Exponentiable {
    public static func sinch(_ x: Self) -> Self {
        return x == Self.zero ? Self.one : Self.sinh(x) / x
    }
}

extension Complex: Exponentiable where RealType: LANumeric {
    public static var two: Complex {
        Complex(2, 0)
    }
    
    public static func raiseTwo(to n: Int) -> Complex<RealType> {
        return Complex(RealType.pow(2, n), 0)
    }
    
    public static func raise(_ a: RealType, to b: RealType) -> RealType {
        return RealType.pow(a, b)
    }
    
    public static func greater(of a: RealType, and b: RealType) -> RealType {
        return RealType.maximum(a, b)
    }
}

extension Double: Exponentiable {
    public static var one: Double {
        return 1
    }
    
    public static var two: Double {
        return 2
    }
    
    public static func raiseTwo(to n: Int) -> Double {
        return pow(2, n)
    }
    
    public static func raise(_ a: Double, to b: Double) -> Double {
        return pow(a, b)
    }
    
    public static func greater(of a: Double, and b: Double) -> Double {
        return maximum(a, b)
    }
}

extension Float: Exponentiable {
    public static var one: Float {
        return 1
    }
    
    public static var two: Float {
        return 2
    }
    
    public static func raiseTwo(to n: Int) -> Float {
        return pow(2, n)
    }
    
    public static func raise(_ a: Float, to b: Float) -> Float {
        return pow(a, b)
    }
    
    public static func greater(of a: Float, and b: Float) -> Float {
        return maximum(a, b)
    }
}
