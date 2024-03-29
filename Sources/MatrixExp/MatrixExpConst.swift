//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 10/11/20.
//

import Numerics
import LANumerics

public enum PadeApproximantOrder: Int {
    case three = 3
    case five = 5
    case seven = 7
    case nine = 9
    case thirteen = 13
}

struct MatrixExpConst<T> where T: Exponentiable, T.Magnitude: Real {
    static func padeCoefficients(_ order: PadeApproximantOrder) -> Vector<T> {
        var padeCoefficients: Vector<T>
        
        switch (order) {
        case .three:
            padeCoefficients = [120.0, 60.0, 12.0, 1.0]
        case .five:
            padeCoefficients = [30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0]
        case .seven:
            padeCoefficients = [17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0]
        case .nine:
            padeCoefficients = [17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0,
                     2162160.0, 110880.0, 3960.0, 90.0, 1.0]
        case .thirteen:
            padeCoefficients = [64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
                     1187353796428800.0,  129060195264000.0,   10559470521600.0,
                     670442572800.0,      33522128640.0,       1323241920.0,
                     40840800.0,          960960.0,            16380.0,
                     182.0,  1.0]
        }
        
        return padeCoefficients
    }
    
    static var coefficientsOfBackwardsErrorFunction : Vector<T.Magnitude> {
        return [1.0/100800.0,
                1.0/10059033600.0,
                1.0/4487938430976000.0,
                1.0/5914384781877411840000.0,
                1.0/113250775606021113483283660800000000.0]
    }
    
    static func theta(for order: PadeApproximantOrder) -> T.Magnitude? {
        var theta: T.Magnitude?
        
        switch (order) {
        case .three:
            theta = 1.495585217958292e-002
        case .five:
            theta = 2.539398330063230e-001
        case .seven:
            theta = 9.504178996162932e-001
        case .nine:
            theta = 2.097847961257068e+000
        case .thirteen:
            theta = 5.371920351148152e+000
        }
        
        return theta
    }
    
}
