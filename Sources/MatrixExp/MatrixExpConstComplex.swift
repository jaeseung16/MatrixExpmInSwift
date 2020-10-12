//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 10/12/20.
//

import Foundation
import LANumerics
import ComplexModule

struct MatrixExpConstComplex {
    
    enum PadeApproximantOrder: Int {
        case three = 3
        case five = 5
        case seven = 7
        case nine = 9
        case thirteen = 13
    }
    
    static func padeCoefficients(for order: Int) -> Vector<Complex<Double>>? {
        guard let m = PadeApproximantOrder(rawValue: order) else {
            return nil
        }
        
        var coeff: Vector<Complex<Double>>
        
        switch (m) {
        case .three:
            coeff = [Complex<Double>(120.0, 0.0), Complex<Double>(60.0, 0.0), Complex<Double>(12.0, 0.0), Complex<Double>(1.0, 0.0)]
        case .five:
            coeff = [Complex<Double>(30240.0, 0.0), Complex<Double>(15120.0, 0.0), Complex<Double>(3360.0, 0.0), Complex<Double>(420.0, 0.0), Complex<Double>(30.0, 0.0), Complex<Double>(1.0, 0.0)]
        case .seven:
            coeff = [Complex<Double>(17297280.0, 0.0), Complex<Double>(8648640.0, 0.0), Complex<Double>(1995840.0, 0.0), Complex<Double>(277200.0, 0.0), Complex<Double>(25200.0, 0.0), Complex<Double>(1512.0, 0.0), Complex<Double>(56.0, 0.0), Complex<Double>(1.0, 0.0)]
        case .nine:
            coeff = [Complex<Double>(17643225600.0, 0.0), Complex<Double>(8821612800.0, 0.0), Complex<Double>(2075673600.0, 0.0), Complex<Double>(302702400.0, 0.0), Complex<Double>(30270240.0, 0.0),
                        Complex<Double>(2162160.0, 0.0), Complex<Double>(110880.0, 0.0), Complex<Double>(3960.0, 0.0), Complex<Double>(90.0, 0.0), Complex<Double>(1.0, 0.0)]
        case .thirteen:
            coeff = [Complex<Double>(64764752532480000.0, 0.0), Complex<Double>(32382376266240000.0, 0.0), Complex<Double>(7771770303897600.0, 0.0),
                        Complex<Double>(1187353796428800.0, 0.0),  Complex<Double>(129060195264000.0, 0.0),   Complex<Double>(10559470521600.0, 0.0),
                        Complex<Double>(670442572800.0, 0.0),      Complex<Double>(33522128640.0, 0.0),       Complex<Double>(1323241920.0, 0.0),
                        Complex<Double>(40840800.0, 0.0),          Complex<Double>(960960.0, 0.0),            Complex<Double>(16380.0, 0.0),
                        Complex<Double>(182.0, 0.0),  Complex<Double>(1.0, 0.0)]
        }
        
        return coeff
    }
    
    static let coefficientsOfBackwardsErrorFunction : Vector<Double>
        = [1.0/100800.0, 1.0/10059033600.0, 1.0/4487938430976000.0,
           1.0/5914384781877411840000.0, 1.0/113250775606021113483283660800000000.0]
    
    static func theta(for order: Int) -> Double? {
        guard let m = PadeApproximantOrder(rawValue: order) else {
            return nil
        }
        
        var theta: Double
        
        switch (m) {
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
