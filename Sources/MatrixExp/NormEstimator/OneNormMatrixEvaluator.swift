//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 11/18/23.
//

import Foundation
import Numerics
import LANumerics

protocol OneNormMatrixEvaluator {
    associatedtype T: Exponentiable
    
    var dimension: Int { get }
    var isReal: Bool { get }
    
    func multiply(by: Matrix<T>) -> Matrix<T>
    func conjugateTransposeAndMultiply(by: Matrix<T>) -> Matrix<T>
}
