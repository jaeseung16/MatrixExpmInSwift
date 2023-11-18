//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 12/8/20.
//

enum NormEst1TerminationReason: String {
    case iterationLimitReached
    case estimateNotIncreased
    case repeatedSignMatrix
    case powerMethodConvergenceTest
    case repeatedUnitVectors
    case noIterationNormComputedExactly
    case notApplicable
}
