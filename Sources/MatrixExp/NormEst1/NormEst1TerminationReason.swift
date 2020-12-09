//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 12/8/20.
//

import Foundation

enum NormEst1TerminationReason: String {
    case IterationLimitReachedn
    case EstimateNotIncreased
    case RepeatedSignMatrix
    case PowerMethodConvergenceTest
    case RepeatedUnitVectors
    case NotApplicable
}
