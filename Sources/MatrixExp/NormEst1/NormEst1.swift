//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 11/3/20.
//

import Foundation
import Numerics
import LANumerics

class NormEst1 {
    // MARK: - Properties, Inputs
    let A: Matrix<Double>
    let n: Int
    let isReal: Bool
    var t = 2
    var X: Matrix<Double>
    
    // MARK: - Properties
    //var prnt: Bool
    
    // MARK: - Properties, Results
    var estimate: Double
    var V = Vector<Double>()
    var W = Vector<Double>()
    var numberOfIterations: Int
    var numberOfProducts: Int
    
    init(F: (String, Matrix<Double>) -> Matrix<Double>, t: Int = 2) {
        self.A = Matrix<Double>()
        
        let n = Int(NormEst1.normapp(F: F, flag: NormEst1Flag.dim, X: Matrix<Double>())[0])
        self.n = n
        
        self.isReal = NormEst1.normapp(F: F, flag: NormEst1Flag.real, X: Matrix<Double>())[0] == 0 ? false : true
        
        let prnt = (t < 0)
        self.t = abs(t)
        
        var rpt_S = 0
        var rpt_e = 0
        
        var Y: Matrix<Double>
        if (t == n || n <= 4) {
            X = Matrix<Double>.eye(n)
            Y = NormEst1.normapp(F: F, flag: .notransp, X: X)
            
            let absY = Y.map { abs($0) }
            
            var vals = Vector<Double>()
            for k in 0..<absY.columns {
                var summation = 0.0
                for l in 0..<absY.rows {
                    summation += absY[l, k]
                }
                vals.append(summation)
            }
            
            let sortedVals = vals.enumerated().sorted(by: {$0.element > $1.element})
            
            vals = sortedVals.map { $0.element }
            let m = sortedVals.map { $0.offset }
            
            self.estimate = vals[0]
            
            V = Vector<Double>(repeating: 0.0, count: n)
            V[m[n]] = 1.0
            W = Y[0..<n, (m[n]-1)..<m[n]].vector
            
            numberOfIterations = 0
            numberOfProducts = 1
            
            return
        }
        
        var Xtemp = Matrix<Double>(repeating : 1.0, rows: n, columns: abs(t))
        
        let B1 = 2.0 * NormEst1.randomMatrix(rows: n, columns: 1)
        let B2 = Matrix<Double>(repeating: 1.0, rows: n, columns: 1)
        Xtemp[0..<n, 1..<t] = NormEst1.mysign(A: B1 - B2)
        self.X = Xtemp.map { $0 / Double(n) }
        
        print("self.n = \(self.n)")
        print("self.t = \(self.t)")
        print("self.X = \(self.X)")
        
        let itmax = 5
        
        var it = 0
        var nmv = 0
        
        var ind = [Int](repeating: 0, count: t)
        var ind_hist = [Int](repeating: 0, count: t)
        var S = Matrix<Double>(rows: n, columns: t)
        
        var est_old = 0.0
        var est_j = 0
        var info = 0
        
        var estimate = 0.0
        
        while (true) {
            it += 1
            
            Y = NormEst1.normapp(F: F, flag: .notransp, X: self.X)
            nmv += 1
            
            let absY = Y.map { abs($0) }
            
            var vals = Vector<Double>()
            for k in 0..<absY.columns {
                var summation = 0.0
                for l in 0..<absY.rows {
                    summation += absY[l, k]
                }
                vals.append(summation)
            }
            let sortedVals = vals.enumerated().sorted(by: {$0.element > $1.element})
            vals = sortedVals.map { $0.element }
            let m = sortedVals.map { $0.offset }
            
            var vals_ind = [Int]()
            for k in 0..<t {
                vals_ind.append(ind[m[k]])
            }
            
            estimate = vals[0]
            
            if (estimate > est_old || it == 2) {
                est_j = vals_ind[0]
                W = Y[0..<n, (m[0]-1)..<m[0]].vector
            }
            
            if (prnt) {
                print("\(it): ")
                for k in 0..<t {
                    print("\(vals_ind[k]), \(vals[k])")
                }
            }
            
            if (it >= 2 && estimate <= est_old) {
                estimate = est_old
                info = 2
                break
            }
            est_old = estimate
            
            if (it > itmax) {
                it = itmax
                info = 1
                break
            }
            
            let oldS = S
            S = NormEst1.mysign(A: Y)
            
            // For Double
            let SS = oldS.transpose * S
            let absSS = SS.map { abs($0) }
            var maxVector = Vector<Double>()
            for k in 0..<t {
                let maxValue = absSS[0..<t, (k-1)..<k].reduce { max($0, $1) }
                maxVector.append(maxValue == Double(n) ? 1 : 0)
            }
            let np = maxVector.reduce(0, +)
            
            if (np == Double(t)) {
                info = 3
                break
            }
            
            let r: Int
            (S, r) = NormEst1.undupli(S: S, oldS: oldS, prnt: prnt)
            rpt_S = rpt_S + r
            
            let Z = NormEst1.normapp(F: F, flag: .transp, X: S)
            nmv += 1
            
            let absZ = Z.map { abs($0) }
            var Zvals = Vector<Double>()
            for k in 0..<t {
                let maxValue = absZ[(k-1)..<k, 0..<t].reduce(-1.0 * Double.greatestFiniteMagnitude, {x, y in max(x,y)})
                Zvals.append(maxValue)
            }
            
            if (it >= 2) {
                let maxZvals = Zvals.reduce(Double.greatestFiniteMagnitude, {x, y in max(x,y)})
                if (maxZvals == Zvals[est_j]) {
                    info = 4
                    break
                }
            }
            
            let sortedZVals = Zvals.enumerated().sorted(by: {$0.element > $1.element})
            let m2 = sortedZVals.map { $0.offset }
            
            var imax = t
            
            if (it == 1) {
                ind = m2
                ind_hist = ind
            } else {
                let ismember = m2.map { ind_hist.contains($0) ? 1 : 0 }
                let rep = ismember.reduce(0, +)
                
                rpt_e += rep
                
                if (rep > 0 && prnt) {
                    print("rep e_j = \(rep)")
                }
                
                if (rep == t) {
                    info = 5
                    break
                }
                
                var j = 1
                for i in 1...t {
                    if (j > n) {
                        imax = i - 1
                        break
                    }
                    while (ind_hist.contains(m2[j-1])) {
                        j += 1
                        if (j > n) {
                            imax = i - 1
                            break
                        }
                    }
                    if (j > n) {
                        break
                    }
                    ind[i-1] = m2[j-1]
                    j += 1
                }
                ind_hist.append(contentsOf: ind[0..<imax])
            }
            
            self.X = Matrix<Double>(repeating: 0.0, rows: n, columns: t)
            for j in 0..<imax {
                X[ind[j], j] = 1
            }
        }
        
        if (prnt) {
            switch (info) {
            case 1:
                print("MATLAB:normest1:TerminateIterationLimitReachedn")
            case 2:
                print("MATLAB:normest1:TerminateEstimateNotIncreased")
            case 3:
                print("MATLAB:normest1:TerminateRepeatedSignMatrix")
            case 4:
                print("MATLAB:normest1:TerminatePowerMethodConvergenceTest")
            case 5:
                print("MATLAB:normest1:TerminateRepeatedUnitVectors")
            default:
                print("MATLAB:normest1")
            }
        }
        
        self.estimate = estimate
        numberOfIterations = it
        numberOfProducts = nmv
        
        V = Vector<Double>(repeating: 0.0, count: n)
        V[est_j] = 1
        
        print("MATLAB:normest1:RepeatedUnitVectors \(rpt_S)")
        
    }
    
    init(A: Matrix<Double>, t: Int = 2) {
        /*
        guard A.rows == A.columns else {
            print("Cannot initialize. The matrix is not square.")
            return
        }
        */
        
        self.A = A
        self.n = A.rows
        self.isReal = true
        
        /*
        guard abs(t) >= 1 && abs(t) <= max(A.rows, 2) else {
            print("Cannot initialize. The number of columns in the iteration matrix is out of range.")
            return
        }
        */
        
        let prnt = (t < 0)
        self.t = abs(t)
        
        var rpt_S = 0
        var rpt_e = 0
        
        var Y: Matrix<Double>
        if (t == A.rows || A.rows <= 4) {
            Y = A
            
            let absY = Y.map { abs($0) }
            
            var vals = Vector<Double>()
            for k in 0..<absY.columns {
                var summation = 0.0
                for l in 0..<absY.rows {
                    summation += absY[l, k]
                }
                vals.append(summation)
            }
            
            let sortedVals = vals.enumerated().sorted(by: {$0.element > $1.element})
            
            vals = sortedVals.map { $0.element }
            let m = sortedVals.map { $0.offset }
            
            self.estimate = vals[0]
            
            V = Vector<Double>(repeating: 0.0, count: A.rows)
            V[m[A.rows]] = 1.0
            W = Y[0..<A.rows, (m[A.rows]-1)..<m[A.rows]].vector
            
            numberOfIterations = 0
            numberOfProducts = 1
            
            self.X = Matrix<Double>()
            return
        }
        
        var Xtemp = Matrix<Double>(repeating : 1.0, rows: A.rows, columns: abs(t))
        
        let B1 = 2.0 * NormEst1.randomMatrix(rows: A.rows, columns: 1)
        let B2 = Matrix<Double>(repeating: 1.0, rows: n, columns: 1)
        Xtemp[0..<n, 1..<t] = NormEst1.mysign(A: B1 - B2)
        self.X = Xtemp.map { $0 / Double(A.rows) }
        
        print("self.n = \(self.n)")
        print("self.t = \(self.t)")
        print("self.X = \(self.X)")
        
        print("mysign(A) = \(NormEst1.mysign(A: A))")
        
        let itmax = 5
        
        var it = 0
        var nmv = 0
        
        var ind = [Int](repeating: 0, count: t)
        var ind_hist = [Int](repeating: 0, count: t)
        var S = Matrix<Double>(rows: n, columns: t)
        
        var est_old = 0.0
        var est_j = 0
        var info = 0
        
        var estimate = 0.0
        
        while (true) {
            it += 1
            
            Y = self.A * self.X
            //var Y = normapp(self.A, 'notransp', self.X)
            nmv += 1
            
            let absY = Y.map { abs($0) }
            
            var vals = Vector<Double>()
            for k in 0..<absY.columns {
                var summation = 0.0
                for l in 0..<absY.rows {
                    summation += absY[l, k]
                }
                vals.append(summation)
            }
            let sortedVals = vals.enumerated().sorted(by: {$0.element > $1.element})
            vals = sortedVals.map { $0.element }
            let m = sortedVals.map { $0.offset }
            
            var vals_ind = [Int]()
            for k in 0..<t {
                vals_ind.append(ind[m[k]])
            }
            
            estimate = vals[0]
            
            if (estimate > est_old || it == 2) {
                est_j = vals_ind[0]
                W = Y[0..<n, (m[0]-1)..<m[0]].vector
            }
            
            if (prnt) {
                print("\(it): ")
                for k in 0..<t {
                    print("\(vals_ind[k]), \(vals[k])")
                }
            }
            
            if (it >= 2 && estimate <= est_old) {
                estimate = est_old
                info = 2
                break
            }
            est_old = estimate
            
            if (it > itmax) {
                it = itmax
                info = 1
                break
            }
            
            let oldS = S
            S = NormEst1.mysign(A: Y)
            
            // For Double
            let SS = oldS.transpose * S
            let absSS = SS.map { abs($0) }
            var maxVector = Vector<Double>()
            for k in 0..<t {
                let maxValue = absSS[0..<t, (k-1)..<k].reduce { max($0, $1) }
                maxVector.append(maxValue == Double(n) ? 1 : 0)
            }
            let np = maxVector.reduce(0, +)
            
            if (np == Double(t)) {
                info = 3
                break
            }
            
            let r: Int
            (S, r) = NormEst1.undupli(S: S, oldS: oldS, prnt: prnt)
            rpt_S = rpt_S + r
            
            let Z = A.transpose * S
            nmv += 1
            
            let absZ = Z.map { abs($0) }
            var Zvals = Vector<Double>()
            for k in 0..<t {
                let maxValue = absZ[(k-1)..<k, 0..<t].reduce(-1.0 * Double.greatestFiniteMagnitude, {x, y in max(x,y)})
                Zvals.append(maxValue)
            }
            
            if (it >= 2) {
                let maxZvals = Zvals.reduce(Double.greatestFiniteMagnitude, {x, y in max(x,y)})
                if (maxZvals == Zvals[est_j]) {
                    info = 4
                    break
                }
            }
            
            let sortedZVals = Zvals.enumerated().sorted(by: {$0.element > $1.element})
            let m2 = sortedZVals.map { $0.offset }
            
            var imax = t
            
            if (it == 1) {
                ind = m2
                ind_hist = ind
            } else {
                let ismember = m2.map { ind_hist.contains($0) ? 1 : 0 }
                let rep = ismember.reduce(0, +)
                
                rpt_e += rep
                
                if (rep > 0 && prnt) {
                    print("rep e_j = \(rep)")
                }
                
                if (rep == t) {
                    info = 5
                    break
                }
                
                var j = 1
                for i in 1...t {
                    if (j > n) {
                        imax = i - 1
                        break
                    }
                    while (ind_hist.contains(m2[j-1])) {
                        j += 1
                        if (j > n) {
                            imax = i - 1
                            break
                        }
                    }
                    if (j > n) {
                        break
                    }
                    ind[i-1] = m2[j-1]
                    j += 1
                }
                ind_hist.append(contentsOf: ind[0..<imax])
            }
            
            self.X = Matrix<Double>(repeating: 0.0, rows: n, columns: t)
            for j in 0..<imax {
                X[ind[j], j] = 1
            }
        }
        
        if (prnt) {
            switch (info) {
            case 1:
                print("MATLAB:normest1:TerminateIterationLimitReachedn")
            case 2:
                print("MATLAB:normest1:TerminateEstimateNotIncreased")
            case 3:
                print("MATLAB:normest1:TerminateRepeatedSignMatrix")
            case 4:
                print("MATLAB:normest1:TerminatePowerMethodConvergenceTest")
            case 5:
                print("MATLAB:normest1:TerminateRepeatedUnitVectors")
            default:
                print("MATLAB:normest1")
            }
        }
        
        self.estimate = estimate
        numberOfIterations = it
        numberOfProducts = nmv
        
        V = Vector<Double>(repeating: 0.0, count: n)
        V[est_j] = 1
        
        print("MATLAB:normest1:RepeatedUnitVectors \(rpt_S)")
        print("MATLAB:normest1:RepeatedUnitVectors \(rpt_e)")
    }
    
    static func mysign(A: Matrix<Double>) -> Matrix<Double> {
        let S = A.map { $0 < 0.0 ? -1.0 : 1.0}
        
        return S
    }
    
    static func undupli(S: Matrix<Double>, oldS: Matrix<Double>, prnt: Bool) -> (Matrix<Double>, Int) {
        var Sout = S
        let n = S.rows
        let t = S.columns
        
        var r = 0
        
        if (t == 1) {
            return (S, r)
        }
        
        var W = Matrix<Double>(rows: n, columns: n)
        var jstart: Int
        var last_col: Int
        if (oldS.count == 0) {
            W[0..<n, 0..<1] = S[0..<n, 0..<1]
            jstart = 2
            last_col = 1
        } else {
            W[0..<n, 0..<t] = oldS
            jstart = 1
            last_col = t
        }
        
        for j in jstart...t {
            var rpt = 0
            while (isMaxN(S: Sout[0..<n, (j-1)..<j], W: W, n: n, last_col: last_col)) {
                rpt += 1
                
                let A = 2.0 * NormEst1.randomMatrix(rows: n, columns: 1)
                let B = Matrix<Double>(repeating: 1.0, rows: n, columns: 1)
                Sout[0..<n, (j-1)..<j] = NormEst1.mysign(A: A-B)
                if (rpt > Int(Float(n)/Float(t))) {
                    break
                }
            }
            
            if prnt && rpt > 0 {
                print("Unduplicate rpt = \(rpt)")
            }
            
            r += rpt > 0 ? 1 : 0
            
            if (j < t) {
                last_col += 1
                W[0..<n, (last_col-1)..<last_col] = Sout[0..<n, (j-1)..<j]
            }
        }
        
        return (Sout, r)
    }
    
    static func isMaxN(S: Matrix<Double>, W: Matrix<Double>, n: Int, last_col: Int) -> Bool {
        let temp = S.transpose * W[0..<n, 0..<(last_col-1)]
        let absTemp = temp.map { abs($0) }
        let maxElement = absTemp.reduce(-1.0 * Double.greatestFiniteMagnitude, {max($0, $1)})
        return maxElement == Double(n)
    }
    
    static func randomMatrix(rows: Int, columns: Int) -> Matrix<Double> {
        var elements = Vector<Double>()
        for _ in 0..<(rows * columns) {
            elements.append(drand48())
        }
        return Matrix<Double>(rows: rows, columns: columns, elements: elements)
    }
    
    static func normapp(F: (String, Matrix<Double>) -> Matrix<Double>, flag: NormEst1Flag, X: Matrix<Double>) -> Matrix<Double> {
        let y = F(flag.rawValue, X)
        
        switch (flag) {
        case .dim:
            guard y.count == 1, y[0] == y[0].rounded() else {
                print("Cannot determine the dimension: \(y)")
                return Matrix<Double>()
            }
        case .real:
            guard y.count == 1, y[0] == 0 || y[0] == 1 else {
                print("Cannot determine real or not: \(y)")
                return Matrix<Double>()
            }
        case .notransp, .transp:
            guard y.rows == X.rows && y.columns == X.columns else {
                print("Matrix size mismatch: \(y), \(X)")
                return Matrix<Double>()
            }
        }
        
        return y
    }
}
