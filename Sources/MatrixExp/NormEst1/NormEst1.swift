//
//  File.swift
//  
//
//  Created by Jae Seung Lee on 11/3/20.
//

import Foundation
import Numerics
import LANumerics

class NormEst1<Type> where Type: Exponentiable, Type.Magnitude: Real {
    // MARK: - Properties, Inputs
    let A: Matrix<Type>
    let n: Int
    var t = 2
    var X: Matrix<Type>
    
    // MARK: - Properties
    //var prnt: Bool
    
    // MARK: - Properties, Results
    var estimate: Double
    var V = Vector<Type>()
    var W = Vector<Type>()
    var numberOfIterations: Int
    var numberOfProducts: Int
    
    init(A: Matrix<Type>, t: Int = 2) {
        /*
        guard A.rows == A.columns else {
            print("Cannot initialize. The matrix is not square.")
            return
        }
        */
        
        self.A = A
        self.n = A.rows
        
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
        
        var Y: Matrix<Type>
        if (self.t == A.rows || A.rows <= 4) {
            Y = A
            
            let absY = Y.map { Double(floatLiteral: $0.length as! Double) }
            
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
            
            V = Vector<Type>(repeating: 0.0, count: A.rows)
            V[m[A.rows-1]] = 1.0
            W = Y[0..<A.rows, m[A.rows-1]..<(m[A.rows-1]+1)].vector
            
            numberOfIterations = 0
            numberOfProducts = 1
            
            self.X = Matrix<Type>()
            return
        }
        
        var Xtemp = Matrix<Type>(repeating : 1.0, rows: A.rows, columns: self.t)
        
        let B1 = 2.0 * NormEst1.randomMatrix(rows: A.rows, columns: self.t - 1)
        let B2 = Matrix<Type>(repeating: 1.0, rows: A.rows, columns: self.t - 1)
        Xtemp[0..<n, 1..<self.t] = NormEst1.mysign(A: B1 - B2)
        (Xtemp, _) = NormEst1.undupli(S: Xtemp, oldS: Matrix<Type>(), prnt: prnt)
        self.X = Xtemp.map { $0 / Type(floatLiteral: Double(A.rows) as! Type.FloatLiteralType) }
        
        let itmax = 5
        
        var it = 0
        var nmv = 0
        
        var ind = [Int](repeating: 0, count: self.t)
        var ind_hist = [Int](repeating: 0, count: self.t)
        var S = Matrix<Type>(rows: n, columns: self.t)
        
        var est_old = 0.0
        var est_j = 0
        var info = 0
        
        var est = 0.0
        while (true) {
            it += 1
            
            Y = self.A * self.X
            nmv += 1
            
            let absY = Y.map { Double(floatLiteral: $0.length as! Double) }
            
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
            for k in 0..<self.t {
                vals_ind.append(ind[m[k]])
            }
   
            est = vals[0]
            
            if (est > est_old || it == 2) {
                est_j = vals_ind[0]
                W = Y[0..<n, m[0]..<(m[0]+1)].vector
            }
            
            if (prnt) {
                print("\(it): ")
                for k in 0..<self.t {
                    print("\(vals_ind[k]), \(vals[k])")
                }
            }
            
            if (it >= 2 && est <= est_old) {
                est = est_old
                info = 2
                break
            }
            est_old = est
            
            if (it > itmax) {
                it = itmax
                info = 1
                break
            }
            
            let oldS = S
            S = NormEst1.mysign(A: Y)
            
            // For Double
            if (Type(floatLiteral: 0.0 as! Type.FloatLiteralType) is Double) {
                let SS = oldS.transpose * S
                let absSS = SS.map { $0.length as! Double }
                var maxVector = Vector<Double>()
                for k in 0..<self.t {
                    let maxValue = absSS[0..<self.t, k..<(k+1)].reduce { max($0, $1) }
                    maxVector.append(maxValue == Double(n) ? 1 : 0)
                }
                let np = maxVector.reduce(0, +)
                
                if (np == Double(self.t)) {
                    info = 3
                    break
                }

                let r: Int
                (S, r) = NormEst1.undupli(S: S, oldS: oldS, prnt: prnt)
                rpt_S = rpt_S + r
            }

            //
            let Z = self.A.adjoint * S
            nmv += 1
        
            let absZ = Z.map { $0.length as! Double }
            var Zvals = Vector<Double>()
            for k in 0..<self.n {
                let maxValue = absZ[k..<(k+1), 0..<self.t].reduce(-1.0 * Double.greatestFiniteMagnitude, {x, y in max(x,y)})
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
            var imax = self.t
            if (it == 1) {
                ind = [Int](m2[0..<self.t])
                ind_hist = ind
            } else {
                let ismember = m2.map { ind_hist.contains($0) ? 1 : 0 }
                let rep = ismember.reduce(0, +)
                
                rpt_e += rep
                
                if (rep > 0 && prnt) {
                    print("rep e_j = \(rep)")
                }
                
                if (rep == self.t) {
                    info = 5
                    break
                }
                
                var j = 1
                for i in 1...self.t {
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
            
            self.X = Matrix<Type>(repeating: 0.0, rows: n, columns: self.t)
            for j in 0..<imax {
                self.X[ind[j], j] = 1
            }
        }
        
        if (true) {
            switch (info) {
            case 1:
                print("TerminateIterationLimitReachedn")
            case 2:
                print("TerminateEstimateNotIncreased")
            case 3:
                print("TerminateRepeatedSignMatrix")
            case 4:
                print("TerminatePowerMethodConvergenceTest")
            case 5:
                print("TerminateRepeatedUnitVectors")
            default:
                print("")
            }
        }
        
        self.estimate = est
        numberOfIterations = it
        numberOfProducts = nmv
        
        V = Vector<Type>(repeating: 0.0, count: n)
        V[est_j] = 1
        
        print("ParallelCol \(rpt_S)")
        print("RepeatedUnitVectors \(rpt_e)")
    }
    
    static func mysign(A: Matrix<Type>) -> Matrix<Type> {
        let S = A.map { $0 == 0.0 ? 1.0 : ($0 / Type(floatLiteral: $0.length as! Type.FloatLiteralType))}
        return S
    }
    
    static func undupli(S: Matrix<Type>, oldS: Matrix<Type>, prnt: Bool) -> (Matrix<Type>, Int) {
        var Sout = S
        let n = S.rows
        let t = S.columns
        
        var r = 0
        
        if (t == 1) {
            return (S, r)
        }
        
        var W = Matrix<Type>(rows: n, columns: 2 * t - 1)
        var jstart: Int
        var last_col: Int
        if (oldS.count == 0) {
            W[0..<n, 0..<1] = S[0..<n, 0..<1]
            jstart = 1
            last_col = 0
        } else {
            W[0..<n, 0..<t] = oldS
            jstart = 0
            last_col = t-1
        }
        
        for j in jstart..<t {
            var rpt = 0
            while (isMaxN(S: Sout[0..<n, j..<(j+1)], W: W, n: n, last_col: last_col)) {
                rpt += 1
                
                let A = 2.0 * NormEst1.randomMatrix(rows: n, columns: 1)
                let B = Matrix<Type>(repeating: 1.0, rows: n, columns: 1)
                Sout[0..<n, j..<(j+1)] = NormEst1.mysign(A: A-B)
                if (rpt > Int(Float(n)/Float(t))) {
                    break
                }
            }
            if prnt && rpt > 0 {
                print("Unduplicate rpt = \(rpt)")
            }
            
            r += rpt > 0 ? 1 : 0
            
            if (j < t-1) {
                last_col += 1
                W[0..<n, last_col..<(last_col+1)] = Sout[0..<n, j..<(j+1)]
            }
        }
        return (Sout, r)
    }
    
    static func isMaxN(S: Matrix<Type>, W: Matrix<Type>, n: Int, last_col: Int) -> Bool {
        let temp = S.transpose * W[0..<n, 0..<last_col]
        let absTemp = temp.map { $0.length as! Double}
        let maxElement = absTemp.reduce(-1.0 * Double.greatestFiniteMagnitude, { max($0, $1) })
        return maxElement == Double(n)
    }
    
    static func randomMatrix(rows: Int, columns: Int) -> Matrix<Type> {
        var elements = Vector<Type>()
        for _ in 0..<(rows * columns) {
            elements.append(Type(floatLiteral: drand48() as! Type.FloatLiteralType))
        }
        return Matrix<Type>(rows: rows, columns: columns, elements: elements)
    }
    
}
