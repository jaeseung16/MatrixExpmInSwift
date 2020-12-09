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
    let t: Int
    var X: Matrix<Type>
    let order: Int
    
    // MARK: - Properties
    //var prnt: Bool
    
    // MARK: - Properties, Results
    var estimate: Double
    var V = Vector<Type>()
    var W = Vector<Type>()
    var numberOfIterations: Int
    var numberOfProducts: Int
    var terminationReason: NormEst1TerminationReason
    
    init(A: Matrix<Type>, t: Int = 2, order: Int = 1, toPrint: Bool = false) {
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
        
        self.t = t
        self.order = order
        
        var rpt_S = 0
        var rpt_e = 0
        
        var Y: Matrix<Type>
        if (self.t == A.rows || A.rows <= 4) {
            print("Calling NormEst1.computeExactly \(A)")
            (self.estimate, self.V, self.W) = NormEst1.computeExactly(A, order: order)
            numberOfIterations = 0
            numberOfProducts = 1
            self.X = Matrix<Type>()
            self.terminationReason = .NotApplicable
        } else {
            self.X = NormEst1.initializeX(rows: self.n, columns: self.t, toPrint: toPrint)
            
            let itmax = 5
            
            var it = 0
            var nmv = 0
            
            var ind = [Int](repeating: 0, count: self.t)
            var ind_hist = [Int](repeating: 0, count: self.t)
            var S = Matrix<Type>(rows: n, columns: self.t)
            
            var est_old = 0.0
            var est_j = 0
            var info = NormEst1TerminationReason.NotApplicable
            
            var est = 0.0
            while (true) {
                it += 1
                
                Y = self.X
                for _ in 0..<order {
                    Y = self.A * Y
                }
                
                nmv += 1
                
                let absY = Y.map { Double(floatLiteral: $0.length as! Double) }
                
                let sortedVals = NormEst1.summationAlongRow(absY)
                    .enumerated()
                    .sorted(by: {$0.element > $1.element})
                let vals = sortedVals.map { $0.element }
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
                
                if (toPrint) {
                    print("\(it): ")
                    for k in 0..<self.t {
                        print("\(vals_ind[k]), \(vals[k])")
                    }
                }
                
                if (it >= 2 && est <= est_old) {
                    est = est_old
                    info = .EstimateNotIncreased
                    break
                }
                est_old = est
                
                if (it > itmax) {
                    it = itmax
                    info = .IterationLimitReachedn
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
                        info = .RepeatedSignMatrix
                        break
                    }

                    let r: Int
                    (S, r) = NormEst1.undupli(S: S, oldS: oldS, toPrint: toPrint)
                    rpt_S = rpt_S + r
                }

                //
                var Z = S
                for _ in 0..<order {
                    Z = self.A.adjoint * S
                }
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
                        info = .PowerMethodConvergenceTest
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
                    
                    if (rep > 0 && toPrint) {
                        print("rep e_j = \(rep)")
                    }
                    
                    if (rep == self.t) {
                        info = .RepeatedUnitVectors
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
            
            if (toPrint) {
                print("Termination Reason: \(info)")
            }
            
            self.terminationReason = info
            self.estimate = est
            numberOfIterations = it
            numberOfProducts = nmv
            
            V = Vector<Type>(repeating: 0.0, count: n)
            V[est_j] = 1
            
            print("ParallelCol \(rpt_S)")
            print("RepeatedUnitVectors \(rpt_e)")
        }
    }
    
    static private func computeExactly(_ A: Matrix<Type>, order: Int = 1) -> (Double, Vector<Type>, Vector<Type>) {
        var B = A
        
        for _ in 0..<(order-1) {
            B = A * B
        }
        
        let absB = B.map { Double(floatLiteral: $0.length as! Double) }
        
        let vals = NormEst1.summationAlongRow(absB)
            .enumerated()
            .sorted(by: {$0.element > $1.element})
        
        let elements = vals.map { $0.element }
        let offests = vals.map { $0.offset }
        
        var V = Vector<Type>(repeating: 0.0, count: A.rows)
        V[offests[A.rows-1]] = 1.0
        
        let W = A[0..<A.rows, offests[A.rows-1]..<(offests[A.rows-1]+1)].vector
        
        return (elements[0], V, W)
    }
    
    static private func initializeX(rows: Int, columns: Int, toPrint: Bool) -> Matrix<Type> {
        var Xtemp = Matrix<Type>(repeating : 1.0, rows: rows, columns: columns)
        
        let B1 = 2.0 * NormEst1.randomMatrix(rows: rows, columns: columns - 1)
        let B2 = Matrix<Type>(repeating: 1.0, rows: rows, columns: columns - 1)
        Xtemp[0..<rows, 1..<columns] = NormEst1.mysign(A: B1 - B2)
        
        (Xtemp, _) = NormEst1.undupli(S: Xtemp, oldS: Matrix<Type>(), toPrint: toPrint)
        
        return Xtemp.map { $0 / Type(floatLiteral: Double(rows) as! Type.FloatLiteralType) }
    }
    
    static func mysign(A: Matrix<Type>) -> Matrix<Type> {
        let S = A.map { $0 == 0.0 ? 1.0 : ($0 / Type(floatLiteral: $0.length as! Type.FloatLiteralType))}
        return S
    }
    
    // Look for and replace columns of S parallel to other columns of S or to columns of oldS
    static func undupli(S: Matrix<Type>, oldS: Matrix<Type>, toPrint: Bool) -> (Matrix<Type>, Int) {
        if (S.columns == 1) {
            return (S, 0)
        }
        
        let nRows = S.rows
        let nColumns = S.columns
        var Sout = S
        var r = 0
        
        var startColumn: Int
        var lastColumn: Int
        var W = Matrix<Type>(rows: nRows, columns: 2 * nColumns - 1)
        if oldS.count == 0 {
            W[0..<nRows, 0..<1] = S[0..<nRows, 0..<1]
            startColumn = 1
            lastColumn = 0
        } else {
            W[0..<nRows, 0..<nColumns] = oldS
            startColumn = 0
            lastColumn = nColumns-1
        }
        
        for col in startColumn..<nColumns {
            var rpt = 0
            while (isMaxN(S: Sout[0..<nRows, col..<(col+1)], W: W, n: nRows, lastColumn: lastColumn)) {
                rpt += 1
                
                let A = 2.0 * NormEst1.randomMatrix(rows: nRows, columns: 1)
                let B = Matrix<Type>(repeating: 1.0, rows: nRows, columns: 1)
                Sout[0..<nRows, col..<(col+1)] = NormEst1.mysign(A: A-B)
                if (rpt > Int(Float(nRows)/Float(nColumns))) {
                    break
                }
            }
            
            if toPrint && rpt > 0 {
                print("Unduplicate rpt = \(rpt)")
            }
            
            r += rpt > 0 ? 1 : 0
            
            if (col < nColumns-1) {
                lastColumn += 1
                W[0..<nRows, lastColumn..<(lastColumn+1)] = Sout[0..<nRows, col..<(col+1)]
            }
        }
        return (Sout, r)
    }
    
    static func summationAlongRow(_ matrix: Matrix<Double>) -> Vector<Double> {
        var vals = Vector<Double>()
        for k in 0..<matrix.columns {
            var summation = 0.0
            for l in 0..<matrix.rows {
                summation += matrix[l, k]
            }
            vals.append(summation)
        }
        return vals
    }
    
    static func isMaxN(S: Matrix<Type>, W: Matrix<Type>, n: Int, lastColumn: Int) -> Bool {
        let temp = S.transpose * W[0..<n, 0..<lastColumn]
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
