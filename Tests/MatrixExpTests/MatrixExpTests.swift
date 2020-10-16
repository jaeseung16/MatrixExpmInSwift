import XCTest
@testable import MatrixExp
@testable import Numerics
@testable import LANumerics

final class MatrixExpTests: XCTestCase {
    func testExample() {
        // This is an example of a functional test case.
        // Use XCTAssert and related functions to verify your tests produce the correct
        // results.
        XCTAssertEqual(MatrixExp<Double>().text, "Hello, World!")
    }
    
    func testIsDiag() {
        let diagMatrix = Matrix<Double>(diagonal: [1.0, 0.5, 0.0, -1.0])
        XCTAssertTrue(MatrixExp.isDiag(diagMatrix))
    }
    
    func testDiag() {
        let diag = Vector<Double>([1.0, 0.5, 0.0, -1.0])
        let diagMatrix = Matrix<Double>(diagonal: diag)
        let result = MatrixExp<Double>.evaluate(for: diagMatrix)!
        let expected = Matrix<Double>(diagonal: diag.map { exp($0) })
        
        for k in 0..<expected.elements.count {
            XCTAssertEqual(expected[k], result[k], accuracy: Double.ulpOfOne)
        }
    }
    
    func testPadeCoefficients() {
        let order = 5
        let coeff: Vector<Double> = [30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0]
        let expected = MatrixExpConst<Double>.padeCoefficients(for: order)

        XCTAssertEqual(coeff, expected)
    }
    
    func testTheta() {
        let order = 7
        let theta = 9.504178996162932e-001
        let expected = MatrixExpConst<Double>.theta(for: order)

        XCTAssertEqual(theta, expected)
    }
    
    func testScalingAndOrder() {
        let M = Matrix<Double>(rows: [[0, 1], [1, 0]])
        let order = 9
        let scaling = 0
        
        let (s, o, _) = MatrixExp<Double>.expmParams(for: M)

        XCTAssertEqual(scaling, s)
        XCTAssertEqual(order, o)
    }
    
    func testEvaluate() {
        let M = Matrix<Double>(rows: [[0, 1], [1, 0]])
        let expected = Matrix<Double>(rows: [[1.5430806348152437, 1.1752011936438014], [1.1752011936438014, 1.5430806348152437]])
        let result = MatrixExp<Double>.evaluate(for: M)!
        
        for k in 0..<expected.elements.count {
            XCTAssertEqual(expected[k], result[k], accuracy: Double.ulpOfOne)
        }
    }
    
    func testEvaluate2() {
        let M = Matrix<Double>(rows: [[0, 1.0 / 1024.0], [1.0 / 2048.0, 0]])
        let expected = Matrix<Double>(rows: [[1.0000002384185886, 0.0009765625776102275], [0.0004882812888051137, 1.0000002384185886]])
        let result = MatrixExp<Double>.evaluate(for: M)!
        
        for k in 0..<expected.elements.count {
            XCTAssertEqual(expected[k], result[k], accuracy: Double.ulpOfOne)
        }
    }
    
    func testEvaluate3() {
        let M = Matrix<Double>(rows: [[0.0, 8.0], [16.0, 0.0]])
        let expected = Matrix<Double>(rows: [[40968.60491465862, 28969.178342277763], [57938.356684555525, 40968.60491465862]])
        let result = MatrixExp<Double>.evaluate(for: M)!
        
        for k in 0..<expected.elements.count {
            XCTAssertEqual(expected[k], result[k], accuracy: Double.ulpOfOne)
        }
    }
    
    func testEvaluate4() {
        let M = Matrix<Double>(rows: [[1.0, 10000.0], [0.0, -1.0]])
        let expected = Matrix<Double>(rows: [[exp(1.0), 5000.0 * (exp(1.0) - exp(-1.0))], [0.0, exp(-1.0)]])
        let result = MatrixExp<Double>.evaluate(for: M)!
        
        for k in 0..<expected.elements.count {
            XCTAssertEqual(expected[k], result[k], accuracy: Double.ulpOfOne)
        }
    }
    
    func testEvaluateComplex() {
        let M = Matrix<Complex<Double>>(rows: [[Complex<Double>(0,0), Complex<Double>(0, -1)],
                                               [Complex<Double>(0, 1), Complex<Double>(0, 0)]])
        let expected = Matrix<Complex<Double>>(rows: [[Complex<Double>(1.5430806348152422, 0), Complex<Double>(0, -1.1752011936438003)],
                                                      [Complex<Double>(0, 1.1752011936438003), Complex<Double>(1.543080634815243, 0)]])
        let result = MatrixExp<Complex>.evaluate(for: M)!
        
        for k in 0..<expected.elements.count {
            XCTAssertEqual(expected[k].real, result[k].real, accuracy: Double.ulpOfOne)
            XCTAssertEqual(expected[k].imaginary, result[k].imaginary, accuracy: Double.ulpOfOne)
        }
    }
    
    func testEll() {
        let M = Matrix<Complex<Double>>(rows: [[Complex<Double>(0,0), Complex<Double>(0, -1)],
                                               [Complex<Double>(0, 1), Complex<Double>(0, 0)]])
        let expected = 10.0
        let result = MatrixExp<Complex>.ell(M, coeff: 1.0/100800.0, order: 2)
        
        XCTAssertEqual(expected, result)
    }

    func testQuasiTrianglularStructure() {
        let M = Matrix<Double>(rows: [[1, -1, 0, 0],
                                      [0, 1, 1, 0],
                                      [0, 0, 1, 1],
                                      [0, 0, 1, 1]])
        
        let expected = [1, 0, 2]
        let result = MatrixExp<Double>.quasiTrianglularStructure(M)
        
        XCTAssertEqual(expected, result)
    }

    static var allTests = [
        ("testExample", testExample),
        ("testIsDiag", testIsDiag),
        ("testDiag", testDiag),
        ("testPadeCoefficients", testPadeCoefficients),
        ("testScalingAndOrder", testScalingAndOrder),
        ("testEvaluate", testEvaluate),
        ("testEvaluate2", testEvaluate2),
        ("testEvaluate3", testEvaluate3),
        ("testEvaluate4", testEvaluate4),
        ("testEvaluateComplex", testEvaluateComplex),
        ("testQuasiTrianglularStructure", testQuasiTrianglularStructure)
    ]
}
