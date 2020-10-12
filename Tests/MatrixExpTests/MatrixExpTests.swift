import XCTest
@testable import MatrixExp
@testable import Numerics
@testable import LANumerics

final class MatrixExpTests: XCTestCase {
    func testExample() {
        // This is an example of a functional test case.
        // Use XCTAssert and related functions to verify your tests produce the correct
        // results.
        XCTAssertEqual(MatrixExp().text, "Hello, World!")
    }
    
    func testIsDiag() {
        let diagMatrix = Matrix<Double>(diagonal: [1.0, 0.5, 0.0, -1.0])
        XCTAssertTrue(MatrixExp.isDiag(diagMatrix))
    }
    
    func testDiag() {
        let diag = Vector<Double>([1.0, 0.5, 0.0, -1.0])
        let diagMatrix = Matrix<Double>(diagonal: diag)
        let result = MatrixExp.evaluate(for: diagMatrix)!
        let expected = Matrix<Double>(diagonal: diag.map { exp($0) })
        
        XCTAssertEqual(result, expected)
    }
    
    func testPadeCoefficients() {
        let order = 5
        let coeff: Vector<Double> = [30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0]
        let expected = MatrixExpConst.padeCoefficients(for: order)

        XCTAssertEqual(coeff, expected)
    }
    
    func testTheta() {
        let order = 7
        let theta = 9.504178996162932e-001
        let expected = MatrixExpConst.theta(for: order)

        XCTAssertEqual(theta, expected)
    }
    
    func testScalingAndOrder() {
        let M = Matrix<Double>(rows: [[0, 1], [1, 0]])
        let order = 13
        let scaling = -2
        
        let (s, o, _) = MatrixExp.expmParams(for: M)

        XCTAssertEqual(scaling, s)
        XCTAssertEqual(order, o)
    }
    
    func testEvaluate() {
        let M = Matrix<Double>(rows: [[0, 1], [1, 0]])
        let expected = Matrix<Double>(rows: [[1.5430806348152437, 1.1752011936438014], [1.1752011936438014, 1.5430806348152437]])
        let result = MatrixExp.evaluate(for: M)
        
        XCTAssertEqual(expected, result)
    }
    
    func testEvaluateComplex() {
        let M = Matrix<Complex<Double>>(rows: [[Complex<Double>(0,0), Complex<Double>(0, -1)],
                                               [Complex<Double>(0, 1), Complex<Double>(0, 0)]])
        let expected = Matrix<Complex<Double>>(rows: [[Complex<Double>(1.5430806348152437, 0), Complex<Double>(0, -1.1752011936438016)],
                                                      [Complex<Double>(0, 1.1752011936438014), Complex<Double>(1.5430806348152437, 0)]])
        let result = MatrixExpComplex.evaluate(for: M)
        
        XCTAssertEqual(expected, result)
    }

    static var allTests = [
        ("testExample", testExample),
        ("testIsDiag", testIsDiag),
        ("testDiag", testDiag),
        ("testPadeCoefficients", testPadeCoefficients),
        ("testScalingAndOrder", testScalingAndOrder),
        ("testEvaluate", testEvaluate),
        ("testEvaluateComplex", testEvaluateComplex)
    ]
}
