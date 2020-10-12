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

    static var allTests = [
        ("testExample", testExample),
        ("testIsDiag", testIsDiag),
        ("testDiag", testDiag),
        ("testPadeCoefficients", testPadeCoefficients)
    ]
}
