// swift-tools-version:5.3
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "MatrixExp",
    products: [
        // Products define the executables and libraries a package produces, and make them visible to other packages.
        .library(
            name: "MatrixExp",
            targets: ["MatrixExp"]),
    ],
    dependencies: [
        // Dependencies declare other packages that this package depends on.
        // .package(url: /* package url */, from: "1.0.0"),
        .package(url: "https://github.com/apple/swift-numerics", from: "0.1.0"),
        .package(url: "https://github.com/phlegmaticprogrammer/LANumerics", from: "0.1.12")
    ],
    targets: [
        // Targets are the basic building blocks of a package. A target can define a module or a test suite.
        // Targets can depend on other targets in this package, and on products in packages this package depends on.
        .target(
            name: "MatrixExp",
            dependencies: [
                .product(name: "Numerics", package: "swift-numerics"),
                .product(name: "LANumerics", package: "LANumerics")
            ]),
        .testTarget(
            name: "MatrixExpTests",
            dependencies: [
                "MatrixExp",
                .product(name: "Numerics", package: "swift-numerics"),
                .product(name: "LANumerics", package: "LANumerics")
            ]),
    ]
)
