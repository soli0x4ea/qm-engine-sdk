// swift-tools-version: 6.0
import PackageDescription

let package = Package(
    name: "EngineKit",
    defaultLocalization: "zh-CN",
    platforms: [.iOS(.v17), .macOS(.v14)],
    products: [
        .library(name: "EngineKit", targets: ["EngineKit"]),
        .library(name: "QMModules", targets: ["QMModules"]),
    ],
    dependencies: [
        // W6：swift-numerics（Apple 官方，Apache-2.0）——QFT 模块复数通道
        // （8×8 酉矩阵构造 + 酉性断言），计划 §W6 唯一新增三方依赖。
        .package(url: "https://github.com/apple/swift-numerics", from: "1.1.1"),
    ],
    targets: [
        .target(
            name: "EngineKit",
            resources: [.process("Resources")]
        ),
        // W3：真值模块层（SOP 移植目标）。
        // W3 包体裁决：移除 mlx-swift（26 MB 包体 + 无 Metal 模拟器必崩），
        // 数值线代改由 EngineKit 直调 Accelerate LAPACK（SymEigh），
        // 模块保持 swift test 可直测 fixtures。MLX 版见 git 历史。
        .target(
            name: "QMModules",
            dependencies: [
                "EngineKit",
                .product(name: "Numerics", package: "swift-numerics"),
            ]
        ),
        .testTarget(
            name: "EngineKitTests",
            dependencies: ["EngineKit"],
            resources: [.copy("Fixtures")]
        ),
        .testTarget(
            name: "QMModulesTests",
            dependencies: ["QMModules"],
            resources: [.copy("Fixtures")]
        ),
    ]
)
