import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// 笔记 11《希尔伯特空间与狄拉克符号_Bloch球》模块测试。
/// fixture 为 3D 线框球（零线条数据）→ 结构级对拍；物理律以不变量网格断言扛主力。
@Suite(.serialized)
struct BlochSphereModuleTests {

    private let module = BlochSphereModule()

    private func constants() throws -> ConstantsSet { try ConstantsSet.load(.v2022) }

    private func values(theta: Double = BlochMath.anchorTheta,
                        phi: Double = BlochMath.anchorPhi) -> ParamValues {
        var v = ParamValues.defaults(for: module.params)
        v.sliders["theta"] = theta
        v.sliders["phi"] = phi
        return v
    }

    // MARK: fixture 结构对拍（3D 图零线条）

    @Test("fixture：结构级对拍（3D 线框球，单图零线条）")
    func fixtureStructure() throws {
        let fixture = try ModuleFixture.load("希尔伯特空间与狄拉克符号_Bloch球__v2022")
        #expect(fixture.document["status"] as? String == "ok")
        let figures = try #require(fixture.document["figures"] as? [[String: Any]])
        #expect(figures.count == 1)
        let axes = try #require(figures[0]["axes"] as? [[String: Any]])
        #expect(axes.count == 1)
        #expect((axes[0]["title"] as? String)?.contains("Bloch sphere") == true)
        let lines = axes[0]["lines"] as? [[String: Any]] ?? []
        #expect(lines.isEmpty, "脚本是 3D 线框球，不导出 lines——曲线对拍由物理律不变量承担")
    }

    // MARK: 物理律 1：单位球 + 双路径一致

    @Test("物理律1：θ/φ 网格上 ‖r‖=1 且密度矩阵 Pauli 迹路径 == 参数化路径")
    func unitSphereAndDualPath() {
        let n = 11
        for i in 0..<n {
            let theta = Double(i) * .pi / Double(n - 1)
            for j in 0..<n {
                let phi = Double(j) * 2.0 * .pi / Double(n - 1)
                let sv = BlochMath.stateVector(theta: theta, phi: phi)
                let rState = BlochMath.blochVectorFromState(alphaRe: sv.alphaRe,
                                                            betaRe: sv.betaRe, betaIm: sv.betaIm)
                let rParam = BlochMath.blochVector(theta: theta, phi: phi)
                let norm = sqrt(rParam.rx * rParam.rx + rParam.ry * rParam.ry + rParam.rz * rParam.rz)
                #expect(abs(norm - 1) < 1e-12,
                        "θ=\(theta) φ=\(phi): |r|=\(norm) 应为单位球面")
                #expect(abs(rState.rx - rParam.rx) < 1e-12
                        && abs(rState.ry - rParam.ry) < 1e-12
                        && abs(rState.rz - rParam.rz) < 1e-12,
                        "θ=\(theta) φ=\(phi): 双路径不一致")
            }
        }
    }

    // MARK: 物理律 2：对跖点正交

    @Test("物理律2：对跖点 (π−θ, φ+π) 重叠 |⟨ψ|ψ⊥⟩|² = (1+r₁·r₂)/2 = 0")
    func antipodalOrthogonality() {
        for i in 1..<9 {
            let theta = Double(i) * .pi / 10
            for j in 0..<8 {
                let phi = Double(j) * .pi / 4
                let r = BlochMath.blochVector(theta: theta, phi: phi)
                let anti = BlochMath.antipode(theta: theta, phi: phi)
                #expect(abs(anti.theta - (.pi - theta)) < 1e-12)
                #expect(abs(anti.phi - (phi + .pi).truncatingRemainder(dividingBy: 2 * .pi)) < 1e-12)
                let rAnti = BlochMath.blochVector(theta: anti.theta, phi: anti.phi)
                #expect(BlochMath.overlapSquared(r, rAnti) < 1e-12,
                        "θ=\(theta) φ=\(phi): 对跖态应严格正交")
                #expect(abs(BlochMath.purity(r) - 1) < 1e-12, "纯态纯度 = 1")
            }
        }
    }

    // MARK: 物理律 3：参数化一致（|α|²+|β|²=1，r_z = |α|²−|β|² = cosθ）

    @Test("物理律3：归一化 + r_z = |α|²−|β|² = cosθ")
    func parameterizationConsistency() {
        for i in 0..<13 {
            let theta = Double(i) * .pi / 12
            let phi = 1.234
            let sv = BlochMath.stateVector(theta: theta, phi: phi)
            let p0 = sv.alphaRe * sv.alphaRe
            let p1 = sv.betaRe * sv.betaRe + sv.betaIm * sv.betaIm
            #expect(abs(p0 + p1 - 1) < 1e-12)
            #expect(abs(p0 - pow(cos(theta / 2), 2)) < 1e-12)
            let r = BlochMath.blochVector(theta: theta, phi: phi)
            #expect(abs(r.rz - (p0 - p1)) < 1e-12)
            #expect(abs(r.rz - cos(theta)) < 1e-12)
            #expect(abs(BlochMath.vonNeumannEntropyBits(r)) < 1e-12, "纯态 S = 0 bit")
        }
    }

    // MARK: 物理律 4：最大混态

    @Test("物理律4：最大混态 r=0、纯度 1/2、S=1 bit")
    func maximallyMixed() {
        let r = (rx: 0.0, ry: 0.0, rz: 0.0)
        #expect(BlochMath.purity(r) == 0.5)
        #expect(abs(BlochMath.vonNeumannEntropyBits(r) - 1) < 1e-12)
        // 中间混态：纯度 (1+|r|²)/2、熵介于 0 与 1
        for radius in stride(from: 0.1, through: 0.9, by: 0.2) {
            let rm = (rx: radius, ry: 0.0, rz: 0.0)
            #expect(abs(BlochMath.purity(rm) - (1 + radius * radius) / 2) < 1e-12)
            let s = BlochMath.vonNeumannEntropyBits(rm)
            #expect(s > 0 && s < 1)
        }
    }

    // MARK: 物理律 5：脚本锚点 θ=π/3, φ=π/4

    @Test("物理律5：脚本锚点 p0=0.75、r=(0.61237, 0.35355, 0.5)")
    func scriptAnchor() {
        let sv = BlochMath.stateVector(theta: BlochMath.anchorTheta, phi: BlochMath.anchorPhi)
        let p0 = sv.alphaRe * sv.alphaRe
        let p1 = sv.betaRe * sv.betaRe + sv.betaIm * sv.betaIm
        #expect(abs(p0 - 0.75) < 1e-12, "脚本打印 |alpha|^2 = 0.75")
        #expect(abs(p1 - 0.25) < 1e-12, "脚本打印 |beta|^2 = 0.25")
        let r = BlochMath.blochVector(theta: BlochMath.anchorTheta, phi: BlochMath.anchorPhi)
        // sin(π/3)·cos(π/4) = sin(π/3)·sin(π/4)（cos(π/4)=sin(π/4)）→ rx 与 ry 相等
        #expect(abs(r.rx - 0.6123724356957945) < 1e-12)
        #expect(abs(r.ry - 0.6123724356957945) < 1e-12)
        #expect(abs(r.rz - 0.5) < 1e-12)
        #expect(abs(sqrt(r.rx * r.rx + r.ry * r.ry + r.rz * r.rz) - 1) < 1e-12,
                "脚本打印 |b| = 1.0000")
    }

    // MARK: compute 输出与预算

    @Test("compute：曲线 p0/p1 网格正确 + 实时档 16 ms 预算")
    func computeAndBudget() async throws {
        let cs = try constants()
        let result = try await module.compute(values(), constants: cs)
        let chart = try #require(chartLineSeries(result, 1))  // W13：曲线图移至 index 1
        #expect(chart.series.count == 2)
        // |α|² = cos²(θ/2) 在 θ=π 处为 0，|β|² 为 1；θ=0 反之
        let s0 = chart.series[0].points
        let s1 = chart.series[1].points
        #expect(abs(s0.first!.y - 1) < 1e-12 && abs(s0.last!.y) < 1e-12)
        #expect(abs(s1.first!.y) < 1e-12 && abs(s1.last!.y - 1) < 1e-12)
        for i in s0.indices {
            #expect(abs(s0[i].y + s1[i].y - 1) < 1e-12, "逐点归一化")
        }
        try await expectComputeUnderBudget(module: module, values: values(),
                                           constants: cs, budgetMillis: 16)
    }
}
