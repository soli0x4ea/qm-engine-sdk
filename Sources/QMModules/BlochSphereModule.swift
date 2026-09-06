import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/希尔伯特空间与狄拉克符号_Bloch球.py，式 6.1）

/// 两能级纯态 Bloch 球纯函数核：θ/φ 参数化 ↔ 密度矩阵 Pauli 迹双路径、
/// 对跖点正交、纯度与 von Neumann 熵。**3D 数据模型（SDK 交付物）**：
/// SceneKit 渲染层（剩余工作）只需消费 `BlochSphereModule` summary/series 与
/// 本枚举的纯函数，即可重建单位球 + 态矢量场景；图表契约未追加 blochSphere 图种
/// （穷尽 switch 需连带渲染器，超出本简化验收线）。
enum BlochMath {

    /// 态矢 |ψ⟩ = cos(θ/2)|0⟩ + e^{iφ} sin(θ/2)|1⟩ → (α 实部, β 实部, β 虚部)。
    /// α 恒为实（整体相位规范固定），β = e^{iφ} sin(θ/2)。
    static func stateVector(theta: Double, phi: Double) -> (alphaRe: Double, betaRe: Double, betaIm: Double) {
        (cos(theta / 2), sin(theta / 2) * cos(phi), sin(theta / 2) * sin(phi))
    }

    /// 密度矩阵 Pauli 迹路径（脚本 bloch_vector 的闭式等价）：
    /// r = Tr(ρσ)，ρ = |ψ⟩⟨ψ|。α 实 → rx = 2αRe(β)，ry = −2Im(αβ*) = 2αIm(β)，rz = |α|²−|β|²。
    static func blochVectorFromState(alphaRe: Double, betaRe: Double, betaIm: Double)
        -> (rx: Double, ry: Double, rz: Double) {
        (2 * alphaRe * betaRe, 2 * alphaRe * betaIm,
         alphaRe * alphaRe - betaRe * betaRe - betaIm * betaIm)
    }

    /// 参数化路径（脚本正文）：r = (sinθcosφ, sinθsinφ, cosθ)。
    static func blochVector(theta: Double, phi: Double) -> (rx: Double, ry: Double, rz: Double) {
        (sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
    }

    /// 对跖点（正交补态）：(θ, φ) → (π−θ, φ+π mod 2π)。
    static func antipode(theta: Double, phi: Double) -> (theta: Double, phi: Double) {
        ((.pi - theta).truncatingRemainder(dividingBy: .pi),
         (phi + .pi).truncatingRemainder(dividingBy: 2 * .pi))
    }

    /// 重叠平方 |⟨ψ₁|ψ₂⟩|² = (1 + r₁·r₂)/2（Bloch 几何：夹角 γ 的 cos²(γ/2)）。
    static func overlapSquared(_ r1: (rx: Double, ry: Double, rz: Double),
                               _ r2: (rx: Double, ry: Double, rz: Double)) -> Double {
        let dot = r1.rx * r2.rx + r1.ry * r2.ry + r1.rz * r2.rz
        return (1 + dot) / 2
    }

    /// 纯度 Tr(ρ²) = (1 + |r|²)/2：纯态 1、最大混态 1/2。
    static func purity(_ r: (rx: Double, ry: Double, rz: Double)) -> Double {
        (1 + r.rx * r.rx + r.ry * r.ry + r.rz * r.rz) / 2
    }

    /// von Neumann 熵 S = −Tr(ρ log₂ρ)（bit）。二能级本征值 (1±|r|)/2 → 二元熵；
    /// 与脚本 eigvalsh + clip(1e-12) 数值路径等价（本征值 0 处 0·log0 := 0）。
    static func vonNeumannEntropyBits(_ r: (rx: Double, ry: Double, rz: Double)) -> Double {
        let norm = sqrt(r.rx * r.rx + r.ry * r.ry + r.rz * r.rz)
        let p = min(max((1 + norm) / 2, 0), 1)
        let q = 1 - p
        func xlog2(_ x: Double) -> Double { x <= 1e-12 ? 0 : x * log2(x) }
        return -(xlog2(p) + xlog2(q))
    }

    /// 脚本示例锚点：θ = π/3、φ = π/4。
    static let anchorTheta = Double.pi / 3
    static let anchorPhi = Double.pi / 4
}

// MARK: - 模块：Bloch 球（笔记 11 · 实时档）

/// 笔记 11《希尔伯特空间与狄拉克符号》：两能级纯态 θ/φ 参数化 → 态矢量 / Bloch 矢量，
/// 纯态族 |α|²、|β|² 随 θ 的演化。计算核 + 数据模型为 SDK 交付物（SceneKit 线框球渲染属 UI 层）。
struct BlochSphereModule: SimModule {

    let meta = ModuleMeta(
        id: "希尔伯特空间与狄拉克符号_Bloch球", title: "Bloch 球 · 两能级纯态参数化",
        subtitle: "|ψ⟩=cos(θ/2)|0⟩+e^{iφ}sin(θ/2)|1⟩：单位球面上的态矢量与对跖正交",
        category: .quantumStates, noteNumber: 11, tier: .realtime, difficulty: .basic,
        keywords: ["Bloch球", "Bloch", "希尔伯特", "狄拉克符号", "态矢量", "极角", "方位角",
                   "纯态", "布洛赫", "qubit", "Bloch sphere", "antipode", "对跖"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "theta", title: "极角 θ", symbol: "θ", unit: "rad",
                              range: 0...(2 * .pi), defaultValue: BlochMath.anchorTheta,
                              scale: .linear, decimalPlaces: 3)),
            .slider(SliderSpec(key: "phi", title: "方位角 φ", symbol: "φ", unit: "rad",
                              range: 0...(2 * .pi), defaultValue: BlochMath.anchorPhi,
                              scale: .linear, decimalPlaces: 3)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .bloch(BlochSpec(title: "Bloch 球：当前态矢 |ψ⟩ 与对跖点 |ψ⊥⟩（随 θ/φ 实时转动）")),
            .lineSeries(LineSeriesSpec(
                title: "纯态族：|α|² 与 |β|² 随极角 θ（φ 固定）",
                xAxis: .init(label: "极角 θ (rad)"),
                yAxis: .init(label: "概率幅平方"),
                seriesNames: ["|α|² = cos²(θ/2)", "|β|² = sin²(θ/2)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let theta = input.slider("theta")
        let phi = input.slider("phi")

        // 纯态族曲线（φ 固定）：极角全域 θ∈[0,π]，|α|² 1→0、|β|² 0→1
        let thGrid = Num.linspace(0.0, Double.pi, count: 240)
        let p0 = thGrid.map { pow(cos($0 / 2), 2) }
        let p1 = thGrid.map { pow(sin($0 / 2), 2) }

        // 当前态：双路径态矢 / Bloch 矢量
        let sv = BlochMath.stateVector(theta: theta, phi: phi)
        let rFromState = BlochMath.blochVectorFromState(alphaRe: sv.alphaRe,
                                                        betaRe: sv.betaRe, betaIm: sv.betaIm)
        let rFromParam = BlochMath.blochVector(theta: theta, phi: phi)
        let dualPathDev = max(abs(rFromState.rx - rFromParam.rx),
                              abs(rFromState.ry - rFromParam.ry),
                              abs(rFromState.rz - rFromParam.rz))
        let norm = sqrt(rFromParam.rx * rFromParam.rx + rFromParam.ry * rFromParam.ry
                        + rFromParam.rz * rFromParam.rz)
        let anti = BlochMath.antipode(theta: theta, phi: phi)
        let rAnti = BlochMath.blochVector(theta: anti.theta, phi: anti.phi)
        let overlapAnti = BlochMath.overlapSquared(rFromParam, rAnti)

        func fmt(_ v: Double) -> String { String(format: "%.5f", v) }

        return SimResult(
            charts: [
                .bloch(BlochData(
                    spec: charts.requireBloch(0),
                    state: Point3D(x: rFromParam.rx, y: rFromParam.ry, z: rFromParam.rz),
                    stateLabel: "|ψ⟩",
                    antipode: Point3D(x: rAnti.rx, y: rAnti.ry, z: rAnti.rz),
                    antipodeLabel: "|ψ⊥⟩",
                    northLabel: "|0⟩", southLabel: "|1⟩")),
                .lineSeries(LineSeriesData(
                    spec: charts.requireLineSeries(1),
                    series: [
                        .init(name: "|α|² = cos²(θ/2)",
                              points: Num.strided(thGrid, p0, stride: 1)),
                        .init(name: "|β|² = sin²(θ/2)",
                              points: Num.strided(thGrid, p1, stride: 1)),
                    ],
                    referenceLines: [
                        ReferenceLine(label: "θ = \(fmt(theta))", axis: .x, value: theta, style: .threshold),
                    ],
                    pointMarkers: [
                        PointMarker(
                            x: theta,
                            y: pow(cos(theta / 2), 2),
                            label: String(format: "θ=%.3f", theta),
                            colorIndex: 0),
                    ])),
            ],
            summary: [
                .init(id: "state", title: "态矢量（α 实规范）",
                      value: String(format: "α=%.5f, β=%.5f+%.5fi", sv.alphaRe, sv.betaRe, sv.betaIm),
                      note: String(format: "θ=%.5f, φ=%.5f rad", theta, phi)),
                .init(id: "p", title: "|α|² / |β|²",
                      value: "\(fmt(pow(cos(theta / 2), 2))) / \(fmt(pow(sin(theta / 2), 2)))",
                      note: "和恒为 1（归一化）"),
                .init(id: "bloch", title: "Bloch 矢量 (rx, ry, rz)",
                      value: "(\(fmt(rFromParam.rx)), \(fmt(rFromParam.ry)), \(fmt(rFromParam.rz)))",
                      note: "sinθcosφ, sinθsinφ, cosθ（脚本正文参数化）"),
                .init(id: "norm", title: "|b|（纯态恒为 1）",
                      value: fmt(norm),
                      note: String(format: "双路径核对：Pauli 迹 vs 参数化，最大偏差 %.1e", dualPathDev)),
                .init(id: "purity", title: "纯度 Tr(ρ²)",
                      value: fmt(BlochMath.purity(rFromParam)), note: "(1+|r|²)/2，纯态 = 1"),
                .init(id: "anti", title: "对跖点 (π−θ, φ+π) 重叠 |⟨ψ|ψ⊥⟩|²",
                      value: String(format: "%.3e", overlapAnti),
                      note: "恒为 0：对跖态正交（(1+r₁·r₂)/2）"),
                .init(id: "entropy", title: "von Neumann 熵 S",
                      value: fmt(BlochMath.vonNeumannEntropyBits(rFromParam)) + " bit",
                      note: "纯态 S=0；最大混态（r=0）S=1 bit"),
            ],
            theory: TheoryCard(
                title: "Bloch 球与两能级纯态（笔记 11 式 6.1）",
                formulas: [
                    "|ψ⟩ = cos(θ/2)|0⟩ + e^{iφ} sin(θ/2)|1⟩（归一化，整体相位规范固定 α 实）",
                    "Bloch 矢量 r = (sinθcosφ, sinθsinφ, cosθ) = Tr(ρσ)，|r| = 1（纯态）",
                    "对跖点 (π−θ, φ+π) ↔ 正交补态：|⟨ψ|ψ⊥⟩|² = (1 + r₁·r₂)/2 = 0",
                    "纯度 Tr(ρ²) = (1+|r|²)/2；S = −Tr(ρ log₂ρ)：纯态 0 bit，最大混态 1 bit",
                ],
                reading: "二维复 Hilbert 空间的归一化态（扣除整体相位）一一对应单位球面上的点："
                    + "极角 θ 给出基矢概率幅 cos²(θ/2)/sin²(θ/2)，方位角 φ 是相对相位。"
                    + "球心是最大混态 I/2（r=0，S=1 bit），球面是纯态（|r|=1，S=0）。"
                    + "对跖两点对应一对正交态——测量基的选取就是在球面上选一条直径。"))
    }
}
