import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/光场量子化_二阶相干度.py）

/// 二阶相干度纯函数核。
/// g⁽²⁾(0) = ⟨n(n−1)⟩/⟨n⟩²（式 12）：相干态 = 1（泊松）、热态 = 2（玻色-爱因斯坦）、
/// Fock |n⟩ = 1 − 1/n（反聚束）；压缩真空 Wigner W(x,p) = (1/π)exp(−x²e^{2r} − p²e^{−2r})。
enum SecondOrderCoherenceMath {

    /// 矩求和域上限：脚本固定 120（μ=5 时热态尾部质量 ~3×10⁻¹⁰）；
    /// μ 可调后按尾部覆盖准则扩展（与 PhotonStatsMath 同准则：μ + 15σ，σ = √(μ(1+μ))）。
    static func nCompute(forMu mu: Double) -> Int {
        let sigma = (mu * (1.0 + mu)).squareRoot()
        return max(120, Int((mu + 15.0 * sigma).rounded(.up)))
    }

    /// 相干态（泊松）g⁽²⁾(0)：截断矩比（对数域求 P(n) 防大 n 溢出）。
    static func g2Coherent(mu: Double, nmax: Int = 120) -> Double {
        let n = mu <= 0 ? 1 : max(nmax, nCompute(forMu: mu))
        var mean = 0.0, mom2 = 0.0
        for k in 0...n {
            let p = exp(-mu + Double(k) * log(mu) - Num.lnFactorial(k))
            let dk = Double(k)
            mean += dk * p
            mom2 += dk * (dk - 1.0) * p
        }
        return mom2 / (mean * mean)
    }

    /// 热态（玻色-爱因斯坦）g⁽²⁾(0)：截断矩比（脚本 q = μ/(1+μ) 递推）。
    static func g2Thermal(mu: Double, nmax: Int = 120) -> Double {
        let n = mu <= 0 ? 1 : max(nmax, nCompute(forMu: mu))
        let q = mu / (1.0 + mu)
        let norm = 1.0 / (1.0 + mu)
        var mean = 0.0, mom2 = 0.0, qk = 1.0
        for k in 0...n {
            let p = norm * qk
            let dk = Double(k)
            mean += dk * p
            mom2 += dk * (dk - 1.0) * p
            qk *= q
        }
        return mom2 / (mean * mean)
    }

    /// Fock 态 |n⟩ 的 g⁽²⁾(0) = 1 − 1/n（n ≥ 1；n = 1 反聚束极限为 0）。
    static func g2Fock(n: Int) -> Double {
        guard n >= 1 else { return .nan }
        return 1.0 - 1.0 / Double(n)
    }

    /// 压缩真空 Wigner 函数 W(x,p) = (1/π)·exp(−x²e^{2r} − p²e^{−2r})。
    static func squeezedVacuumWigner(x: Double, p: Double, r: Double) -> Double {
        exp(-x * x * exp(2.0 * r) - p * p * exp(-2.0 * r)) / Double.pi
    }

    /// 1σ 等值线高度（指数 = −1）：W₁σ = e⁻¹/π。
    static func oneSigmaLevel() -> Double { exp(-1.0) / Double.pi }

    /// 压缩/反压缩正交分量方差（真空 1/2）：⟨x²⟩ = e^{−2r}/2，⟨p²⟩ = e^{+2r}/2。
    static func quadratureVariances(r: Double) -> (x: Double, p: Double) {
        (exp(-2.0 * r) / 2.0, exp(2.0 * r) / 2.0)
    }
}

// MARK: - 模块

/// 笔记 35《光场量子化与相干态》：二阶相干度 g⁽²⁾(0) 三态对比
/// （相干 = 1 / 热态 = 2 / Fock = 反聚束）与压缩真空 Wigner 函数等高线（320×320）。
struct SecondOrderCoherenceModule: SimModule {

    let meta = ModuleMeta(
        id: "光场量子化_二阶相干度", title: "二阶相干度 · g⁽²⁾(0)",
        subtitle: "相干态 = 1、热态 = 2、Fock 反聚束 < 1——光子统计的「指纹」",
        category: .quantumOptics, noteNumber: 35, tier: .seconds, difficulty: .advanced,
        keywords: ["二阶相干度", "g2", "光子统计", "反聚束", "聚束", "压缩真空", "Wigner 函数",
                   "second-order coherence", "antibunching", "squeezed vacuum", "Wigner function"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "mu", title: "平均光子数", symbol: "μ", unit: "",
                               range: 0.5...50, defaultValue: 5,
                               scale: .log, decimalPlaces: 2)),
            .slider(SliderSpec(key: "r", title: "压缩参量", symbol: "r", unit: "",
                               range: 0...2, defaultValue: 1,
                               decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .bars(BarSpec(
                title: "零延迟二阶相干度 g⁽²⁾(0) 三态对比",
                xAxis: .init(label: "光场态"),
                yAxis: .init(label: "g⁽²⁾(0)"))),
            .contour(ContourSpec(
                title: "压缩真空 Wigner 函数（320×320 等高线）",
                xAxis: .init(label: "x（压缩正交分量）"),
                yAxis: .init(label: "p（反压缩正交分量）"),
                valueLabel: "W(x, p)",
                levels: 30)),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let mu = input.slider("mu")
        let r = input.slider("r")

        // 图 1：三态 g⁽²⁾(0)（脚本 μ=5 基线；μ 可调演示 g⁽²⁾ 与光强的无关性）
        let g2Coh = SecondOrderCoherenceMath.g2Coherent(mu: mu)
        let g2Th = SecondOrderCoherenceMath.g2Thermal(mu: mu)
        let g2Fock1 = SecondOrderCoherenceMath.g2Fock(n: 1)

        // 图 2：压缩真空 Wigner（320×320；窗口随 r 自适应保 3σ 覆盖 + 等比纵横）
        let nGrid = 320
        let xMax = 3.0 * exp(-r)
        let pMax = 3.0 * exp(r)
        let xs = Num.linspace(-xMax, xMax, count: nGrid)
        let ps = Num.linspace(-pMax, pMax, count: nGrid)
        var w: [[Double]] = []
        w.reserveCapacity(nGrid)
        for p in ps {
            var row = [Double]()
            row.reserveCapacity(nGrid)
            for x in xs {
                row.append(SecondOrderCoherenceMath.squeezedVacuumWigner(x: x, p: p, r: r))
            }
            w.append(row)
        }

        let (varX, varP) = SecondOrderCoherenceMath.quadratureVariances(r: r)

        return SimResult(
            charts: [
                .bars(BarData(
                    spec: charts[0].barSpec!,
                    bars: [
                        .init(label: "相干态（泊松）", value: g2Coh),
                        .init(label: "热态（玻色-爱因斯坦）", value: g2Th),
                        .init(label: "Fock |1⟩（反聚束）", value: g2Fock1),
                    ])),
                .contour(ContourData(
                    spec: charts[1].contourSpec!,
                    xGrid: xs, yGrid: ps, values: w,
                    highlightLevels: [SecondOrderCoherenceMath.oneSigmaLevel()])),
            ],
            summary: [
                .init(id: "coh", title: "相干态 g⁽²⁾(0)",
                      value: String(format: "%.6f", g2Coh),
                      note: "= 1：泊松统计（散粒噪声极限）"),
                .init(id: "th", title: "热态 g⁽²⁾(0)",
                      value: String(format: "%.6f", g2Th),
                      note: "= 2：超泊松（混沌光聚束）"),
                .init(id: "fock", title: "Fock |1⟩ g⁽²⁾(0)",
                      value: String(format: "%.6f", g2Fock1),
                      note: "1 − 1/n = 0：单光子反聚束（理想单光子源）"),
                .init(id: "mu", title: "平均光子数 μ",
                      value: String(format: "%.2f", mu),
                      note: "三态 g⁽²⁾(0) 均与 μ 无关——态的「指纹」而非强度"),
                .init(id: "squeeze", title: "压缩正交方差",
                      value: String(format: "⟨x²⟩ = %.4f，⟨p²⟩ = %.4f", varX, varP),
                      note: "真空各 1/2；乘积 ⟨x²⟩⟨p²⟩ = 1/4（最小不确定度）"),
            ],
            theory: TheoryCard(
                title: "二阶相干度与压缩真空（笔记 35 式 12）",
                formulas: [
                    "g⁽²⁾(0) = ⟨n(n−1)⟩/⟨n⟩²——零延迟强度关联",
                    "相干态（泊松）：g⁽²⁾ = 1；热态（BE）：g⁽²⁾ = 2；Fock |n⟩：g⁽²⁾ = 1 − 1/n",
                    "压缩真空 Wigner：W(x,p) = (1/π)exp(−x²e^{2r} − p²e^{−2r})",
                    "1σ 等值线：x²e^{2r} + p²e^{−2r} = 1；⟨x²⟩ = e^{−2r}/2 < 1/2（低于真空）",
                ],
                reading: "左图：g⁽²⁾(0) 是光场量子态的指纹——热态超泊松（2，聚束），"
                    + "相干态恰为 1（散粒噪声），Fock 单光子态为 0（反聚束，Hanbury Brown-Twiss"
                    + "零符合）。右图：压缩真空的 Wigner 椭圆沿 x 压扁、沿 p 拉长（白色虚线为 1σ 圈），"
                    + "两方向方差乘积仍是最小不确定度 1/4——压缩不违反海森伯关系，只是重新分配噪声。"))
    }
}