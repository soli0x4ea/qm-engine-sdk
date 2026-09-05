import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/光场量子化_相干态与光子统计.py）

/// 光场量子态光子数分布纯函数核。
/// 相干态（式 8）：P(n) = e^(−|α|²)|α|^(2n)/n!（泊松）；
/// 热态（式 9）：P(n) = (1+μ)⁻¹(μ/(1+μ))ⁿ（玻色-爱因斯坦）。
enum PhotonStatsMath {

    /// 相干态泊松分布 P(n)（对数域求值防大 n 溢出；与 Python 逐位差 < 1e-15）。
    static func coherent(_ n: Int, mu: Double) -> Double {
        exp(-mu + Double(n) * log(mu) - Num.lnFactorial(n))
    }

    /// 热态玻色-爱因斯坦分布 P(n)。
    static func thermal(_ n: Int, mu: Double) -> Double {
        exp(-log(1 + mu) + Double(n) * (log(mu) - log(1 + mu)))
    }

    /// 统计求和域上限：Python 固定 120（μ=5 时热态方差 30 已收敛）；
    /// μ 可调后按尾部覆盖准则扩展：μ + 15σ（σ = √(μ(1+μ))），保 ΣP ≈ 1。
    static func nCompute(forMu mu: Double) -> Int {
        let sigma = sqrt(mu * (1 + mu))
        return max(120, Int(ceil(mu + 15 * sigma)))
    }

    /// 前 n+1 项分布数组。
    static func distribution(_ f: (Int, Double) -> Double, mu: Double, n: Int) -> [Double] {
        (0...n).map { f($0, mu) }
    }

    /// 一阶/二阶矩：(⟨n⟩, ⟨n²⟩)。
    static func moments(_ probs: [Double]) -> (mean: Double, meanSq: Double) {
        var m1 = 0.0, m2 = 0.0
        for (n, p) in probs.enumerated() {
            let dn = Double(n)
            m1 += dn * p
            m2 += dn * dn * p
        }
        return (m1, m2)
    }
}

// MARK: - 模块

/// 笔记 35《光场量子化与相干态》：同均值 μ 下相干态（泊松）与热态（玻色-爱因斯坦）
/// 的光子数分布对比——亚泊松/超泊松统计与散粒噪声极限。
struct PhotonStatisticsModule: SimModule {

    let meta = ModuleMeta(
        id: "光场量子化_相干态与光子统计", title: "相干态与光子统计",
        subtitle: "P(n) = e^(−μ)μⁿ/n!（泊松）vs 玻色-爱因斯坦：Δn 的两个极限",
        category: .quantumOptics, noteNumber: 35, tier: .realtime, difficulty: .basic,
        keywords: ["相干态", "光子统计", "泊松分布", "玻色-爱因斯坦", "热态", "散粒噪声",
                   "光场量子化", "Poisson", "coherent", "thermal", "Mandel Q"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "mu", title: "平均光子数", symbol: "μ = |α|²", unit: "",
                               range: 0.1...50, defaultValue: 5,
                               scale: .log, decimalPlaces: 2)),
            .slider(SliderSpec(key: "nshow", title: "显示光子数上限", symbol: "n", unit: "",
                               range: 10...60, defaultValue: 25,
                               step: 5, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .bars(BarSpec(
                title: "相干态 P(n)（泊松）",
                xAxis: .init(label: "光子数 n"),
                yAxis: .init(label: "概率 P(n)"))),
            .bars(BarSpec(
                title: "热态 P(n)（玻色-爱因斯坦）",
                xAxis: .init(label: "光子数 n"),
                yAxis: .init(label: "概率 P(n)"))),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let mu = input.slider("mu")
        let nShow = Int(input.slider("nshow"))
        let nMax = PhotonStatsMath.nCompute(forMu: mu)

        // ln(n!) 增量数组（对齐 Python fact 数组的逐句构造，O(n)；单点函数仅供测试对拍）
        var lnFact = [Double](repeating: 0, count: nMax + 1)
        if nMax >= 2 {
            for i in 2...nMax { lnFact[i] = lnFact[i - 1] + log(Double(i)) }
        }
        let lnMu = log(mu), ln1Mu = log(1 + mu), lnQ = lnMu - ln1Mu
        var coh = [Double](repeating: 0, count: nMax + 1)
        var th = [Double](repeating: 0, count: nMax + 1)
        for n in 0...nMax {
            let dn = Double(n)
            coh[n] = exp(-mu + dn * lnMu - lnFact[n])
            th[n] = exp(-ln1Mu + dn * lnQ)
        }

        let sumCoh = coh.reduce(0, +), sumTh = th.reduce(0, +)
        let (mCoh, m2Coh) = PhotonStatsMath.moments(coh)
        let (mTh, m2Th) = PhotonStatsMath.moments(th)
        let dCoh = sqrt(max(0, m2Coh - mCoh * mCoh))
        let dTh = sqrt(max(0, m2Th - mTh * mTh))

        let cohShow = coh[0...min(nShow, nMax)]
        let thShow = th[0...min(nShow, nMax)]
        let cohBars = cohShow.enumerated().map {
            BarItem(label: "\($0.offset)", value: $0.element)
        }
        let thBars = thShow.enumerated().map {
            BarItem(label: "\($0.offset)", value: $0.element)
        }

        return SimResult(
            charts: [
                .bars(BarData(spec: charts[0].barSpec!, bars: cohBars)),
                .bars(BarData(spec: charts[1].barSpec!, bars: thBars)),
            ],
            summary: [
                .init(id: "mean_coh", title: "⟨n⟩ 相干", value: String(format: "%.6g", mCoh),
                      note: String(format: "ΣP = %.10f", sumCoh)),
                .init(id: "delta_coh", title: "Δn 相干", value: String(format: "%.6g", dCoh),
                      note: "泊松极限 √μ（散粒噪声）"),
                .init(id: "mean_th", title: "⟨n⟩ 热态", value: String(format: "%.6g", mTh),
                      note: String(format: "ΣP = %.10f", sumTh)),
                .init(id: "delta_th", title: "Δn 热态", value: String(format: "%.6g", dTh),
                      note: "超泊松 √(μ(1+μ))"),
                .init(id: "excess", title: "Mandel Q(热态)",
                      value: String(format: "%.4g", (m2Th - mTh * mTh) / mTh - 1),
                      note: "超泊松度 = (Δn)²/⟨n⟩ − 1 = μ"),
            ],
            theory: TheoryCard(
                title: "光子数分布（笔记 35 式 8/9）",
                formulas: [
                    "相干态 |α⟩：P(n) = e^(−|α|²)|α|^(2n)/n!，μ = |α|²",
                    "热态：P(n) = (1+μ)⁻¹ · (μ/(1+μ))ⁿ",
                    "相干 Δn = √μ；热态 Δn = √(μ(1+μ))",
                ],
                reading: "同均值下热态长尾、双峰性消失：热噪声是散粒噪声的 √(1+μ) 倍——这决定光电探测的量子噪声极限。"))
    }
}
