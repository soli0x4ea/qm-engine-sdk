import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/不确定性关系_等号与势阱.py）

/// 高斯波包（等号）与无限深势阱（超等号）的不确定度纯函数核。
/// 全部无量纲（动量以 p/ħ 计，与脚本打印口径一致）；
/// 高斯动量分布走 EngineKit FFT 数据层（vDSP，W3 裁决的 Accelerate 路线）。
enum UncertaintyMath {

    /// 脚本 trapz：均匀网格 dx·(Σy − (y₀+y_{n−1})/2)。
    static func trapz(_ y: [Double], dx: Double) -> Double {
        var s = 0.0
        for v in y { s += v }
        return dx * (s - 0.5 * (y[0] + y[y.count - 1]))
    }

    /// 高斯波包：ψ = (2πσ₀²)^(−1/4)·exp(−x²/4σ₀²)，x ∈ [−12σ₀, 12σ₀] × 4000。
    static func gaussianState(sigma0: Double, count: Int = 4000)
        -> (x: [Double], psi: [Double]) {
        let x = Num.linspace(-12 * sigma0, 12 * sigma0, count: count)
        let norm = pow(2 * .pi * sigma0 * sigma0, -0.25)
        return (x, x.map { norm * exp(-$0 * $0 / (4 * sigma0 * sigma0)) })
    }

    /// FFT 动量分布与矩（脚本口径：psif = fftshift(fft(ifftshift(ψ)))·dx）。
    /// k = fftshift(fftfreq(N,d))·2π；Δx、Δp/ħ 由 trapz 加权矩给出。
    static func gaussianMoments(
        sigma0: Double, count: Int = 4000
    ) -> (k: [Double], density: [Double], dx: Double, dpOverHbar: Double) {
        let (x, psi) = gaussianState(sigma0: sigma0, count: count)
        let dxStep = x[1] - x[0]

        var sq = [Double](repeating: 0, count: count)
        var x2w = [Double](repeating: 0, count: count)
        for i in 0..<count {
            sq[i] = psi[i] * psi[i]
            x2w[i] = x[i] * x[i] * sq[i]
        }
        let norm = trapz(sq, dx: dxStep)
        let dxStd = (trapz(x2w, dx: dxStep) / norm).squareRoot()

        let shifted = FFT.ifftshift(psi)
        let f = FFT.fft(shifted)
        var psifRe = FFT.fftshift(f.re)
        var psifIm = FFT.fftshift(f.im)
        for i in 0..<count {
            psifRe[i] *= dxStep
            psifIm[i] *= dxStep
        }
        let k = FFT.fftshift(FFT.fftfreq(count, d: dxStep)).map { $0 * 2 * .pi }

        var pd = [Double](repeating: 0, count: count)
        var k2w = [Double](repeating: 0, count: count)
        for i in 0..<count {
            pd[i] = psifRe[i] * psifRe[i] + psifIm[i] * psifIm[i]
            k2w[i] = k[i] * k[i] * pd[i]
        }
        let dk = abs(k[1] - k[0])
        let knorm = trapz(pd, dx: dk)
        let dpStd = (trapz(k2w, dx: dk) / knorm).squareRoot()
        return (k, pd, dxStd, dpStd)
    }

    /// 势阱基态：ψ₁ = √(2/L)·sin(πx/L)，x ∈ [0, L] × 4000。
    static func wellGround(L: Double, count: Int = 4000)
        -> (x: [Double], psi: [Double]) {
        let x = Num.linspace(0, L, count: count)
        let norm = (2 / L).squareRoot()
        return (x, x.map { norm * sin(.pi * $0 / L) })
    }

    /// 势阱矩：位置直接积分；动量用有限差分二阶导（<p²>/ħ² = −∫ψψ″dx/norm）。
    static func wellMoments(
        L: Double, count: Int = 4000
    ) -> (x: [Double], density: [Double], dx: Double, dpOverHbar: Double) {
        let (x, psi) = wellGround(L: L, count: count)
        let dxStep = x[1] - x[0]

        var sq = [Double](repeating: 0, count: count)
        var xw = [Double](repeating: 0, count: count)
        var x2w = [Double](repeating: 0, count: count)
        for i in 0..<count {
            sq[i] = psi[i] * psi[i]
            xw[i] = x[i] * sq[i]
            x2w[i] = x[i] * x[i] * sq[i]
        }
        let norm = trapz(sq, dx: dxStep)
        let xmean = trapz(xw, dx: dxStep) / norm
        let x2 = trapz(x2w, dx: dxStep) / norm
        let dxStd = (x2 - xmean * xmean).squareRoot()

        var second = [Double](repeating: 0, count: count)
        for i in 1..<(count - 1) {
            second[i] = (psi[i + 1] - 2 * psi[i] + psi[i - 1]) / (dxStep * dxStep)
        }
        var prod = [Double](repeating: 0, count: count)
        for i in 0..<count { prod[i] = psi[i] * second[i] }
        let p2 = -trapz(prod, dx: dxStep) / norm
        return (x, sq, dxStd, p2.squareRoot())
    }

    /// 势阱基态解析参考：Δx = L·√((π²−6)/(12π²))，Δp/ħ = π/L。
    static func wellAnalytic(L: Double) -> (dx: Double, dpOverHbar: Double) {
        let dx = L * ((.pi * .pi - 6) / (12 * .pi * .pi)).squareRoot()
        return (dx, .pi / L)
    }
}

// MARK: - 模块

/// 笔记 08《不确定性关系》：高斯最小不确定态（等号）vs 无限深势阱基态（超等号）。
/// 秒级档：4000 点 FFT + 有限差分（轻量数值；动量无量纲 p/ħ 口径与脚本一致）。
struct UncertaintyModule: SimModule {

    let meta = ModuleMeta(
        id: "不确定性关系_等号与势阱", title: "不确定性关系 · 等号与势阱",
        subtitle: "高斯 ΔxΔp = ħ/2 精确取等；势阱基态 ΔxΔp ≈ 0.568ħ",
        category: .quantumStates, noteNumber: 8, tier: .seconds, difficulty: .basic,
        keywords: ["不确定性", "海森伯", "Heisenberg", "高斯波包", "势阱",
                   "最小不确定态", "傅里叶", "FFT", "方差"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "sigma0", title: "高斯宽度", symbol: "σ₀", unit: "arb.",
                               range: 0.5...2.0, defaultValue: 1.0,
                               scale: .linear, decimalPlaces: 2)),
            .slider(SliderSpec(key: "L", title: "势阱宽度", symbol: "L", unit: "arb.",
                               range: 0.5...2.0, defaultValue: 1.0,
                               scale: .linear, decimalPlaces: 2)),
            .constant(ConstantSpec(key: "hbar", title: "ħ", note: "约化普朗克常量（矩以 p/ħ 无量纲计）")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "高斯波包：Δx·Δp = ħ/2（等号）",
                xAxis: .init(label: "x（左）/ p/3ħ（右，任意单位）"),
                yAxis: .init(label: "归一化密度"),
                seriesNames: ["|ψ(x)|²", "|φ(p)|²"])),
            .lineSeries(LineSeriesSpec(
                title: "无限深势阱基态：Δx·Δp ≈ 0.568 ħ > ħ/2",
                xAxis: .init(label: "x / L"),
                yAxis: .init(label: "归一化 |ψ₁|²"),
                seriesNames: ["|ψ₁(x)|²"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let sigma0 = input.slider("sigma0")
        let L = input.slider("L")

        // 高斯：位置 + 动量分布（脚本口径：动量按 |p| 排序后画在 k/3 轴对照）
        let (k, pd, dxG, dpG) = UncertaintyMath.gaussianMoments(sigma0: sigma0)
        let (gx, gpsi) = UncertaintyMath.gaussianState(sigma0: sigma0)
        var gPeak = 0.0
        for v in gpsi where v > gPeak { gPeak = v }
        var positionPts: [Point] = []
        positionPts.reserveCapacity(2000)
        for i in stride(from: 0, to: gx.count, by: 2) {
            positionPts.append(Point(x: gx[i], y: gpsi[i] * gpsi[i] / gPeak))
        }
        let order = k.indices.sorted {
            abs(k[$0]) != abs(k[$1]) ? abs(k[$0]) < abs(k[$1]) : k[$0] < k[$1]
        }
        var pdPeak = 0.0
        for v in pd where v > pdPeak { pdPeak = v }
        var momentumPts: [Point] = []
        momentumPts.reserveCapacity(2000)
        for (n, i) in order.enumerated() where n % 2 == 0 {
            momentumPts.append(Point(x: k[i] / 3, y: pd[i] / pdPeak))
        }

        // 势阱基态（脚本归一口径：y = ψ²/ψ_max）
        let (wx, wd, dxW, dpW) = UncertaintyMath.wellMoments(L: L)
        var wPeak = 0.0
        for v in wd where v > wPeak { wPeak = v }
        let wellPeak = wPeak.squareRoot()
        let wellPts = stride(from: 0, to: wx.count, by: 8).map {
            Point(x: wx[$0] / L, y: wd[$0] / wellPeak)
        }

        let prodG = dxG * dpG
        let prodW = dxW * dpW
        let ana = UncertaintyMath.wellAnalytic(L: L)

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "x（左）/ p/3ħ（右，任意单位）"),
                        yAxis: .init(label: "归一化密度"),
                        seriesNames: ["|ψ(x)|²", "|φ(p)|²"]),
                    series: [
                        .init(name: "|ψ(x)|²", points: positionPts),
                        .init(name: "|φ(p)|²", points: momentumPts, colorIndex: 1),
                    ])),
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "x / L"),
                        yAxis: .init(label: "归一化 |ψ₁|²"),
                        seriesNames: ["|ψ₁(x)|²"]),
                    series: [.init(name: "|ψ₁(x)|²", points: wellPts)])),
            ],
            summary: [
                .init(id: "gdx", title: "高斯 Δx",
                      value: String(format: "%.4f", dxG), note: "σ₀"),
                .init(id: "gdp", title: "高斯 Δp/ħ",
                      value: String(format: "%.4f", dpG), note: "1/(2σ₀)"),
                .init(id: "gprod", title: "高斯 Δx·Δp",
                      value: String(format: "%.4f ħ", prodG), note: "下界 ħ/2 精确取等"),
                .init(id: "wprod", title: "势阱 Δx·Δp",
                      value: String(format: "%.4f ħ", prodW),
                      note: String(format: "解析 %.4f ħ", ana.dx * ana.dpOverHbar)),
            ],
            theory: TheoryCard(
                title: "不确定性关系",
                formulas: [
                    "Δx·Δp ≥ ħ/2（海森伯不等式）",
                    "高斯波包取等号：Δx = σ₀，Δp = ħ/(2σ₀)",
                    "势阱基态 Δx·Δp = 0.568ħ：空间受限 → 动量弥散",
                    "动量分布 = 波函数的傅里叶变换（FFT 数值验证）",
                ],
                reading: "左图高斯的位置/动量分布同步收窄（乘积不变）；右图势阱越窄动量越弥散——受限与弥散的对价。"))
    }
}
