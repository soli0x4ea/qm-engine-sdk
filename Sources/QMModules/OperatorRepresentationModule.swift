import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/算符表象_高斯波包傅里叶对偶.py）

/// 位置/动量表象傅里叶对偶的纯函数核（ħ = 1 自然单位）。
/// N = 4096（2 幂，vDSP 直连路径），L = 40；trapz 为脚本版 Σy·dx（无端点修正）。
enum OperatorRepresentationMath {

    /// 脚本 trapz：均匀网格 Σy·dx。
    static func trapz(_ y: [Double], dx: Double) -> Double {
        var s = 0.0
        for v in y { s += v }
        return s * dx
    }

    /// 位置表象高斯：ψ = (πσ²)^{−1/4}·e^{−x²/2σ²}，数值归一化。
    static func gaussianState(sigmaX: Double, count: Int = 4096, L: Double = 40.0)
        -> (x: [Double], psi: [Double], dx: Double) {
        let x = Num.linspace(-L / 2, L / 2, count: count)
        let dx = x[1] - x[0]
        let pref = pow(1.0 / (.pi * sigmaX * sigmaX), 0.25)
        var psi = x.map { pref * exp(-$0 * $0 / (2 * sigmaX * sigmaX)) }
        var sq = [Double](repeating: 0, count: count)
        for i in 0..<count { sq[i] = psi[i] * psi[i] }
        let nrm = trapz(sq, dx: dx).squareRoot()
        for i in 0..<count { psi[i] /= nrm }
        return (x, psi, dx)
    }

    /// 动量表象：φ = fftshift(fft(ifftshift(ψ)))·dx/√(2π)，数值归一化后取 |φ|²。
    static func momentumDensity(psi: [Double], dx: Double)
        -> (p: [Double], phi2: [Double]) {
        let n = psi.count
        let f = FFT.fft(FFT.ifftshift(psi))
        var re = FFT.fftshift(f.re)
        var im = FFT.fftshift(f.im)
        let scale = dx / (2 * .pi).squareRoot()
        var phi2 = [Double](repeating: 0, count: n)
        for i in 0..<n {
            re[i] *= scale
            im[i] *= scale
            phi2[i] = re[i] * re[i] + im[i] * im[i]
        }
        let p = FFT.fftshift(FFT.fftfreq(n, d: dx)).map { $0 * 2 * .pi }
        let dp = abs(p[1] - p[0])
        let nrm = trapz(phi2, dx: dp)
        for i in 0..<n { phi2[i] /= nrm }
        return (p, phi2)
    }

    /// 脚本口径矩：mean = trapz(x·ρ)，var = trapz((x−mean)²·ρ)（先减后平方）。
    static func moments(grid: [Double], density: [Double]) -> (mean: Double, std: Double) {
        let dx = abs(grid[1] - grid[0])
        var m0 = 0.0, m1 = 0.0
        for i in 0..<grid.count {
            m0 += density[i]
            m1 += grid[i] * density[i]
        }
        let norm = m0 * dx
        let mean = m1 * dx / norm
        var m2 = 0.0
        for i in 0..<grid.count {
            let d = grid[i] - mean
            m2 += d * d * density[i]
        }
        return (mean, (m2 * dx / norm).squareRoot())
    }
}

// MARK: - 模块

/// 笔记 12《算符与表象》：同一量子态在位置/动量两个表象中的傅里叶对偶。
/// 秒级档：4096 点 FFT（vDSP 2 幂路径），σ 拖动即重算。
struct OperatorRepresentationModule: SimModule {

    let meta = ModuleMeta(
        id: "算符表象_高斯波包傅里叶对偶", title: "算符表象 · 傅里叶对偶",
        subtitle: "ψ(x) 与 φ(p) 是同一态的不同基展开——酉等价",
        category: .formalTheory, noteNumber: 12, tier: .seconds, difficulty: .basic,
        keywords: ["表象", "算符", "傅里叶变换", "representation", "位置表象",
                   "动量表象", "酉变换", "对偶", "最小不确定态"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "sigmaX", title: "波包宽度", symbol: "σₓ", unit: "arb.",
                               range: 0.5...4.0, defaultValue: 2.0,
                               scale: .linear, decimalPlaces: 2)),
            .constant(ConstantSpec(key: "hbar", title: "ħ",
                                   note: "自然单位 ħ = 1（无量纲动量）")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "位置表象 |ψ(x)|²",
                xAxis: .init(label: "位置 x"), yAxis: .init(label: "概率密度"),
                seriesNames: ["|ψ(x)|²"])),
            .lineSeries(LineSeriesSpec(
                title: "动量表象 |φ(p)|²",
                xAxis: .init(label: "动量 p (ħ)"), yAxis: .init(label: "概率密度"),
                seriesNames: ["|φ(p)|²"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let sigmaX = input.slider("sigmaX")

        let (x, psi, dx) = OperatorRepresentationMath.gaussianState(sigmaX: sigmaX)
        let (p, phi2) = OperatorRepresentationMath.momentumDensity(psi: psi, dx: dx)

        var psi2 = [Double](repeating: 0, count: x.count)
        var psi2Max = 0.0
        for i in 0..<x.count {
            psi2[i] = psi[i] * psi[i]
            if psi2[i] > psi2Max { psi2Max = psi2[i] }
        }
        var phi2Max = 0.0
        for v in phi2 where v > phi2Max { phi2Max = v }

        // 位置图：全网格抽稀（stride 8 → 512 点）
        let positionPts = Num.strided(x, psi2.map { $0 / psi2Max }, stride: 8)

        // 动量图：信号区 |p| ≤ 4 全分辨率（频段 ±321 绝大部分为 0，截窗显示）
        var momentumPts: [Point] = []
        momentumPts.reserveCapacity(64)
        for i in 0..<p.count where abs(p[i]) <= 4.0 {
            momentumPts.append(Point(x: p[i], y: phi2[i] / phi2Max))
        }

        let (_, dxStd) = OperatorRepresentationMath.moments(grid: x, density: psi2)
        let (_, dpStd) = OperatorRepresentationMath.moments(grid: p, density: phi2)
        let product = dxStd * dpStd

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "位置 x"), yAxis: .init(label: "概率密度"),
                        seriesNames: ["|ψ(x)|²"]),
                    series: [.init(name: "|ψ(x)|²", points: positionPts)])),
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "动量 p (ħ)"), yAxis: .init(label: "概率密度"),
                        seriesNames: ["|φ(p)|²"]),
                    series: [.init(name: "|φ(p)|²", points: momentumPts, colorIndex: 1)])),
            ],
            summary: [
                .init(id: "dx", title: "Δx",
                      value: String(format: "%.4f", dxStd), note: "σₓ/√2"),
                .init(id: "dp", title: "Δp",
                      value: String(format: "%.4f", dpStd), note: "1/(√2σₓ)"),
                .init(id: "prod", title: "Δx·Δp",
                      value: String(format: "%.4f", product), note: "下界 ħ/2 = 0.5"),
                .init(id: "ratio", title: "比值 ΔxΔp/(ħ/2)",
                      value: String(format: "%.4f", product / 0.5), note: "高斯精确取等"),
            ],
            theory: TheoryCard(
                title: "表象与傅里叶对偶",
                formulas: [
                    "⟨x|ψ⟩ = ψ(x)，⟨p|ψ⟩ = φ(p)：同一态在不同基下的展开",
                    "φ(p) = (2πħ)^{−1/2}·∫e^{−ipx/ħ}ψ(x)dx（酉傅里叶变换）",
                    "高斯最小不确定态：Δx = σₓ/√2，Δp = ħ/(√2σₓ)，乘积恰为 ħ/2",
                    "表象变换不改变物理：谱与期望值在两表象中一致",
                ],
                reading: "拖 σₓ 看双图联动——位置图变宽时动量图同步变窄（乘积恒为 ħ/2），这就是对偶。"))
    }
}
