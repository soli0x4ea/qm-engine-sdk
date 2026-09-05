import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子谐振子_有限差分与相干态.py，自然单位 ħ=m=ω=1）

/// 量子谐振子有限差分纯函数核：三对角哈密顿全谱对角化（Accelerate dsyevr float64）+
/// Hermite 解析波函数 + 相干态叠加演化。
enum OscillatorMath {

    /// 差分网格与三对角哈密顿元素：T = −½d²/dx²，V = ½x²。
    static func grid(N: Int, L: Double) -> (x: [Double], dx: Double) {
        let x = Num.linspace(-L, L, count: N)
        return (x, x[1] - x[0])
    }

    /// 行主序扁平三对角矩阵。
    static func tridiagonalFlat(N: Int, L: Double) -> (a: [Double], x: [Double]) {
        let (x, dx) = grid(N: N, L: L)
        let dx2 = dx * dx
        var a = [Double](repeating: 0, count: N * N)
        for i in 0..<N {
            a[i * N + i] = 1.0 / dx2 + 0.5 * x[i] * x[i]
            if i + 1 < N {
                a[i * N + i + 1] = -0.5 / dx2
                a[(i + 1) * N + i] = -0.5 / dx2
            }
        }
        return (a, x)
    }

    /// 全谱对角化（秒级通道）：返回前 nstates 个 (E_n, ψₙ(x_grid))。
    /// eigh 特征向量符号任意——与解析 ψₙ 对齐后返回。
    static func solve(N: Int, L: Double, nstates: Int) -> (x: [Double], energies: [Double], vectors: [[Double]]) {
        precondition(N >= 20 && nstates >= 1 && nstates <= 12)
        let (a, x) = tridiagonalFlat(N: N, L: L)
        let (w, v) = SymEigh.eighSymmetric(a, n: N)
        var energies = [Double]()
        var vectors = [[Double]]()
        for s in 0..<min(nstates, N) {
            energies.append(w[s])
            let raw = (0..<N).map { v[$0 * N + s] }
            vectors.append(aligned(raw, x: x, n: s))
        }
        return (x, energies, vectors)
    }

    /// 数值特征向量符号与解析 ψₙ 对齐（逐点内积取号）。
    static func aligned(_ vec: [Double], x: [Double], n: Int) -> [Double] {
        let analytic = Hermite.normalizedWavefunction(n, xs: x)
        var dot = 0.0
        for i in x.indices { dot += vec[i] * analytic[i] }
        return dot >= 0 ? vec : vec.map { -$0 }
    }

    /// 节点计数（Python count_nodes：阈值 1e-3×峰值之下的噪声不计入）。
    static func nodeCount(_ vec: [Double]) -> Int {
        guard let mx = vec.map({ abs($0) }).max(), mx > 0 else { return 0 }
        let threshold = 1e-3 * mx
        var lastSign = 0
        var flips = 0
        for value in vec {
            guard abs(value) > threshold else { continue }
            let sign = value > 0 ? 1 : -1
            if lastSign != 0 && sign != lastSign { flips += 1 }
            lastSign = sign
        }
        return flips
    }

    /// 线性插值到目标网格（网格单调递增）。
    static func interpolate(xs: [Double], ys: [Double], onto grid: [Double]) -> [Double] {
        precondition(xs.count == ys.count && xs.count > 1)
        var out = [Double](repeating: 0, count: grid.count)
        var j = 0
        for (i, g) in grid.enumerated() {
            while j + 2 < xs.count && xs[j + 1] < g { j += 1 }
            let x0 = xs[j], x1 = xs[min(j + 1, xs.count - 1)]
            let y0 = ys[j], y1 = ys[min(j + 1, ys.count - 1)]
            if x1 == x0 {
                out[i] = y0
            } else {
                let t = min(1, max(0, (g - x0) / (x1 - x0)))
                out[i] = y0 + t * (y1 - y0)
            }
        }
        return out
    }

    /// 相干态叠加 |ψ(x,t)|²：ψ = Σₙ cₙψₙ(x)，cₙ = e^(−α²/2)αⁿe^(−it(n+½))/√(n!)。
    /// Nmax = 44 与 Python 一致（α ≤ 3 时尾部已收敛）；α = 0 退化为基态。
    ///
    /// 收尾包 1 性能专项重组（bit 级无损）：系数 amp/cr/ci、归一化系数 pre、
    /// Hermite 三项递推均改为每 (n) / 每 (x) 只算一次——旧实现每点每 n 各自
    /// 从 H₀ 递推并重算 exp（44 阶 × 2400 点实测 176 ms → 重组后 ~5 ms）。
    /// 浮点语义保持：每点按 n 升序累加 Σ cₙψₙ；ψₙ = (pre·H)·exp 乘法结合序、
    /// H 递推中间值与旧逐点 `Hermite.value(n, x)` 完全一致——逐元素 bit 级相等
    /// （一致性由 HarmonicOscillatorModuleTests 的 naive 参考实现锁死）。
    static func coherentDensity(alpha: Double, t: Double, xs: [Double], nmax: Int = 44) -> [Double] {
        if alpha <= 0 {
            let psi0 = Hermite.normalizedWavefunction(0, xs: xs)
            return psi0.map { $0 * $0 }
        }
        // 每 n 一次：lnAmp → amp、cr、ci（与旧实现同公式同序）
        let lnAlpha = log(alpha)
        var crs = [Double](repeating: 0, count: nmax)
        var cis = [Double](repeating: 0, count: nmax)
        for n in 0..<nmax {
            let lnAmp = -alpha * alpha / 2 + Double(n) * lnAlpha - 0.5 * Num.lnFactorial(n)
            let amp = exp(lnAmp)
            let phase = -t * (Double(n) + 0.5)
            crs[n] = amp * cos(phase)
            cis[n] = amp * sin(phase)
        }
        // 每 n 一次：归一化系数 pre（与 Hermite.normalizedWavefunction 同公式）
        var pres = [Double](repeating: 0, count: nmax)
        for n in 0..<nmax {
            var lnCoef = Double(n) * log(2.0)
            if n > 1 { for k in 2...n { lnCoef += log(Double(k)) } }
            pres[n] = exp(-0.5 * lnCoef - 0.25 * log(.pi))
        }
        // 逐点：一次递推出 H₀…H_{nmax−1}（中间值与逐点 value(n, x) bit 级一致），
        // 一次 exp(−x²/2)，按 n 升序累加 Σ cₙψₙ（与旧 re[i] += cr·ψ[i] 顺序一致）
        var re = [Double](repeating: 0, count: xs.count)
        var im = [Double](repeating: 0, count: xs.count)
        for i in xs.indices {
            let x = xs[i]
            let gauss = exp(-x * x / 2)
            var hm2 = 1.0            // H₀
            var hm1 = 2 * x          // H₁
            var sumR = 0.0
            var sumI = 0.0
            for n in 0..<nmax {
                let h: Double
                switch n {
                case 0: h = hm2
                case 1: h = hm1
                default:
                    // 与 Hermite.value 递推体一致：current = 2x·Hₙ₋₁ − 2(n−1)·Hₙ₋₂
                    let current = 2 * x * hm1 - 2 * Double(n - 1) * hm2
                    hm2 = hm1
                    hm1 = current
                    h = current
                }
                let psi = pres[n] * h * gauss
                sumR += crs[n] * psi
                sumI += cis[n] * psi
            }
            re[i] = sumR
            im[i] = sumI
        }
        return zip(re, im).map { $0 * $0 + $1 * $1 }
    }

    /// 概率密度的（均值, 方差）——梯形积分。
    static func meanAndVariance(ofRho rho: [Double], xs: [Double]) -> (mean: Double, variance: Double) {
        let z = Num.trapezoid(xs, rho)
        var mx = [Double](repeating: 0, count: xs.count)
        var dx2 = [Double](repeating: 0, count: xs.count)
        for i in xs.indices {
            mx[i] = xs[i] * rho[i]
        }
        let mean = Num.trapezoid(xs, mx) / z
        for i in xs.indices {
            let d = xs[i] - mean
            dx2[i] = d * d * rho[i]
        }
        return (mean, Num.trapezoid(xs, dx2) / z)
    }
}

// MARK: - 模块

/// 笔记 15《量子谐振子》：有限差分本征值问题（数值 vs 解析 n+½）+ 相干态不扩散演化。
/// 秒级档：N×N 三对角全谱对角化走 Accelerate dsyevr float64 通道
/// （W3 前为 MLX eigh：N=800 ≈ 1.9 s；dsyevr 实测见 SymEighTests，包体裁决后替换）。
struct HarmonicOscillatorModule: SimModule {

    /// 相干态展示的四个相位（与 Python times 一致）。
    private static let phases: [(id: String, label: String, t: Double)] = [
        ("p0", "t = 0", 0),
        ("p1", "t = π/2", .pi / 2),
        ("p2", "t = π", .pi),
        ("p3", "t = 3π/2", 3 * .pi / 2),
    ]

    let meta = ModuleMeta(
        id: "量子谐振子_有限差分与相干态", title: "量子谐振子 · 有限差分与相干态",
        subtitle: "三对角哈密顿 eigh 全谱 vs 解析 E_n = n+½；相干态质心振荡不扩散",
        category: .schrodinger1D, noteNumber: 15, tier: .seconds, difficulty: .basic,
        keywords: ["谐振子", "有限差分", "本征值", "对角化", "eigh", "Hermite", "厄米多项式",
                   "相干态", "波包", "不扩散", "oscillator", "coherent state"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "N", title: "差分网格点数", symbol: "N", unit: "点",
                               range: 200...800, defaultValue: 500,
                               step: 100, decimalPlaces: 0)),
            .slider(SliderSpec(key: "L", title: "计算箱半宽", symbol: "L", unit: "√(ħ/mω)",
                               range: 8...18, defaultValue: 12,
                               step: 0.5, decimalPlaces: 1)),
            .slider(SliderSpec(key: "nstates", title: "本征态数", symbol: "n", unit: "态",
                               range: 4...10, defaultValue: 6,
                               step: 1, decimalPlaces: 0)),
            .slider(SliderSpec(key: "alpha", title: "相干态振幅", symbol: "α", unit: "",
                               range: 0...3, defaultValue: 2,
                               step: 0.1, decimalPlaces: 1)),
            .multiCompare(MultiCompareSpec(
                key: "phases", title: "演化时刻",
                candidates: Self.phases.map { .init(id: $0.id, label: $0.label, value: $0.t) },
                defaultSelectionIDs: ["p0", "p1", "p2", "p3"], maxSelection: 4)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .bars(BarSpec(
                title: "本征能量 E_n：数值 vs 解析 n+½",
                xAxis: .init(label: "量子数 n"),
                yAxis: .init(label: "E (ħω)"),
                seriesNames: ["数值", "解析"])),
            .lineSeries(LineSeriesSpec(
                title: "波函数 ψn(x)：数值 vs 解析（n=0…2）",
                xAxis: .init(label: "x √(mω/ħ)"),
                yAxis: .init(label: "ψ(x)"),
                seriesNames: ["n=0 数值", "n=0 解析", "n=1 数值", "n=1 解析", "n=2 数值", "n=2 解析"])),
            .lineSeries(LineSeriesSpec(
                title: "相干态 |ψ(x,t)|²：质心振荡 · 宽度恒定",
                xAxis: .init(label: "x √(mω/ħ)"),
                yAxis: .init(label: "|ψ(x,t)|²"),
                seriesNames: Self.phases.map(\.label))),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let N = Int(input.slider("N"))
        let L = input.slider("L")
        let nstates = Int(input.slider("nstates"))
        let alpha = input.slider("alpha")
        let selectedPhases = input.multiCompare("phases")

        try Task.checkCancellation()
        let (x, energies, vectors) = OscillatorMath.solve(N: N, L: L, nstates: nstates)

        // 图 1：本征能量数值 vs 解析（分组柱）
        var energyBars: [BarItem] = []
        var maxRelErr = 0.0
        var nodesAll = 0
        for n in 0..<energies.count {
            let analytic = Double(n) + 0.5
            let rel = abs(energies[n] - analytic) / analytic
            maxRelErr = max(maxRelErr, rel)
            nodesAll += OscillatorMath.nodeCount(vectors[n]) == n ? 1 : 0
            energyBars.append(BarItem(label: "n=\(n)", value: energies[n], series: "数值"))
            energyBars.append(BarItem(label: "n=\(n)", value: analytic, series: "解析"))
        }

        // 图 2：波函数数值 vs 解析（zoom 网格，数值插值后重归一化——与 Python 同法）
        try Task.checkCancellation()
        let zoom = Num.linspace(-5, 5, count: 400)
        var waveSeries: [SeriesPoints] = []
        for n in 0..<min(3, vectors.count) {
            var num = OscillatorMath.interpolate(xs: x, ys: vectors[n], onto: zoom)
            let norm = Num.trapezoid(zoom, num.map { $0 * $0 })
            if norm > 0 {
                let inv = 1.0 / sqrt(norm)
                num = num.map { $0 * inv }
            }
            let ana = Hermite.normalizedWavefunction(n, xs: zoom)
            waveSeries.append(.init(name: "n=\(n) 数值", points: Num.strided(zoom, num, stride: 2),
                                    colorIndex: n * 2))
            waveSeries.append(.init(name: "n=\(n) 解析", points: Num.strided(zoom, ana, stride: 2),
                                    colorIndex: n * 2 + 1))
        }

        // 图 3：相干态 |ψ(x,t)|²（2400 网格归一化后抽稀显示）
        try Task.checkCancellation()
        let xg = Num.linspace(-9, 9, count: 2400)
        var phaseSeries: [SeriesPoints] = []
        var varAtT0: Double?
        var centroidAtT0: Double?
        for p in Self.phases where selectedPhases.contains(p.id) {
            let rho = OscillatorMath.coherentDensity(alpha: alpha, t: p.t, xs: xg)
            let norm = Num.trapezoid(xg, rho)
            let rhoN = rho.map { $0 / norm }
            let (mean, variance) = OscillatorMath.meanAndVariance(ofRho: rhoN, xs: xg)
            if p.id == "p0" { varAtT0 = variance; centroidAtT0 = mean }
            phaseSeries.append(.init(name: p.label, points: Num.strided(xg, rhoN, stride: 6)))
        }

        var summary: [SummaryItem] = [
            .init(id: "e0", title: "E₀ 数值", value: String(format: "%.8f", energies.first ?? 0),
                  note: String(format: "解析 0.5，相对误差 %.1e",
                               abs((energies.first ?? 0) - 0.5) / 0.5)),
            .init(id: "maxerr", title: "最大相对误差", value: String(format: "%.1e", maxRelErr),
                  note: "前 \(energies.count) 态 vs E_n = n+½"),
            .init(id: "nodes", title: "节点校验",
                  value: "\(nodesAll)/\(energies.count) 态通过",
                  note: "第 n 本征态恰有 n 个节点"),
        ]
        if let v = varAtT0 {
            summary.append(.init(id: "var0", title: "Var(t=0)", value: String(format: "%.6f", v),
                                 note: "解析 0.5（基态宽度，不随 t 变化）"))
        }
        if let c = centroidAtT0 {
            summary.append(.init(id: "xc0", title: "⟨x⟩(t=0)", value: String(format: "%.3f", c),
                                 note: "√2·Re(αe^(−it)) = √2·α"))
            // 注：Python 源脚本标签 x_c=2cos(t) 系系数笔误，正确质心为 √2·α·cos(t)
        }

        let barSpec = BarSpec(
            xAxis: .init(label: "量子数 n"),
            yAxis: .init(label: "E (ħω)"),
            seriesNames: ["数值", "解析"])
        let waveSpec = LineSeriesSpec(
            xAxis: .init(label: "x √(mω/ħ)"),
            yAxis: .init(label: "ψ(x)"),
            seriesNames: ["n=0 数值", "n=0 解析", "n=1 数值", "n=1 解析", "n=2 数值", "n=2 解析"])
        let phaseSpec = LineSeriesSpec(
            xAxis: .init(label: "x √(mω/ħ)"),
            yAxis: .init(label: "|ψ(x,t)|²"),
            seriesNames: Self.phases.map(\.label))

        return SimResult(
            charts: [
                .bars(BarData(spec: barSpec, bars: energyBars)),
                .lineSeries(LineSeriesData(spec: waveSpec, series: waveSeries)),
                .lineSeries(LineSeriesData(spec: phaseSpec, series: phaseSeries)),
            ],
            summary: summary,
            theory: TheoryCard(
                title: "谐振子：差分本征值与相干态",
                formulas: [
                    "H = −½d²/dx² + ½x²（自然单位 ħ = m = ω = 1）",
                    "解析解：E_n = n + ½；ψₙ ∝ Hₙ(x)e^(−x²/2)",
                    "相干态：|α⟩ = e^(−|α|²/2)Σ αⁿ/√(n!) |n⟩",
                    "演化：⟨x⟩(t) = √2·Re(αe^(−it))；Var 恒为 0.5",
                ],
                reading: "数值能级与解析 n+½ 的偏差随网格加密按 dx² 收敛；相干态质心沿经典轨迹振荡而波包宽度不变——这是谐振子谱等间距的标志性结果。"))
    }
}
