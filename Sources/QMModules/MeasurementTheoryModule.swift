import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/测量理论_指针耦合与芝诺.py）

/// von Neumann 指针耦合与量子芝诺的纯函数核。
/// 指针平移按脚本口径复刻 np.interp（分段线性插值 + 边缘 clamp），
/// 保证与 Python 导出的 fixtures 曲线逐点对拍。
enum MeasurementTheoryMath {

    /// 复刻 np.interp(xq, x − λ, φ₀v)：键为左移 λ 的均匀网格、值取 φ₀ 数组。
    /// 数值上即 φ₀(xq + λ) 的分段线性插值近似；两侧超界 clamp 到首末值。
    @inline(__always)
    static func interpShiftedPhi(_ xq: Double, grid: [Double], phi0v: [Double],
                                 lambda: Double) -> Double {
        let n = grid.count
        let dx = grid[1] - grid[0]
        let kFirst = grid[0] - lambda
        let kLast = grid[n - 1] - lambda
        if xq <= kFirst { return phi0v[0] }
        if xq >= kLast { return phi0v[n - 1] }
        let t = (xq - kFirst) / dx
        let j = Int(t)
        if j >= n - 1 { return phi0v[n - 1] }
        let frac = t - Double(j)
        return phi0v[j] * (1 - frac) + phi0v[j + 1] * frac
    }

    /// 脚本 trapz：均匀网格 dx·(Σy − (y₀+y_{n−1})/2)。
    static func trapz(_ y: [Double], dx: Double) -> Double {
        var s = 0.0
        for v in y { s += v }
        return dx * (s - 0.5 * (y[0] + y[y.count - 1]))
    }

    /// 指针初态网格：x ∈ [−10σ, 10σ] × 4000，φ₀ = exp(−x²/2σ²)/(√π·σ)^½。
    static func pointerGrid(sigma: Double, count: Int = 4000)
        -> (x: [Double], phi0: [Double]) {
        let x = Num.linspace(-10 * sigma, 10 * sigma, count: count)
        let norm = ( .pi.squareRoot() * sigma).squareRoot()
        return (x, x.map { exp(-$0 * $0 / (2 * sigma * sigma)) / norm })
    }

    /// 单个 λ 的联合概率曲线（归一化到自身峰值）与指针重叠积分。
    /// 脚本口径：prob = α²·φ₊² + β²·φ₋²，其中
    /// φ₊ = interp(x, x−λ, φ₀)（≈φ₀(x+λ)，峰移向 −λ），
    /// φ₋ = interp(x, x+λ, φ₀)（≈φ₀(x−λ)，峰移向 +λ）。
    /// 单循环原地累计（trapz 端点修正），避免中间数组分配。
    ///
    /// 收尾包 1 性能专项（bit 级无损）：λ=0（含 −0.0）时 φ₊ 与 φ₋ 是同一
    /// 插值的重复计算——`grid[0] − (−0.0) = grid[0] + 0.0` 与 `grid[0] − 0.0`
    /// bit 级相等，故复用一次插值结果（4λ 曲线族中 λ=0 恒在）。
    static func jointProbability(
        lambda: Double, sigma: Double, alpha2: Double, grid: [Double], phi0v: [Double]
    ) -> (y: [Double], overlap: Double) {
        let n = grid.count
        let dx = grid[1] - grid[0]
        let beta2 = 1 - alpha2
        var prob = [Double](repeating: 0, count: n)
        var sum = 0.0
        var peak = 0.0
        let lambdaZero = lambda == 0
        for i in 0..<n {
            let xq = grid[i]
            let plus = interpShiftedPhi(xq, grid: grid, phi0v: phi0v, lambda: lambda)
            let minus = lambdaZero
                ? plus
                : interpShiftedPhi(xq, grid: grid, phi0v: phi0v, lambda: -lambda)
            sum += plus * minus
            let p = alpha2 * plus * plus + beta2 * minus * minus
            prob[i] = p
            if p > peak { peak = p }
        }
        let overlap = trapzTailCorrect(sum, grid: grid, phi0v: phi0v,
                                       lambda: lambda, dx: dx)
        for i in 0..<n { prob[i] /= peak }
        return (prob, overlap)
    }

    /// trapz 的首尾半权重修正项：dx·(Σ − (y₀+y_{n−1})/2) 中扣除首尾的另一半。
    private static func trapzTailCorrect(
        _ sum: Double, grid: [Double], phi0v: [Double], lambda: Double, dx: Double
    ) -> Double {
        let first = interpShiftedPhi(grid[0], grid: grid, phi0v: phi0v, lambda: lambda)
            * interpShiftedPhi(grid[0], grid: grid, phi0v: phi0v, lambda: -lambda)
        let n = grid.count
        let last = interpShiftedPhi(grid[n - 1], grid: grid, phi0v: phi0v, lambda: lambda)
            * interpShiftedPhi(grid[n - 1], grid: grid, phi0v: phi0v, lambda: -lambda)
        return dx * (sum - 0.5 * (first + last))
    }

    /// 芝诺存活概率 P(N) = |cos(Ω·τ/2)|^{2N}，τ = t_total/N。
    static func zenoSurvival(N: Int, omega: Double, tTotal: Double) -> Double {
        let tau = tTotal / Double(N)
        let c = cos(omega * tau / 2)
        return pow(abs(c), 2 * Double(N))
    }

    /// 芝诺测量次数序列（脚本：1 起倍增到 512）。
    static let zenoCounts = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
}

// MARK: - 模块

/// 笔记 10《测量理论与波函数坍缩》：指针耦合退相干 + 量子芝诺效应。
/// 实时档：指针曲线族（4 个 λ）+ 芝诺存活概率（semilogx）。
struct MeasurementTheoryModule: SimModule {

    private static let lambdas = [0.0, 1.0, 2.0, 4.0]

    let meta = ModuleMeta(
        id: "测量理论_指针耦合与芝诺", title: "测量理论 · 指针耦合与芝诺",
        subtitle: "纠缠使相干性随 λ 指数衰减；频繁测量冻结演化",
        category: .quantumStates, noteNumber: 10, tier: .realtime, difficulty: .basic,
        keywords: ["测量", "坍缩", "指针", "von Neumann", "芝诺", "Zeno",
                   "退相干", "纠缠", "存活概率", "投影"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "alpha2", title: "初态 |0⟩ 权重", symbol: "α²", unit: "",
                               range: 0.1...0.9, defaultValue: 0.6,
                               scale: .linear, decimalPlaces: 2)),
            .slider(SliderSpec(key: "sigma", title: "指针宽度", symbol: "σ", unit: "σ",
                               range: 0.5...2.0, defaultValue: 1.0,
                               scale: .linear, decimalPlaces: 2)),
            .slider(SliderSpec(key: "omega", title: "芝诺驱动强度", symbol: "Ω", unit: "rad/s",
                               range: 0.5...2.0, defaultValue: 1.0,
                               scale: .linear, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "von Neumann 指针耦合：系统-指针纠缠",
                xAxis: .init(label: "指针坐标 x (σ)"),
                yAxis: .init(label: "联合概率 |Ψ|²（归一）"),
                seriesNames: Self.lambdas.map { String(format: "λ = %.1f σ", $0) })),
            .lineSeries(LineSeriesSpec(
                title: "量子芝诺：频繁投影冻结演化",
                xAxis: .init(label: "测量次数 N", scale: .log),
                yAxis: .init(label: "存活概率 P(N)"),
                seriesNames: ["P(N)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let alpha2 = input.slider("alpha2")
        let sigma = input.slider("sigma")
        let omega = input.slider("omega")
        let tTotal = Double.pi

        // 图 1：指针耦合曲线族（x/σ 轴，2000 点抽稀显示）
        let (grid, phi0v) = MeasurementTheoryMath.pointerGrid(sigma: sigma)
        let xsScaled = grid.map { $0 / sigma }
        var pointerSeries: [SeriesPoints] = []
        var overlapNotes: [String] = []
        for lam in Self.lambdas {
            let (y, overlap) = MeasurementTheoryMath.jointProbability(
                lambda: lam, sigma: sigma, alpha2: alpha2, grid: grid, phi0v: phi0v)
            let name = String(format: "λ = %.1f σ", lam)
            pointerSeries.append(.init(name: name, points: Num.strided(xsScaled, y, stride: 2)))
            overlapNotes.append(String(format: "λ=%.0f: %.4f", lam, overlap))
        }

        // 图 2：芝诺存活概率（semilogx）
        let Ns = MeasurementTheoryMath.zenoCounts
        let Ps = Ns.map { MeasurementTheoryMath.zenoSurvival(N: $0, omega: omega, tTotal: tTotal) }

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "指针坐标 x (σ)"),
                        yAxis: .init(label: "联合概率 |Ψ|²（归一）"),
                        seriesNames: Self.lambdas.map { String(format: "λ = %.1f σ", $0) }),
                    series: pointerSeries)),
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "测量次数 N", scale: .log),
                        yAxis: .init(label: "存活概率 P(N)"),
                        seriesNames: ["P(N)"]),
                    series: [.init(name: "P(N)",
                                   points: zip(Ns, Ps).map { Point(x: Double($0.0), y: $0.1) })],
                    referenceLines: [ReferenceLine(id: "freeze", label: "芝诺极限 P → 1",
                                                   axis: .y, value: 1.0)])),
            ],
            summary: [
                .init(id: "ov", title: "指针重叠积分",
                      value: overlapNotes.joined(separator: "，"),
                      note: "∫φ₊φ₋dx ~ e^(−λ²/σ²)"),
                .init(id: "coherence", title: "相干性衰减",
                      value: "e^(−λ²/2σ²)", note: "非对角元 ~ 指数衰减"),
                .init(id: "zeno1", title: "P(N=1)",
                      value: String(format: "%.2e", Ps[0]), note: "单次半周期测量"),
                .init(id: "zeno512", title: "P(N=512)",
                      value: String(format: "%.4f", Ps.last ?? 0), note: "频繁测量 → 冻结"),
            ],
            theory: TheoryCard(
                title: "测量理论两支柱",
                formulas: [
                    "指针耦合：|ψ⟩|φ₀⟩ → α|0⟩|φ₀(x−λ)⟩ + β|1⟩|φ₀(x+λ)⟩",
                    "非对角元 ~ exp(−λ²/2σ²)：指针可分辨 → 相干性指数衰减",
                    "芝诺存活概率 P(N) = |cos(Ω·τ/2)|^(2N)，τ = t/N",
                    "N → ∞ 时 P → 1：演化被频繁投影冻结",
                ],
                reading: "左图 λ 增大两峰分离、重叠积分塌向 0（测量完成）；右图 N 越大存活率越高——盯着看的量子不衰变。"))
    }
}
