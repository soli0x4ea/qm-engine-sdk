import Foundation
import Accelerate
import EngineKit

// MARK: - 计算核（逐句移植 code/波粒二象性_波长与累积模拟.py）

/// 波粒二象性的纯函数核：模型 A 波长标定曲线 + 模型 B 单电子累积模拟。
enum WaveParticleMath {

    // MARK: 模型 A：德布罗意波长

    /// 非相对论 λ = h/√(2·m_e·e·V)。
    static func lambdaNonrel(_ V: Double, h: Double, me: Double, e: Double) -> Double {
        h / (2 * me * e * V).squareRoot()
    }

    /// 相对论：pc = √(E_k(E_k + 2m_ec²))，λ = h/p。
    static func lambdaRel(_ V: Double, h: Double, me: Double, e: Double, c: Double) -> Double {
        let ek = e * V
        let p = (ek * (ek + 2 * me * c * c)).squareRoot() / c
        return h / p
    }

    /// np.logspace(0, 5, count)：10 的幂等比网格（V 扫描用）。
    static func logspace(_ fromExp: Double, _ toExp: Double, count: Int) -> [Double] {
        Num.linspace(fromExp, toExp, count: count).map { pow(10, $0) }
    }

    // MARK: 模型 B：双缝累积（自然单位 L=d=1, a=0.4）

    /// 理论分布 P₁₂(x) = 4·sinc²(a·x/L)·cos²(π·d·x/L)（归一化）。
    static func theoreticalP12(a: Double = 0.4, d: Double = 1.0, L: Double = 1.0,
                               count: Int = 2001) -> (x: [Double], p: [Double]) {
        let x = Num.linspace(-2.5, 2.5, count: count)
        var p = x.map { xi -> Double in
            let s = xi == 0 ? 1.0 : sin(.pi * a * xi / L) / (.pi * a * xi / L)
            return 4 * s * s * pow(cos(.pi * d * xi / L), 2)
        }
        let dx = x[1] - x[0]
        let norm = p.reduce(0, +) * dx
        p = p.map { $0 / norm }
        return (x, p)
    }

    /// CDF 逆变换抽样：hits = x[searchsorted(cdf, r)]（np.searchsorted 'left' 语义）。
    static func sampleHits(x: [Double], p: [Double], count n: Int,
                           seed: UInt64) -> [Double] {
        let dx = x[1] - x[0]
        var cdf = [Double](repeating: 0, count: p.count)
        var acc = 0.0
        for i in p.indices { acc += p[i] * dx; cdf[i] = acc }
        var rng = SplitMix64(seed: seed)
        return (0..<n).map { _ in
            let r = rng.nextDouble()
            var lo = 0, hi = cdf.count
            while lo < hi {           // 首个 cdf[i] >= r
                let mid = (lo + hi) / 2
                if cdf[mid] < r { lo = mid + 1 } else { hi = mid }
            }
            return x[min(lo, x.count - 1)]
        }
    }

    /// 显示用抽样：从 P₁₂ 分布抽 min(n, maxPoints) 个命中点。
    ///
    /// 工程权衡（实时档 16 ms 预算）：N=60000 全量抽样在 Debug 下约 19 ms
    /// （RNG 5.6 + sort 3.7 + walk 4.0 + 图表构建），超预算。显示层抽 2000 个
    /// i.i.d. 样本与「60000 个再步进取 2000」统计等价（散点渲染与次序无关），
    /// 而 N=60000 的统计一致性断言（卡方 |z|<3）在单测中按全量验证。
    static func displayHits(x: [Double], p: [Double], count n: Int,
                            seed: UInt64, maxPoints: Int = 2000) -> [Point] {
        let m = min(n, maxPoints)
        let dx = x[1] - x[0]
        var cdf = [Double](repeating: 0, count: p.count)
        var acc = 0.0
        for i in p.indices { acc += p[i] * dx; cdf[i] = acc }
        var rng = SplitMix64(seed: seed)
        var rs = [Double](repeating: 0, count: m)
        for i in 0..<m { rs[i] = rng.nextDouble() }
        vDSP.sort(&rs, sortOrder: .ascending)
        var out: [Point] = []
        out.reserveCapacity(m)
        var j = 0
        for (i, r) in rs.enumerated() {
            while j < cdf.count - 1 && cdf[j] < r { j += 1 }
            out.append(Point(x: x[j], y: 1.0 + Double(i % 3) * 0.02))
        }
        return out
    }

    /// 卡方检验（与脚本同口径：100 bin、同权期望、低计数合并）。
    /// 返回 (chi2, dof, z)；z = (χ²−dof)/√(2·dof)，|z| < 3 ≈ p > 1e-3。
    static func chiSquare(hits: [Double], x: [Double], p: [Double],
                          bins: Int = 100) -> (chi2: Double, dof: Int, z: Double) {
        let lo = -2.5, hi = 2.5, w = (hi - lo) / Double(bins)
        var hist = [Double](repeating: 0, count: bins)
        for v in hits {
            let b = min(bins - 1, max(0, Int((v - lo) / w)))
            hist[b] += 1
        }
        var expected = [Double](repeating: 0, count: bins)
        for (xi, pi) in zip(x, p) {
            let b = min(bins - 1, max(0, Int((xi - lo) / w)))
            expected[b] += pi
        }
        let n = Double(hits.count)
        let esum = expected.reduce(0, +)
        expected = expected.map { $0 / esum * n }
        // 合并 expected < 5 的 bin
        var hOk = [Double](), eOk = [Double]()
        var hTail = 0.0, eTail = 0.0
        for i in 0..<bins {
            if expected[i] >= 5 { hOk.append(hist[i]); eOk.append(expected[i]) }
            else { hTail += hist[i]; eTail += expected[i] }
        }
        hOk.append(hTail); eOk.append(eTail)
        var chi2 = 0.0
        for i in hOk.indices { chi2 += pow(hOk[i] - eOk[i], 2) / eOk[i] }
        let dof = eOk.count - 1
        let z = (chi2 - Double(dof)) / (2 * Double(dof)).squareRoot()
        return (chi2, dof, z)
    }
}

/// 定种子均匀 RNG（SplitMix64）——替代 np.random.default_rng(20260901)。
/// 断言为统计检验（|z| < 3），不逐位对拍 PCG64 流；种子固定保证可复现。
struct SplitMix64 {
    private var state: UInt64
    init(seed: UInt64) { state = seed }
    mutating func next() -> UInt64 {
        state &+= 0x9E3779B97F4A7C15
        var z = state
        z = (z ^ (z >> 30)) &* 0xBF58476D1CE4E5B9
        z = (z ^ (z >> 27)) &* 0x94D049BB133111EB
        return z ^ (z >> 31)
    }
    /// [0, 1) 均匀双精度（高 53 位）。
    mutating func nextDouble() -> Double {
        Double(next() >> 11) * (1.0 / 9007199254740992.0)
    }
}

// MARK: - 模块

/// 笔记 02《波粒二象性与德布罗意关系》：波长标定曲线 + 单电子累积（Tonomura 1989 序列）。
/// 实时档：模型 A 解析曲线；模型 B 定种子 CDF 抽样（N 四档切换）。
struct WaveParticleModule: SimModule {

    let meta = ModuleMeta(
        id: "波粒二象性_波长与累积模拟", title: "波粒二象性 · 波长与累积模拟",
        subtitle: "德布罗意波长标定 + 单电子双缝累积（Tonomura 序列）",
        category: .oldQuantum, noteNumber: 2, tier: .realtime, difficulty: .basic,
        keywords: ["德布罗意", "de Broglie", "波长", "累积", "Tonomura",
                   "Davisson", "Germer", "双缝", "统计", "波粒二象性"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "V", title: "加速电压", symbol: "V", unit: "V",
                               range: 1...100000, defaultValue: 54,
                               scale: .log, decimalPlaces: 1)),
            .discrete(DiscreteSpec(key: "N", title: "累积电子数",
                                   options: [8, 270, 2000, 60000].map {
                                       .init(id: "\($0)", title: "N = \($0)",
                                             subtitle: $0 == 60000 ? "Tonomura 终图" : nil)
                                   }, defaultOptionID: "2000")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "模型 A：电子德布罗意波长 vs 加速电压（log-log）",
                xAxis: .init(label: "V (Volt)", scale: .log),
                yAxis: .init(label: "λ (nm)", scale: .log),
                seriesNames: ["非相对论", "相对论"])),
            .scatter(ScatterSpec(
                title: "模型 B：单电子累积（红点击中位置 + 理论包络）",
                xAxis: .init(label: "屏位置 x（自然单位）"),
                yAxis: .init(label: "击中"),
                seriesNames: ["击中位置"])),
            .lineSeries(LineSeriesSpec(
                title: "理论分布 P₁₂(x)（归一化峰值）",
                xAxis: .init(label: "屏位置 x（自然单位）"),
                yAxis: .init(label: "P₁₂ / max"),
                seriesNames: ["理论分布"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let V = input.slider("V")
        let n = Int(input.discrete("N")) ?? 2000
        let h = try constants.value("h")
        let me = try constants.value("m_e")
        let e = try constants.value("e")
        let c = try constants.value("c")

        // 模型 A：log-log 双曲线
        let Vs = WaveParticleMath.logspace(0, 5, count: 500)
        let nr = Vs.map { WaveParticleMath.lambdaNonrel($0, h: h, me: me, e: e) * 1e9 }
        let rl = Vs.map { WaveParticleMath.lambdaRel($0, h: h, me: me, e: e, c: c) * 1e9 }
        let lamNr = WaveParticleMath.lambdaNonrel(V, h: h, me: me, e: e) * 1e9
        let lamRl = WaveParticleMath.lambdaRel(V, h: h, me: me, e: e, c: c) * 1e9
        let corr = (lamRl / lamNr - 1) * 100

        // 模型 B：定种子抽样（种子 20260901，与脚本一致；显示侧抽样 ≤2000 点）
        let (x, p12) = WaveParticleMath.theoreticalP12()
        let pmax = p12.max() ?? 1
        let shownHits = WaveParticleMath.displayHits(x: x, p: p12, count: n,
                                                     seed: 20260901)

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "V (Volt)", scale: .log),
                        yAxis: .init(label: "λ (nm)", scale: .log),
                        seriesNames: ["非相对论", "相对论"]),
                    series: [
                        .init(name: "非相对论", points: Num.strided(Vs, nr, stride: 2)),
                        .init(name: "相对论", points: Num.strided(Vs, rl, stride: 2)),
                    ],
                    referenceLines: [ReferenceLine(
                        id: "vsel", label: String(format: "V = %.1f V", V),
                        axis: .x, value: V)])),
                .scatter(ScatterData(
                    spec: ScatterSpec(
                        xAxis: .init(label: "屏位置 x（自然单位）"),
                        yAxis: .init(label: "击中"),
                        seriesNames: ["击中位置"]),
                    series: [.init(name: "击中位置", points: shownHits)])),
                .lineSeries(LineSeriesData(
                    spec: LineSeriesSpec(
                        xAxis: .init(label: "屏位置 x（自然单位）"),
                        yAxis: .init(label: "P₁₂ / max"),
                        seriesNames: ["理论分布"]),
                    series: [.init(name: "理论分布",
                                   points: Num.strided(x, p12.map { $0 / pmax }, stride: 2))])),
            ],
            summary: [
                .init(id: "lnr", title: "λ（非相对论）",
                      value: String(format: "%.4f nm", lamNr), note: "h/√(2m_e·eV)"),
                .init(id: "lrl", title: "λ（相对论）",
                      value: String(format: "%.4f nm", lamRl), note: "pc=√(E_k(E_k+2mc²))"),
                .init(id: "corr", title: "相对论修正",
                      value: String(format: "%+.2f%%", corr), note: "50 kV 约 −2.4%"),
                .init(id: "hits", title: "累积电子", value: "\(n)",
                      note: "种子 20260901；显示抽样 ≤2000 点"),
            ],
            theory: TheoryCard(
                title: "德布罗意关系",
                formulas: [
                    "λ = h/p；非相对论 λ = h/√(2m_e·eV)",
                    "相对论：pc = √(E_k(E_k + 2m_ec²))",
                    "单电子累积：每个电子落点随机，N → ∞ 显形干涉分布",
                    "P₁₂ = 4·sinc²(πax/L)·cos²(πdx/L)",
                ],
                reading: "模型 A 看波长随电压下降及相对论偏离；模型 B 切换 N 观察条纹如何从随机点云中「长」出来。"))
    }
}
