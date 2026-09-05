import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/一维势场_方势垒与有限深势阱.py）

/// 方势垒透射 + 有限深方势阱束缚态的纯函数核。
///
/// Python 口径：
/// - 势垒：T(E) = 1/(1 + V₀²sinh²(κa)/(4E(V₀−E)))，κ = √(2m(V₀−E))/ħ；
///   厚垒极限 T ≈ 16E(V₀−E)/V₀²·e^(−2κa)。
/// - 势阱：200000 点 k 网格扫描偶/奇超越方程符号变化 → 每括号二分 80 次迭代。
enum BarrierWellMath {

    /// k 网格点数与二分迭代次数（脚本固定）。
    static let scanPoints = 200_000
    static let bisectionIterations = 80

    // MARK: 方势垒（E < V₀ 隧穿支；脚本 transmission 的 NaN 域即 E ≤ 0 / E ≥ V₀）

    /// 矩形势垒精确透射系数（E < V₀）。E 越域返回 NaN（与脚本一致）。
    static func transmission(E: Double, v0: Double, width a: Double,
                             mass: Double, hbar: Double) -> Double {
        guard E > 0, E < v0 else { return .nan }
        let kap = sqrt(2.0 * mass * (v0 - E)) / hbar
        let sinh2 = sinh(kap * a) * sinh(kap * a)
        return 1.0 / (1.0 + v0 * v0 * sinh2 / (4.0 * E * (v0 - E)))
    }

    /// 厚垒极限 T ≈ 16E(V₀−E)/V₀²·e^(−2κa)（脚本 thick_limit）。
    static func thickLimit(E: Double, v0: Double, width a: Double,
                           mass: Double, hbar: Double) -> Double {
        let kap = sqrt(2.0 * mass * (v0 - E)) / hbar
        return 16.0 * E * (v0 - E) / (v0 * v0) * exp(-2.0 * kap * a)
    }

    /// E > V₀ 越垒支（物理律验证扩展：脚本未覆盖；T = 1/(1 + V₀²sin²(kL)/(4E(E−V₀)))）。
    /// 共振条件 kL = nπ 处 T = 1。
    static func transmissionAbove(E: Double, v0: Double, width a: Double,
                                  mass: Double, hbar: Double) -> Double {
        guard E > v0 else { return .nan }
        let k = sqrt(2.0 * mass * E) / hbar
        let s = sin(k * a)
        return 1.0 / (1.0 + v0 * v0 * s * s / (4.0 * E * (E - v0)))
    }

    // MARK: 有限深势阱（−L/2 ≤ x ≤ L/2，深度 V₀）

    /// 偶宇称超越方程 f_e(k) = k·tan(kL/2) − √(k²ₘₐₓ−k²)。
    static func evenF(_ k: Double, kmax: Double, width L: Double) -> Double {
        k * tan(k * L / 2.0) - sqrt(max(0, kmax * kmax - k * k))
    }

    /// 奇宇称超越方程 f_o(k) = −k·cot(kL/2) − √(k²ₘₐₓ−k²)。
    static func oddF(_ k: Double, kmax: Double, width L: Double) -> Double {
        -k / tan(k * L / 2.0) - sqrt(max(0, kmax * kmax - k * k))
    }

    /// 单个符号变化括号内的二分求根（脚本口径：80 次迭代，保 f[i] 符号端）。
    /// - Parameters:
    ///   - lo0, hi0: 括号端点（f(lo0) 与扫描端符号一致）
    ///   - anchorSign: 扫描左端点 f 值符号
    static func bisect(even: Bool, lo0: Double, hi0: Double, anchorSign: Double,
                       kmax: Double, width L: Double) -> Double {
        var lo = lo0, hi = hi0
        for _ in 0..<bisectionIterations {
            let mid = 0.5 * (lo + hi)
            let fmid = even ? evenF(mid, kmax: kmax, width: L)
                            : oddF(mid, kmax: kmax, width: L)
            if (fmid.sign == anchorSign) || fmid == 0 {
                lo = mid
            } else {
                hi = mid
            }
        }
        return 0.5 * (lo + hi)
    }

    /// 扫描 + 二分求根（脚本 well_energies/wavefunction 重构路径共用），
    /// **返回原始根（含 tan 极点伪根）**——与脚本逐句对应。
    static func rawRoots(v0: Double, width L: Double, mass: Double, hbar: Double)
        -> [(k: Double, even: Bool)] {
        let kmax = sqrt(2.0 * mass * v0) / hbar
        let ks = Num.linspace(1.0e-6, kmax * (1 - 1.0e-6), count: scanPoints)
        var out: [(k: Double, even: Bool)] = []
        for even in [true, false] {
            var prevSign = 0.0
            for i in 0..<ks.count {
                let f = even ? evenF(ks[i], kmax: kmax, width: L)
                             : oddF(ks[i], kmax: kmax, width: L)
                let s = f.sign
                if i > 0 && (prevSign == 0 || s != prevSign) {
                    let kroot = bisect(even: even, lo0: ks[i - 1], hi0: ks[i],
                                       anchorSign: prevSign, kmax: kmax, width: L)
                    out.append((kroot, even))
                }
                prevSign = s
            }
        }
        return out
    }

    /// 伪根剔除（**基准脚本疑点，见周报偏差说明**）：
    /// tan/cot 极点（kL/2 = π/2 + mπ）处超越方程符号跳变被脚本误判为根——
    /// 二分收敛到极点而非零点，残差 |f(k_root)| 巨大。真根残差 ~0（机器精度收敛），
    /// 伪根残差 ≥ k·|tan| ≫ 1e-6·k_max，以残差阈值甄别。
    static func isValidRoot(_ k: Double, even: Bool, kmax: Double, width L: Double) -> Bool {
        let f = even ? evenF(k, kmax: kmax, width: L)
                     : oddF(k, kmax: kmax, width: L)
        return abs(f) < 1.0e-6 * kmax
    }

    /// 束缚态：能量（J，负值）升序 + 宇称（按能量序交替 even/odd，与脚本 idx%2 约定一致）。
    /// 已剔除 tan 极点伪根。
    static func boundStates(v0: Double, width L: Double, mass: Double, hbar: Double)
        -> [(E: Double, k: Double, even: Bool)] {
        let kmax = sqrt(2.0 * mass * v0) / hbar
        let roots = rawRoots(v0: v0, width: L, mass: mass, hbar: hbar)
            .filter { isValidRoot($0.k, even: $0.even, kmax: kmax, width: L) }
        return roots.map {
            let e = hbar * hbar * $0.k * $0.k / (2.0 * mass) - v0
            return (e, $0.k, $0.even)
        }.sorted { $0.E < $1.E }
    }

    /// 束缚态个数的解析基准：N = ⌈2z₀/π⌉，z₀ = k_max·L/2（对称有限深方势阱标准结果）。
    static func analyticBoundCount(v0: Double, width L: Double, mass: Double, hbar: Double) -> Int {
        let z0 = sqrt(2.0 * mass * v0) / hbar * L / 2.0
        return Int(ceil(2.0 * z0 / .pi))
    }

    /// 束缚态波函数（脚本 well_wavefunction）：阱内 cos/sin + 阱外指数衰减，网格归一化。
    static func wavefunction(k: Double, v0: Double, width L: Double, even: Bool,
                             xs: [Double], mass: Double, hbar: Double) -> [Double] {
        let kap = sqrt(max(0, 2.0 * mass * v0 / (hbar * hbar) - k * k))
        let half = L / 2.0
        let psi = xs.map { x -> Double in
            let ax = abs(x)
            if ax <= half {
                return even ? cos(k * x) : sin(k * x)
            }
            let edge = even ? cos(k * half) : sin(k * half)
            return edge * exp(-kap * (ax - half))
        }
        // 梯形法归一（按实际区间宽）：脚本在自身均匀 2000 点网格上用 Σψ²·dx，
        // 与梯形法逐位一致（~1e-12），本实现口径不变；但 fixture 导出管线的 x 经
        // 2000 点网格索引重映射而抖动（首间距≠均值间距，黎曼和偏差达 1e-2 量级），
        // 梯形法对任意网格稳健，故取之。
        var integral = 0.0
        for i in xs.indices.dropLast() {
            integral += 0.5 * (psi[i] * psi[i] + psi[i + 1] * psi[i + 1]) * (xs[i + 1] - xs[i])
        }
        let norm = sqrt(integral)
        return psi.map { $0 / norm }
    }
}

private extension Double {
    /// np.sign 语义（−1/0/+1）
    var sign: Double { self > 0 ? 1 : (self < 0 ? -1 : 0) }
}

// MARK: - 模块

/// 笔记 14《一维势场》：矩形势垒精确透射 T(E)/T(a) + 有限深方势阱束缚态
/// （超越方程 200000 点扫描 + 二分 80 次；**剔除 tan 极点伪根，脚本疑点见周报**）。秒级档。
struct BarrierWellModule: SimModule {

    let meta = ModuleMeta(
        id: "一维势场_方势垒与有限深势阱", title: "一维势场 · 方势垒与有限深势阱",
        subtitle: "精确 T(E) 隧穿/共振 + 超越方程二分求根与束缚态谱",
        category: .schrodinger1D, noteNumber: 14, tier: .seconds, difficulty: .basic,
        keywords: ["方势垒", "隧穿", "透射系数", "有限深势阱", "超越方程", "二分法",
                   "束缚态", "共振", "tunneling", "transmission"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "V0b", title: "势垒高度", symbol: "V₀", unit: "eV",
                               range: 0.1...5, defaultValue: 1.0, decimalPlaces: 2)),
            .slider(SliderSpec(key: "ab", title: "势垒宽度", symbol: "a", unit: "nm",
                               range: 0.1...3, defaultValue: 1.0, decimalPlaces: 2)),
            .slider(SliderSpec(key: "V0w", title: "阱深", symbol: "V₀", unit: "eV",
                               range: 0.1...5, defaultValue: 0.5, decimalPlaces: 2)),
            .slider(SliderSpec(key: "Lw", title: "阱宽", symbol: "L", unit: "nm",
                               range: 0.5...5, defaultValue: 2.0, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "矩形势垒：透射系数 T 随入射能量（a 按滑杆）",
                xAxis: .init(label: "E / V₀（无量纲）"),
                yAxis: .init(label: "透射系数 T"),
                seriesNames: ["精确 T(E)"])),
            .lineSeries(LineSeriesSpec(
                title: "透射系数随垒宽（E = 0.3V₀）：精确 vs 厚垒极限",
                xAxis: .init(label: "垒宽 a (nm)"),
                yAxis: .init(label: "透射系数 T", scale: .log),
                seriesNames: ["精确 T(a)", "厚垒极限 e^(−2κa)"])),
            .lineSeries(LineSeriesSpec(
                title: "有限深势阱：最低束缚态波函数（叠加 0.6·(n−1)，阱区阴影见参考线）",
                xAxis: .init(label: "x (nm)"),
                yAxis: .init(label: "ψ(x) + 0.6(n−1) (m⁻¹ᐟ²)"),
                seriesNames: (1...5).map { "n=\($0)" })),
            .lineSeries(LineSeriesSpec(
                title: "有限深势阱：束缚态能级",
                xAxis: .init(label: "（示意横轴）"),
                yAxis: .init(label: "能量 E (eV)"),
                seriesNames: ["束缚态能级"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let hbar = try constants.value("hbar")
        let mElectron = try constants.value("m_e")
        let eCharge = try constants.value("e")

        let v0bEv = input.slider("V0b")
        let aNm = input.slider("ab")
        let v0wEv = input.slider("V0w")
        let lNm = input.slider("Lw")
        let v0b = v0bEv * eCharge
        let aMeters = aNm * 1.0e-9
        let v0w = v0wEv * eCharge
        let lMeters = lNm * 1.0e-9

        let progress = ComputeProgress()
        let payload = try await SecondsChannel.run(progress: progress)
            { () -> (tE: [Point], tA: [Point], tlA: [Point], states: [(E: Double, k: Double, even: Bool)]) in
            progress.update(0.15, phase: "势垒 T(E) 扫描")
            // 图 1：E ∈ [0.02, 0.98]·V₀，400 点（脚本口径）
            var tE: [Point] = []
            for i in 0..<400 {
                let ratio = 0.02 + 0.96 * Double(i) / 399.0
                let t = BarrierWellMath.transmission(
                    E: ratio * v0b, v0: v0b, width: aMeters, mass: mElectron, hbar: hbar)
                tE.append(Point(x: ratio, y: t.isNaN ? 0 : t))
            }
            try Task.checkCancellation()
            progress.update(0.4, phase: "T(a) 垒宽扫描（含厚垒极限）")
            // 图 2：E = 0.3·V₀，a ∈ [0.1, 3] nm，400 点（脚本口径，默认参数下即 E=0.3 eV）
            let eFix = 0.3 * v0b
            var tA: [Point] = [], tlA: [Point] = []
            for i in 0..<400 {
                let aw = (0.1 + 2.9 * Double(i) / 399.0) * 1.0e-9
                tA.append(Point(x: aw * 1.0e9, y: BarrierWellMath.transmission(
                    E: eFix, v0: v0b, width: aw, mass: mElectron, hbar: hbar)))
                tlA.append(Point(x: aw * 1.0e9, y: BarrierWellMath.thickLimit(
                    E: eFix, v0: v0b, width: aw, mass: mElectron, hbar: hbar)))
            }
            try Task.checkCancellation()
            progress.update(0.6, phase: "超越方程扫描 + 二分求根（\(BarrierWellMath.scanPoints) 点 ×80 次）")
            let states = BarrierWellMath.boundStates(
                v0: v0w, width: lMeters, mass: mElectron, hbar: hbar)
            try Task.checkCancellation()
            progress.update(0.9, phase: "波函数重构")
            return (tE, tA, tlA, states)
        }

        // 图 3：波函数（最低 min(5, N) 态，x ∈ [−3L, 3L]，2000 点，叠加 0.6(n−1)）
        let xs = Num.linspace(-3.0 * lMeters, 3.0 * lMeters, count: 2000).map { $0 * 1.0e9 }
        var psiSeries: [SeriesPoints] = []
        var nodeCounts: [Int] = []
        let nshow = min(5, payload.states.count)
        for idx in 0..<nshow {
            let st = payload.states[idx]
            let psi = BarrierWellMath.wavefunction(
                k: st.k, v0: v0w, width: lMeters, even: st.even,
                xs: xs.map { $0 * 1.0e-9 }, mass: mElectron, hbar: hbar)
            nodeCounts.append(SchrodingerFDWellMath.countNodes(psi))
            let offset = 0.6 * Double(idx)
            psiSeries.append(.init(name: "n=\(idx + 1)",
                                   points: zip(xs, psi).map { Point(x: $0, y: $1 + offset) }))
        }

        // 图 4：能级横线（束缚态 + 阱底/渐近线参考线）
        let energySeries = payload.states.enumerated().map { idx, st -> SeriesPoints in
            .init(name: idx == 0 ? "束缚态能级" : "E_\(idx + 1)",
                  points: [Point(x: 0, y: st.E / eCharge), Point(x: 1, y: st.E / eCharge)])
        }

        // 摘要：脚本打印口径 T(E=0.3V₀, a) 与束缚态计数
        let tAnchor = BarrierWellMath.transmission(
            E: 0.3 * v0b, v0: v0b, width: aMeters, mass: mElectron, hbar: hbar)
        let tlAnchor = BarrierWellMath.thickLimit(
            E: 0.3 * v0b, v0: v0b, width: aMeters, mass: mElectron, hbar: hbar)
        let analyticCount = BarrierWellMath.analyticBoundCount(
            v0: v0w, width: lMeters, mass: mElectron, hbar: hbar)
        let z0 = sqrt(2.0 * mElectron * v0w) / hbar * lMeters / 2.0

        let chart0 = LineSeriesData(spec: charts[0].lineSeriesSpec!,
                                    series: [.init(name: "精确 T(E)", points: payload.tE)])
        let chart1 = LineSeriesData(spec: charts[1].lineSeriesSpec!,
                                    series: [
                                        .init(name: "精确 T(a)", points: payload.tA),
                                        .init(name: "厚垒极限 e^(−2κa)", points: payload.tlA),
                                    ])
        let chart2 = LineSeriesData(spec: charts[2].lineSeriesSpec!,
                                    series: psiSeries)
        let chart3 = LineSeriesData(
            spec: charts[3].lineSeriesSpec!,
            series: energySeries,
            referenceLines: [
                ReferenceLine(label: "V = 0（渐近线）", axis: .y, value: 0, style: .subtle),
                ReferenceLine(label: String(format: "阱底 = −%.1f eV", v0wEv), axis: .y,
                              value: -v0wEv, style: .subtle),
            ])

        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chart1),
                     .lineSeries(chart2), .lineSeries(chart3)],
            summary: [
                .init(id: "t", title: "T(E = 0.3V₀, a)",
                      value: String(format: "%.4e", tAnchor.isNaN ? 0 : tAnchor),
                      note: String(format: "厚垒极限 %.4e（a → 大时两者重合）",
                                   tlAnchor.isNaN ? 0 : tlAnchor)),
                .init(id: "count", title: "束缚态个数",
                      value: "\(payload.states.count)（解析 ⌈2z₀/π⌉ = \(analyticCount)）",
                      note: String(format: "z₀ = k_max·L/2 = %.3f", z0)),
                .init(id: "e1", title: "基态能量",
                      value: payload.states.isEmpty ? "—"
                          : String(format: "%.4f eV", payload.states[0].E / eCharge),
                      note: "深阱极限 → 阱底 + ħ²π²/2mL²"),
                .init(id: "nodes", title: "节点数（n=1…）",
                      value: nodeCounts.isEmpty ? "—" : nodeCounts.map(String.init).joined(separator: ", "),
                      note: "已剔除 tan 极点伪根（脚本疑点，见周报）"),
            ],
            theory: TheoryCard(
                title: "方势垒隧穿与有限深势阱",
                formulas: [
                    "T(E) = 1/[1 + V₀²sinh²(κa)/(4E(V₀−E))]，κ = √(2m(V₀−E))/ħ（E < V₀）",
                    "厚垒极限：T ≈ 16E(V₀−E)/V₀²·e^(−2κa)（指数主导）",
                    "E > V₀ 共振：T = 1 当且仅当 kL = nπ（越垒全透）",
                    "偶/奇宇称束缚态：k·tan(kL/2) = κ / −k·cot(kL/2) = κ，个数 ⌈2z₀/π⌉",
                ],
                reading: "左上 T(E) 在深垒区指数压低；右上精确曲线随 a 增大贴上厚垒极限，"
                    + "薄垒区明显偏离。势阱部分：能级数随阱深/阱宽增加（⌈2z₀/π⌉），"
                    + "基态恒为偶宇称、能级按宇称交替。"))
    }
}
