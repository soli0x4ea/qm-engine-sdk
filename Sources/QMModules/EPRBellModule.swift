import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/EPR与Bell不等式_CHSH.py）

/// EPR 佯谬与 Bell-CHSH 不等式的纯函数核。
enum EPRBellMath {

    static let tsirelson = 2.0 * sqrt(2.0)

    /// LHV 枚举：A(a), A(a'), B(b), B(b') 全部 16 种 ±1 取值组合，
    /// 验证 |S| ≤ 2（脚本 lhv_enumerate）。
    static func lhvMaxS() -> Double {
        var maxAbs = 0.0
        for aa in [-1.0, 1.0] {
            for aap in [-1.0, 1.0] {
                for bb in [-1.0, 1.0] {
                    for bbp in [-1.0, 1.0] {
                        let s = aa * bb - aa * bbp + aap * bb + aap * bbp
                        maxAbs = max(maxAbs, abs(s))
                    }
                }
            }
        }
        return maxAbs
    }

    /// CHSH 组合量 S = E(a,b) − E(a,b') + E(a',b) + E(a',b')（笔记式 2.9）。
    static func chsh(_ eab: Double, _ eabp: Double, _ eapb: Double, _ eapbp: Double) -> Double {
        eab - eabp + eapb + eapbp
    }

    /// 自旋单态关联 E(a,b) = −a·b，共面方向即 −cos(θa − θb)（笔记式 2.12）。
    static func E_singlet(_ thetaA: Double, _ thetaB: Double) -> Double {
        -cos(thetaA - thetaB)
    }

    /// 单态扫描族（a=0, b=φ, a'=2φ, b'=3φ）：S(φ) = cos3φ − 3cosφ，
    /// φ = π/4 处取极值 −2√2（脚本 singlet_scan）。
    static func singletS(_ phi: Double) -> Double {
        let eab = E_singlet(0, phi)
        let eabp = E_singlet(0, 3 * phi)
        let eapb = E_singlet(2 * phi, phi)
        let eapbp = E_singlet(2 * phi, 3 * phi)
        return chsh(eab, eabp, eapb, eapbp)
    }

    /// 光子偏振纠缠 |Φ⁺⟩ 关联 E = cos 2(θa − θb)（笔记式 2.13）。
    static func E_phipol(_ thetaA: Double, _ thetaB: Double) -> Double {
        cos(2.0 * (thetaA - thetaB))
    }

    /// 偏振扫描族（a=0, a'=2φ, b=φ, b'=3φ）：S(φ) = 3cos2φ − cos6φ，
    /// φ = 22.5°（Bell 角配置）处达 +2√2（脚本 photon_scan）。
    static func photonS(_ phi: Double) -> Double {
        let eab = E_phipol(0, phi)
        let eabp = E_phipol(0, 3 * phi)
        let eapb = E_phipol(2 * phi, phi)
        let eapbp = E_phipol(2 * phi, 3 * phi)
        return chsh(eab, eabp, eapb, eapbp)
    }

    /// Bell (1964) 原始不等式违反样例（笔记式 2.10）：
    /// a⊥b，c 与两者各成 45°；返回 (|P(a,b)−P(a,c)|, 1+P(b,c))。
    static func bell1964() -> (lhs: Double, rhs: Double) {
        let pab = 0.0
        let pac = -cos(Double.pi / 4)
        let pbc = -cos(Double.pi / 4)
        return (abs(pab - pac), 1.0 + pbc)
    }

    /// Werner（退极化）信道：|S|(p) = p·2√2，违反 CHSH 需 p > 1/√2。
    /// 返回 (p 网格, S 网格, 首个 S > 2 的 p)（脚本 werner_threshold）。
    static func wernerThreshold() -> (p: [Double], s: [Double], pThr: Double) {
        let p = Num.linspace(0, 1, count: 10001)
        let s = p.map { $0 * tsirelson }
        let idx = s.firstIndex { $0 > 2.0 } ?? p.count - 1
        return (p, s, p[idx])
    }

    /// Tsirelson 界随机采样（脚本 tsirelson_sampling）：SplitMix64 定种子
    /// （种子 = Python 脚本原值 42，序列不同分布相同，见 docs/RNG_SEED_POLICY.md）
    /// + Box-Muller 三维正态 → 单位矢量，S = −a·b + a·b' − a'·b − a'·b'。
    /// |S| ≤ 2√2 对自旋单态是数学必然（Tsirelson 定理），任意采样不得超出。
    /// W13（B2）：keepSamples > 0 时沿样本序列等步长保留至多该数量的 S 值
    /// （UI 采样分布散点用；N 越大覆盖越广，滑杆真实可见）。
    static func tsirelsonSampling(n: Int, seed: UInt64 = 42, keepSamples: Int = 0)
        -> (maxS: Double, minS: Double, violations: Int, samples: [Double]) {
        var rng = SplitMix64(seed: seed)
        var spare: Double? = nil
        // Box-Muller：缓存正态对的第二个样本
        func nextNormal() -> Double {
            if let z = spare { spare = nil; return z }
            let u1 = max(rng.nextDouble(), Double.leastNonzeroMagnitude)
            let u2 = rng.nextDouble()
            let r = sqrt(-2 * log(u1))
            let ang = 2 * .pi * u2
            spare = r * sin(ang)
            return r * cos(ang)
        }
        func unitVectorDot(_ u: (Double, Double, Double), _ v: (Double, Double, Double)) -> Double {
            u.0 * v.0 + u.1 * v.1 + u.2 * v.2
        }
        func randUnit() -> (Double, Double, Double) {
            var v = (nextNormal(), nextNormal(), nextNormal())
            let norm = sqrt(unitVectorDot(v, v))
            v = (v.0 / norm, v.1 / norm, v.2 / norm)
            return v
        }

        var vMax = -Double.infinity
        var vMin = Double.infinity
        var violations = 0
        var samples: [Double] = []
        let stride = keepSamples > 0 ? max(1, n / keepSamples) : 0
        for i in 0..<n {
            let a = randUnit(), ap = randUnit(), b = randUnit(), bp = randUnit()
            let s = -unitVectorDot(a, b) + unitVectorDot(a, bp)
                - unitVectorDot(ap, b) - unitVectorDot(ap, bp)
            vMax = max(vMax, s)
            vMin = min(vMin, s)
            if abs(s) > tsirelson + 1e-9 { violations += 1 }
            if stride > 0, i % stride == 0 { samples.append(s) }
        }
        return (vMax, vMin, violations, samples)
    }
}

// MARK: - 实验数据（脚本 EXPERIMENTS：S 实测值 + 不确定度）

/// 历代 Bell 实验选列（S 为 CHSH/Bell 参数实测值）。
enum EPRBellExperiments {
    /// (标签, S, ±σ, 年份, 体系)
    static let all: [(label: String, s: Double, err: Double, year: Int, system: String)] = [
        ("Aspect et al. 1982", 2.697, 0.015, 1982, "Ca 级联"),
        ("Rowe et al. 2001", 2.25, 0.03, 2001, "Be⁺ 离子"),
        ("Hensen et al. 2015", 2.42, 0.20, 2015, "NV 色心"),
    ]
}

// MARK: - 模块

/// 笔记 30《EPR 佯谬与贝尔不等式》：CHSH 组合量数值验证——
/// LHV 枚举界 2 / 量子极值 2√2（自旋单态与光子偏振两族扫描）+
/// 20 万定种子随机采样 + 实验时间线。秒级档。
struct EPRBellModule: SimModule {

    /// 采样网格口径与脚本一致：单态 1000 点、光子 2000 点。
    static let singletGridN = 1000
    static let photonGridN = 2000

    let meta = ModuleMeta(
        id: "EPR佯谬与Bell不等式_CHSH", title: "EPR 佯谬与 Bell 不等式",
        subtitle: "CHSH 组合量：LHV 界 2 与量子极值 2√2",
        category: .quantumInfo, noteNumber: 30, tier: .seconds, difficulty: .basic,
        keywords: ["EPR", "Bell", "CHSH", "不等式", "定域隐变量", "LHV",
                   "Tsirelson", "纠缠", "自旋单态", "偏振纠缠", "关联函数"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "N", title: "随机采样数", symbol: "N", unit: "次",
                               range: 20000...200000, defaultValue: 200000,
                               decimalPlaces: 0)),
            // W13h #30：N 滑杆只影响统计收敛（视觉不可见），补一个直接进入图形的
            // Werner 退极化强度滑杆——拖 p 沿 |S|(p) 直线滑动、跨 1/√2 违反阈。
            .slider(SliderSpec(key: "wernerP", title: "Werner 态权重", symbol: "p", unit: "",
                               range: 0...1, defaultValue: 0.9,
                               scale: .linear, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "自旋单态 CHSH 扫描：E(a,b) = −a·b",
                xAxis: .init(label: "共面角间隔 φ (°)"),
                yAxis: .init(label: "CHSH 组合量 S"),
                seriesNames: ["自旋单态 S(φ)"])),
            .lineSeries(LineSeriesSpec(
                title: "光子偏振纠缠 |Φ⁺⟩ CHSH 扫描",
                xAxis: .init(label: "偏振片角间隔 θ (°)"),
                yAxis: .init(label: "CHSH 组合量 S"),
                seriesNames: ["偏振纠缠 S(θ)"])),
            .scatter(ScatterSpec(
                title: "历代实验实测 CHSH 值（选列）",
                xAxis: .init(label: "实测 CHSH 值 S"),
                yAxis: .init(label: "年份"),
                seriesNames: EPRBellExperiments.all.map(\.label))),
            .scatter(ScatterSpec(
                title: "Tsirelson 采样分布（N 次随机方向的 S 值，均匀抽样 ≤600 点）",
                xAxis: .init(label: "样本序号"),
                yAxis: .init(label: "CHSH 值 S"),
                seriesNames: ["随机方向采样 S"])),
            // W13h #30：Werner 信道随 p 联动
            .lineSeries(LineSeriesSpec(
                title: "Werner 信道：|S|(p) = p·2√2（CHSH 违反阈 p > 1/√2）",
                xAxis: .init(label: "Werner 态权重 p"),
                yAxis: .init(label: "CHSH 组合量 |S|"),
                seriesNames: ["|S|(p) = p·2√2"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let n = Int(input.slider("N"))
        let pWerner = input.slider("wernerP")

        // --- [2] 自旋单态扫描：φ ∈ [0, π/2]，1000 点网格（脚本口径） ---
        let phiGrid = Num.linspace(0, .pi / 2, count: Self.singletGridN)
        let sSinglet = phiGrid.map(EPRBellMath.singletS)
        var idxSinglet = 0
        for i in sSinglet.indices where abs(sSinglet[i]) > abs(sSinglet[idxSinglet]) {
            idxSinglet = i
        }
        let phiBest = phiGrid[idxSinglet]
        let sBest = sSinglet[idxSinglet]

        // --- [3] 光子偏振扫描：θ ∈ [0, π/4]，2000 点网格（脚本口径） ---
        let thetaGrid = Num.linspace(0, .pi / 4, count: Self.photonGridN)
        let sPhoton = thetaGrid.map(EPRBellMath.photonS)
        var idxPhoton = 0
        for i in sPhoton.indices where abs(sPhoton[i]) > abs(sPhoton[idxPhoton]) {
            idxPhoton = i
        }
        let thetaBest = thetaGrid[idxPhoton]
        let sphBest = sPhoton[idxPhoton]

        // 标准 Bell 角配置 (0°, 45° | 22.5°, 67.5°)
        let d2r = Double.pi / 180
        let sBellAngles = EPRBellMath.chsh(
            EPRBellMath.E_phipol(0, 22.5 * d2r),
            EPRBellMath.E_phipol(0, 67.5 * d2r),
            EPRBellMath.E_phipol(45 * d2r, 22.5 * d2r),
            EPRBellMath.E_phipol(45 * d2r, 67.5 * d2r))

        // --- [4] Tsirelson 界随机采样（SplitMix64 定种子 42；W13 附带均匀抽样序列） ---
        // W13h #30：显示点数随 N 缩放（N/333，60→600）——拖 N 散点密度可见变化
        let sample = EPRBellMath.tsirelsonSampling(
            n: n, keepSamples: max(60, min(600, n / 333)))

        // --- [4b] Bell 1964 原始不等式违反样例 ---
        let b1964 = EPRBellMath.bell1964()

        // --- [4c] Werner 阈值 ---
        let werner = EPRBellMath.wernerThreshold()

        // 图 1：自旋单态扫描（参考线：LHV ±2、Tsirelson ±2√2、φ*）
        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "共面角间隔 φ (°)"),
                        yAxis: .init(label: "CHSH 组合量 S"),
                        seriesNames: ["自旋单态 S(φ)"]),
            series: [.init(name: "自旋单态 S(φ)",
                           points: Num.strided(phiGrid.map { $0 / d2r }, sSinglet, stride: 2))],
            referenceLines: [
                ReferenceLine(label: "LHV 界 |S| = 2", axis: .y, value: 2),
                ReferenceLine(label: "LHV 界 |S| = −2", axis: .y, value: -2, style: .subtle),
                ReferenceLine(label: "Tsirelson 界 2√2", axis: .y, value: EPRBellMath.tsirelson),
                ReferenceLine(label: "Tsirelson 界 −2√2", axis: .y, value: -EPRBellMath.tsirelson, style: .subtle),
                ReferenceLine(label: String(format: "φ* = %.1f°", phiBest / d2r),
                              axis: .x, value: phiBest / d2r, style: .subtle),
            ],
            pointMarkers: [
                PointMarker(x: phiBest / d2r, y: sBest,
                            label: String(format: "φ* → %.3f", sBest)),
            ])

        // 图 2：光子偏振扫描（参考线：LHV 2、Tsirelson 2√2、θ*）
        let chart2 = LineSeriesData(
            spec: .init(xAxis: .init(label: "偏振片角间隔 θ (°)"),
                        yAxis: .init(label: "CHSH 组合量 S"),
                        seriesNames: ["偏振纠缠 S(θ)"]),
            series: [.init(name: "偏振纠缠 S(θ)",
                           points: Num.strided(thetaGrid.map { $0 / d2r }, sPhoton, stride: 4))],
            referenceLines: [
                ReferenceLine(label: "LHV 界 S = 2", axis: .y, value: 2),
                ReferenceLine(label: "Tsirelson 界 2√2", axis: .y, value: EPRBellMath.tsirelson),
                ReferenceLine(label: String(format: "θ* = %.2f°", thetaBest / d2r),
                              axis: .x, value: thetaBest / d2r, style: .subtle),
            ])

        // 图 3：实验时间线（散点 + LHV/Tsirelson 竖线）
        let chart3 = ScatterData(
            spec: .init(xAxis: .init(label: "实测 CHSH 值 S"),
                        yAxis: .init(label: "年份"),
                        seriesNames: EPRBellExperiments.all.map(\.label)),
            series: EPRBellExperiments.all.map {
                .init(name: $0.label, points: [.init(x: $0.s, y: Double($0.year))])
            },
            referenceLines: [
                ReferenceLine(label: "LHV 界 S = 2", axis: .x, value: 2),
                ReferenceLine(label: "Tsirelson 界 2√2", axis: .x, value: EPRBellMath.tsirelson),
            ])

        // 图 4（W13/B2）：采样分布散点——W13h #30：点数随 N 缩放（密度可见），
        // 极值 max/min 以横参考线标注（N 越大越逼近 ±2√2）
        let chart4 = ScatterData(
            spec: .init(xAxis: .init(label: "样本序号"),
                        yAxis: .init(label: "CHSH 值 S"),
                        seriesNames: ["随机方向采样 S"]),
            series: [.init(name: "随机方向采样 S",
                           points: sample.samples.enumerated().map {
                               Point(x: Double($0.offset), y: $0.element)
                           })],
            referenceLines: [
                ReferenceLine(label: "LHV 界 ±2", axis: .y, value: 2),
                ReferenceLine(label: String(format: "max S = %.4f", sample.maxS),
                              axis: .y, value: sample.maxS, style: .subtle),
                ReferenceLine(label: String(format: "min S = %.4f", sample.minS),
                              axis: .y, value: sample.minS, style: .subtle),
            ])

        // 图 5（W13h #30）：Werner 信道 |S|(p)——当前 p 大圆点沿线滑动，
        // 跨 p = 1/√2 违反阈即进可 violating 区
        let pGrid = Num.linspace(0, 1, count: 200)
        let chart5 = LineSeriesData(
            spec: .init(xAxis: .init(label: "Werner 态权重 p"),
                        yAxis: .init(label: "CHSH 组合量 |S|"),
                        seriesNames: ["|S|(p) = p·2√2"]),
            series: [.init(name: "|S|(p) = p·2√2",
                           points: pGrid.map { Point(x: $0, y: $0 * EPRBellMath.tsirelson) })],
            referenceLines: [
                ReferenceLine(label: "LHV 界 |S| = 2", axis: .y, value: 2),
                ReferenceLine(label: "违反阈 p = 1/√2", axis: .x, value: 1 / sqrt(2), style: .subtle),
            ],
            pointMarkers: [
                PointMarker(x: pWerner, y: pWerner * EPRBellMath.tsirelson,
                            label: String(format: "p = %.2f → |S| = %.4f",
                                          pWerner, pWerner * EPRBellMath.tsirelson),
                            colorIndex: 1),
            ])

        let lhvMax = EPRBellMath.lhvMaxS()
        return SimResult(
            charts: [.lineSeries(chart1), .lineSeries(chart2),
                     .scatter(chart3), .scatter(chart4), .lineSeries(chart5)],
            summary: [
                .init(id: "lhv", title: "LHV 枚举界",
                      value: String(format: "max|S| = %.4f", lhvMax),
                      note: "16 种 ±1 取值组合穷举，定域实在论上界恰为 2"),
                .init(id: "singlet", title: "自旋单态极值",
                      value: String(format: "|S(φ*)| = %.4f", abs(sBest)),
                      note: String(format: "φ* = %.2f° 处取 −2√2", phiBest / d2r)),
                .init(id: "photon", title: "偏振纠缠极值",
                      value: String(format: "S(θ*) = %.4f", sphBest),
                      note: String(format: "θ* = %.2f°（Bell 角 22.5°），配置值 %.4f",
                                   thetaBest / d2r, sBellAngles)),
                .init(id: "sampling", title: "Tsirelson 采样",
                      value: String(format: "n = %d，越界 %d", n, sample.violations),
                      note: String(format: "max S = %.4f，min S = %.4f（种子 42 固定）",
                                   sample.maxS, sample.minS)),
                .init(id: "bell1964", title: "Bell (1964) 违反",
                      value: String(format: "%.4f > %.4f", b1964.lhs, b1964.rhs),
                      note: "|P(a,b)−P(a,c)| 超出 1+P(b,c)（a⊥b，c 居中 45°）"),
                .init(id: "werner", title: "Werner 阈值",
                      value: String(format: "p* = %.4f", werner.pThr),
                      note: String(format: "解析值 1/√2 = %.6f", 1 / sqrt(2))),
                .init(id: "wernerCur", title: "Werner 当前 |S|(p)",
                      value: String(format: "p = %.2f → %.4f", pWerner,
                                    pWerner * EPRBellMath.tsirelson),
                      note: pWerner > 1 / sqrt(2)
                          ? "p > 1/√2：CHSH 违反区"
                          : "p ≤ 1/√2：不违反（图 5 大圆点在阈左侧）"),
                .init(id: "experiments", title: "实验选列",
                      value: "2.697 / 2.25 / 2.42",
                      note: "Aspect 1982 ±0.015 · Rowe 2001 ±0.03 · Hensen 2015 ±0.20"),
            ],
            theory: TheoryCard(
                title: "CHSH 不等式与 Tsirelson 界",
                formulas: [
                    "S = E(a,b) − E(a,b′) + E(a′,b) + E(a′,b′)，LHV：|S| ≤ 2",
                    "自旋单态 E(a,b) = −a·b，光子偏振 E = cos 2(θa−θb)",
                    "量子极值 |S|max = 2√2 ≈ 2.8284（Tsirelson 界，任意纠缠态不得超出）",
                    "Werner 态 ρ = p|Ψ⁻⟩⟨Ψ⁻| + (1−p)I/4：违反 CHSH 需 p > 1/√2",
                ],
                reading: "左图自旋单态在 φ=45° 处到 −2√2（绿色点线为 Tsirelson 界，"
                    + "灰色虚线为 LHV 界 ±2）——量子关联在经典界之外、又在 Tsirelson 界之内；"
                    + "右图偏振族在 Bell 角配置 (0°, 45° | 22.5°, 67.5°) 恰达 +2√2。"
                    + "20 万次随机方向采样（SplitMix64 种子 42）无一越界：2√2 是数学必然，"
                    + "而非统计涨落。实验点全部越过 LHV 界——定域实在论被否定。"))
    }
}
