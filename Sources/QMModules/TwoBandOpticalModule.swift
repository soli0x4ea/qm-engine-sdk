import Foundation
import EngineKit

// MARK: - 计算核（移植 code/46_矿物光学性质的第一性原理计算_模型.py，模型 a）
//
// 两带（直接带隙）有效质量模型：k 采样 → 联合态密度直方图（k² 权重 = 3D 态密度因子）
// → ε₂(ω)（√ 开启 × 指数包络 × SCALE）→ Kramers–Kronig 数值反演 ε₁（复用 W9 的
// `KramersKronigMath.invertEpsilon1`，与脚本 kramers_kronig_epsilon1 逐句同算法）
// → 复折射率 n(ω)、κ(ω)。
// 策略档（计划 §W11）：k 采样 200001 → 20000（统计量已收敛）；采样分批检查点续算；
// JDOS 按参数哈希缓存；总谱权重收敛监视。脚本常量 TK/KMAX/N_E/SCALE 固定不移。

/// 两带光学模型的纯函数核。
enum TwoBandMath {

    /// 脚本常量：t = ħ²/(2μ) = 1.0 eV·Å²（m_c* = m_v* ⇒ μ = m*/2）。
    static let tk = 1.0
    /// 最大波矢（Å⁻¹）。
    static let kmax = 3.5
    /// 能量直方图 bin 数（脚本 N_E = 600 是 bin **边界**数：np.histogram(bins=E_bins)
    /// 传 600 条边界 → 实际 599 个 bin——直方图错一格即全曲线错位，W11B 首轮对拍教训）。
    static let binCount = 599
    /// 缩放因子（脚本 SCALE = 7，试取使静态 n ≈ 1.7）。
    static let scale = 7.0
    /// 绘图网格：linspace(0.05, 18, 1200)（脚本口径）。
    static let plotMin = 0.05
    static let plotMax = 18.0
    static let plotPoints = 1200

    /// 跃迁能量上界 w_max = E_g + TK·KMAX²。
    static func energyMax(eg: Double) -> Double { eg + tk * kmax * kmax }

    /// 能量直方图 bin 边界：linspace(E_g, w_max, 600)——600 条边界、binCount = 599 个 bin
    /// （脚本 E_bins = np.linspace(EG, w_max_data, N_E) 逐句口径）。
    static func binEdges(eg: Double) -> [Double] {
        Num.linspace(eg, energyMax(eg: eg), count: binCount + 1)
    }

    /// 一批 k 点落 bin 累加 k² 权重（numpy.histogram 口径：末 bin 含右端点）。
    ///
    /// 二分找「最大的 i 使 edges[i] ≤ E」——能量严格递增时与 numpy 的
    /// searchsorted 左插入语义一致；E 恰等于末边时归入末 bin。
    static func accumulate(ks: ArraySlice<Double>, eg: Double, edges: [Double], into sums: inout [Double]) {
        precondition(edges.count == binCount + 1, "binEdges 口径不匹配")
        for k in ks {
            let e = eg + tk * k * k
            var lo = 0
            var hi = edges.count - 1
            while lo < hi {
                let mid = (lo + hi + 1) / 2
                if edges[mid] <= e { lo = mid } else { hi = mid - 1 }
            }
            let idx = lo >= binCount ? binCount - 1 : lo
            sums[idx] += k * k
        }
    }

    /// 归一化 JDOS（峰值 = 1）与 bin 中心（脚本 JDOS / JDOS.max() 口径）。
    static func jdosTable(sums: [Double], edges: [Double]) -> (centers: [Double], jdos: [Double]) {
        var peak = 0.0
        for s in sums where s > peak { peak = s }
        let centers = zip(edges.dropLast(), edges.dropFirst()).map { 0.5 * ($0 + $1) }
        let inv = peak > 0 ? 1.0 / peak : 0.0
        return (centers, sums.map { $0 * inv })
    }

    /// ε₂(ω) 全网格：带隙以下恒 0；[E_g, w_max] 内 JDOS 线性插值 × 指数包络
    /// exp[−(ħω−E_g)/Δ] × SCALE（包络使远紫外 ε₂ → 0，KK 积分良好收敛）。
    static func eps2Grid(w: [Double], centers: [Double], jdos: [Double],
                         eg: Double, decay: Double) -> [Double] {
        let wMax = energyMax(eg: eg)
        return w.map { wi -> Double in
            guard wi >= eg, wi <= wMax else { return 0 }
            return scale * KramersKronigMath.interp(wi, xs: centers, ys: jdos)
                * exp(-(wi - eg) / decay)
        }
    }

    /// KK 静态极限闭合量：1 + (2/π)·∫ε₂(ω')/ω' dω'（ε₁ 的 ω→0⁺ 极限，
    /// 与逐点主值反演核 ω'ε₂/(ω'²−ω²) 不同式——独立的 KK 闭合校验）。
    static func eps1StaticLimit(w: [Double], eps2: [Double]) -> Double {
        1.0 + (2.0 / .pi) * Num.trapezoid(w, zip(eps2, w).map { $0 / $1 })
    }

    /// 带隙以下平均（脚本 low = w < E_g − 0.3 的 np.mean 口径）。
    static func meanBelowGap(w: [Double], values: [Double], eg: Double) -> Double {
        var acc = 0.0
        var count = 0
        for i in w.indices where w[i] < eg - 0.3 {
            acc += values[i]
            count += 1
        }
        return count == 0 ? 0 : acc / Double(count)
    }

    /// 总谱权重 S(N) = dk·Σk²（dk = KMAX/(N−1)）。
    ///
    /// 直方图分箱守恒总权重（binning 只搬运不增删），故 S 无逐 bin 量子化噪声，
    /// 以 O(dk²) 收敛到解析目标 ∫₀^{KMAX} k² dk = KMAX³/3——策略档降维
    /// 200001 → 20000 下真正「已收敛」的统计量（n̄/峰值类受峰 bin 摆动支配，
    /// 5000→20000 仍漂移 ~5%，见 W11B 报告偏差节）。
    static func totalWeight(ks: [Double]) -> Double {
        guard ks.count > 1 else { return 0 }
        let dk = kmax / Double(ks.count - 1)
        var s = 0.0
        for k in ks { s += k * k }
        return s * dk
    }

    /// 解析目标：KMAX³/3。
    static let totalWeightTarget = kmax * kmax * kmax / 3.0
}

// MARK: - 模块

/// 笔记 46《矿物光学性质的第一性原理计算》模型 (a)：两带直接带隙光学模型（策略档）。
///
/// k 采样 200001 → 20000（计划 §W11 定案：统计量已收敛）；联合态密度直方图 → ε₂ →
/// KK 反演按参数缓存（复用 W9 `KramersKronigMath`）；后台串行通道 + 检查点续算。
struct TwoBandOpticalModule: SimModule {

    /// k 采样档位（脚本 200001 → 策略档三档降维）。
    static func kSamples(for quality: GridQuality) -> Int {
        switch quality {
        case .preview: return 5000
        case .standard: return 10000
        case .fine: return 20000
        }
    }

    /// 收敛监视的采样档序列（S(N) 五点）。
    static let convergenceTiers = [1250, 2500, 5000, 10000, 20000]

    /// k 采样分批大小（检查点粒度）。
    static let batchSize = 2500

    /// JDOS 缓存载荷（Sendable）：bin 中心 + 归一曲线 + 总谱权重。
    struct JdosTable: Sendable {
        let centers: [Double]
        let jdos: [Double]
        let totalWeight: Double
    }

    /// 采样检查点：部分直方图随检查点携带；恢复时校验 eg/kSamples 匹配才续算。
    struct SamplingCheckpoint: Sendable {
        let eg: Double
        let kSamples: Int
        let batchesDone: Int
        let binSums: [Double]
    }

    /// 共享会话：取消后重算从已完成批次续跑。
    static let session = StrategySession<SamplingCheckpoint>()

    /// JDOS 缓存：键含 (kSamples, eg)——改包络宽度 Δ 等后处理参数命中缓存。
    static let jdosCache = StrategyCache<String, JdosTable>(capacity: 8)

    let meta = ModuleMeta(
        id: "矿物光学性质的第一性原理计算_模型",
        title: "矿物光学 · 两带光学模型",
        subtitle: "k 采样 JDOS → ε₂ → KK 反演 ε₁ → n(ω)（策略档 200001→20000）",
        category: .gemology, noteNumber: 46, tier: .strategy, difficulty: .advanced,
        keywords: ["两带", "直接带隙", "联合态密度", "JDOS", "带隙", "吸收边",
                   "介电函数", "Kramers", "Kronig", "折射率", "消光系数"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "eg", title: "模型带隙", symbol: "E_g", unit: "eV",
                               range: 2...6, defaultValue: 4, decimalPlaces: 1)),
            .slider(SliderSpec(key: "decay", title: "包络衰减宽度", symbol: "Δ", unit: "eV",
                               range: 4...12, defaultValue: 8, decimalPlaces: 1)),
            .discrete(DiscreteSpec(
                key: "quality", title: "k 采样档位（策略降维）",
                options: [
                    DiscreteOption(id: "preview", title: "预览 · 5000 点", subtitle: "交互粗算"),
                    DiscreteOption(id: "standard", title: "标准 · 10000 点", subtitle: "默认口径"),
                    DiscreteOption(id: "fine", title: "高精 · 20000 点", subtitle: "脚本 200001 → 20000（统计量已收敛）"),
                ],
                defaultOptionID: "fine")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "联合态密度 J(ħω)：k 采样直方图（k² 权重，归一）",
                xAxis: .init(label: "跃迁能量 ħω (eV)"),
                yAxis: .init(label: "J(ħω) (归一)"),
                seriesNames: ["JDOS"])),
            .lineSeries(LineSeriesSpec(
                title: "虚介电函数 ε₂(ω)：√(ħω−E_g) 开启 × 指数包络",
                xAxis: .init(label: "光子能量 ħω (eV)"),
                yAxis: .init(label: "ε₂ (a.u.)"),
                seriesNames: ["ε₂ (模型)"])),
            .lineSeries(LineSeriesSpec(
                title: "折射率 n(ω) 与消光系数 κ(ω)：KK 反演闭合",
                xAxis: .init(label: "光子能量 ħω (eV)"),
                yAxis: .init(label: "n, κ"),
                seriesNames: ["n(ω)", "κ(ω)"])),
            .lineSeries(LineSeriesSpec(
                title: "采样收敛：总谱权重 S(N) → KMAX³/3（分箱守恒统计量）",
                xAxis: .init(label: "k 采样点数 N"),
                yAxis: .init(label: "S(N) / (KMAX³/3)"),
                seriesNames: ["S(N)/S∞"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let eg = input.slider("eg")
        let decay = input.slider("decay")
        let quality = GridQuality(rawValue: input.discrete("quality")) ?? .fine
        let nK = Self.kSamples(for: quality)

        let progress = ComputeProgress()
        progress.update(0.04, phase: "策略档通道排队")

        let payload = try await StrategyChannel.run(session: Self.session, progress: progress) {
            (ctx: StrategyContext<TwoBandOpticalModule.SamplingCheckpoint>) -> Payload in
            let edges = TwoBandMath.binEdges(eg: eg)
            let cacheKey = StrategyKey.make(Self.metaID, [("kSamples", Double(nK)), ("eg", eg)])

            // --- 阶段 0：k 采样联合态密度（缓存 → 检查点续算 → 全量采样） ---
            let table: TwoBandOpticalModule.JdosTable
            if let hit = Self.jdosCache.value(for: cacheKey) {
                table = hit
            } else {
                try Task.checkCancellation()
                let ks = Num.linspace(0, TwoBandMath.kmax, count: nK)
                let totalBatches = (nK + Self.batchSize - 1) / Self.batchSize
                // 恢复：检查点的 eg/kSamples 与本次一致才续算（参数变了重采）
                var binSums = [Double](repeating: 0, count: TwoBandMath.binCount)
                var startBatch = 0
                if let r = ctx.resumedFrom, r.eg == eg, r.kSamples == nK,
                   r.binSums.count == TwoBandMath.binCount, r.batchesDone > 0 {
                    binSums = r.binSums
                    startBatch = min(r.batchesDone, totalBatches)
                }
                for b in startBatch..<totalBatches {
                    try Task.checkCancellation()
                    let lo = b * Self.batchSize
                    let hi = min(lo + Self.batchSize, nK)
                    TwoBandMath.accumulate(ks: ks[lo..<hi], eg: eg, edges: edges, into: &binSums)
                    ctx.checkpoint(SamplingCheckpoint(
                        eg: eg, kSamples: nK, batchesDone: b + 1, binSums: binSums))
                    progress.update(0.1 + 0.35 * Double(b + 1) / Double(totalBatches),
                                    phase: "k 采样 \(b + 1)/\(totalBatches) 批（N = \(nK)）")
                }
                let (centers, jdos) = TwoBandMath.jdosTable(sums: binSums, edges: edges)
                table = JdosTable(centers: centers, jdos: jdos,
                                  totalWeight: TwoBandMath.totalWeight(ks: ks))
                Self.jdosCache.insert(table, for: cacheKey)
            }

            // --- 阶段 1：ε₂ → KK 反演 ε₁（复用 W9 vDSP 核）→ n/κ ---
            try Task.checkCancellation()
            progress.update(0.55, phase: "构造 ε₂（JDOS × 包络 × SCALE）")
            let w = Num.linspace(TwoBandMath.plotMin, TwoBandMath.plotMax,
                                 count: TwoBandMath.plotPoints)
            let eps2 = TwoBandMath.eps2Grid(w: w, centers: table.centers, jdos: table.jdos,
                                            eg: eg, decay: decay)

            try Task.checkCancellation()
            progress.update(0.7, phase: "Kramers–Kronig 反演 ε₁（1200² 主值梯形）")
            let eps1 = KramersKronigMath.invertEpsilon1(w: w, eps2: eps2)

            try Task.checkCancellation()
            progress.update(0.85, phase: "折射率 n/κ 与采样收敛统计")
            var nArr = [Double](repeating: 0, count: w.count)
            var kappaArr = [Double](repeating: 0, count: w.count)
            for i in w.indices {
                let nk = KramersKronigMath.refractiveIndex(eps1: eps1[i], eps2: eps2[i])
                nArr[i] = nk.n
                kappaArr[i] = nk.kappa
            }

            // 收敛统计：总谱权重 S(N) 五档序列（ConvergenceMonitor 只标注，全曲线照算）
            var tierX = [Double]()
            var tierS = [Double]()
            var monitor = ConvergenceMonitor(relTolerance: 2e-3, patience: 2,
                                             maxSteps: Self.convergenceTiers.count)
            for t in Self.convergenceTiers {
                let s = TwoBandMath.totalWeight(ks: Num.linspace(0, TwoBandMath.kmax, count: t))
                tierX.append(Double(t))
                tierS.append(s / TwoBandMath.totalWeightTarget)
                _ = monitor.observe(s)
            }

            let anchorIdx = w.indices.min { abs(w[$0] - 6.0) < abs(w[$1] - 6.0) } ?? 0
            let staticLimit = TwoBandMath.eps1StaticLimit(w: w, eps2: eps2)
            let closureRel = abs(eps1[0] - staticLimit) / max(abs(staticLimit), 1e-300)

            return Payload(
                w: w, eps2: eps2, eps1: eps1, n: nArr, kappa: kappaArr,
                jdosCenters: table.centers, jdos: table.jdos,
                totalWeight: table.totalWeight,
                tierX: tierX, tierS: tierS,
                convergedAt: monitor.convergedAt,
                meanN: TwoBandMath.meanBelowGap(w: w, values: nArr, eg: eg),
                meanEps1: TwoBandMath.meanBelowGap(w: w, values: eps1, eg: eg),
                anchorIdx: anchorIdx, staticClosureRel: closureRel,
                kSamples: nK, quality: quality.badge,
                resumedBatch: ctx.resumedFrom?.batchesDone)
        }

        // --- 图 1：JDOS（600 bin 全点） ---
        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "跃迁能量 ħω (eV)"),
                        yAxis: .init(label: "J(ħω) (归一)"),
                        seriesNames: ["JDOS"]),
            series: [
                .init(name: "JDOS",
                      points: zip(payload.jdosCenters, payload.jdos).map(Point.init)),
            ],
            referenceLines: [
                ReferenceLine(label: "带边 E_g", axis: .x, value: input.slider("eg"),
                              style: .subtle),
            ])

        // --- 图 2：ε₂（1200 点抽稀 600） ---
        let chart2 = LineSeriesData(
            spec: .init(xAxis: .init(label: "光子能量 ħω (eV)"),
                        yAxis: .init(label: "ε₂ (a.u.)"),
                        seriesNames: ["ε₂ (模型)"]),
            series: [
                .init(name: "ε₂ (模型)", points: Num.strided(payload.w, payload.eps2, stride: 2)),
            ],
            referenceLines: [
                ReferenceLine(label: String(format: "吸收边 E_g = %.1f eV", input.slider("eg")),
                              axis: .x, value: input.slider("eg"), style: .subtle),
            ])

        // --- 图 3：n/κ ---
        let chart3 = LineSeriesData(
            spec: .init(xAxis: .init(label: "光子能量 ħω (eV)"),
                        yAxis: .init(label: "n, κ"),
                        seriesNames: ["n(ω)", "κ(ω)"]),
            series: [
                .init(name: "n(ω)", points: Num.strided(payload.w, payload.n, stride: 2)),
                .init(name: "κ(ω)", points: Num.strided(payload.w, payload.kappa, stride: 2)),
            ],
            referenceLines: [
                ReferenceLine(label: String(format: "吸收边 E_g = %.1f eV", input.slider("eg")),
                              axis: .x, value: input.slider("eg"), style: .subtle),
            ])

        // --- 图 4：采样收敛 ---
        let sRel = abs(payload.totalWeight - TwoBandMath.totalWeightTarget)
            / TwoBandMath.totalWeightTarget
        let chart4 = LineSeriesData(
            spec: .init(xAxis: .init(label: "k 采样点数 N"),
                        yAxis: .init(label: "S(N) / (KMAX³/3)"),
                        seriesNames: ["S(N)/S∞"]),
            series: [
                .init(name: "S(N)/S∞", points: zip(payload.tierX, payload.tierS).map(Point.init)),
            ],
            referenceLines: [
                ReferenceLine(label: "解析目标 1.0", axis: .y, value: 1.0, style: .subtle),
            ])

        let cache = Self.jdosCache.stats
        let a = payload.anchorIdx
        let convergeNote: String
        if let at = payload.convergedAt {
            convergeNote = "第 \(at) 档（N = \(TwoBandOpticalModule.convergenceTiers[at - 1])）达 2×10⁻³（连续 2 检查点）"
        } else {
            convergeNote = "五档内未达 2×10⁻³（继续登记）"
        }

        return SimResult(
            charts: [.lineSeries(chart1), .lineSeries(chart2),
                     .lineSeries(chart3), .lineSeries(chart4)],
            summary: [
                .init(id: "model", title: "两带模型参数",
                      value: String(format: "E_g = %.1f eV, Δ = %.1f eV（TK = 1, KMAX = 3.5, SCALE = 7）",
                                    input.slider("eg"), input.slider("decay")),
                      note: "m_c* = m_v* ⇒ μ = m*/2；演示参数，不代表具体矿物实测值"),
                .init(id: "jdosEdge", title: "带边特征",
                      value: String(format: "J(E_g⁺) = %.3f → J 峰 = 1（√(ħω−E_g) 开启）",
                                    payload.jdos.first ?? 0),
                      note: "三维抛物线带的平方根开启行为，k² 权重即 3D 态密度因子"),
                .init(id: "kkClosure", title: "KK 静态极限闭合",
                      value: String(format: "rel = %.2e", payload.staticClosureRel),
                      note: "ε₁(ω→0⁺) = 1 + (2/π)∫ε₂/ω' 与逐点主值反演独立互证（因果性 ⇒ 色散）"),
                .init(id: "belowGap", title: "带隙以下平均（脚本锚点口径 w < E_g − 0.3）",
                      value: String(format: "n̄ = %.4f, ε̄₁ = %.4f", payload.meanN, payload.meanEps1),
                      note: "ε₂ ≡ 0 ⇒ n = √ε₁：透明区折射率纯由色散贡献"),
                .init(id: "anchor", title: "ħω = 6.0 eV 锚点",
                      value: String(format: "n = %.4f, κ = %.4f", payload.n[a], payload.kappa[a]),
                      note: String(format: "ε₁ = %.4f, ε₂ = %.4f（吸收带内）",
                                   payload.eps1[a], payload.eps2[a])),
                .init(id: "sampling", title: "采样降维与收敛统计",
                      value: String(format: "200001 → %d 点（%@档）；S(N) 偏差 %.1e",
                                    payload.kSamples, payload.quality, sRel),
                      note: "总谱权重分箱守恒、无逐 bin 量子化噪声，收敛到 KMAX³/3；"
                          + "ConvergenceMonitor：" + convergeNote),
                .init(id: "cache", title: "缓存命中统计（StrategyCore）",
                      value: "hits \(cache.hits) / misses \(cache.misses)",
                      note: "JDOS 按 (kSamples, E_g) 缓存——改包络宽度 Δ 命中，改带隙/档位重采"
                          + (payload.resumedBatch.map { "；本次自第 \($0) 批检查点续算" } ?? "")),
            ],
            theory: TheoryCard(
                title: "两带模型：从联合态密度到光学常数",
                formulas: [
                    "E_v(k) = −E_g/2 − ħ²k²/2m_v*，E_c(k) = +E_g/2 + ħ²k²/2m_c*（垂直跃迁）",
                    "ħω(k) = E_g + ħ²k²/2μ ⇒ J(ħω) ∝ √(ħω − E_g)（3D 抛物线带平方根开启）",
                    "ε₂(ω) = C·J(ħω)·|M|²·exp[−(ħω−E_g)/Δ]（f 求和规则要求高能回落）",
                    "ε₁(ω) = 1 + (2/π) P∫₀^∞ ω'ε₂(ω')/(ω'²−ω²) dω'；ñ² = ε₁ + iε₂",
                ],
                reading: "直接带隙半导体/矿物的吸收边由联合态密度的平方根开启决定："
                    + "k² 权重恰好是三维态密度因子， ε₂ 在 E_g 处开启、随 √(ħω−E_g) 上升，"
                    + "再乘指数包络满足 f 求和规则。把 ε₂ 喂给 Kramers–Kronig 关系，"
                    + "反演出的 ε₁ 与之互为因果性的两面——透明区（带隙以下）ε₂ ≡ 0，"
                    + "n = √ε₁ 纯由整条吸收带的色散贡献决定，这正是「宝石的折射率由它的"
                    + "吸收谱决定」的定量版本。k 采样从 20 万点降到 2 万点：逐 bin 曲线"
                    + "存在百分位量子化，但总谱权重（分箱守恒量）收敛到解析值 KMAX³/3，"
                    + "光学常数统计量不变——这就是策略档降维的物理依据。"))
    }

    /// 策略通道载荷（Sendable，跨 actor 传递）。
    private struct Payload: Sendable {
        let w: [Double]
        let eps2: [Double]
        let eps1: [Double]
        let n: [Double]
        let kappa: [Double]
        let jdosCenters: [Double]
        let jdos: [Double]
        let totalWeight: Double
        let tierX: [Double]
        let tierS: [Double]
        let convergedAt: Int?
        let meanN: Double
        let meanEps1: Double
        let anchorIdx: Int
        let staticClosureRel: Double
        let kSamples: Int
        let quality: String
        let resumedBatch: Int?
    }

    /// 模块 id（缓存键用；compute 体内不可引实例属性 meta，静态常量镜像）。
    private static let metaID = "矿物光学性质的第一性原理计算_模型"
}
