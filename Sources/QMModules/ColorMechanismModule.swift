import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/43_过渡金属离子致色的量子机制_振动耦合.py / _模型.py）

/// 过渡金属离子致色机制纯函数核：振动耦合饱和强度、晶体场分裂→吸收波长反比律。
enum ColorMechanismMath {
    /// 无量纲振动耦合强度 ξ 下的相对 d-d 跃迁强度（归一化饱和形式）。
    /// I(ξ) = ξ² / (ξ² + ξ0²)，ξ0 为半强度耦合常数。
    static func relativeIntensity(xi: Double, xi0: Double) -> Double {
        xi * xi / (xi * xi + xi0 * xi0)
    }

    /// 由晶体场分裂 Δ (cm⁻¹) 反算吸收波长 λ (nm)：λ = 10⁷ / Δ。
    static func deltaToLambdaNm(deltaCm: Double) -> Double { 1.0e7 / deltaCm }

    /// h·c（eV·nm），CODATA 精确派生量（脚本固定 1239.841984）。
    static let hcEVnm: Double = 1239.841984

    /// 八面体 O_h 场 d 轨道分裂（Dq 单位）：E(eg)=+6Dq, E(t2g)=−4Dq, Δ_o=10Dq。
    static func octahedralSplitting(dqCm: Double) -> (eEg: Double, eT2g: Double, delta: Double, barycenter: Double) {
        let eEg = 6.0 * dqCm
        let eT2g = -4.0 * dqCm
        let delta = eEg - eT2g
        let bary = 2.0 * eEg + 3.0 * eT2g
        return (eEg: eEg, eT2g: eT2g, delta: delta, barycenter: bary)
    }
}

// MARK: - 模块一：振动耦合饱和强度

/// 笔记 43《过渡金属离子致色的量子机制》：奇宇称振动耦合弛豫 Laporte 禁戒
/// I(ξ) = ξ²/(ξ²+ξ0²)——小 ξ 线性、大 ξ 饱和。
struct ColorVibrationalCouplingModule: SimModule {

    let meta = ModuleMeta(
        id: "过渡金属离子致色的量子机制_振动耦合", title: "振动耦合 · d-d 强度饱和",
        subtitle: "I(ξ)=ξ²/(ξ²+ξ₀²)：宇称禁戒靠奇宇称振动弛豫，小ξ线性、大ξ饱和",
        category: .gemology, noteNumber: 43, tier: .realtime, difficulty: .advanced,
        keywords: ["振动耦合", "Laporte", "宇称禁戒", "d-d跃迁", "强度饱和", "vibronic",
                   "vibronic coupling", "Laporte forbidden", "intensity", "selection rule"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "xi0", title: "半强度常数 ξ₀", symbol: "ξ₀", unit: "",
                              range: 0.1...1.5, defaultValue: 0.5,
                              scale: .linear, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "d-d 跃迁相对强度 vs 奇宇称振动耦合强度 ξ",
                xAxis: .init(label: "振动耦合强度 ξ (无量纲，示意)"),
                yAxis: .init(label: "相对强度 I(ξ)"),   // W13d：回线性——I(ξ)∈[0,1] 饱和型，log 下形态与笔记 43 图2 不符
                seriesNames: ["I(ξ) = ξ²/(ξ²+ξ₀²)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let xi0 = input.slider("xi0")
        let xi = Num.linspace(1e-3, 3.0, count: 600)  // 网格保持 1e-3 起（fixture 对拍口径）
        let I = xi.map { ColorMechanismMath.relativeIntensity(xi: $0, xi0: xi0) }
        let area = Num.trapezoid(xi, I)
        // W13h #43：显示窗随 ξ₀ 收缩——x ≤ min(3, max(0.6, 4ξ₀))。原恒画 [0,3]：
        // ξ₀ 小时右侧 2/3 是平坦饱和死区。裁剪后饱和膝点恒居画面中部，
        // 拖 ξ₀ 时窗口与 ξ₀ 竖线联动。
        let xWin = min(3.0, max(0.6, 4.0 * xi0))

        let refLines = [
            ReferenceLine(label: "I = 0.5（半强度，ξ = ξ₀）", axis: .y, value: 0.5, style: .subtle),
            ReferenceLine(label: "ξ₀ = \(String(format: "%.2f", xi0))", axis: .x, value: xi0, style: .threshold),
        ]

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: [
                        .init(name: "I(ξ) = ξ²/(ξ²+ξ₀²)",
                              points: Num.strided(xi, I, stride: 1)
                                  .filter { $0.x <= xWin }),
                    ],
                    referenceLines: refLines,
                    // W13（B3）：ξ₀ 半强度交点大圆点（阈值线视觉强化）
                    pointMarkers: [
                        PointMarker(x: xi0, y: 0.5,
                                    label: String(format: "ξ₀ = %.2f, I = 0.5", xi0),
                                    colorIndex: 1),
                    ])),
            ],
            summary: [
                .init(id: "half", title: "半强度点 I(ξ₀)=0.5",
                      value: String(format: "%.4f", ColorMechanismMath.relativeIntensity(xi: xi0, xi0: xi0)),
                      note: "ξ = ξ₀ 处强度恰为 1/2"),
                .init(id: "small", title: "小 ξ 近似 I≈ξ²/ξ₀²",
                      value: String(format: "%.4f", ColorMechanismMath.relativeIntensity(xi: 0.1, xi0: xi0)),
                      note: "ξ=0.1：≈ ξ²/ξ₀²（线性于 ξ²）"),
                .init(id: "sat", title: "饱和值 I(∞)=1",
                      value: String(format: "%.4f", ColorMechanismMath.relativeIntensity(xi: 3.0, xi0: xi0)),
                      note: "ξ=3：趋于饱和上限 1"),
                .init(id: "area", title: "∫ I(ξ) dξ (ξ∈[0,3])",
                      value: String(format: "%.3f", area), note: "np.trapezoid 数值积分"),
            ],
            theory: TheoryCard(
                title: "振动耦合弛豫 Laporte 禁戒（笔记 43 第六节）",
                formulas: [
                    "宇称禁戒 d-d 跃迁本征偶极矩 ≈ 0；奇宇称振动（eg/T1u）瞬时破缺反演对称",
                    "有效跃迁矩 ∝ ξ → 相对强度 I ∝ |μ_eff|² ∝ ξ²（小 ξ）",
                    "归一化饱和形式：I(ξ) = ξ²/(ξ² + ξ₀²)，ξ₀ 为半强度耦合常数",
                ],
                reading: "d-d 跃迁本身宇称禁戒（强度仅 ε~1–100），靠奇宇称振动混入 u 成分获得微弱强度。"
                    + "强度随振动耦合 ξ 先按 ξ² 上升、后饱和于 1——ξ₀ 是半强度点。电荷转移跃迁（强耦合）"
                    + "可达 ε~10³–10⁴，对应右上方饱和区。"))
    }
}

// MARK: - 模块二：致色模型（晶体场分裂 → 吸收波长）

/// 笔记 43《过渡金属离子致色的量子机制》：λ = 10⁷/Δ 反比律 + 矿物数据表
/// 代表矿物（红宝石/祖母绿/橄榄石）的晶体场分裂 Δ 来自 MaterialDB。
struct ColorMechanismModelModule: SimModule {

    let meta = ModuleMeta(
        id: "过渡金属离子致色的量子机制_模型", title: "致色模型 · λ = 10⁷/Δ",
        subtitle: "吸收波长 λ(nm)=10⁷/Δ(cm⁻¹) 反比律；红宝石/祖母绿/橄榄石对照 MaterialDB",
        category: .gemology, noteNumber: 43, tier: .realtime, difficulty: .advanced,
        keywords: ["致色", "晶体场分裂", "吸收波长", "红宝石", "祖母绿", "橄榄石", "颜色",
                   "color", "crystal field splitting", "ruby", "emerald", "peridot", "gemstone"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "Dq", title: "晶体场参数 Dq", symbol: "Dq", unit: "cm⁻¹",
                              range: 500...3000, defaultValue: 1830,
                              scale: .linear, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "吸收波长 vs 晶体场分裂（宝石致色离子）",
                xAxis: .init(label: "晶体场分裂 Δ (cm⁻¹)"),
                yAxis: .init(label: "吸收波长 λ (nm)"),
                seriesNames: ["λ = 10⁷ / Δ (nm)"])),
            .scatter(ScatterSpec(
                title: "代表矿物（MaterialDB）",
                xAxis: .init(label: "Δ (cm⁻¹)"),
                yAxis: .init(label: "λ (nm)"),
                seriesNames: ["Ruby (Cr³⁺)", "Emerald (Cr³⁺)", "Peridot (Fe²⁺)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let db = try MaterialDB.load()
        // W13h #43：Dq 滑杆原只进摘要——补 λ(Δ) 曲线上的 10Dq 标记点 + 竖参考线，
        // 拖 Dq 沿反比律曲线滑动（Δ_o = 10Dq 处的吸收波长即时可见）
        let dq = input.slider("Dq")
        let dqDelta = 10.0 * dq
        let dqLambda = ColorMechanismMath.deltaToLambdaNm(deltaCm: dqDelta)

        let delta = Num.linspace(5000.0, 25000.0, count: 400)
        let lam = delta.map { ColorMechanismMath.deltaToLambdaNm(deltaCm: $0) }

        // 可见光对应 Δ ∈ [12820, 25000] cm⁻¹（λ 780–400 nm）
        let refLines = [
            ReferenceLine(label: "可见光下界 λ=780 nm → Δ≈12820", axis: .x, value: 12820, style: .subtle),
            ReferenceLine(label: "可见光上界 λ=400 nm → Δ≈25000", axis: .x, value: 25000, style: .subtle),
            ReferenceLine(label: String(format: "10Dq = %.0f cm⁻¹ (Dq = %.0f)", dqDelta, dq),
                          axis: .x, value: dqDelta, style: .threshold),
        ]

        // 矿物散点（来自 MaterialDB）
        let order = ["ruby", "emerald", "peridot"]
        var scatterSeries: [SeriesPoints] = []
        var mineralSummary: [SummaryItem] = []
        for (i, mid) in order.enumerated() {
            guard let m = db.mineral(id: mid) else { continue }
            let l = ColorMechanismMath.deltaToLambdaNm(deltaCm: m.deltaCm)
            let ev = ColorMechanismMath.hcEVnm / l
            scatterSeries.append(.init(name: "\(m.short) (\(m.ion))",
                                       points: [Point(x: m.deltaCm, y: l)], colorIndex: i))
            mineralSummary.append(.init(
                id: "min_\(mid)", title: "\(m.short) · Δ = \(Int(m.deltaCm)) cm⁻¹",
                value: String(format: "%.0f nm (%.2f eV)", l, ev),
                note: "\(m.ion) in \(m.host)"))
        }

        // O_h 八面体 d 轨道分裂核验（Dq 滑块）
        let sp = ColorMechanismMath.octahedralSplitting(dqCm: dq)

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: [
                        .init(name: "λ = 10⁷ / Δ (nm)",
                              points: Num.strided(delta, lam, stride: 1)),
                    ],
                    referenceLines: refLines,
                    pointMarkers: [
                        PointMarker(x: min(dqDelta, 25000), y: dqLambda,
                                    label: String(format: "10Dq = %.0f cm⁻¹ → λ = %.0f nm",
                                                  dqDelta, dqLambda),
                                    colorIndex: 1),
                    ])),
                .scatter(ScatterData(
                    spec: charts[1].scatterSpec!,
                    series: scatterSeries)),
            ],
            summary: mineralSummary + [
                .init(id: "inv", title: "反比律 λ = 10⁷/Δ",
                      value: "λ·Δ = 10⁷ nm·cm⁻¹", note: "Δ 越大 → λ 越短（越偏蓝/紫）"),
                .init(id: "split", title: "O_h 分裂（Dq=\(Int(dq)) cm⁻¹）",
                      value: String(format: "E(eg)=+%.0f, E(t2g)=−%.0f", sp.eEg, sp.eT2g),
                      note: "Δ_o = E(eg)−E(t2g) = 10Dq = \(Int(sp.delta)) cm⁻¹"),
                .init(id: "bary", title: "质心守恒检验 2E(eg)+3E(t2g)",
                      value: String(format: "%.0f", sp.barycenter), note: "应 = 0（配位场无迹）"),
            ],
            theory: TheoryCard(
                title: "晶体场致色：λ = 10⁷/Δ（笔记 43 第八节）",
                formulas: [
                    "吸收波长 λ(nm) = 10⁷ / Δ(cm⁻¹)（波数↔波长换算）",
                    "O_h 场 d 轨道：E(eg)=+6Dq, E(t2g)=−4Dq, Δ_o=10Dq（质心守恒）",
                    "红宝石 Cr³⁺ Δ≈18000 cm⁻¹ → λ≈556 nm（绿光吸收，呈红）",
                ],
                reading: "跃迁能量 = 晶体场分裂 Δ，吸收波长与 Δ 成反比。同一 Cr³⁺ 在刚玉（红宝石）与"
                    + "绿柱石（祖母绿）中因配体场强不同（Δ 不同）呈现红/绿两色——这就是晶体场致色。"
                    + "橄榄石 Fe²⁺ 近红外主带（~1050 nm）对应 Δ≈9524 cm⁻¹。下方散点为 MaterialDB 收录的代表矿物。"))
    }
}
