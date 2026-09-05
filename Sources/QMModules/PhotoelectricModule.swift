import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/光电效应_截止电压线性.py）

/// 光电效应纯函数核：eV₀ = hν − W 的截止电压直线族。
/// 单位约定与 Python 一致：W 用 eV、ν 用 Hz、V₀ 用 V；h/e 为普适斜率。
enum PhotoelectricMath {

    /// 普适斜率 h/e（V·s）。
    static func hOverE(h: Double, e: Double) -> Double { h / e }

    /// hc（eV·nm）= 1239.84…，截止波长换算用。
    static func hcEVNm(h: Double, e: Double, c: Double) -> Double { h * c / e * 1e9 }

    /// 截止电压 V₀(ν) = (h/e)ν − W。
    static func stoppingVoltage(_ nu: Double, W: Double, hOverE: Double) -> Double {
        hOverE * nu - W
    }

    /// 截止频率 ν₀ = W/(h/e)；截止波长 λ₀ = hc/W（nm）。
    static func cutoffs(W: Double, hOverE: Double, hcEVNm: Double) -> (nu0: Double, lambda0Nm: Double) {
        (W / hOverE, hcEVNm / W)
    }
}

// MARK: - 模块

/// 笔记 04《光电效应与光子》：四种标准金属的截止电压-频率直线（理论关系，非实验数据点）。
/// 实时档：400 点 × 4 金属，斜率恒为 h/e（直线族的普适性）。
struct PhotoelectricModule: SimModule {

    /// 标准逸出功（文献常见值，eV）——与 Python WORK 字典一致。
    static let metals: [(id: String, symbol: String, work: Double)] = [
        ("Na", "钠", 2.28), ("K", "钾", 2.30), ("Cs", "铯", 2.10), ("W", "钨", 4.50),
    ]

    let meta = ModuleMeta(
        id: "光电效应_截止电压线性", title: "光电效应 · 截止电压线性",
        subtitle: "eV₀ = hν − W：直线斜率 h/e 的普适性与截止波长",
        category: .oldQuantum, noteNumber: 4, tier: .realtime, difficulty: .basic,
        keywords: ["光电效应", "爱因斯坦", "Einstein", "逸出功", "功函数", "截止电压",
                   "截止频率", "截止波长", "密立根", "Millikan", "光子"])

    var params: [ParamSpec] {
        [
            .multiCompare(MultiCompareSpec(
                key: "metals", title: "金属（逸出功 W）",
                candidates: Self.metals.map {
                    .init(id: $0.id, label: "\($0.symbol) \($0.id) · W=\($0.work) eV", value: $0.work)
                },
                defaultSelectionIDs: ["Na", "K", "Cs", "W"], maxSelection: 4)),
            .slider(SliderSpec(key: "numax", title: "最大频率", symbol: "ν", unit: "×10¹⁵ Hz",
                               range: 0.6...3.0, defaultValue: 1.6,
                               step: 0.1, decimalPlaces: 1)),
            .constant(ConstantSpec(key: "h", title: "h", note: "普朗克常量")),
            .constant(ConstantSpec(key: "e", title: "e", note: "元电荷")),
        ]
    }

    var charts: [ChartSpec] {
        [.lineSeries(LineSeriesSpec(
            title: "截止电压 V₀ = (h/e)ν − W",
            xAxis: .init(label: "ν (10¹⁵ Hz)"),
            yAxis: .init(label: "V₀ (V)"),
            seriesNames: Self.metals.map(\.id)))]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let selected = input.multiCompare("metals")
        let nuMax = input.slider("numax") * 1e15
        let h = try constants.value("h")
        let e = try constants.value("e")
        let c = try constants.value("c")
        let hE = PhotoelectricMath.hOverE(h: h, e: e)
        let hcE = PhotoelectricMath.hcEVNm(h: h, e: e, c: c)

        // 与 Python nu = linspace(0, 1.6e15, 400) 同构（上限可调）
        let grid = Num.linspace(0, nuMax, count: 400)

        var series: [SeriesPoints] = []
        var refLines: [ReferenceLine] = []
        var summary: [SummaryItem] = [
            .init(id: "slope", title: "斜率 h/e",
                  value: String(format: "%.7g V·s", hE),
                  note: "所有金属直线共用的普适斜率"),
            .init(id: "hc", title: "hc",
                  value: String(format: "%.2f eV·nm", hcE),
                  note: "光子能量-波长换算"),
        ]

        for m in Self.metals where selected.contains(m.id) {
            // 物理区 V₀ ≥ 0（与 Python mask 一致）
            var pts: [Point] = []
            for nu in grid {
                let v0 = PhotoelectricMath.stoppingVoltage(nu, W: m.work, hOverE: hE)
                if v0 >= 0 { pts.append(Point(x: nu / 1e15, y: v0)) }
            }
            if !pts.isEmpty {
                series.append(.init(name: m.id, points: pts, colorIndex: nil))
            }
            let cut = PhotoelectricMath.cutoffs(W: m.work, hOverE: hE, hcEVNm: hcE)
            refLines.append(ReferenceLine(
                id: "nu0_\(m.id)",
                label: String(format: "ν₀(%@) = %.2f", m.id as NSString, cut.nu0 / 1e15),
                axis: .x, value: cut.nu0 / 1e15))
            summary.append(.init(
                id: "lam0_\(m.id)", title: "λ₀ (\(m.symbol))",
                value: String(format: "%.0f nm", cut.lambda0Nm),
                note: "截止波长 hc/W"))
        }

        let lineSpec = LineSeriesSpec(
            xAxis: .init(label: "ν (10¹⁵ Hz)"),
            yAxis: .init(label: "V₀ (V)"),
            seriesNames: Self.metals.map(\.id))

        return SimResult(
            charts: [.lineSeries(LineSeriesData(spec: lineSpec, series: series,
                                                referenceLines: refLines))],
            summary: summary,
            theory: TheoryCard(
                title: "爱因斯坦光电方程",
                formulas: [
                    "eV₀ = hν − W（W：金属逸出功）",
                    "截止频率 ν₀ = W/(h/e)；截止波长 λ₀ = hc/W",
                    "直线斜率 = h/e = 4.1357×10⁻¹⁵ V·s（与金属无关）",
                ],
                reading: "四条直线平行：斜率普适即「光量子」证据；V₀ < 0 的区段无论光强多强都无光电流（经典波动论无法解释）。"))
    }
}
