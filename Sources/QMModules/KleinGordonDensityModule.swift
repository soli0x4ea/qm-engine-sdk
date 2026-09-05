import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/克莱因-戈登方程_概率密度振荡.py）

/// KG 守恒概率密度振荡的纯函数核。自然单位 ħ = c = m = 1。
enum KleinGordonDensityMath {

    /// ω_j = √(k_j² + m²)，m = 1。
    static func omega(_ k: Double) -> Double { sqrt(k * k + 1.0) }

    /// 场：φ(x,t) = exp(i(k₁x − ω₁t)) + β·exp(−i(k₂x + ω₂t))。
    /// 返回 (re, im)。
    static func field(x: Double, t: Double, k1: Double, k2: Double, beta: Double)
        -> (re: Double, im: Double) {
        let w1 = omega(k1), w2 = omega(k2)
        let a1 = k1 * x - w1 * t          // e^{+i a1}
        let a2 = k2 * x + w2 * t          // e^{−i a2}
        let re = cos(a1) + beta * cos(a2)
        let im = sin(a1) - beta * sin(a2)
        return (re, im)
    }

    /// 密度按定义计算（脚本 density，避免手推符号错误）：
    /// ρ = (i/2)(φ* ∂_t φ − φ ∂_t φ*)，∂_t φ = (−iω₁)e^{i(k₁x−ω₁t)} + (iβω₂)e^{−i(k₂x+ω₂t)}
    /// （脚本 dphi_dt 逐句：第二项 +iβω₂ 与解析式自洽，见 fixture 恒等断言）。
    static func density(x: Double, t: Double, k1: Double, k2: Double, beta: Double) -> Double {
        let w1 = omega(k1), w2 = omega(k2)
        let a1 = k1 * x - w1 * t
        let a2 = k2 * x + w2 * t
        // φ 与 ∂_t φ 的 (re, im)
        let phiRe = cos(a1) + beta * cos(a2)
        let phiIm = sin(a1) - beta * sin(a2)
        // (−iω₁)e^{ia1}：Re = ω₁ sin a1, Im = −ω₁ cos a1；(iβω₂)e^{−ia2}：Re = βω₂ sin a2, Im = βω₂ cos a2
        let dphiRe = w1 * sin(a1) + beta * w2 * sin(a2)
        let dphiIm = -w1 * cos(a1) + beta * w2 * cos(a2)
        // z = φ*·dφ：Im(z) = φRe·dφIm − φIm·dφRe
        let cross = phiRe * dphiIm - phiIm * dphiRe
        // φ·dφ* = conj(z) → ρ = (i/2)(z − z̄) = −Im(z)
        return -cross
    }

    /// 解析形式（验证用，脚本 density_analytic）：
    /// ρ = ω₁ − ω₂β² − β(ω₂ − ω₁)cos((k₁+k₂)x + (ω₂ − ω₁)t)。
    static func densityAnalytic(x: Double, t: Double, k1: Double, k2: Double, beta: Double)
        -> Double {
        let w1 = omega(k1), w2 = omega(k2)
        return w1 - w2 * beta * beta
            - beta * (w2 - w1) * cos((k1 + k2) * x + (w2 - w1) * t)
    }

    /// 空间周期 Lx = 2π/(k₁+k₂)。
    static func spatialPeriod(k1: Double, k2: Double) -> Double {
        2.0 * .pi / (k1 + k2)
    }

    /// 时间周期（拍频）Lt = 2π/|ω₂ − ω₁|。
    static func temporalPeriod(k1: Double, k2: Double) -> Double {
        2.0 * .pi / abs(omega(k2) - omega(k1))
    }
}

// MARK: - 模块

/// 笔记 37《克莱因-戈登方程》：双动量叠加的守恒概率密度 ρ(x,t)——
/// 400×600 等高线（W9 大网格降采样基建首用）+ 双时刻空间剖面 + ρ(0,t) 时序，
/// 演示 KG 概率密度非正定困难。秒级档。
struct KleinGordonDensityModule: SimModule {

    /// 网格口径与脚本一致：x 400 点 × t 600 点。
    static let nx = 400
    static let nt = 600

    let meta = ModuleMeta(
        id: "克莱因-戈登方程_概率密度振荡",
        title: "克莱因-戈登方程 · 概率密度振荡",
        subtitle: "双动量叠加 ρ(x,t)：KG 概率密度非正定困难",
        category: .relativisticQFT, noteNumber: 37, tier: .seconds, difficulty: .basic,
        keywords: ["概率密度振荡", "Klein-Gordon", "KG方程", "负能", "拍频",
                   "非正定", "守恒流", "相对论量子"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "beta", title: "负能分量振幅", symbol: "β", unit: "",
                               range: 0.1...1.2, defaultValue: 0.9, decimalPlaces: 2)),
            .slider(SliderSpec(key: "k2", title: "负能动量", symbol: "k₂", unit: "（自然单位）",
                               range: 1...5, defaultValue: 3, decimalPlaces: 1)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .contour(ContourSpec(
                title: "KG 概率密度 ρ(x,t)（一个时间周期）",
                xAxis: .init(label: "x（自然单位，一个空间周期）"),
                yAxis: .init(label: "t（自然单位，一个时间周期）"),
                valueLabel: "ρ(x,t)", levels: 12, diverging: true)),
            .lineSeries(LineSeriesSpec(
                title: "KG 概率密度：双时刻空间剖面",
                xAxis: .init(label: "x（自然单位）"),
                yAxis: .init(label: "ρ(x,t)"),
                seriesNames: ["t = 0", "t = T/2"])),
            .lineSeries(LineSeriesSpec(
                title: "KG 概率密度 x = 0：振荡进入负值",
                xAxis: .init(label: "t（自然单位，一个时间周期）"),
                yAxis: .init(label: "ρ(0,t)"),
                seriesNames: ["ρ(0,t)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let beta = input.slider("beta")
        let k2 = input.slider("k2")
        let k1 = 1.0   // 脚本口径：正能动量固定 k₁ = 1

        let w1 = KleinGordonDensityMath.omega(k1)
        let w2 = KleinGordonDensityMath.omega(k2)
        let lx = KleinGordonDensityMath.spatialPeriod(k1: k1, k2: k2)
        let lt = KleinGordonDensityMath.temporalPeriod(k1: k1, k2: k2)

        // --- 秒级档调度通道（400×600 = 24 万点解析求值） ---
        let progress = ComputeProgress()
        progress.update(0.05, phase: "生成网格与解析场")

        let payload = try await SecondsChannel.run(progress: progress) {
            () -> (rho: [[Double]], rhoMin: Double, rhoMax: Double, fracPos: Double,
                   maxDefDiff: Double, spatialMeans: [Double]) in
            try Task.checkCancellation()
            progress.update(0.3, phase: "ρ(x,t) 定义式 vs 解析式全网格互验")

            let x = Num.linspace(0, lx, count: Self.nx)
            let t = Num.linspace(0, lt, count: Self.nt)

            // 网格 ρ（定义式）+ 与解析式的最大偏差（脚本 assert < 1e-12 口径）
            var rho: [[Double]] = []
            rho.reserveCapacity(Self.nt)
            var maxDefDiff = 0.0
            var fracPos = 0
            var rhoMin = Double.infinity
            var rhoMax = -Double.infinity
            for ti in 0..<Self.nt {
                var row: [Double] = []
                row.reserveCapacity(Self.nx)
                for xi in 0..<Self.nx {
                    let byDef = KleinGordonDensityMath.density(
                        x: x[xi], t: t[ti], k1: k1, k2: k2, beta: beta)
                    let analytic = KleinGordonDensityMath.densityAnalytic(
                        x: x[xi], t: t[ti], k1: k1, k2: k2, beta: beta)
                    maxDefDiff = max(maxDefDiff, abs(byDef - analytic))
                    row.append(byDef)
                    if byDef > 0 { fracPos += 1 }
                    rhoMin = min(rhoMin, byDef)
                    rhoMax = max(rhoMax, byDef)
                }
                rho.append(row)
            }
            let fracPosD = Double(fracPos) / Double(Self.nx * Self.nt)

            // 物理律：空间平均 ∫ρdx 不随 t（梯形积分，取 5 个时刻）
            var spatialMeans: [Double] = []
            for frac in [0.0, 0.25, 0.5, 0.75, 1.0] {
                let ti = min(Int(frac * Double(Self.nt - 1)), Self.nt - 1)
                spatialMeans.append(Num.trapezoid(x, rho[ti]))
            }
            progress.update(0.85, phase: "剖面与时序曲线")
            return (rho, rhoMin, rhoMax, fracPosD, maxDefDiff, spatialMeans)
        }

        // --- 图 1：400×600 等高线（max(nx,ny)=600 > 384 → 渲染层自动降采样） ---
        let xGrid = Num.linspace(0, lx, count: Self.nx)
        let tGrid = Num.linspace(0, lt, count: Self.nt)
        let contour = ContourData(spec: charts[0].contourSpec ?? ContourSpec(
            xAxis: .init(label: "x"), yAxis: .init(label: "t"), valueLabel: "ρ"),
            xGrid: xGrid, yGrid: tGrid, values: payload.rho)

        // --- 图 2：双时刻空间剖面（脚本 fig ax1：t=0 与 t=T/2，400 点） ---
        let prof0 = payload.rho[0]
        let profHalf = payload.rho[Self.nt / 2]
        let chart2 = LineSeriesData(
            spec: .init(xAxis: .init(label: "x（自然单位）"),
                        yAxis: .init(label: "ρ(x,t)"),
                        seriesNames: ["t = 0", "t = T/2"]),
            series: [
                .init(name: "t = 0", points: Num.strided(xGrid, prof0, stride: 1)),
                .init(name: "t = T/2", points: Num.strided(xGrid, profHalf, stride: 1)),
            ],
            referenceLines: [ReferenceLine(label: "ρ = 0", axis: .y, value: 0, style: .subtle)])

        // --- 图 3：ρ(0,t) 时序（脚本 fig ax2：400 点） ---
        let tt = Num.linspace(0, lt, count: 400)
        let rho0 = tt.map {
            KleinGordonDensityMath.density(x: 0, t: $0, k1: k1, k2: k2, beta: beta)
        }
        let chart3 = LineSeriesData(
            spec: .init(xAxis: .init(label: "t（自然单位，一个时间周期）"),
                        yAxis: .init(label: "ρ(0,t)"),
                        seriesNames: ["ρ(0,t)"]),
            series: [.init(name: "ρ(0,t)", points: Num.strided(tt, rho0, stride: 1))],
            referenceLines: [ReferenceLine(label: "ρ = 0", axis: .y, value: 0, style: .subtle)])

        return SimResult(
            charts: [.contour(contour), .lineSeries(chart2), .lineSeries(chart3)],
            summary: [
                .init(id: "params", title: "模式参数",
                      value: String(format: "k₁=%.1f, k₂=%.1f, β=%.2f", k1, k2, beta),
                      note: String(format: "ω₁ = %.4f, ω₂ = %.4f（自然单位 m = 1）", w1, w2)),
                .init(id: "period", title: "拍频周期",
                      value: String(format: "T = 2π/|ω₂−ω₁| = %.4f", lt),
                      note: String(format: "空间周期 Lx = 2π/(k₁+k₂) = %.4f", lx)),
                .init(id: "nonpositive", title: "非正定困难",
                      value: String(format: "ρ_min = %.4f, ρ_max = %.4f",
                                    payload.rhoMin, payload.rhoMax),
                      note: String(format: "ρ>0 网格占比 %.3f：密度可正可负（单粒子诠释失效）",
                                   payload.fracPos)),
                .init(id: "identity", title: "定义式 ⇔ 解析式",
                      value: String(format: "max|Δ| = %.1e", payload.maxDefDiff),
                      note: "ρ = (i/2)(φ*∂_tφ − φ∂_tφ*) 与闭式逐点一致（全网格）"),
                .init(id: "conservation", title: "空间平均守恒",
                      value: String(format: "∫ρdx = %.4f ± %.1e",
                                    payload.spatialMeans[2],
                                    payload.spatialMeans.map { abs($0 - payload.spatialMeans[2]) }
                                        .max() ?? 0),
                      note: "全时段不变：KG 守恒流的全域电荷（正负抵消）"),
            ],
            theory: TheoryCard(
                title: "KG 守恒概率密度与非正定困难",
                formulas: [
                    "ρ = (iħ/2mc²)(φ*∂_tφ − φ∂_tφ*)，自然单位 ρ = (i/2)(φ*∂_tφ − φ∂_tφ*)",
                    "双动量叠加：ρ = ω₁ − ω₂β² − β(ω₂−ω₁)cos((k₁+k₂)x + (ω₂−ω₁)t)",
                    "拍频周期 T = 2π/|ω₂ − ω₁|：正负能分量的干涉拍",
                    "ρ 可正可负 → 单粒子概率诠释失效 → 多粒子场论的入口",
                ],
                reading: "负能分量（振幅 β）与正能分量的干涉给出空间拍频振荡：等高线图里"
                    + "正负区带以 T = 2π/|ω₂−ω₁| 为周期整体翻转。密度穿过零线进入负值——"
                    + "KG 的『概率密度』不保号，这是把 φ 重新诠释为场算符、走向"
                    + "量子场论的历史动因。空间平均 ∫ρdx 全时段恒定：守恒的是"
                    + "电荷（正负贡献相消后的净值），不是概率。"))
    }
}
