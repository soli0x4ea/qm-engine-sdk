import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子测量与退相干_退相干动力学.py）

/// Caldeira-Leggett 退相干时间尺度 + 两能级/指针基/猫态衰减的纯函数核。
enum DecoherenceDynamicsMath {

    /// 退相干时间 τ_D = ℏ²/(12π·η·k_BT·a³)（布朗模型，Δx = a）。
    static func tauD(_ a: Double, eta: Double, T: Double,
                     hbar: Double, kB: Double) -> Double {
        hbar * hbar / (12 * .pi * eta * kB * T * pow(a, 3))
    }

    /// 两能级：ρ_ee = e^(−t/T1)（式 6.1）。
    static func rhoEE(_ t: Double, gamma: Double) -> Double {
        exp(-gamma * t)
    }

    /// 两能级：|ρ_eg| = e^(−t/2T1)（式 6.2）——T2 = 2T1。
    static func rhoEG(_ t: Double, gamma: Double) -> Double {
        exp(-gamma * t / 2)
    }

    /// 指针基：位置叠加相干衰减核 exp(−(2x₀)²t/(2τ₀))。
    static func pointerCoherence(_ t: Double, x0: Double, tau0: Double) -> Double {
        exp(-pow(2 * x0, 2) * t / (2 * tau0))
    }

    /// 腔 QED 猫态 Wigner 负性 N(t) = e^(−Γ_cat·t)，Γ_cat = 2κ|α|²。
    static func catNegativity(_ t: Double, kappa: Double, alpha: Double) -> Double {
        exp(-2 * kappa * alpha * alpha * t)
    }
}

// MARK: - 模块

/// 笔记 31《量子测量与退相干》：Caldeira-Leggett 时间尺度（loglog）
/// + 两能级 T2=2T1 / 指针基 einselection / 腔 QED 猫态负性三子图。实时档。
struct DecoherenceDynamicsModule: SimModule {

    let meta = ModuleMeta(
        id: "量子测量与退相干_退相干动力学", title: "量子测量 · 退相干动力学",
        subtitle: "τ_D∝a⁻³ 时间尺度 + T2=2T1 / 指针基 / 猫态负性",
        category: .quantumInfo, noteNumber: 31, tier: .realtime, difficulty: .basic,
        keywords: ["退相干", "decoherence", "Caldeira-Leggett", "指针基", "einselection",
                   "猫态", "Schrodinger", "Wigner", "负性", "T1", "T2", "布朗"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "eta", title: "环境黏度", symbol: "η", unit: "Pa·s",
                               range: 1e-6...1e-3, defaultValue: 1.8e-5,
                               scale: .log, decimalPlaces: 2)),
            .slider(SliderSpec(key: "Tenv", title: "环境温度", symbol: "T", unit: "K",
                               range: 100...600, defaultValue: 300, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "退相干时间尺度 vs 物体尺寸（布朗模型）",
                xAxis: .init(label: "物体半径 a (m)", scale: .log),
                yAxis: .init(label: "退相干时间 τ_D (s)", scale: .log),
                seriesNames: ["τ_D(a) = ℏ²/(12πηk_BT·a³)"])),
            .lineSeries(LineSeriesSpec(
                title: "(a) 两能级衰减：T2 = 2T1",
                xAxis: .init(label: "时间 t (T1 = 1)"),
                yAxis: .init(label: "密度矩阵元"),
                seriesNames: ["ρ_ee(t)", "|ρ_eg(t)|"])),
            .lineSeries(LineSeriesSpec(
                title: "(b) 指针基（einselection）：环境耦合位置",
                xAxis: .init(label: "时间 t"),
                yAxis: .init(label: "非对角相干 |⟨x₀|ρ|−x₀⟩|"),
                seriesNames: ["位置叠加", "位置对角混态"])),
            .lineSeries(LineSeriesSpec(
                title: "(c) 腔 QED 猫态负性衰减",
                xAxis: .init(label: "时间 t (1/κ = 1)"),
                yAxis: .init(label: "Wigner 负性 N(t)"),
                seriesNames: ["|α|=1", "|α|=2", "|α|=3"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let eta = input.slider("eta")
        let tenv = input.slider("Tenv")
        let hbar = try constants.value("hbar")
        let kB = try constants.value("kB")

        // ---- 图 1：τ_D(a)，a = logspace(−9, −3, 400)（与 Python 一致）----
        let lnA = Num.linspace(-9, -3, count: 400)
        let aGrid = lnA.map { pow(10, $0) }
        let tau = aGrid.map {
            DecoherenceDynamicsMath.tauD($0, eta: eta, T: tenv, hbar: hbar, kB: kB)
        }

        // ---- 图 2(a)：两能级，t = linspace(0, 5, 400)，Γ = 1 ----
        let t1 = Num.linspace(0, 5, count: 400)
        let rhoEE = t1.map { DecoherenceDynamicsMath.rhoEE($0, gamma: 1) }
        let rhoEG = t1.map { DecoherenceDynamicsMath.rhoEG($0, gamma: 1) }

        // ---- 图 2(b)：指针基，tt = linspace(0, 4, 400)，x₀=1，τ₀=0.6 ----
        let t2 = Num.linspace(0, 4, count: 400)
        let cohSuper = t2.map {
            DecoherenceDynamicsMath.pointerCoherence($0, x0: 1, tau0: 0.6)
        }
        let cohMix = [Double](repeating: 0, count: 400)

        // ---- 图 2(c)：猫态负性，t = linspace(0, 5, 400)，κ = 1 ----
        let catSeries: [SeriesPoints] = [1.0, 2.0, 3.0].map { alpha in
            let neg = t1.map {
                DecoherenceDynamicsMath.catNegativity($0, kappa: 1, alpha: alpha)
            }
            return SeriesPoints(name: String(format: "|α|=%.0f, Γ_cat=%.0f",
                                             alpha, 2 * alpha * alpha),
                                points: zip(t1, neg).map { Point(x: $0, y: $1) })
        }

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "物体半径 a (m)", scale: .log),
                        yAxis: .init(label: "退相干时间 τ_D (s)", scale: .log),
                        seriesNames: ["τ_D(a) = ℏ²/(12πηk_BT·a³)"]),
            series: [.init(name: "τ_D(a) = ℏ²/(12πηk_BT·a³)",
                           points: zip(aGrid, tau).map { Point(x: $0, y: $1) })],
            referenceLines: [
                ReferenceLine(label: "1 s（人类）", axis: .y, value: 1, style: .subtle),
                ReferenceLine(label: "1 ms", axis: .y, value: 1e-3, style: .subtle),
                ReferenceLine(label: "1 ps", axis: .y, value: 1e-12, style: .subtle),
            ])

        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "时间 t (T1 = 1)"),
                        yAxis: .init(label: "密度矩阵元"),
                        seriesNames: ["ρ_ee(t)", "|ρ_eg(t)|"]),
            series: [
                .init(name: "ρ_ee(t) = e^(−t/T1)", points: zip(t1, rhoEE).map { Point(x: $0, y: $1) }),
                .init(name: "|ρ_eg(t)| = e^(−t/2T1)",
                      points: zip(t1, rhoEG).map { Point(x: $0, y: $1) }),
            ])

        let chart2 = LineSeriesData(
            spec: .init(xAxis: .init(label: "时间 t"),
                        yAxis: .init(label: "非对角相干 |⟨x₀|ρ|−x₀⟩|"),
                        seriesNames: ["位置叠加", "位置对角混态"]),
            series: [
                .init(name: "位置叠加 (|x₀⟩+|−x₀⟩)/√2",
                      points: zip(t2, cohSuper).map { Point(x: $0, y: $1) }),
                .init(name: "位置对角混态（经典记录）",
                      points: zip(t2, cohMix).map { Point(x: $0, y: $1) }),
            ])

        let chart3 = LineSeriesData(
            spec: .init(xAxis: .init(label: "时间 t (1/κ = 1)"),
                        yAxis: .init(label: "Wigner 负性 N(t)"),
                        seriesNames: catSeries.map(\.name)),
            series: catSeries)

        let tauDust = DecoherenceDynamicsMath.tauD(1e-6, eta: eta, T: tenv,
                                                   hbar: hbar, kB: kB)
        let tauNano = DecoherenceDynamicsMath.tauD(1e-9, eta: eta, T: tenv,
                                                  hbar: hbar, kB: kB)
        let tauMacro = DecoherenceDynamicsMath.tauD(1e-3, eta: eta, T: tenv,
                                                    hbar: hbar, kB: kB)
        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chart1),
                     .lineSeries(chart2), .lineSeries(chart3)],
            summary: [
                .init(id: "nano", title: "τ_D（1 nm）",
                      value: String(format: "%.2e s", tauNano),
                      note: "Δx = a 口径下已 < 飞秒"),
                .init(id: "dust", title: "τ_D（1 µm 尘埃）",
                      value: String(format: "%.2e s", tauDust),
                      note: "比 ns 级能量弛豫快 18 个量级——退相干远快于耗散"),
                .init(id: "macro", title: "τ_D（1 mm）",
                      value: String(format: "%.2e s", tauMacro),
                      note: "宏观叠加瞬间湮灭"),
            ],
            theory: TheoryCard(
                title: "退相干动力学：环境选择指针基",
                formulas: [
                    "τ_D = ℏ²/(12π·η·k_BT·a³)：Δx≈a 时退相干比耗散快 (a/λ_th)² 倍",
                    "两能级：ρ_ee = e^(−t/T1)，|ρ_eg| = e^(−t/2T1)——T2 = 2T1",
                    "指针基：退相干核 exp(−(x−x′)²t/2τ₀)——位置叠加湮灭，对角混态幸存",
                    "猫态：N(t) = e^(−Γ_cat·t)，Γ_cat = 2κ|α|²（Brune 1996 腔 QED）",
                ],
                reading: "τ_D 随 a⁻³ 暴跌：1 nm 粒子 ~1e-18 s、尘埃 ~1e-27 s、毫米物"
                    + " ~1e-36 s——「宏观叠加态从未有机会存在」；而猫态负性衰减率"
                    + " 随 |α|² 增大正是「越宏观越难看见猫」的定量表述。"))
    }
}
