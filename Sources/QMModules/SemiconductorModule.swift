import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/半导体_载流子与pn结.py）

/// Si 本征载流子 + 理想二极管的纯函数核。
enum SemiconductorMath {

    /// 有效态密度 N_c/N_v 随 T^1.5 标度（300 K 参考值），cm⁻³。
    static func effectiveDOS(_ ref300: Double, T: Double) -> Double {
        ref300 * pow(T / 300, 1.5)
    }

    /// 本征载流子浓度 n_i = √(N_c N_v)·exp(−E_g/2k_BT)，cm⁻³（质量作用定律）。
    static func ni(_ T: Double, eg: Double, kB: Double) -> Double {
        let nc = effectiveDOS(2.8e19, T: T)
        let nv = effectiveDOS(1.04e19, T: T)
        return sqrt(nc * nv) * exp(-eg / (2 * kB * T))
    }

    /// Shockley 二极管电流 I = I_s(e^(V/kT) − 1)，A。
    static func diodeCurrent(_ V: Double, isat: Double, vT: Double) -> Double {
        isat * (exp(V / vT) - 1)
    }

    /// 热电压 k_BT/e（V）。
    static func thermalVoltage(T: Double, kB: Double, e: Double) -> Double {
        kB * T / e
    }
}

// MARK: - 模块

/// 笔记 25《半导体与 pn 结》：Si 的 n_i(T)（semilogy）+ 理想二极管
/// Shockley I-V（semilogy）。实时档。
struct SemiconductorModule: SimModule {

    let meta = ModuleMeta(
        id: "半导体_载流子与pn结", title: "半导体 · 载流子与 pn 结",
        subtitle: "本征载流子 n_i(T) + Shockley 二极管 I-V",
        category: .manyBody, noteNumber: 25, tier: .realtime, difficulty: .basic,
        keywords: ["半导体", "pn结", "二极管", "本征", "载流子", "Shockley",
                   "肖克利", "禁带", "带隙", "质量作用定律", "热电压"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "Eg", title: "禁带宽度", symbol: "E_g", unit: "eV",
                               range: 0.6...1.5, defaultValue: 1.12, decimalPlaces: 2)),
            .slider(SliderSpec(key: "T", title: "标记温度", symbol: "T", unit: "K",
                               range: 100...600, defaultValue: 300, decimalPlaces: 0)),
            .slider(SliderSpec(key: "Is", title: "反向饱和电流", symbol: "I_s", unit: "A",
                               range: 1e-13...1e-9, defaultValue: 1e-12,
                               scale: .log, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "Si 本征载流子浓度随温度（质量作用定律）",
                xAxis: .init(label: "温度 T (K)"),
                yAxis: .init(label: "n_i (cm⁻³)", scale: .log),
                seriesNames: ["n_i(T)"])),
            .lineSeries(LineSeriesSpec(
                title: "理想 pn 结 I-V（Shockley 方程）",
                xAxis: .init(label: "外加电压 V (V)"),
                yAxis: .init(label: "|I| (A)", scale: .log),
                seriesNames: ["|I(V)|"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let eg_eV = input.slider("Eg")
        let tMark = input.slider("T")
        let iSat = input.slider("Is")
        let kB = try constants.value("kB")
        let eV = try constants.value("eV")
        let eg = eg_eV * eV

        // ---- 段 1：n_i(T)，T = linspace(100, 600, 400)（与 Python 一致）----
        let Tgrid = Num.linspace(100, 600, count: 400)
        let niVals = Tgrid.map { SemiconductorMath.ni($0, eg: eg, kB: kB) }

        // ---- 段 2：二极管 I-V，V = linspace(−1, 0.85, 900)，T = 300 K ----
        let vT = SemiconductorMath.thermalVoltage(T: 300, kB: kB, e: eV)
        let V = Num.linspace(-1, 0.85, count: 900)
        let I = V.map { SemiconductorMath.diodeCurrent($0, isat: iSat, vT: vT) }
        let iAt07 = SemiconductorMath.diodeCurrent(0.7, isat: iSat, vT: vT)

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "温度 T (K)"),
                        yAxis: .init(label: "n_i (cm⁻³)", scale: .log),
                        seriesNames: ["n_i(T)"]),
            series: [.init(name: "n_i(T)",
                           points: zip(Tgrid, niVals).map { Point(x: $0, y: $1) })],
            referenceLines: [ReferenceLine(
                label: String(format: "T = %.0f K", tMark), axis: .x, value: tMark)])

        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "外加电压 V (V)"),
                        yAxis: .init(label: "|I| (A)", scale: .log),
                        seriesNames: ["|I(V)|"]),
            series: [.init(name: "|I(V)|",
                           points: Num.strided(V, I.map { abs($0) }, stride: 2))],
            referenceLines: [
                ReferenceLine(label: "开启 ~0.7 V", axis: .x, value: 0.7),
                ReferenceLine(label: String(format: "I_s = %.0e A", iSat),
                              axis: .y, value: iSat, style: .subtle),
            ])

        let niMark = SemiconductorMath.ni(tMark, eg: eg, kB: kB)
        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chart1)],
            summary: [
                .init(id: "ni", title: String(format: "n_i(%.0f K)", tMark),
                      value: String(format: "%.3e cm⁻³", niMark),
                      note: "教科书 Si 值 ~1e10 cm⁻³ @300 K"),
                .init(id: "vt", title: "热电压 kT/e",
                      value: String(format: "%.2f mV", vT * 1000),
                      note: "300 K"),
                .init(id: "ion", title: "开启电流 I(0.7 V)",
                      value: String(format: "%.3e A", iAt07),
                      note: "e^(0.7/kT) ≈ 5.7×10¹¹"),
            ],
            theory: TheoryCard(
                title: "质量作用定律与 Shockley 方程",
                formulas: [
                    "n_i(T) = √(N_cN_v)·e^(−E_g/2k_BT)，N_c/N_v ∝ T^(3/2)",
                    "E_g 指数主导：每升温 50 K，n_i 约翻数倍",
                    "I(V) = I_s(e^(V/k_BT/e) − 1)：正向指数、反向饱和 −I_s",
                    "热电压 k_BT/e ≈ 25.85 mV @300 K（开启 ~0.7 V ≈ 27 倍）",
                ],
                reading: "n_i 曲线在 semilog 下近直线——斜率即 E_g；二极管每升高"
                    + "热电压电流翻 e 倍，0.7 V 处电流暴涨 11 个量级即整流本质。"))
    }
}
