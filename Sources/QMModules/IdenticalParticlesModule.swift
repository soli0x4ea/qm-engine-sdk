import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/全同粒子_费米气体与HOM.py）

/// 自由电子费米气 + HOM 双光子干涉的纯函数核。
enum IdenticalParticlesMath {

    /// 费米波矢 k_F = (3π²n)^(1/3)，m⁻¹。
    static func kF(_ n: Double) -> Double {
        pow(3 * .pi * .pi * n, 1.0 / 3)
    }

    /// 费米能量 E_F = ℏ²k_F²/2m_e，J。
    static func eF(_ n: Double, hbar: Double, me: Double) -> Double {
        let k = kF(n)
        return hbar * hbar * k * k / (2 * me)
    }

    /// 态密度 D(E) = (2m)^(3/2)/(2π²ℏ³) · E^(1/2)，J⁻¹m⁻³。
    static func dos(_ E: Double, hbar: Double, me: Double) -> Double {
        let pref = pow(2 * me, 1.5) / (2 * .pi * .pi * pow(hbar, 3))
        return pref * sqrt(E)
    }

    /// HOM 符合率 R_c(τ) = ½(1 − |g₁(τ)|²)，g₁ = exp(−σ_ω²τ²/2)。
    static func homRate(_ tau: Double, sigmaW: Double) -> Double {
        let g1 = exp(-sigmaW * sigmaW * tau * tau / 2)
        return 0.5 * (1 - g1 * g1)
    }

    /// 可见度 V = (max−min)/(max+min)。
    static func visibility(_ rc: [Double]) -> Double {
        guard let mx = rc.max(), let mn = rc.min(), mx + mn > 0 else { return 0 }
        return (mx - mn) / (mx + mn)
    }
}

// MARK: - 模块

/// 笔记 21《全同粒子与泡利不相容原理》：铜的费米气体（D(E)∝E^½）
/// + Hong-Ou-Mandel 双光子干涉（全同 V=1 vs 非全同 V=0）。实时档。
struct IdenticalParticlesModule: SimModule {

    let meta = ModuleMeta(
        id: "全同粒子_费米气体与HOM", title: "全同粒子 · 费米气体与 HOM",
        subtitle: "态密度 D(E)∝E^(1/2) + HOM 干涉凹陷",
        category: .manyBody, noteNumber: 21, tier: .realtime, difficulty: .basic,
        keywords: ["全同粒子", "泡利", "Pauli", "费米", "Fermi", "费米能", "态密度",
                   "HOM", "Hong-Ou-Mandel", "双光子干涉", "可见度"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "n", title: "电子数密度", symbol: "n", unit: "m⁻³",
                               range: 1e27...1e30, defaultValue: 8.49e28,
                               scale: .log, decimalPlaces: 2)),
            .slider(SliderSpec(key: "sigma_w", title: "光谱宽度", symbol: "σ_ω",
                               unit: "rad/s", range: 1e10...1e13, defaultValue: 1e12,
                               scale: .log, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "自由电子态密度与费米能",
                xAxis: .init(label: "E (eV)"),
                yAxis: .init(label: "D(E) (态/J/m³)"),   // W13d：回线性——D(E) 仅跨 1 个量级，log 收益小且与笔记 21 图不符
                seriesNames: ["D(E)∝E^(1/2)"])),
            .lineSeries(LineSeriesSpec(
                title: "HOM 符合率：全同 vs 非全同",
                xAxis: .init(label: "延时 τ (ps)"),
                yAxis: .init(label: "符合率 R_c"),
                seriesNames: ["全同光子", "非全同光子"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let n = input.slider("n")
        let sigmaW = input.slider("sigma_w")
        let hbar = try constants.value("hbar")
        let me = try constants.value("m_e")
        let eV = try constants.value("eV")

        // ---- 段 1：费米气体（与 Python E = linspace(1e-3, 3EF, 600) 一致）----
        let kF = IdenticalParticlesMath.kF(n)
        let eF = IdenticalParticlesMath.eF(n, hbar: hbar, me: me)
        let tF = eF / (try constants.value("kB"))
        let vF = hbar * kF / me

        let E = Num.linspace(1e-3, 3 * eF, count: 600)
        let D = E.map { IdenticalParticlesMath.dos($0, hbar: hbar, me: me) }
        let eFAbscissa = eF / eV

        // ---- 段 2：HOM 干涉（τ = linspace(−4e-12, 4e-12, 800)，ps 表示）----
        let tau = Num.linspace(-4e-12, 4e-12, count: 800)
        let rcIdent = tau.map { IdenticalParticlesMath.homRate($0, sigmaW: sigmaW) }
        let rcNon = tau.map { _ in 0.5 }
        let vIdent = IdenticalParticlesMath.visibility(rcIdent)

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "E (eV)"),
                        yAxis: .init(label: "D(E) (态/J/m³)"),
                        seriesNames: ["D(E)∝E^(1/2)"]),
            series: [.init(name: "D(E)∝E^(1/2)",
                           points: Num.strided(E.map { $0 / eV }, D.map { $0 * eV },
                                               stride: 2))],
            referenceLines: [ReferenceLine(
                id: "ef", label: String(format: "E_F = %.2f eV", eFAbscissa),
                axis: .x, value: eFAbscissa)])

        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "延时 τ (ps)"),
                        yAxis: .init(label: "符合率 R_c"),
                        seriesNames: ["全同光子", "非全同光子"]),
            series: [
                .init(name: String(format: "全同光子 (V=%.3f)", vIdent),
                      points: Num.strided(tau.map { $0 * 1e12 }, rcIdent, stride: 2)),
                .init(name: "非全同光子 (V=0)",
                      points: Num.strided(tau.map { $0 * 1e12 }, rcNon, stride: 2)),
            ],
            referenceLines: [ReferenceLine(label: "τ = 0", axis: .x, value: 0, style: .subtle)])

        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chart1)],
            summary: [
                .init(id: "kf", title: "k_F", value: String(format: "%.4e m⁻¹", kF),
                      note: "(3π²n)^(1/3)"),
                .init(id: "ef", title: "E_F", value: String(format: "%.4f eV", eF / eV),
                      note: "费米能量"),
                .init(id: "tf", title: "T_F", value: String(format: "%.4e K", tF),
                      note: "简并温度 E_F/k_B"),
                .init(id: "vf", title: "v_F", value: String(format: "%.4e m/s", vF),
                      note: "费米速度 ℏk_F/m"),
                .init(id: "vhom", title: "HOM 可见度",
                      value: String(format: "%.6f", vIdent),
                      note: "全同光子 V→1，非全同 V=0"),
            ],
            theory: TheoryCard(
                title: "全同粒子：费米统计与双光子干涉",
                formulas: [
                    "k_F = (3π²n)^(1/3)；E_F = ℏ²k_F²/2m；T_F = E_F/k_B",
                    "D(E) = (2m)^(3/2)/(2π²ℏ³) · E^(1/2)（三维自由电子气）",
                    "R_c(τ) = ½(1 − |g₁(τ)|²)，g₁(τ) = exp(−σ_ω²τ²/2)",
                    "HOM：全同光子在 τ=0 处符合率塌缩（聚束），V=1；非全同恒 ½",
                ],
                reading: "铜中 E_F≈7 eV 远高于室温 k_BT≈0.026 eV——电子气始终简并；"
                    + "HOM 凹陷宽度 ~1/σ_ω 即单光子相干时间，全同性是干涉的门票。"))
    }
}
