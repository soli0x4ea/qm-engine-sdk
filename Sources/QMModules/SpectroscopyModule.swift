import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/45_光谱学仪器的量子基础_荧光与寿命.py / _拉曼谱.py）

/// 光谱学仪器纯函数核：爱因斯坦 A 系数、荧光/磷光衰减、拉曼 Stokes/anti-Stokes 谱。
enum SpectroscopyMath {
    /// 1 Debye = 3.33564095×10⁻³⁰ C·m
    static let debyeCm: Double = 3.33564095e-30
    /// h·c（eV·nm），CODATA 精确派生量
    static let hcEVnm = Num.hcEVnm

    /// 爱因斯坦 A 系数（SI）：A = ω³ |d|² / (3π ε₀ ħ c³)
    static func einsteinA(wavelength_m: Double, d_Cm: Double,
                          h: Double, eps0: Double, hbar: Double, c: Double) -> Double {
        let omega = 2.0 * .pi * c / wavelength_m
        return omega * omega * omega * d_Cm * d_Cm / (3.0 * .pi * eps0 * hbar * c * c * c)
    }

    /// 振幅型洛伦兹（峰值幅度 = amp，参数 gamma = FWHM）：L = amp·(γ/2)²/((x−x0)²+(γ/2)²)
    static func lorentzAmp(_ x: Double, _ x0: Double, _ gamma: Double, _ amp: Double) -> Double {
        let h = gamma / 2.0
        return amp * (h * h) / ((x - x0) * (x - x0) + h * h)
    }

    /// 拉曼谱：横轴为散射光波数 (cm⁻¹)，含各模式 Stokes(ν0−νv) 与 anti-Stokes(ν0+νv) 洛伦兹线。
    /// 强度标度 (ν_s)⁴·(n_v+1)（Stokes）/ (ν_s)⁴·n_v（anti-Stokes），n_v 为玻色占据。
    static func ramanSpectrum(_ xs: [Double], nu0: Double, c_cm: Double,
                              hbar: Double, kB: Double, T: Double,
                              modes: [(nuV: Double, fwhm: Double)]) -> [Double] {
        var spec = [Double](repeating: 0.0, count: xs.count)
        for i in xs.indices {
            var s = 0.0
            for m in modes {
                let omegaV = 2.0 * .pi * c_cm * m.nuV
                let nv = 1.0 / (exp(hbar * omegaV / (kB * T)) - 1.0)
                let nuS = nu0 - m.nuV
                let nuAS = nu0 + m.nuV
                let aS = pow(nuS, 4) * (nv + 1.0)
                let aAS = pow(nuAS, 4) * nv
                s += lorentzAmp(xs[i], nuS, m.fwhm, aS) + lorentzAmp(xs[i], nuAS, m.fwhm, aAS)
            }
            spec[i] = s
        }
        return spec
    }
}

// MARK: - 模块一：荧光与寿命

/// 笔记 45《光谱学仪器的量子基础》：爱因斯坦 A 系数、荧光/磷光双指数衰减、Jablonski 图。
struct FluorescenceLifetimeModule: SimModule {

    let meta = ModuleMeta(
        id: "光谱学仪器的量子基础_荧光与寿命", title: "荧光与寿命 · A 系数 + 双指数衰减",
        subtitle: "A = ω³d²/(3πε₀ħc³)；τ = 1/A；荧光(ns) 与磷光(ms) 双指数衰减 + Jablonski 图",
        category: .gemology, noteNumber: 45, tier: .realtime, difficulty: .advanced,
        keywords: ["荧光", "磷光", "爱因斯坦A系数", "寿命", "Jablonski", "自发辐射",
                   "fluorescence", "phosphorescence", "Einstein A", "lifetime", "Jablonski diagram"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "dipole", title: "跃迁偶极矩", symbol: "|d|", unit: "D",
                              range: 0.1...3.0, defaultValue: 1.0,
                              scale: .linear, decimalPlaces: 2)),
            .slider(SliderSpec(key: "tauPhos", title: "磷光寿命 τ_phos", symbol: "τ_phos", unit: "ms",
                              range: 0.1...20.0, defaultValue: 3.0,
                              scale: .linear, decimalPlaces: 1)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .dualAxisLineSeries(DualAxisLineSeriesSpec(
                title: "爱因斯坦 A 系数与寿命随发射波长变化（对数轴）",
                xAxis: .init(label: "发射波长 λ (nm)"),
                primaryAxis: .init(label: "自发辐射速率 A (s⁻¹，log)", scale: .log),
                secondaryAxis: .init(label: "寿命 τ (ns，log)", scale: .log),
                primaryNames: ["A(λ)"], secondaryNames: ["τ(λ)"])),
            .lineSeries(LineSeriesSpec(
                title: "荧光与磷光：指数衰减（双指数分量）",
                xAxis: .init(label: "时间 t (s)", scale: .log),
                yAxis: .init(label: "归一化强度"),
                seriesNames: [
                    "fluorescence (tau~1.6 ns, allowed E1)",
                    "phosphorescence (tau~3 ms, spin-forbidden e.g. ruby Cr³⁺)",
                ])),
            .schematic(SchematicSpec(
                title: "Jablonski 图",
                xAxis: .init(label: ""),
                yAxis: .init(label: "能量 E (eV)"))),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let h = try constants.value("h")
        let hbar = try constants.value("hbar")
        let c = try constants.value("c")
        let eps0 = try constants.value("epsilon_0")
        let e = try constants.value("e")
        let a0 = try constants.value("a0")

        // 标定：H 2p→1s（Lyman-α），d_H = 0.7449 e a0
        let dH = 0.7449 * e * a0
        let lamH = 121.567e-9
        let aH = SpectroscopyMath.einsteinA(wavelength_m: lamH, d_Cm: dH,
                                            h: h, eps0: eps0, hbar: hbar, c: c)
        let tauH = 1.0 / aH

        // 图 1：A(λ) 与 τ(λ)（双轴）
        let d = input.slider("dipole") * SpectroscopyMath.debyeCm
        let wl = Num.linspace(300.0e-9, 800.0e-9, count: 300)
        let aWl = wl.map { SpectroscopyMath.einsteinA(wavelength_m: $0, d_Cm: d,
                                                     h: h, eps0: eps0, hbar: hbar, c: c) }
        let tauWl = aWl.map { 1.0 / $0 }

        let wlNm = wl.map { $0 * 1.0e9 }

        // 图 2：荧光 vs 磷光衰减
        let tauPhos = input.slider("tauPhos") * 1.0e-3
        let t = Num.linspace(-11.0, -2.0, count: 400).map { pow(10.0, $0) }
        let iFluor = t.map { exp(-$0 / tauH) }
        let iPhos = t.map { exp(-$0 / tauPhos) }

        // 图 3：Jablonski 图（示意）
        let levels: [(String, Double)] = [("S0", 0.0), ("S1", 3.0), ("T1", 1.6)]
        let refLines = levels.map {
            ReferenceLine(label: $0.0, axis: .y, value: $0.1, style: .subtle)
        }
        let callouts: [SchematicCallout] = [
            .init(text: "absorption S0→S1", anchor: .init(x: 1.0, y: 1.5),
                 arrowEnd: .init(x: 1.0, y: 3.0)),
            .init(text: "fluorescence S1→S0 (ns)", anchor: .init(x: 3.0, y: 1.5),
                 arrowEnd: .init(x: 3.0, y: 0.1)),
            .init(text: "ISC S1→T1", anchor: .init(x: 2.0, y: 2.3),
                 arrowEnd: .init(x: 2.0, y: 1.7)),
            .init(text: "phosphorescence T1→S0 (ms)", anchor: .init(x: 0.5, y: 0.8),
                 arrowEnd: .init(x: 0.5, y: 0.1)),
        ]

        return SimResult(
            charts: [
                .dualAxisLineSeries(DualAxisLineSeriesData(
                    spec: charts.requireDualAxisLineSeries(0),
                    primary: [.init(name: "A(λ)",
                                   points: Num.strided(wlNm, aWl, stride: 1))],
                    secondary: [.init(name: "τ(λ)",
                                     points: Num.strided(wlNm, tauWl.map { $0 * 1e9 }, stride: 1),
                                     colorIndex: 1)])),
                .lineSeries(LineSeriesData(
                    spec: charts.requireLineSeries(1),
                    series: [
                        .init(name: "fluorescence (tau~1.6 ns, allowed E1)",
                              points: Num.strided(t, iFluor, stride: 1)),
                        .init(name: "phosphorescence (tau~3 ms, spin-forbidden e.g. ruby Cr³⁺)",
                              points: Num.strided(t, iPhos, stride: 1), colorIndex: 1),
                    ])),
                .schematic(SchematicData(
                    spec: charts.requireSchematic(2),
                    callouts: callouts, referenceLines: refLines)),
            ],
            summary: [
                .init(id: "AH", title: "H 2p→1s 标定 A",
                      value: String(format: "%.3e s⁻¹", aH),
                      note: "τ = 1/A ≈ \(String(format: "%.2f", tauH * 1e9)) ns"),
                .init(id: "tauFluor", title: "荧光寿命 τ_fluor",
                      value: String(format: "%.2f ns", tauH * 1e9),
                      note: "= 1/A_H（允许电偶极跃迁，ns 量级）"),
                .init(id: "tauPhos", title: "磷光寿命 τ_phos",
                      value: String(format: "%.1f ms", tauPhos * 1e3),
                      note: "自旋禁戒（如红宝石 Cr³⁺），ms 量级"),
                .init(id: "ratio", title: "τ_phos/τ_fluor",
                      value: String(format: "%.1e", tauPhos / tauH),
                      note: "磷光比荧光慢 ~6 个数量级"),
            ],
            theory: TheoryCard(
                title: "荧光与寿命：爱因斯坦 A 系数（笔记 45 第二节 2.4、第六节 6.2）",
                formulas: [
                    "A = ω³ |d|² / (3π ε₀ ħ c³)（SI 单位制自发辐射速率）",
                    "寿命 τ = 1/A；固定偶极下 A ∝ 1/λ³（短波辐射更快）",
                    "荧光（ns，允许 E1）与磷光（ms，自旋禁戒）为双指数衰减分量",
                ],
                reading: "左上为双轴图：蓝线 A(λ) 随波长变短（能量升高）急剧上升（A∝λ⁻³），"
                    + "红虚线为对应寿命 τ=1/A（右下轴，ns）。左下为荧光与磷光两条指数衰减——"
                    + "荧光 ns 级、磷光 ms 级，相差约 10⁶ 倍。右下 Jablonski 图给出吸收、内转换、"
                    + "荧光、系间窜越(ISC)、磷光的能级与跃迁关系。"))
    }
}

// MARK: - 模块二：拉曼谱

/// 笔记 45《光谱学仪器的量子基础》：拉曼谱（785 nm 激发）
/// Stokes(损失 ħω_v) 与 anti-Stokes(获得 ħω_v) 带；anti-Stokes 含玻尔兹曼布居 n_v。
struct RamanSpectrumModule: SimModule {

    let meta = ModuleMeta(
        id: "光谱学仪器的量子基础_拉曼谱", title: "拉曼谱 · Stokes/anti-Stokes",
        subtitle: "785 nm 激发：散射光子损失/获得 ħω_v；anti-Stokes 含玻尔兹曼布居 n_v",
        category: .gemology, noteNumber: 45, tier: .seconds, difficulty: .advanced,
        keywords: ["拉曼", "拉曼谱", "Stokes", "anti-Stokes", "玻尔兹曼", "布居", "金刚石",
                   "Raman", "Raman spectrum", "Stokes", "anti-Stokes", "Boltzmann", "diamond"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "lambda0", title: "激发波长", symbol: "λ₀", unit: "nm",
                              range: 400...1064, defaultValue: 785,
                              scale: .linear, decimalPlaces: 0)),
            .slider(SliderSpec(key: "T", title: "温度", symbol: "T", unit: "K",
                              range: 50...600, defaultValue: 300,
                              scale: .linear, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "拉曼谱：Stokes 与 anti-Stokes 带（785 nm 激发，对数强度轴）",
                xAxis: .init(label: "散射光波数 / cm⁻¹"),
                yAxis: .init(label: "相对强度 (a.u.，log)", scale: .log),
                seriesNames: ["Raman (Stokes + anti-Stokes)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let c = try constants.value("c")
        let hbar = try constants.value("hbar")
        let kB = try constants.value("kB")
        let db = try MaterialDB.load()

        let lambda0 = input.slider("lambda0") * 1.0e-9
        let T = input.slider("T")
        let c_cm = c * 100.0                         // cm/s
        let nu0 = 1.0 / lambda0 * 1.0e-2            // 入射激光波数 cm⁻¹

        let modes = db.ramanModes.map { (nuV: $0.nuVcm, fwhm: $0.fwhmCm) }
        let shiftMax = 3000.0
        // W13j #45：60001 → 8001 点（0.75 cm⁻¹ 步长，最窄 FWHM=6 仍有 8 采样），
        // 6 万 LineMark 在 iOS 真机渲染过重；峰形由 ChartKit 显式域保证满幅。
        // W13h #45：x 轴改「拉曼位移 Δν = ν散射 − ν₀」——绝对波数轴随 λ₀ 平移、
        // 刻度值无物理含义；位移轴下 Rayleigh 恒在 0、Stokes/anti-Stokes 恒在 ∓ν_v。
        // λ₀ 的影响保留在 ν⁴ 权重（anti/Stokes 强度比随 λ₀ 微变）。
        let shifts = Num.linspace(-shiftMax, shiftMax, count: 8001)
        let xsAbs = shifts.map { nu0 + $0 }
        let spec = SpectroscopyMath.ramanSpectrum(xsAbs, nu0: nu0, c_cm: c_cm,
                                                 hbar: hbar, kB: kB, T: T, modes: modes)

        // 金刚石 1332 cm⁻¹ 模式的玻色占据与 Stokes/anti-Stokes 强度比
        // （W8 收尾：消除 force-unwrap。行存在性由 MaterialDBTests /
        //   SpectroscopyModuleTests 双侧断言锚定，此处 guard 兜底防数据回退）
        guard let diamond = db.ramanMode(id: "diamond_raman") else {
            throw ComputeError.invalidParameter("MaterialDB 缺少 diamond_raman 行（数据回退）")
        }
        let omegaV = 2.0 * .pi * c_cm * diamond.nuVcm
        let nv = 1.0 / (exp(hbar * omegaV / (kB * T)) - 1.0)
        let boltzmann = exp(-hbar * omegaV / (kB * T))

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: .init(xAxis: .init(label: "拉曼位移 Δν = ν散射 − ν₀ (cm⁻¹)"),
                                yAxis: .init(label: "相对强度 (a.u.，log)", scale: .log),
                                seriesNames: ["Raman (Stokes + anti-Stokes)"]),
                    series: [
                        .init(name: "Raman (Stokes + anti-Stokes)",
                              points: Num.strided(shifts, spec, stride: 1)),
                    ],
                    referenceLines: [
                        ReferenceLine(label: "Rayleigh（弹性，Δν = 0）", axis: .x, value: 0, style: .subtle),
                        ReferenceLine(label: "金刚石 1332 cm⁻¹ Stokes", axis: .x,
                                     value: -diamond.nuVcm, style: .threshold),
                        ReferenceLine(label: "金刚石 1332 cm⁻¹ anti-Stokes", axis: .x,
                                     value: diamond.nuVcm, style: .threshold),
                    ])),
            ],
            summary: [
                .init(id: "nu0", title: "入射激光波数 ν₀",
                      value: String(format: "%.1f cm⁻¹", nu0),
                      note: "λ₀ = \(Int(input.slider("lambda0"))) nm"),
                .init(id: "nv", title: "金刚石 1332 cm⁻¹ 布居 n_v",
                      value: String(format: "%.4f", nv),
                      note: "Bose-Einstein：1/(e^{ħω/kT}−1)，T=\(Int(T)) K"),
                .init(id: "boltz", title: "anti-Stokes 布居因子 e^{−ħω/kT}",
                      value: String(format: "%.4f", boltzmann),
                      note: "≈ n_v/(n_v+1)：anti-Stokes 被热布居压低"),
                .init(id: "Tcheck", title: "玻尔兹曼一致性",
                      value: abs(nv/(nv+1) - boltzmann) < 1e-9 ? "✓" : "✗",
                      note: "n_v/(n_v+1) = e^{−ħω/kT}（Stokes/anti-Stokes 布居比）"),
            ],
            theory: TheoryCard(
                title: "拉曼光谱（笔记 45 第二节、第六节 6.1）",
                formulas: [
                    "Kramers-Heisenberg / Placzek 极化率选择定则",
                    "Stokes 散射：光子损失 ħω_v → 散射波数 ν₀−ν_v",
                    "anti-Stokes 散射：光子获得 ħω_v → 散射波数 ν₀+ν_v（含布居 n_v）",
                    "强度比（布居部分）= n_v/(n_v+1) = e^{−ħω/kT}（玻尔兹曼）",
                ],
                reading: "以 785 nm 激光激发，在瑞利线（弹性散射，ν₀）两侧出现拉曼位移带："
                    + "低频侧为 Stokes（分子获得振动能），高频侧为 anti-Stokes（分子损失振动能）。"
                    + "anti-Stokes 强度受玻尔兹曼布居 n_v 压制——常温下 n_v ≪ 1，故 anti-Stokes 远弱于 Stokes。"
                    + "金刚石一阶拉曼位移 1332 cm⁻¹ 是宝石学鉴别钻石与仿制品的公认标志（图中标注两条 1332 带）。"))
    }
}
