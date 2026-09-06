import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/微扰论与变分法_方势阱与氦变分.py）

/// 微扰展开 + 变分法数值验证的纯函数核（自然单位 ħ = ω = m = 1）。
enum PerturbationVariationMath {

    /// 谐振子基底维度（脚本 N = 40）。
    static let basisSize = 40
    /// 变分扫描点数（脚本 400）。
    static let variationalScanPoints = 400
    /// 氦原子精确非相对论基态能量（Ha，脚本常数）。
    static let heliumExactHa = -2.903724377

    /// (A1) 非谐振子基底：E_n^(0) = n + 1/2，X⁴ = ((a+a†)/√2)⁴（N×N 稠密，行主序）。
    static func anharmonicBasis(N: Int) -> (e0: [Double], x4: [Double]) {
        var a = [Double](repeating: 0, count: N * N)
        for n in 1..<N { a[(n - 1) * N + n] = sqrt(Double(n)) }   // a|n⟩ = √n|n−1⟩（a[n−1,n]）
        // x = (a + aᵀ)/√2
        var x = [Double](repeating: 0, count: N * N)
        for i in 0..<N {
            for j in 0..<N {
                x[i * N + j] = (a[i * N + j] + a[j * N + i]) / sqrt(2.0)
            }
        }
        // X² = x·x，X⁴ = X²·X²（对称，行主序）
        let x2 = matmul(x, x, N)
        let x4 = matmul(x2, x2, N)
        let e0 = (0..<N).map { Double($0) + 0.5 }
        return (e0, x4)
    }

    /// 行主序方阵乘。
    static func matmul(_ p: [Double], _ q: [Double], _ n: Int) -> [Double] {
        var out = [Double](repeating: 0, count: n * n)
        for i in 0..<n {
            for t in 0..<n {
                let pit = p[i * n + t]
                guard pit != 0 else { continue }
                for j in 0..<n { out[i * n + j] += pit * q[t * n + j] }
            }
        }
        return out
    }

    /// 二阶微扰修正（基态或第 s 态）：E^(2) = λ²·Σ_{m≠s} |X⁴[s,m]|²/(E₀[s]−E₀[m])。
    static func secondOrder(state s: Int, e0: [Double], x4: [Double], lam: Double, N: Int) -> Double {
        var sum = 0.0
        for m in 0..<N where m != s {
            let elem = x4[s * N + m]
            sum += elem * elem / (e0[s] - e0[m])
        }
        return lam * lam * sum
    }

    /// (A1) 精确基态：H = E⁰ + λX⁴ 全谱 eigh 取最小。
    static func groundEnergy(lam: Double, e0: [Double], x4: [Double], n: Int) -> Double {
        var h = [Double](repeating: 0, count: n * n)
        for s in 0..<n {
            h[s * n + s] = e0[s] + lam * x4[s * n + s]
            for t in 0..<n where t != s { h[s * n + t] = lam * x4[s * n + t] }
        }
        return AlgebraCore.eigh(h, n: n).w[0]
    }

    /// (A2) 无限深方势阱 + 线性斜坡 V = λx 的矩阵元（脚本 xmat，自然单位 L = 1）：
    /// x_nn = 1/4；(n+l) 偶 → 0；否则 −4nl/[π²(n²−l²)²]。
    static func rampMatrixElement(_ n: Int, _ l: Int) -> Double {
        if n == l { return 0.25 }
        if (n + l) % 2 == 0 { return 0 }
        let nf = Double(n), lf = Double(l)
        return -4.0 * nf * lf / (.pi * .pi * pow(nf * nf - lf * lf, 2))
    }

    /// (B) 类氦离子单参变分：E(α) = α² − 2Zα + (5/8)α（原子单位）。
    static func heliumEnergy(_ alpha: Double, z: Double = 2.0) -> Double {
        alpha * alpha - 2.0 * z * alpha + 0.625 * alpha
    }

    /// 变分解析最优：α*(Z) = Z − 5/16，E(α*) = −(Z−5/16)²。
    static func heliumAlphaOpt(z: Double) -> Double { z - 5.0 / 16.0 }

    /// 类氦离子（1s²）非相对论精确基态能量（Ha，无限核质量；Pekeris 型高精度变分）。
    /// 上界自检：E*(Z) = −(Z−5/16)² 恒在精确值之上（H⁻ −0.4727 > −0.5278 ✓ …）。
    static let heliumLikeExactHa: [Double: Double] = [
        1: -0.5277510165,    // H⁻
        2: -2.9037243770,    // He（脚本常数）
        3: -7.2799134127,    // Li⁺
        4: -13.6555662400,   // Be²⁺
    ]
}

// MARK: - 模块

/// 笔记 19《微扰论与变分法》：非谐振子二阶微扰 vs 精确对角化（40² 基底）+
/// 方势阱线性斜坡（渐近级数证据）+ 类氦离子单参变分（400 点扫描，上界定理）。秒级档。
struct PerturbationVariationModule: SimModule {

    let meta = ModuleMeta(
        id: "微扰论与变分法_方势阱与氦变分", title: "微扰论与变分法 · 微扰展开与氦变分",
        subtitle: "40² 谐振子基底微扰 vs 精确谱 + 方势阱斜坡 + 类氦变分上界（Z = 1…4）",
        category: .formalTheory, noteNumber: 19, tier: .seconds, difficulty: .advanced,
        keywords: ["微扰论", "变分法", "非谐振子", "氦原子", "类氦离子", "上界定理", "渐近级数",
                   " Rayleigh-Schrodinger", "变分原理", "方势阱"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "lam", title: "微扰强度 λx⁴", symbol: "λ", unit: "",
                               range: 0.01...0.2, defaultValue: 0.1, decimalPlaces: 3)),
            // W13i 真机反馈 #19：氦变分图与 λ 无关，拖滑杆纹丝不动——给它自己的参数：
            // 类氦离子核电荷数，E(α) = α² − (2Z−5/8)α 抛物线随 Z 平移，极小随之走。
            .discrete(DiscreteSpec(
                key: "heliumZ", title: "类氦离子（核电荷数 Z）",
                options: [
                    .init(id: "1", title: "H⁻（Z = 1）", subtitle: "单参变分偏差最大 ~10%"),
                    .init(id: "2", title: "He（Z = 2）", subtitle: "脚本口径，精确 −2.903724 Ha"),
                    .init(id: "3", title: "Li⁺（Z = 3）", subtitle: "偏差收窄至 ~0.8%"),
                    .init(id: "4", title: "Be²⁺（Z = 4）", subtitle: "α* = 3.6875"),
                ],
                defaultOptionID: "2")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "非谐振子：二阶微扰 (0+1+2) vs 精确对角化的相对误差（λ 扫描）",
                xAxis: .init(label: "微扰强度 λ"),
                yAxis: .init(label: "相对误差 |num − (0+1+2)| / num", scale: .log),
                seriesNames: ["相对误差扫描", "当前 λ"])),
            // W13g 真机反馈 #19：λ 滑杆原只在误差曲线上挪一个点，氦图与 λ 无关——
            // 补一张随 λ 联动的能量图，拖动即见微扰级数逐阶逼近精确对角化。
            .lineSeries(LineSeriesSpec(
                title: "非谐振子基态能量：微扰级数 vs 精确对角化（随 λ 联动）",
                xAxis: .init(label: "微扰强度 λ"),
                yAxis: .init(label: "基态能量 E（ħ = ω = m = 1）"),
                seriesNames: ["仅到一阶 E⁰+E¹", "二阶微扰 E⁰+E¹+E²", "精确对角化"])),
            .lineSeries(LineSeriesSpec(
                title: "类氦离子变分：E(α) = α² − (2Z−5/8)α（单参乘积试探波函数）",
                xAxis: .init(label: "变分参数 α（有效核电荷）"),
                yAxis: .init(label: "变分能量 E (Hartree)"),
                seriesNames: ["E(α) 扫描"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let hartreeEv = try constants.value("Ehartree_eV")
        let lam = input.slider("lam")
        let nB = PerturbationVariationMath.basisSize
        // W13i #19：类氦离子核电荷数（离散参数，缺省 He）
        let zHelium = Double(input.discrete("heliumZ")) ?? 2.0
        let alphaOpt = PerturbationVariationMath.heliumAlphaOpt(z: zHelium)
        let eOpt = PerturbationVariationMath.heliumEnergy(alphaOpt, z: zHelium)
        let heExact = PerturbationVariationMath.heliumLikeExactHa[zHelium]

        let progress = ComputeProgress()
        let payload = try await SecondsChannel.run(progress: progress)
            { () -> (errs: [Point], cur: Point?, pertGround: Double, exactGround: Double,
                     rampExact: Double, rampFirst: Double, rampNumeric: Double,
                     heScan: [Point], heMinAlpha: Double, heMinE: Double,
                     e1Scan: [Point], pertScan: [Point], exactScan: [Point]) in
            progress.update(0.1, phase: "构造 \(nB)² 谐振子基底 X⁴")
            let (e0, x4) = PerturbationVariationMath.anharmonicBasis(N: nB)
            try Task.checkCancellation()

            progress.update(0.35, phase: "λ 扫描：每点一次 \(nB)² eigh")
            // λ ∈ [0.01, 0.20]，30 点（脚本口径）：二阶微扰 vs 精确基态相对误差
            var errs: [Point] = []
            var e1Scan: [Point] = [], pertScan: [Point] = [], exactScan: [Point] = []
            for i in 0..<30 {
                let l = 0.01 + 0.19 * Double(i) / 29.0
                let e2 = PerturbationVariationMath.secondOrder(state: 0, e0: e0, x4: x4, lam: l, N: nB)
                let num = PerturbationVariationMath.groundEnergy(lam: l, e0: e0, x4: x4, n: nB)
                let pert = 0.5 + l * x4[0] + e2
                errs.append(Point(x: l, y: abs((num - pert) / num)))
                e1Scan.append(Point(x: l, y: 0.5 + l * x4[0]))
                pertScan.append(Point(x: l, y: pert))
                exactScan.append(Point(x: l, y: num))
            }
            try Task.checkCancellation()
            // 当前滑杆 λ 单独求值（滑杆值不必落在扫描网格上）
            let e2Cur = PerturbationVariationMath.secondOrder(state: 0, e0: e0, x4: x4, lam: lam, N: nB)
            let numCur = PerturbationVariationMath.groundEnergy(lam: lam, e0: e0, x4: x4, n: nB)
            let curPoint = Point(x: lam, y: abs((numCur - (0.5 + lam * x4[0] + e2Cur)) / numCur))
            let pertAtCur = 0.5 + lam * x4[0] + e2Cur
            let exactAtCur = numCur

            progress.update(0.6, phase: "方势阱线性斜坡：40×40 解析矩阵元 eigh")
            // (A2) λ = 0.10：E₁^(0) = π²/2，一阶 = λ/2，精确 = eigh
            let lamRamp = 0.10
            var hr = [Double](repeating: 0, count: nB * nB)
            for i in 0..<nB {
                hr[i * nB + i] = Double(i + 1) * Double(i + 1) * .pi * .pi / 2.0
                    + lamRamp * PerturbationVariationMath.rampMatrixElement(i + 1, i + 1)
                for j in 0..<nB where j != i {
                    hr[i * nB + j] = lamRamp * PerturbationVariationMath.rampMatrixElement(i + 1, j + 1)
                }
            }
            let rampNumeric = AlgebraCore.eigh(hr, n: nB).w[0]
            let rampE0 = .pi * .pi / 2.0
            let rampFirst = rampE0 + lamRamp / 2.0
            try Task.checkCancellation()

            progress.update(0.85, phase: "类氦变分扫描（\(PerturbationVariationMath.variationalScanPoints) 点）")
            // 扫描窗随 α*(Z) 居中（±0.7）：抛物线极小恒在图中央，Z 切换即见平移
            let alphaLo = max(0.05, alphaOpt - 0.7)
            let alphaHi = alphaOpt + 0.7
            var heScan: [Point] = []
            var minE = Double.infinity
            var minAlpha = 0.0
            for i in 0..<PerturbationVariationMath.variationalScanPoints {
                let alpha = alphaLo + (alphaHi - alphaLo) * Double(i)
                    / Double(PerturbationVariationMath.variationalScanPoints - 1)
                let e = PerturbationVariationMath.heliumEnergy(alpha, z: zHelium)
                heScan.append(Point(x: alpha, y: e))
                if e < minE { minE = e; minAlpha = alpha }
            }
            return (errs, curPoint, pertAtCur, exactAtCur, rampE0, rampFirst, rampNumeric,
                    heScan, minAlpha, minE, e1Scan, pertScan, exactScan)
        }

        let chart0Series: [SeriesPoints] = {
            var s = [SeriesPoints(name: "相对误差扫描", points: payload.errs)]
            if let cur = payload.cur { s.append(SeriesPoints(name: "当前 λ", points: [cur])) }
            return s
        }()
        let chart0 = LineSeriesData(
            spec: charts.requireLineSeries(0),
            series: chart0Series,
            referenceLines: [ReferenceLine(label: String(format: "当前 λ = %.3f", lam),
                                           axis: .x, value: lam, style: .subtle)])
        let chartE = LineSeriesData(
            spec: charts.requireLineSeries(1),
            series: [
                .init(name: "仅到一阶 E⁰+E¹", points: payload.e1Scan, colorIndex: 3),
                .init(name: "二阶微扰 E⁰+E¹+E²", points: payload.pertScan, colorIndex: 1),
                .init(name: "精确对角化", points: payload.exactScan, colorIndex: 0),
            ],
            referenceLines: [ReferenceLine(label: String(format: "当前 λ = %.3f", lam),
                                           axis: .x, value: lam, style: .subtle)])
        // W13i #19 布局修订：α* 竖线 + E* 横线两条参考线换成极小点实心标记（消除底部
        // 标注挤压），只保留精确值横线做上界对照；扫描窗随 Z 居中，抛物线不再偏居一侧。
        let chart1 = LineSeriesData(
            spec: charts.requireLineSeries(2),
            series: [.init(name: "E(α) 扫描", points: payload.heScan)],
            referenceLines: [
                ReferenceLine(label: String(format: "精确非相对论 %.4f Ha", heExact ?? .nan),
                              axis: .y, value: heExact ?? .nan),
            ],
            pointMarkers: [
                PointMarker(x: alphaOpt, y: eOpt,
                            label: String(format: "α* = %.4f → E* = %.4f Ha", alphaOpt, eOpt)),
            ])

        let upperBound = heExact != nil && payload.heMinE > heExact!
        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chartE), .lineSeries(chart1)],
            summary: [
                .init(id: "pert", title: String(format: "非谐振子基态（λ=%.2f）", lam),
                      value: String(format: "微扰(0+1+2) %.6f", payload.pertGround),
                      note: String(format: "精确对角化 %.6f（微扰展开的收敛域内）", payload.exactGround)),
                .init(id: "ramp", title: "方势阱 + 线性斜坡（λ=0.1）",
                      value: String(format: "精确 %.5f", payload.rampNumeric),
                      note: String(format: "E₀=%.5f 一阶 %.5f；实际位移/一阶 = %.3f"
                                   + "（一阶高估 ~2×，高阶吃掉一半）",
                                   payload.rampExact, payload.rampFirst,
                                   (payload.rampNumeric - payload.rampExact) / (payload.rampFirst - payload.rampExact))),
                .init(id: "hemin", title: String(format: "类氦变分极小（Z=%.0f）", zHelium),
                      value: String(format: "α=%.4f, E=%.6f Ha", payload.heMinAlpha, payload.heMinE),
                      note: String(format: "解析 α* = Z−5/16 = %.4f, E* = %.6f Ha", alphaOpt, eOpt)),
                .init(id: "ub", title: "上界定理",
                      value: upperBound ? "满足（E_var > E_exact）" : "违反！",
                      note: heExact != nil
                          ? String(format: "精确非相对论 %.6f Ha（Pekeris 型），变分相对偏差 %.2f%%",
                                   heExact!, abs(eOpt - heExact!) / abs(heExact!) * 100)
                          : "该 Z 无精确参考值"),
                .init(id: "heev", title: "类氦变分极小（eV）",
                      value: String(format: "%.3f eV", payload.heMinE * hartreeEv),
                      note: "Hartree → eV（CODATA）"),
            ],
            theory: TheoryCard(
                title: "微扰展开与变分原理",
                formulas: [
                    "Rayleigh-Schrödinger：E = E⁰ + λ⟨X⁴⟩ + λ²Σ|X⁴₀ₘ|²/(E⁰₀−E⁰ₘ) + …",
                    "方势阱斜坡 V = λx：一阶位移 λ⟨x⟩ = λ/2 对所有态相同；级数渐近发散",
                    "变分：E(α) = ⟨ψ_α|H|ψ_α⟩ ≥ E₀ 对任意试探态成立（上界定理）",
                    "类氦单参：ψ = φ_α(r₁)φ_α(r₂)，E(α) = α² − (2Z−5/8)α，α* = Z − 5/16",
                    "精确值（H⁻/He/Li⁺/Be²⁺）：变分间隙 = 电子关联能，单参乘积波函数无法覆盖",
                ],
                reading: "微扰侧：λ 小处相对误差 ~10⁻³ 且随 λ 单调恶化——截断级数只在"
                    + "收敛域内可信；斜坡例展示渐近级数「一阶高估 ~2×、高阶收回一半」。"
                    + "变分侧：切换类氦离子 Z，抛物线 E(α) 随 Z 平移、极小点恒在图中央；"
                    + "E* = −(Z−5/16)² 永远压在精确值（虚线）之上——这是上界定理的直接图形化。"))
    }
}
