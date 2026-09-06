import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/中心力场与氢原子_径向方程.py）

/// 氢原子径向方程有限差分求解的纯函数核。
///
/// Python 口径：r_max 内 N = 3000 步，h = r_max/N，内点 M = N−1（r_i = i·h, i=1…M），
///   d = ħ²/(μh²) + l(l+1)ħ²/(2μr²) − e²/(4πε₀r)，off = −ħ²/(2μh²)，
/// `np.linalg.eigh` 升序取前 n_states（Swift 侧经 AlgebraCore.eighLowest 部分谱等价）。
enum HydrogenRadialMath {

    /// 每步长分段数（脚本 N = 3000，矩阵 M = N−1 = 2999）。
    static let steps = 3000
    /// 每个 l 求解的最低本征对个数（脚本 n_states = 5）。
    static let eigenpairsPerL = 5
    /// 径向密度图展示的每个 l 最低态数（脚本 min(3, …) → l=0,1,2 各 3 态）。
    static let plottedStatesPerL = 3
    /// 参与求解的角量子数（脚本 l = 0, 1, 2）。
    static let angularList = [0, 1, 2]

    /// 径向三对角哈密顿（行主序稠密）+ 内点网格 r。
    /// 约化质量 μ = m_e·m_p/(m_e+m_p) 由调用方传入。
    static func radialHamiltonian(l: Int, rMax: Double, points N: Int,
                                  mu: Double, hbar: Double,
                                  eCharge: Double, eps0: Double, z: Double = 1)
        -> (r: [Double], a: [Double]) {
        let h = rMax / Double(N)
        let m = N - 1                          // 内点数
        let r = (1...m).map { Double($0) * h }
        let kinDiag = hbar * hbar / (mu * h * h)
        let off = -hbar * hbar / (2.0 * mu * h * h)
        var a = [Double](repeating: 0, count: m * m)
        for i in 0..<m {
            let cent = Double(l * (l + 1)) * hbar * hbar / (2.0 * mu * r[i] * r[i])
            let coul = -z * eCharge * eCharge / (4.0 * .pi * eps0 * r[i])
            a[i * m + i] = kinDiag + cent + coul
        }
        for i in 0..<(m - 1) {
            a[i * m + i + 1] = off
            a[(i + 1) * m + i] = off
        }
        return (r, a)
    }

    /// 径向三对角哈密顿的 (d, e) 表示：d = 对角，e = 次对角（常数 off）。
    /// 数值与 `radialHamiltonian` 稠密矩阵的对应元素逐一相同（同公式同序）。
    static func tridiagonalElements(l: Int, rMax: Double, points N: Int,
                                    mu: Double, hbar: Double,
                                    eCharge: Double, eps0: Double, z: Double = 1)
        -> (d: [Double], e: [Double]) {
        let h = rMax / Double(N)
        let m = N - 1                          // 内点数
        var d = [Double](repeating: 0, count: m)
        for i in 1...m {
            let r = Double(i) * h
            let cent = Double(l * (l + 1)) * hbar * hbar / (2.0 * mu * r * r)
            let coul = -z * eCharge * eCharge / (4.0 * .pi * eps0 * r)
            d[i - 1] = hbar * hbar / (mu * h * h) + cent + coul
        }
        let off = -hbar * hbar / (2.0 * mu * h * h)
        return (d, [Double](repeating: off, count: m - 1))
    }

    /// 给定 l 的前 k 个本征对：w 升序（J），v 行主序 M×k。
    ///
    /// 收尾包 1 性能专项：FD 哈密顿天然三对角，改走 `eighLowestTridiagonal`
    /// （dstevr 直达，1441→4.8 ms @3000）。本征值与稠密 dsyevr 路径 bit 级一致；
    /// 特征向量逐元相对差 ~1e-9（MRRR 回转置路径差异，密度对拍容差 1e-6 之内）。
    static func solveLowest(l: Int, rMax: Double, points N: Int, count k: Int,
                            mu: Double, hbar: Double, eCharge: Double, eps0: Double,
                            z: Double = 1)
        -> (r: [Double], w: [Double], v: [Double]) {
        let (d, e) = tridiagonalElements(l: l, rMax: rMax, points: N,
                                         mu: mu, hbar: hbar, eCharge: eCharge, eps0: eps0,
                                         z: z)
        let m = N - 1
        let s = AlgebraCore.eighLowestTridiagonal(d, e, count: k)
        let r = (1...m).map { Double($0) * (rMax / Double(N)) }
        return (r, s.w, s.v)
    }

    /// 里德伯能量（约化质量版）Ry = μe⁴/(2(4πε₀)²ħ²)（J）。
    static func rydberg(mu: Double, eCharge: Double, eps0: Double, hbar: Double) -> Double {
        mu * eCharge * eCharge * eCharge * eCharge
            / (2.0 * pow(4.0 * .pi * eps0, 2) * hbar * hbar)
    }

    /// 玻尔半径（约化质量版）a₀ = 4πε₀ħ²/(μe²)（m）。
    static func bohrRadius(mu: Double, eCharge: Double, eps0: Double, hbar: Double) -> Double {
        4.0 * .pi * eps0 * hbar * hbar / (mu * eCharge * eCharge)
    }

    /// u(r) 网格归一化（∫u²dr = 1，脚本口径 u/√(Σu²·dr)）。
    static func normalizeU(_ u: [Double], dr: Double) -> [Double] {
        let norm = sqrt(u.reduce(0) { $0 + $1 * $1 } * dr)
        return u.map { $0 / norm }
    }

    /// 径向概率密度 |R|²r² = u²/r（r → 0 处钳 0，脚本 np.where(r > 1e-12, …, 0)）。
    static func radialDensity(u: [Double], r: [Double]) -> [Double] {
        zip(u, r).map { ui, ri in ri > 1.0e-12 ? ui * ui / ri : 0 }
    }

    /// 符号变化计数（节点定理：u 的节点数 = n−l−1）。
    static func countNodes(_ u: [Double]) -> Int {
        guard !u.isEmpty else { return 0 }
        func sgn(_ x: Double) -> Double { x > 0 ? 1 : (x < 0 ? -1 : 0) }
        var count = 0
        var prev = sgn(u[0])
        for i in 1..<u.count {
            let cur = sgn(u[i])
            if cur != prev { count += 1 }
            prev = cur
        }
        return count
    }
}

// MARK: - 模块

/// 笔记 17《中心力场与氢原子》：径向方程有限差分（类氢离子 Z = 1…8 × l = 0,1,2 × 前 5 态），
/// 3000² 稠密 eigh 部分谱（dsyevr RANGE='I'）+ Rydberg 解析对照 + r_max 自适应。秒级档。
struct HydrogenRadialModule: SimModule {

    let meta = ModuleMeta(
        id: "中心力场与氢原子_径向方程", title: "中心力场与氢原子 · 径向方程",
        subtitle: "类氢离子 Z² 标度 × 3000² 径向 FD 部分谱 ×3 角动量 + Rydberg 级数与 n² 简并",
        category: .angularCentral, noteNumber: 17, tier: .seconds, difficulty: .basic,
        keywords: ["氢原子", "类氢离子", "径向方程", "径向波函数", "里德伯", "Rydberg", "玻尔半径",
                   "简并", "拉盖尔", "Laguerre", "束缚态", "r_max 自适应"])

    var params: [ParamSpec] {
        [
            // W13i 真机反馈 #17：原唯一参数 r_max 对 0–25 a₀ 窗内低激发态几乎无影响
            // （束缚态振幅在箱壁处 ~e^{−r/a₀}，1–6 nm 全程肉眼不可辨）——补核电荷数：
            // 密度按 a₀/Z 收缩、能级按 Z² 加深，拖动即见；显示窗 25/Z 随之收缩。
            .slider(SliderSpec(key: "z", title: "核电荷数（类氢离子）", symbol: "Z", unit: "",
                               range: 1...8, defaultValue: 1.0, decimalPlaces: 0)),
            .slider(SliderSpec(key: "rmax", title: "径向箱宽", symbol: "r_max", unit: "nm",
                               range: 1...6, defaultValue: 3.0, decimalPlaces: 2)),
            .discrete(DiscreteSpec(
                key: "rmaxMode", title: "箱宽模式",
                options: [
                    .init(id: "fixed", title: "固定（脚本口径）", subtitle: "r_max = 滑杆值，可观察箱效应"),
                    .init(id: "adaptive", title: "自适应", subtitle: "r_max ≥ 2n²ₘₐₓa₀/Z，保护高激发态"),
                ],
                defaultOptionID: "fixed")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "类氢离子径向概率密度 |R(r)|²r²（最低 9 态；Z 越大按 a₀/Z 收缩）",
                xAxis: .init(label: "r / a₀"),
                yAxis: .init(label: "径向概率密度 |R(r)|²r² (m⁻¹)"),
                seriesNames: ["n=1, l=0", "n=2, l=0", "n=3, l=0",
                              "n=2, l=1", "n=3, l=1", "n=4, l=1",
                              "n=3, l=2", "n=4, l=2", "n=5, l=2"])),
            .lineSeries(LineSeriesSpec(
                title: "束缚态能量：数值 vs Rydberg 解析 E_n = −Z²Ry/n²",
                xAxis: .init(label: "主量子数 n"),
                yAxis: .init(label: "能量 E (eV)"),
                seriesNames: ["解析 −Z²Ry/n²", "有限差分数值 (l=0)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let hbar = try constants.value("hbar")
        let mElectron = try constants.value("m_e")
        let mProton = try constants.value("m_p")
        let eCharge = try constants.value("e")
        let eps0 = try constants.value("epsilon_0")

        let mu = mElectron * mProton / (mElectron + mProton)
        let a0 = HydrogenRadialMath.bohrRadius(mu: mu, eCharge: eCharge, eps0: eps0, hbar: hbar)
        let ry = HydrogenRadialMath.rydberg(mu: mu, eCharge: eCharge, eps0: eps0, hbar: hbar)

        let z = max(1.0, input.slider("z"))
        let rmaxSliderNm = input.slider("rmax")
        let mode = input.discrete("rmaxMode")
        // 自适应：箱宽 ≥ 2n²a₀/Z（n_max = l_max + plotted = 2 + 3 = 5），保护最高激发态尾部
        let nMax = 5
        let adaptiveFloorNm = 2.0 * Double(nMax * nMax) * a0 / z * 1.0e9
        let rmaxNm = mode == "adaptive" ? max(rmaxSliderNm, adaptiveFloorNm) : rmaxSliderNm
        let rMaxMeters = rmaxNm * 1.0e-9

        let progress = ComputeProgress()
        let solved = try await SecondsChannel.run(progress: progress) {
            () -> [(l: Int, r: [Double], w: [Double], v: [Double])] in
            var out: [(l: Int, r: [Double], w: [Double], v: [Double])] = []
            for (idx, l) in HydrogenRadialMath.angularList.enumerated() {
                progress.update(Double(idx) / 3.0, phase: "l=\(l)：组装 \(HydrogenRadialMath.steps)² 哈密顿量")
                let s = HydrogenRadialMath.solveLowest(
                    l: l, rMax: rMaxMeters, points: HydrogenRadialMath.steps,
                    count: HydrogenRadialMath.eigenpairsPerL,
                    mu: mu, hbar: hbar, eCharge: eCharge, eps0: eps0, z: z)
                try Task.checkCancellation()
                out.append((l, s.r, s.w, s.v))
                progress.update(Double(idx + 1) / 3.0, phase: "l=\(l)：dsyevr 部分谱完成")
            }
            return out
        }

        // 图 1：9 条径向密度（l=0: n=1,2,3；l=1: n=2,3,4；l=2: n=3,4,5）
        // W13i #17：显示窗 25/Z a₀（Z=1 即脚本 xlim(0, 25)）；窗内采样步长随窗宽自适应，
        // 窄窗（大 Z）加密，保证窗内仍有 ~300 点的平滑曲线。
        let xWin = 25.0 / z
        let rMaxA0 = rMaxMeters / a0
        let stride = max(1, min(6, Int((Double(HydrogenRadialMath.steps - 1) * xWin / rMaxA0 / 300.0)
            .rounded(.down))))
        var densitySeries: [SeriesPoints] = []
        var nodeReport: [String] = []
        for sol in solved {
            let dr = sol.r[1] - sol.r[0]
            for j in 0..<HydrogenRadialMath.plottedStatesPerL {
                var u = (0..<sol.r.count).map { sol.v[$0 * HydrogenRadialMath.eigenpairsPerL + j] }
                u = HydrogenRadialMath.normalizeU(u, dr: dr)
                let n = j + 1 + sol.l
                let nodes = HydrogenRadialMath.countNodes(u)
                nodeReport.append("n=\(n),l=\(sol.l):\(nodes)")
                let density = HydrogenRadialMath.radialDensity(u: u, r: sol.r)
                let pts = Num.strided(sol.r.map { $0 / a0 }, density, stride: stride)
                    .filter { $0.x <= xWin }
                densitySeries.append(.init(name: "n=\(n), l=\(sol.l)", points: pts))
            }
        }
        // 箱壁参考线：固定模式且箱壁落在显示窗内时画出——r_max 滑杆拖动即见箱壁移动，
        // 高激发态在壁处的抬升也随之可见（W13i #17）。
        let xWall = rMaxA0
        var densityRefs: [ReferenceLine] = []
        if mode == "fixed", xWall < xWin {
            densityRefs.append(ReferenceLine(
                label: String(format: "箱壁 r_max = %.1f a₀", xWall),
                axis: .x, value: xWall, style: .threshold))
        }

        // 图 2：能量对照（l=0 前 5 态，n = 1...5），解析侧 Z² 标度
        let z2RyEv = ry * z * z / eCharge
        let w0 = solved[0].w.map { $0 / eCharge }
        let numericPts = (1...w0.count).map { Point(x: Double($0), y: w0[$0 - 1]) }
        let analyticPts = (1...w0.count).map {
            Point(x: Double($0), y: -z2RyEv / (Double($0) * Double($0)))
        }
        let relErrs = (1...w0.count).map { n -> Double in
            let exact = -z2RyEv / (Double(n) * Double(n))
            return abs(w0[n - 1] - exact) / abs(exact)
        }
        let worstLowN = relErrs.prefix(4).max() ?? 0

        let chart0 = LineSeriesData(
            spec: charts[0].lineSeriesSpec!,
            series: densitySeries,
            referenceLines: densityRefs)
        let chart1 = LineSeriesData(
            spec: charts[1].lineSeriesSpec!,
            series: [
                .init(name: "解析 −Z²Ry/n²", points: analyticPts),
                .init(name: "有限差分数值 (l=0)", points: numericPts),
            ],
            referenceLines: [ReferenceLine(label: "电离极限 E = 0", axis: .y, value: 0,
                                           style: .subtle)])

        let boxRatio = rmaxNm * 1.0e-9 / (2.0 * Double(nMax * nMax) * a0 / z)
        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chart1)],
            summary: [
                .init(id: "a0", title: "玻尔半径 a₀",
                      value: String(format: "%.6f Å", a0 * 1.0e10),
                      note: String(format: "约化质量版（CODATA 0.529177 Å）；类氢有效玻尔半径 a₀/Z = %.4f Å",
                                   a0 * 1.0e10 / z)),
                .init(id: "ry", title: String(format: "里德伯能量 Z²Ry（Z=%.0f）", z),
                      value: String(format: "%.4f eV", z2RyEv),
                      note: "E_n = −Z²Ry/n²（约化质量）"),
                .init(id: "rel", title: "最大相对误差（n≤4）",
                      value: String(format: "%.2e", worstLowN),
                      note: "高激发态受箱宽污染：n=5 起误差增大"),
                .init(id: "box", title: "径向箱宽",
                      value: String(format: "%.2f nm（%.1f a₀）", rmaxNm, rMaxA0),
                      note: String(format: "箱宽/2n²ₘₐₓa₀/Z = %.2f", boxRatio)),
                .init(id: "nodes", title: "u 节点数（n−l−1）",
                      value: nodeReport.joined(separator: " "),
                      note: "节点定理：n=1,l=0 → 0；n=3,l=1 → 1 …"),
            ],
            theory: TheoryCard(
                title: "氢原子径向方程的有限差分求解",
                formulas: [
                    "u(r) = r·R(r)：[−ħ²/2m·d²/dr² + l(l+1)ħ²/2μr² − Ze²/4πε₀r]u = Eu",
                    "解析谱：E_n = −Z²μe⁴/(2(4πε₀)²ħ²)/n² = −Z²·13.6 eV/n²（约化质量）",
                    "Z² 标度：库仑势 ze²/r → 玻尔半径 a₀/Z（密度收缩 Z 倍）、能级加深 Z² 倍",
                    "简并度 n²：E 只依赖 n，与 l 无关（库仑势的隐藏对称性）",
                    "节点定理：u 的节点数 = n−l−1；箱宽需 ≥ 2n²a₀/Z 才能容纳第 n 态",
                ],
                reading: "拖 Z（1→3→8）：九条密度整体向原点收缩 a₀/Z、能级按 Z² 加深——"
                    + "氦离子 He⁺（Z=2）的 E₁ = −54.4 eV 恰是氢的 4 倍。数值能级落 Rydberg 级数上"
                    + "（n≤4 相对误差 ~10⁻⁴，箱效应主导）；固定箱宽拖到 1 nm 可见箱壁竖线压进窗口、"
                    + "高激发态被抬升——切「自适应」后箱宽自动 ≥ 2n²ₘₐₓa₀/Z，低能级恢复。"
                    + "九条密度的节点数恰为 n−l−1。"))
    }
}
