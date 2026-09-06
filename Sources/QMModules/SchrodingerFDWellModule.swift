import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/薛定谔方程_有限差分势阱.py）

/// 一维无限深势阱有限差分求解的纯函数核。
///
/// Python 口径：N = 800 内点，dx = L/(N+1)，三对角哈密顿
///   diag = ħ²/(m·dx²)，off = −ħ²/(2m·dx²)（H = −ħ²/2m · d²/dx²），
/// `np.linalg.eigh` 全谱 → 升序取前 6（Swift 侧经 AlgebraCore.eighLowest 部分谱等价）。
enum SchrodingerFDWellMath {

    /// 内点数（脚本固定 N = 800）。
    static let gridPoints = 800
    /// 参与图表/摘要的本征对个数（脚本取前 6）。
    static let eigenpairs = 6
    /// 波函数图叠加的本征态个数（W13d 对齐笔记 13 图：前 6 个本征函数）。
    static let plottedStates = 6

    /// 三对角 FD 哈密顿（行主序稠密存放，只填三对角，其余为 0）。
    /// - Parameters:
    ///   - widthMeters: 阱宽 L（m）
    ///   - points: 内点数 N
    static func hamiltonian(widthMeters: Double, points: Int,
                            hbar: Double, mass: Double) -> [Double] {
        let dx = widthMeters / Double(points + 1)
        let diag = hbar * hbar / (mass * dx * dx)
        let off = -hbar * hbar / (2.0 * mass * dx * dx)
        var a = [Double](repeating: 0, count: points * points)
        for i in 0..<points { a[i * points + i] = diag }
        for i in 0..<(points - 1) {
            a[i * points + i + 1] = off
            a[(i + 1) * points + i] = off
        }
        return a
    }

    /// 前 k 个本征对（部分谱）：w 升序（J），v 行主序 n×k（第 s 列 = 第 s 个特征向量）。
    static func solveLowest(widthMeters: Double, points: Int, count k: Int,
                            hbar: Double, mass: Double) -> (w: [Double], v: [Double]) {
        let a = hamiltonian(widthMeters: widthMeters, points: points, hbar: hbar, mass: mass)
        return AlgebraCore.eighLowest(a, n: points, count: k)
    }

    /// 解析本征值 E_n = n²π²ħ²/(2mL²)（J），n = 1...count。
    static func analyticEnergies(widthMeters: Double, count: Int,
                                 hbar: Double, mass: Double) -> [Double] {
        let hbar2Over2m = hbar * hbar / (2.0 * mass)
        let base = .pi * .pi * hbar2Over2m / (widthMeters * widthMeters)
        var out = [Double](repeating: 0, count: count)
        for n in 1...count { out[n - 1] = base * Double(n) * Double(n) }
        return out
    }

    /// 网格归一化（脚本口径 ψ/√(Σψ²·dx)）。
    static func normalizeOnGrid(_ psi: [Double], dx: Double) -> [Double] {
        let norm = sqrt(psi.reduce(0) { $0 + $1 * $1 } * dx)
        return psi.map { $0 / norm }
    }

    /// 符号变化计数 = Σ[diff(sign(ψ)) ≠ 0]（np.sign 语义：零点计两次过渡）。
    static func countNodes(_ psi: [Double]) -> Int {
        guard !psi.isEmpty else { return 0 }
        func sgn(_ x: Double) -> Double { x > 0 ? 1 : (x < 0 ? -1 : 0) }
        var count = 0
        var prev = sgn(psi[0])
        for i in 1..<psi.count {
            let cur = sgn(psi[i])
            if cur != prev { count += 1 }
            prev = cur
        }
        return count
    }
}

// MARK: - 模块

/// 笔记 13《薛定谔方程》：一维无限深势阱的有限差分本征问题。
/// 800² 三对角 eigh（部分本征对）+ E_n ∝ n² 解析对照 + 节点定理验证。秒级档。
struct SchrodingerFDWellModule: SimModule {

    let meta = ModuleMeta(
        id: "薛定谔方程_有限差分势阱", title: "薛定谔方程 · 有限差分势阱",
        subtitle: "800² 三对角 eigh 部分谱 + E_n ∝ n² 解析对照与节点定理",
        category: .schrodinger1D, noteNumber: 13, tier: .seconds, difficulty: .basic,
        keywords: ["有限差分", "势阱", "本征值", "本征函数", "节点定理", "eigh",
                   "定态薛定谔方程", "无限深势阱", "eigenvalue"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "L", title: "阱宽", symbol: "L", unit: "nm",
                               range: 0.5...5, defaultValue: 1.0, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "数值本征函数（节点数 = n−1，曲线叠加 1.2·(n−1) 便于分辨）",
                xAxis: .init(label: "位置 x (nm)"),
                yAxis: .init(label: "ψ(x) + 1.2(n−1) (m⁻¹)"),
                seriesNames: (1...SchrodingerFDWellMath.plottedStates).map { "n=\($0)" })),
            .lineSeries(LineSeriesSpec(
                title: "本征值：有限差分数值 vs 解析 E_n = n²π²ħ²/(2mL²)",
                xAxis: .init(label: "量子数 n"),
                yAxis: .init(label: "能量 E (eV)"),
                seriesNames: ["解析 E_n ∝ n²", "有限差分数值"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let hbar = try constants.value("hbar")
        let mElectron = try constants.value("m_e")
        let eCharge = try constants.value("e")
        let widthNm = input.slider("L")
        let widthMeters = widthNm * 1.0e-9
        let n = SchrodingerFDWellMath.gridPoints
        let k = SchrodingerFDWellMath.eigenpairs
        let dx = widthMeters / Double(n + 1)

        let progress = ComputeProgress()
        let solved = try await SecondsChannel.run(progress: progress) { () -> (w: [Double], v: [Double]) in
            progress.update(0.2, phase: "组装 \(n)² 三对角哈密顿量")
            try Task.checkCancellation()
            progress.update(0.45, phase: "dsyevr 部分谱特征分解（前 \(k) 对）")
            let s = SchrodingerFDWellMath.solveLowest(
                widthMeters: widthMeters, points: n, count: k, hbar: hbar, mass: mElectron)
            try Task.checkCancellation()
            progress.update(0.85, phase: "归一化与节点计数")
            return s
        }

        let valsEv = solved.w.map { $0 / eCharge }
        let analyticEv = SchrodingerFDWellMath.analyticEnergies(
            widthMeters: widthMeters, count: k, hbar: hbar, mass: mElectron)
            .map { $0 / eCharge }
        let relErrs = zip(valsEv, analyticEv).map { abs($0 - $1) / $1 }

        // 图 1：前 6 个本征函数（网格归一化 + 1.2(n−1) 叠加偏置，x 用 nm；W13d 对齐笔记 13）
        let xNm = (1...n).map { Double($0) * dx * 1.0e9 }
        var psiSeries: [SeriesPoints] = []
        var nodeCounts: [Int] = []
        for s in 0..<SchrodingerFDWellMath.plottedStates {
            var psi = (0..<n).map { solved.v[$0 * k + s] }   // 第 s 列 = 第 s 个特征向量
            psi = SchrodingerFDWellMath.normalizeOnGrid(psi, dx: dx)
            nodeCounts.append(SchrodingerFDWellMath.countNodes(psi))
            let offset = 1.2 * Double(s)
            psiSeries.append(.init(name: "n=\(s + 1)",
                                   points: zip(xNm, psi).map { Point(x: $0, y: $1 + offset) }))
        }

        // 图 2：数值 vs 解析（n = 1...6）
        let numericPts = (1...k).map { Point(x: Double($0), y: valsEv[$0 - 1]) }
        let analyticPts = (1...k).map { Point(x: Double($0), y: analyticEv[$0 - 1]) }

        let chart0 = LineSeriesData(
            spec: charts.requireLineSeries(0),
            series: psiSeries)
        let chart1 = LineSeriesData(
            spec: charts.requireLineSeries(1),
            series: [
                .init(name: "解析 E_n ∝ n²", points: analyticPts),
                .init(name: "有限差分数值", points: numericPts),
            ])

        let nodesText = nodeCounts.map(String.init).joined(separator: ", ")
        let worstRel = relErrs.max() ?? 0
        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chart1)],
            summary: [
                .init(id: "e1", title: "基态 E₁",
                      value: String(format: "%.6f eV", valsEv[0]),
                      note: String(format: "解析 %.6f eV", analyticEv[0])),
                .init(id: "rel", title: "最大相对误差（n≤6）",
                      value: String(format: "%.2e", worstRel),
                      note: "FD 误差 ~ (dx)²，随 N 增大按平方收敛"),
                .init(id: "nodes", title: "节点数（前 6 态）",
                      value: nodesText,
                      note: "节点定理：第 n 态恰有 n−1 个节点"),
                .init(id: "grid", title: "网格",
                      value: "\(n)²（部分谱取 \(k) 对）",
                      note: String(format: "dx = %.4f nm", dx * 1.0e9)),
            ],
            theory: TheoryCard(
                title: "有限差分求解定态薛定谔方程",
                formulas: [
                    "−(ħ²/2m)ψ″ = Eψ，硬墙边界 ψ(0) = ψ(L) = 0",
                    "三对角化：diag = ħ²/(m·dx²)，off = −ħ²/(2m·dx²)",
                    "解析谱：E_n = n²π²ħ²/(2mL²)（无限深方势阱）",
                    "节点定理：一维束缚态第 n 本征函数恰有 n−1 个节点",
                ],
                reading: "数值本征值从下方逼近解析谱（FD 误差 ∝ dx²）；本征函数节点数"
                    + "随 n 线性增长是一维束缚态的拓扑指纹——两条线索互相印证"
                    + "「离散哈密顿量正确逼近了连续问题」。"))
    }
}
