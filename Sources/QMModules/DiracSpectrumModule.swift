import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/狄拉克方程_能谱.py）

/// 自由狄拉克哈密顿量能谱的纯函数核。
/// H = c·α_z·p + β·mc²，Dirac 标准表象，自然单位 c = 1：
/// α_z = [[0, σ_z], [σ_z, 0]]，β = diag(I, −I)——**实对称** 4×4，走 dsyevr（AlgebraCore.eigh）。
enum DiracSpectrumMath {

    /// 自由狄拉克哈密顿量 H(p)（实对称 4×4，行主序）。
    /// p 单位 MeV/c（c = 1 → p·c 数值上等于 p），meC2 单位 MeV。
    static func hamiltonian(_ p: Double, meC2: Double) -> [Double] {
        // α_z = [[0,0,1,0],[0,0,0,-1],[1,0,0,0],[0,-1,0,0]]（σ_z = diag(1,−1)）
        // β = diag(1,1,−1,−1)·mc²
        var h = [Double](repeating: 0, count: 16)
        h[0 * 4 + 0] = meC2
        h[1 * 4 + 1] = meC2
        h[2 * 4 + 2] = -meC2
        h[3 * 4 + 3] = -meC2
        // p·α_z 非对角块
        h[0 * 4 + 2] = p
        h[2 * 4 + 0] = p
        h[1 * 4 + 3] = -p
        h[3 * 4 + 1] = -p
        return h
    }

    /// 解析能谱 |E| = √(p²c² + m²c⁴)（c = 1）。
    static func analyticEnergy(_ p: Double, meC2: Double) -> Double {
        sqrt(p * p + meC2 * meC2)
    }

    /// 质量隙 = 2mc²。
    static func massGap(meC2: Double) -> Double { 2.0 * meC2 }

    /// 动量扫描：每点全谱对角化（脚本 E_num shape (N,4)，本征值升序）。
    static func scan(pGrid: [Double], meC2: Double) -> [[Double]] {
        pGrid.map { AlgebraCore.eigh(hamiltonian($0, meC2: meC2), n: 4).w }
    }
}

// MARK: - 模块

/// 笔记 38《狄拉克方程》：自由粒子能谱——4×4 实对称哈密顿量 400 点动量扫描
/// （W8 dsyevr 通道复用），正负能四支 + 解析包络 ±√(p²+m²)，质量隙 2mc² 标注。秒级档。
struct DiracSpectrumModule: SimModule {

    /// 扫描口径与脚本一致：p = linspace(0, 3, 400) MeV/c。
    static let scanPoints = 400

    let meta = ModuleMeta(
        id: "狄拉克方程_能谱",
        title: "狄拉克方程 · 能谱",
        subtitle: "自由粒子 E = ±√(p²c² + m²c⁴)：四支与质量隙",
        category: .relativisticQFT, noteNumber: 38, tier: .seconds, difficulty: .basic,
        keywords: ["能谱", "Dirac方程", "负能解", "质量隙",
                   "自旋简并", "非相对论极限", "相对论量子"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "pMax", title: "动量上限", symbol: "p_max", unit: "MeV/c",
                               range: 1...6, defaultValue: 3, decimalPlaces: 1)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "自由狄拉克能谱 E(p)：数值四支 + 解析包络",
                xAxis: .init(label: "动量 p (MeV/c)"),
                yAxis: .init(label: "能量 E (MeV)"),
                seriesNames: ["数值 E₁", "数值 E₂", "数值 E₃", "数值 E₄",
                              "解析 |E|", "−解析|E|"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let meC2 = try constants.value("me_c2_MeV")   // 电子静能 m_e c²，MeV
        let pMax = input.slider("pMax")

        // --- 秒级档调度通道：400 点 × 4×4 dsyevr ---
        let progress = ComputeProgress()
        progress.update(0.05, phase: "生成动量网格")

        let payload = try await SecondsChannel.run(progress: progress) {
            () -> (pGrid: [Double], eNum: [[Double]], e0: [Double], e1: [Double]) in
            try Task.checkCancellation()
            progress.update(0.25, phase: "p = 0 与 p = 1 锚点对角化")
            let e0 = AlgebraCore.eigh(DiracSpectrumMath.hamiltonian(0.0, meC2: meC2), n: 4).w
            let e1 = AlgebraCore.eigh(DiracSpectrumMath.hamiltonian(1.0, meC2: meC2), n: 4).w

            progress.update(0.45, phase: "400 点动量扫描（dsyevr 全谱）")
            let pGrid = Num.linspace(0, pMax, count: Self.scanPoints)
            let eNum = DiracSpectrumMath.scan(pGrid: pGrid, meC2: meC2)
            try Task.checkCancellation()
            progress.update(0.9, phase: "解析包络")
            return (pGrid, eNum, e0, e1)
        }

        let pGrid = payload.pGrid
        let eAnalytic = pGrid.map { DiracSpectrumMath.analyticEnergy($0, meC2: meC2) }

        // --- 图 1：四支 + 解析包络（6 条曲线，400 点各） ---
        let branch0 = payload.eNum.map { $0[0] }
        let branch1 = payload.eNum.map { $0[1] }
        let branch2 = payload.eNum.map { $0[2] }
        let branch3 = payload.eNum.map { $0[3] }
        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "动量 p (MeV/c)"),
                        yAxis: .init(label: "能量 E (MeV)"),
                        seriesNames: ["数值 E₁", "数值 E₂", "数值 E₃", "数值 E₄",
                                      "解析 |E|", "−解析|E|"]),
            series: [
                .init(name: "数值 E₁", points: Num.strided(pGrid, branch0, stride: 1)),
                .init(name: "数值 E₂", points: Num.strided(pGrid, branch1, stride: 1)),
                .init(name: "数值 E₃", points: Num.strided(pGrid, branch2, stride: 1)),
                .init(name: "数值 E₄", points: Num.strided(pGrid, branch3, stride: 1)),
                .init(name: "解析 |E|", points: Num.strided(pGrid, eAnalytic, stride: 1)),
                .init(name: "−解析|E|", points: Num.strided(pGrid, eAnalytic.map(-), stride: 1)),
            ],
            referenceLines: [
                ReferenceLine(label: "E = 0", axis: .y, value: 0, style: .subtle),
                ReferenceLine(label: String(format: "+mc² = %.4f", meC2), axis: .y, value: meC2,
                              style: .subtle),
                ReferenceLine(label: String(format: "−mc²", meC2), axis: .y, value: -meC2,
                              style: .subtle),
            ])

        let nonrelP = 0.1   // 非相对论极限检查点
        let nonrelNum = payload.eNum.count > 0
            ? DiracSpectrumMath.scan(pGrid: [nonrelP], meC2: meC2)[0][3] : 0

        return SimResult(
            charts: [.lineSeries(chart1)],
            summary: [
                .init(id: "rest", title: "p = 0 静止能",
                      value: String(format: "±%.6f MeV（各 2 重）", payload.e0[0]),
                      note: String(format: "= ±mc²：正负能各 2 重自旋简并（mc² = %.6f）", meC2)),
                .init(id: "p1", title: "p = 1 MeV/c",
                      value: String(format: "±%.6f MeV", payload.e1[3]),
                      note: String(format: "解析 √2·mc² = %.6f（rel < 1e-12）",
                                   DiracSpectrumMath.analyticEnergy(1.0, meC2: meC2))),
                .init(id: "gap", title: "质量隙",
                      value: String(format: "2mc² = %.4f MeV",
                                    DiracSpectrumMath.massGap(meC2: meC2)),
                      note: "正负能支间距 = 2mc² ≈ 1.022 MeV（狄拉克海 → 正电子）"),
                .init(id: "nonrel", title: "非相对论极限",
                      value: String(format: "p = %.1f：E₊ − mc² ≈ p²/2m", nonrelP),
                      note: String(format: "数值 %.6e vs p²/2m %.6e（低动量回归薛定谔）",
                                   nonrelNum - meC2,
                                   nonrelP * nonrelP / (2.0 * meC2))),
                .init(id: "scan", title: "扫描规模",
                      value: String(format: "%d 点 × 4×4 dsyevr", Self.scanPoints),
                      note: "实对称 → dsyevr 通道（W8 线代基建复用）"),
            ],
            theory: TheoryCard(
                title: "自由狄拉克哈密顿量与负能解",
                formulas: [
                    "H = c·α·p + βmc²，α_z = [[0,σ_z],[σ_z,0]]，β = diag(I,−I)",
                    "E = ±√(p²c² + m²c⁴)：正负能各 2 重自旋简并",
                    "质量隙 2mc² ≈ 1.022 MeV：负能海 → 湮没产生正电子",
                    "非相对论极限 E₊ = mc² + p²/2m − (p²)²/8m³c² + …",
                ],
                reading: "对角化 4×4 实对称 H(p) 得到的四支能谱与解析包络 ±√(p²+m²) 逐点"
                    + "重合（rel < 10⁻¹²）：两支正能（电子两个自旋态）与两支负能在 p = 0"
                    + "处以质量隙 2mc² 相隔。负能解不是灾难而是预言——把海填满，"
                    + "空穴就是正电子。低动量端 E₊ − mc² ≈ p²/2m：薛定谔力学是"
                    + "狄拉克力学的非相对论极限。"))
    }
}
