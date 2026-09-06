import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/狄拉克方程_氢原子精细结构.py）

/// 狄拉克氢原子精确谱纯函数核。
/// E(n, κ) = mc²·[1 + (Zα)²/(n − κ + √(κ² − (Zα)²))²]^(−1/2)（κ = ±(j+1/2)）；
/// 薛定谔能级 E_S(n) = −mc²(Zα)²/(2n²)（相对连续谱 E = 0 的偏移）。
enum DiracHydrogenMath {

    /// 狄拉克总能量（eV，含静能 mc²）。
    static func eDiracTotal(n: Int, kappa: Int, Z: Double, alpha: Double,
                            meC2eV: Double) -> Double {
        let zac = Z * alpha
        let denom = Double(n - kappa) + (Double(kappa * kappa) - zac * zac).squareRoot()
        return meC2eV * pow(1.0 + zac * zac / (denom * denom), -0.5)
    }

    /// 狄拉克束缚能级相对连续谱偏移 E − mc²（eV，负值）。
    static func eDiracBinding(n: Int, kappa: Int, Z: Double, alpha: Double,
                              meC2eV: Double) -> Double {
        eDiracTotal(n: n, kappa: kappa, Z: Z, alpha: alpha, meC2eV: meC2eV) - meC2eV
    }

    /// 薛定谔束缚能级 E_S(n) = −mc²(Zα)²/(2n²)（eV，负值；l、j 全简并）。
    static func eSchrodinger(n: Int, Z: Double, alpha: Double,
                             meC2eV: Double) -> Double {
        -meC2eV * pow(Z * alpha, 2) / (2.0 * Double(n * n))
    }

    /// 2P₃/₂ − 2S₁/₂ 精细结构劈裂（eV）；主阶项 mc²(Zα)⁴/32。
    static func fineStructureSplit(Z: Double, alpha: Double,
                                   meC2eV: Double) -> Double {
        let e2S = eDiracBinding(n: 2, kappa: 1, Z: Z, alpha: alpha, meC2eV: meC2eV)
        let e2P3 = eDiracBinding(n: 2, kappa: 2, Z: Z, alpha: alpha, meC2eV: meC2eV)
        return e2P3 - e2S
    }

    /// 精细结构主阶近似 mc²(Zα)⁴/32（n=2，j=1/2 ↔ 3/2）。
    static func fineStructureLeading(Z: Double, alpha: Double,
                                      meC2eV: Double) -> Double {
        meC2eV * pow(Z * alpha, 4) / 32.0
    }
}

// MARK: - 模块

/// 笔记 38《狄拉克方程》：氢原子狄拉克精确谱 vs 薛定谔能级——
/// 2S₁/₂ = 2P₁/₂ 简并、2P₃/₂ 精细结构劈裂 ≈ mc²(Zα)⁴/32、Lamb 位移标注。
struct DiracHydrogenModule: SimModule {

    let meta = ModuleMeta(
        id: "狄拉克方程_氢原子精细结构", title: "狄拉克方程 · 氢原子精细结构",
        subtitle: "精确谱 vs 薛定谔简并能级——2P₃/₂−2S₁/₂ ≈ mc²(Zα)⁴/32",
        category: .relativisticQFT, noteNumber: 38, tier: .realtime, difficulty: .advanced,
        keywords: ["狄拉克方程", "氢原子", "精细结构", "能级", "简并", "Lamb 位移",
                   "Dirac", "fine structure", "Lamb shift", "degeneracy"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "Z", title: "原子序数", symbol: "Z", unit: "",
                               range: 1...30, defaultValue: 1,
                               step: 1, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .levelDiagram(LevelDiagramSpec(
                title: "氢原子能级：狄拉克 vs 薛定谔",
                energyAxis: "E − m_e c² (eV)",
                showTransitions: false)),
            .levelDiagram(LevelDiagramSpec(
                title: "n = 2 精细结构（放大）",
                energyAxis: "E − m_e c² (eV) [zoom]",
                showTransitions: false)),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let Z = input.slider("Z").rounded()
        let alpha = try 1.0 / constants.value("alpha_inv")
        let meC2eV = try constants.value("me_c2_MeV") * 1e6

        // 能级（相对连续谱 E = 0 的偏移，负值）——与脚本同组
        let e1S = DiracHydrogenMath.eDiracBinding(n: 1, kappa: 1, Z: Z, alpha: alpha, meC2eV: meC2eV)
        let e2S = DiracHydrogenMath.eDiracBinding(n: 2, kappa: 1, Z: Z, alpha: alpha, meC2eV: meC2eV)
        let e2P3 = DiracHydrogenMath.eDiracBinding(n: 2, kappa: 2, Z: Z, alpha: alpha, meC2eV: meC2eV)
        let e2Schr = DiracHydrogenMath.eSchrodinger(n: 2, Z: Z, alpha: alpha, meC2eV: meC2eV)

        let split = e2P3 - e2S
        let splitLeading = DiracHydrogenMath.fineStructureLeading(Z: Z, alpha: alpha, meC2eV: meC2eV)
        let eVHz = try constants.value("eV_Hz")
        let splitGHz = split * eVHz / 1e9

        return SimResult(
            charts: [
                .levelDiagram(LevelDiagramData(
                    spec: charts.requireLevelDiagram(0),
                    levels: [
                        .init(label: "1S₁/₂ (Dirac)", energy: e1S),
                        .init(label: "2S₁/₂ = 2P₁/₂ (Dirac)", energy: e2S),
                        .init(label: "2P₃/₂ (Dirac)", energy: e2P3),
                        .init(label: "n=2 Schrödinger（简并）", energy: e2Schr),
                    ])),
                .levelDiagram(LevelDiagramData(
                    spec: charts.requireLevelDiagram(1),
                    levels: [
                        .init(label: "2S₁/₂ = 2P₁/₂ (Dirac)", energy: e2S),
                        .init(label: "2P₃/₂ (Dirac)", energy: e2P3),
                        .init(label: "n=2 Schrödinger", energy: e2Schr),
                    ],
                    transitions: [
                        .init(fromIndex: 0, toIndex: 1,
                              label: "ΔE_fs = \(String(format: "%.3e", split)) eV"),
                    ])),
            ],
            summary: [
                .init(id: "e1S", title: "1S₁/₂ 狄拉克",
                      value: String(format: "%.6f eV", e1S),
                      note: "薛定谔 −13.605693 eV + 相对论修正"),
                .init(id: "e2S", title: "2S₁/₂ = 2P₁/₂ 狄拉克",
                      value: String(format: "%.6f eV", e2S),
                      note: "点核狄拉克谱的 j 简并（同 |κ|）"),
                .init(id: "e2P3", title: "2P₃/₂ 狄拉克",
                      value: String(format: "%.6f eV", e2P3),
                      note: "结合更松（j 大 → 相对论修正小）"),
                .init(id: "split", title: "精细结构劈裂 2P₃/₂ − 2S₁/₂",
                      value: String(format: "%.6e eV = %.4f GHz", split, splitGHz),
                      note: "主阶 mc²(Zα)⁴/32 = \(String(format: "%.3e", splitLeading)) eV"),
                .init(id: "lamb", title: "Lamb 位移（QED）",
                      value: "~1058 MHz",
                      note: "解除 2S₁/₂/2P₁/₂ 简并的辐射修正（狄拉克谱之外）"),
            ],
            theory: TheoryCard(
                title: "狄拉克氢原子谱（笔记 38 式 16-17）",
                formulas: [
                    "E(n,κ) = mc²·[1 + (Zα)²/(n − κ + √(κ² − (Zα)²))²]^(−1/2)，κ = ±(j+1/2)",
                    "薛定谔：E_S(n) = −mc²(Zα)²/(2n²)（l、j 全简并）",
                    "狄拉克简并：同 n、同 |κ|（j）简并——2S₁/₂ = 2P₁/₂",
                    "精细结构主阶：ΔE(2P₃/₂ − 2S₁/₂) ≈ mc²(Zα)⁴/32",
                ],
                reading: "上图：n=1 与 n=2 群的整体能级——狄拉克（红）与薛定谔（蓝虚线）"
                    + "几乎重合（差 ~5×10⁻⁵ eV）。下图放大 n=2：2P₃/₂ 比 2S₁/₂ 高约"
                    + " 4.5×10⁻⁵ eV（氢），2S₁/₂ = 2P₁/₂ 的狄拉克简并要等 Lamb 位移"
                    + "（QED，~1058 MHz）才解除。增大 Z：劈裂按 (Zα)⁴ 急剧放大。"))
    }
}