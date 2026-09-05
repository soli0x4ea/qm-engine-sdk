import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子隐形传态_保真度.py 与 密集编码容量.py）

/// Werner 资源态 ρ(p) = p|Φ⁺⟩⟨Φ⁺| + (1−p)·I/4 的共享信道指标。
/// 单态（贝尔）分数 f = ⟨Φ⁺|ρ|Φ⁺⟩ = (3p+1)/4；
/// 传态保真度 F = (2f+1)/3 = (p+1)/2；密集编码容量 I = 1 + max(0, S(ρ_B) − S(ρ))。
enum WernerChannelMath {

    /// 平均传态保真度 F(p) = (p+1)/2（脚本第二节式 (7)）。
    static func fidelity(p: Double) -> Double { (p + 1) / 2 }

    /// 单态分数 f(p) = (3p+1)/4。
    static func singletFraction(p: Double) -> Double { (3 * p + 1) / 4 }

    /// ρ(p) 的冯·诺依曼熵（bit）：特征值 λ₀=(3p+1)/4（|Φ⁺⟩，重数 1），
    /// λ₁=(1−p)/4（三个正交贝尔态，各重数 1）。log2(0) 以 0 计。
    static func vonNeumannEntropy(p: Double) -> Double {
        let lam0 = (3 * p + 1) / 4
        let lam1 = (1 - p) / 4
        var s = 0.0
        if lam0 > 0 { s += -lam0 * log2(lam0) }
        if lam1 > 0 { s += -3 * lam1 * log2(lam1) }
        return s
    }

    /// 密集编码总可达互信息 I_total(p) = 1 + max(0, S(ρ_B) − S(ρ))，
    /// S(ρ_B) = S(I/2) = 1 bit（脚本第二节式 (10)）。
    static func denseCodingCapacity(p: Double) -> Double {
        1 + max(0, 1 - vonNeumannEntropy(p: p))
    }

    /// 可分边界：ρ(p) 可分 ⟺ p ≤ 1/3（Werner 判据）。
    static let separableBound = 1.0 / 3.0

    /// 密集编码增益阈值 p*：S(ρ(p*)) = 1 的根（二分求根，≈ 0.7476）。
    /// p ∈ (1/3, p*) 资源态纠缠但 C_extra < 0——纠缠存在与可用性门槛不同。
    static let gainThreshold: Double = {
        var lo = 1.0 / 3.0, hi = 1.0
        for _ in 0..<80 {
            let mid = 0.5 * (lo + hi)
            if vonNeumannEntropy(p: mid) > 1 { lo = mid } else { hi = mid }
        }
        return 0.5 * (lo + hi)
    }()
}

// MARK: - 模块一：保真度

/// 笔记 33《量子隐形传态与密集编码》：Werner 信道的平均传态保真度
/// F(p) = (p+1)/2 与经典极限 F = 2/3——p > 1/3（纠缠）才有量子优势。
struct TeleportationFidelityModule: SimModule {

    let meta = ModuleMeta(
        id: "量子隐形传态_保真度", title: "量子隐形传态 · 保真度",
        subtitle: "共享信道的贝尔权重如何决定传态品质——F(p)=(p+1)/2 与经典极限 2/3",
        category: .quantumInfo, noteNumber: 33, tier: .realtime, difficulty: .basic,
        keywords: ["量子隐形传态", "保真度", "Werner 态", "贝尔态", "经典极限",
                   "teleportation", "fidelity", "Werner state", "entanglement"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "p", title: "贝尔权重", symbol: "p", unit: "",
                               range: 0...1, defaultValue: 0.9, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "平均传态保真度 F(p) vs 信道贝尔权重",
                xAxis: .init(label: "共享信道的贝尔权重 p"),
                yAxis: .init(label: "平均传态保真度 F"),
                seriesNames: ["F(p) = (p+1)/2", "工作点"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let p = input.slider("p")

        // 脚本同口径：p ∈ [0,1] 401 点
        let ps = Num.linspace(0, 1, count: 401)
        let curve = ps.map { WernerChannelMath.fidelity(p: $0) }
        let pClassical = WernerChannelMath.separableBound
        let fClassical = 2.0 / 3.0

        let f0 = WernerChannelMath.fidelity(p: p)
        let f = WernerChannelMath.singletFraction(p: p)
        let quantumAdvantage = p > pClassical

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: [
                        .init(name: "F(p) = (p+1)/2",
                              points: zip(ps, curve).map { Point(x: $0, y: $1) }),
                        .init(name: "工作点", points: [Point(x: p, y: f0)], colorIndex: 2),
                    ],
                    referenceLines: [
                        .init(label: "经典极限 F = 2/3", axis: .y, value: fClassical),
                        .init(label: "完美传态 F = 1", axis: .y, value: 1, style: .subtle),
                        .init(label: "可分边界 p = 1/3", axis: .x, value: pClassical,
                              style: .subtle),
                    ])),
            ],
            summary: [
                .init(id: "F", title: "平均传态保真度 F(p)",
                      value: String(format: "%.6f", f0), note: "(p+1)/2"),
                .init(id: "f", title: "单态分数 f",
                      value: String(format: "%.6f", f),
                      note: "⟨Φ⁺|ρ|Φ⁺⟩ = (3p+1)/4"),
                .init(id: "verdict", title: "判定",
                      value: quantumAdvantage ? "超越经典极限" : "未超经典极限",
                      note: quantumAdvantage
                          ? "p > 1/3：信道纠缠，F > 2/3"
                          : "p ≤ 1/3：信道可分，F ≤ 2/3（经典极限）"),
                .init(id: "F_classical", title: "经典极限",
                      value: String(format: "%.6f", fClassical),
                      note: "p=1/3（可分边界）处的 F"),
            ],
            theory: TheoryCard(
                title: "传态保真度（笔记 33 式 7）",
                formulas: [
                    "共享资源：ρ(p) = p|Φ⁺⟩⟨Φ⁺| + (1−p)·I/4（Werner 形式）",
                    "单态分数 f = ⟨Φ⁺|ρ|Φ⁺⟩ = (3p+1)/4",
                    "平均传态保真度 F = (2f+1)/3 = (p+1)/2",
                    "经典极限：可分信道（p ≤ 1/3）F ≤ 2/3；纯贝尔态（p=1）F = 1",
                ],
                reading: "蓝线为 F(p)=(p+1)/2；红色虚线是任何可分/经典信道不可逾越的 2/3。"
                    + "拖动工作点跨过 p=1/3：只有共享纠缠（p>1/3）时传态才优于经典方案，"
                    + "p=1 的纯贝尔信道实现完美传态 F=1。"))
    }
}

// MARK: - 模块二：密集编码容量

/// 笔记 33《量子隐形传态与密集编码》：superdense coding 容量
/// I_total(p) = 1 + max(0, S(ρ_B) − S(ρ))——纠缠提供的额外经典信道容量。
struct DenseCodingModule: SimModule {

    let meta = ModuleMeta(
        id: "量子隐形传态_密集编码容量", title: "量子隐形传态 · 密集编码容量",
        subtitle: "1 ebit + 1 量子比特 = 2 经典比特——Holevo 界之上的纠缠增益",
        category: .quantumInfo, noteNumber: 33, tier: .realtime, difficulty: .basic,
        keywords: ["密集编码", "超密编码", "Holevo 界", "冯·诺依曼熵", "纠缠增益",
                   "superdense coding", "Holevo bound", "von Neumann entropy"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "p", title: "贝尔权重", symbol: "p", unit: "",
                               range: 0...1, defaultValue: 0.9, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "密集编码容量 I_total(p) vs 信道贝尔权重",
                xAxis: .init(label: "共享信道的贝尔权重 p"),
                yAxis: .init(label: "每发送量子比特可载经典比特 I"),
                seriesNames: ["I_total(p)", "工作点"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let p = input.slider("p")

        // 脚本同口径：p ∈ [0,1] 401 点
        let ps = Num.linspace(0, 1, count: 401)
        let curve = ps.map { WernerChannelMath.denseCodingCapacity(p: $0) }

        let sRho = WernerChannelMath.vonNeumannEntropy(p: p)
        let iTotal = WernerChannelMath.denseCodingCapacity(p: p)
        let cExtra = max(0, 1 - sRho)

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: [
                        .init(name: "I_total(p)",
                              points: zip(ps, curve).map { Point(x: $0, y: $1) }),
                        .init(name: "工作点", points: [Point(x: p, y: iTotal)], colorIndex: 2),
                    ],
                    referenceLines: [
                        .init(label: "Holevo 界 = 1 bit（无纠缠）", axis: .y, value: 1),
                        .init(label: "最大容量 = 2 bit", axis: .y, value: 2, style: .subtle),
                        .init(label: "增益阈值 p* ≈ 0.748（S(ρ) = 1）", axis: .x,
                              value: WernerChannelMath.gainThreshold, style: .subtle),
                    ])),
            ],
            summary: [
                .init(id: "I", title: "总可达容量 I_total",
                      value: String(format: "%.6f bit", iTotal),
                      note: "1 + max(0, S(ρ_B)−S(ρ))"),
                .init(id: "S", title: "资源态熵 S(ρ)",
                      value: String(format: "%.6f bit", sRho),
                      note: "p=1 时 0（纯贝尔），p=0 时 2（最大混合）"),
                .init(id: "Cextra", title: "纠缠增益 C_extra",
                      value: String(format: "%.6f bit", cExtra),
                      note: "S(ρ_B)−S(ρ)，S(ρ_B)=1 bit"),
                .init(id: "verdict", title: "判定",
                      value: cExtra > 0 ? "纠缠带来增益" : "无增益（Holevo 界）",
                      note: cExtra > 0
                          ? "S(ρ)<1（p > p* ≈ 0.748）：贝尔测量读出 > 1 bit"
                          : "S(ρ)≥1：增益被截断——注意 p∈(1/3, p*) 纠缠却不增容量"),
            ],
            theory: TheoryCard(
                title: "密集编码容量（笔记 33 式 10）",
                formulas: [
                    "Holevo 界：独立发送 1 量子比特至多携带 1 经典比特",
                    "先共享 1 ebit 再发送 1 量子比特 ⇒ 联合贝尔测量读 2 bit",
                    "纠缠增益 C_extra = S(ρ_B) − S(ρ_AB)（可负，取 max(0,·)）",
                    "I_total(p) = 1 + max(0, 1 − S(ρ(p)))：p=1 → 2 bit，p=0 → 1 bit",
                    "增益阈值 p* ≈ 0.748：S(ρ) = 1 处（高于可分边界 1/3）",
                ],
                reading: "蓝线从 Holevo 界 1 bit（无纠缠，红色虚线）升到 2 bit（纯贝尔资源，"
                    + "绿色点线）——超出部分正是纠缠换来的信道容量。注意增益阈值 p* ≈ 0.748"
                    + "高于可分边界 1/3：p ∈ (1/3, p*) 的资源态虽纠缠，却不足以提升密集编码"
                    + "容量——纠缠存在与可用性是两个不同的门槛。"))
    }
}