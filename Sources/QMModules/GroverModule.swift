import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子计算_Grover搜索.py）

/// Grover 搜索纯函数核。
/// 旋转半角 sin(θ) = 1/√N；标记态振幅 = sin((2n+1)θ)，
/// 成功概率 P(n) = sin²((2n+1)θ)；最优迭代 n* = round(π/(4θ) − 0.5)。
enum GroverMath {

    /// 旋转半角 θ = arcsin(1/√N)。
    static func theta(N: Double) -> Double {
        asin(1.0 / N.squareRoot())
    }

    /// 第 n 次迭代后命中标记态的成功概率 P(n) = sin²((2n+1)θ)。
    static func probability(n: Int, N: Double) -> Double {
        let t = theta(N: N)
        let a = (2.0 * Double(n) + 1.0) * t
        return sin(a) * sin(a)
    }

    /// P(n) 数组（n = 0…nMax，与脚本 n = arange(0, n_max+1) 同口径）。
    static func probabilities(N: Double, nMax: Int) -> [Double] {
        let t = theta(N: N)
        return (0...nMax).map { n in
            let a = (2.0 * Double(n) + 1.0) * t
            return sin(a) * sin(a)
        }
    }

    /// 最优迭代次数：使 (2n+1)θ 最接近 π/2（脚本 int(round(π/(4θ) − 0.5))）。
    static func optimalIterations(N: Double) -> Int {
        Int((Double.pi / (4.0 * theta(N: N)) - 0.5).rounded())
    }

    /// 峰值概率 P(n*)。
    static func peakProbability(N: Double) -> Double {
        probability(n: optimalIterations(N: N), N: N)
    }

    /// 经典查询次数 ~ N/2。
    static func classicalQueries(N: Double) -> Double { N / 2 }

    /// 量子查询次数 ~ (π/4)√N。
    static func quantumQueries(N: Double) -> Double {
        Double.pi / 4.0 * N.squareRoot()
    }

    /// 二次加速比：(N/2) / ((π/4)√N) = 2√N/π。
    static func speedup(N: Double) -> Double {
        classicalQueries(N: N) / quantumQueries(N: N)
    }
}

// MARK: - 模块

/// 笔记 34《量子计算与量子算法》：Grover 搜索成功概率
/// P(n) = sin²((2n+1)θ) 与最优迭代次数 n* ≈ (π/4)√N——平方根加速。
struct GroverModule: SimModule {

    let meta = ModuleMeta(
        id: "量子计算_Grover搜索", title: "Grover 搜索 · 成功概率",
        subtitle: "P(n)=sin²((2n+1)θ)——n* ≈ (π/4)√N 次迭代平方根加速",
        category: .quantumInfo, noteNumber: 34, tier: .realtime, difficulty: .basic,
        keywords: ["Grover", "量子搜索", "振幅放大", "二次加速", "oracle",
                   "success probability", "amplitude amplification", "quadratic speedup"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "N", title: "搜索空间大小", symbol: "N", unit: "",
                               range: 16...1_048_576, defaultValue: 1024,
                               scale: .log, decimalPlaces: 0)),
            .slider(SliderSpec(key: "nmax", title: "迭代次数上限", symbol: "n_max", unit: "",
                               range: 10...240, defaultValue: 60,
                               step: 10, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "Grover 成功概率 P(n) vs 迭代次数",
                xAxis: .init(label: "迭代次数 n"),
                yAxis: .init(label: "成功概率 P(n)"),
                seriesNames: ["P(n) = sin²((2n+1)θ)", "峰值 n*"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let n = max(2.0, input.slider("N")).rounded()
        let nMax = Int(input.slider("nmax"))

        // 脚本同口径：n = 0…n_max 整数格点（脚本 N=1024、n_max=60）
        let ns = Array(0...nMax)
        let probs = GroverMath.probabilities(N: n, nMax: nMax)

        let nOpt = GroverMath.optimalIterations(N: n)
        let pPeak = GroverMath.peakProbability(N: n)
        let theta = GroverMath.theta(N: n)
        let qQ = GroverMath.quantumQueries(N: n)
        let cQ = GroverMath.classicalQueries(N: n)

        return SimResult(
            charts: [
                .lineSeries(LineSeriesData(
                    spec: charts[0].lineSeriesSpec!,
                    series: [
                        .init(name: "P(n) = sin²((2n+1)θ)",
                              points: zip(ns.map(Double.init), probs).map { Point(x: $0, y: $1) }),
                        .init(name: "峰值 n*", points: [Point(x: Double(nOpt), y: pPeak)],
                              colorIndex: 2),
                    ],
                    referenceLines: [
                        .init(label: "最优迭代 n* = \(nOpt)", axis: .x, value: Double(nOpt)),
                        .init(label: "峰值概率 P(n*) = \(String(format: "%.4f", pPeak))",
                              axis: .y, value: pPeak, style: .subtle),
                    ])),
            ],
            summary: [
                .init(id: "theta", title: "旋转半角 θ",
                      value: String(format: "%.6f rad", theta),
                      note: "sin θ = 1/√N = \(String(format: "%.6f", 1.0 / n.squareRoot()))"),
                .init(id: "nOpt", title: "最优迭代次数 n*",
                      value: String(format: "%d", nOpt),
                      note: "round(π/(4θ) − 0.5) ≈ (π/4)√N"),
                .init(id: "pPeak", title: "峰值成功概率",
                      value: String(format: "%.6f", pPeak),
                      note: "P(n*) = sin²((2n*+1)θ)"),
                .init(id: "speedup", title: "查询次数对比",
                      value: String(format: "%.1f vs %.0f", qQ, cQ),
                      note: "量子 (π/4)√N vs 经典 N/2，加速 \(String(format: "%.1f", GroverMath.speedup(N: n)))×"),
            ],
            theory: TheoryCard(
                title: "Grover 搜索（笔记 34）",
                formulas: [
                    "旋转半角：sin θ = 1/√N（标记态与均匀叠加态张成的二维子空间）",
                    "每次迭代态矢旋转 2θ，标记态振幅 = sin((2n+1)θ)",
                    "成功概率 P(n) = sin²((2n+1)θ)",
                    "最优迭代 n* = round(π/(4θ) − 0.5) ≈ (π/4)√N（N ≫ 1）",
                ],
                reading: "概率从 P(0) = 1/N 出发按 sin² 律攀升，在 n* ≈ (π/4)√N 处达到峰值"
                    + "（N=1024 时 n*=25、P≈0.999）——远小于经典 N/2 次查询。"
                    + "继续迭代会「过冲」：P(n) 周期回落，故迭代数必须停在峰值附近。"))
    }
}
