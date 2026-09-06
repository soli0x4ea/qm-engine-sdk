import Foundation
import EngineKit
import Numerics

// MARK: - 计算核（逐句移植 code/量子计算_QFT.py，swift-numerics 复数通道）

/// 量子傅里叶变换纯函数核：n 比特 QFT 规范矩阵
/// U[y][x] = exp(2πi·x·y/N)/√N（N = 2ⁿ），及其酉性 / 谱性质验证。
enum QFTMath {

    /// 复数矩阵乘法 C = A·B（N×N）。
    static func matmul(_ a: [[Complex<Double>]], _ b: [[Complex<Double>]]) -> [[Complex<Double>]] {
        let n = a.count
        var c = [[Complex<Double>]](repeating: .init(repeating: .zero, count: n), count: n)
        for i in 0..<n {
            for k in 0..<n {
                let aik = a[i][k]
                if aik == .zero { continue }
                for j in 0..<n { c[i][j] += aik * b[k][j] }
            }
        }
        return c
    }

    /// 共轭转置 A†。
    static func dagger(_ a: [[Complex<Double>]]) -> [[Complex<Double>]] {
        let n = a.count
        var c = [[Complex<Double>]](repeating: .init(repeating: .zero, count: n), count: n)
        for i in 0..<n {
            for j in 0..<n { c[j][i] = a[i][j].conjugate }
        }
        return c
    }

    /// 酉性偏差 max|A†A − I|。
    static func unitaryError(_ a: [[Complex<Double>]]) -> Double {
        let n = a.count
        let prod = matmul(dagger(a), a)
        var worst = 0.0
        for i in 0..<n {
            for j in 0..<n {
                let d = (i == j) ? (prod[i][j] - Complex(1.0)).length : prod[i][j].length
                worst = max(worst, d)
            }
        }
        return worst
    }

    /// n 比特 QFT 矩阵：U[y][x] = exp(2πi·x·y/N)/√N。
    static func qftMatrix(n: Int) -> [[Complex<Double>]] {
        let N = 1 << n
        let norm = 1.0 / Double(N).squareRoot()
        var m = [[Complex<Double>]](repeating: .init(repeating: .zero, count: N), count: N)
        for y in 0..<N {
            for x in 0..<N {
                let angle = 2.0 * .pi * Double((x * y) % N) / Double(N)
                m[y][x] = Complex(cos(angle) * norm, sin(angle) * norm)
            }
        }
        return m
    }

    /// 取反置换 P_neg[x → (−x) mod N]（QFT² 的标志性质）。
    static func negationPermutation(n: Int) -> [[Complex<Double>]] {
        let N = 1 << n
        var m = [[Complex<Double>]](repeating: .init(repeating: .zero, count: N), count: N)
        for x in 0..<N { m[((-x) % N + N) % N][x] = Complex(1) }
        return m
    }

    /// QFT² == 取反置换偏差 max|U² − P_neg|。
    static func squareError(n: Int) -> Double {
        let u = qftMatrix(n: n)
        let uu = matmul(u, u)
        let p = negationPermutation(n: n)
        var worst = 0.0
        for i in 0..<uu.count {
            for j in 0..<uu.count { worst = max(worst, (uu[i][j] - p[i][j]).length) }
        }
        return worst
    }

    /// 元素等幅偏差 max||U[y][x]|² − 1/N|（脚本口径：每列各元素模长均 1/√N，
    /// 傅里叶变换的标志）。
    static func columnNormError(n: Int) -> Double {
        let u = qftMatrix(n: n)
        let N = u.count
        var worst = 0.0
        for y in 0..<N {
            for x in 0..<N {
                let m2 = u[y][x].real * u[y][x].real + u[y][x].imaginary * u[y][x].imaginary
                worst = max(worst, abs(m2 - 1.0 / Double(N)))
            }
        }
        return worst
    }

    /// 2 比特受控相位门 R_k = diag(1, 1, 1, e^{2πi/2^k})。
    static func controlledPhase(k: Int) -> [[Complex<Double>]] {
        let phase = Complex(length: 1.0, phase: 2.0 * .pi / Double(1 << k))
        return [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, phase]]
    }

    /// 将 2 比特门 |11⟩ 分量嵌入 n 比特空间（control、target 为线索引，0 = 最低位）。
    static func embed2q(_ gate: [[Complex<Double>]], n: Int, control: Int, target: Int) -> [[Complex<Double>]] {
        let N = 1 << n
        var m = [[Complex<Double>]](repeating: .init(repeating: .zero, count: N), count: N)
        for idx in 0..<N {
            let bControl = (idx >> control) & 1, bTarget = (idx >> target) & 1
            m[idx][idx] = (bControl == 1 && bTarget == 1) ? gate[3][3] : Complex(1)
        }
        return m
    }

    /// 嵌入单比特 H 门到第 q 线（矩阵张量积 I⊗…⊗H⊗…⊗I）。
    static func hadamardEmbedded(n: Int, q: Int) -> [[Complex<Double>]] {
        let h = 1.0 / 2.0.squareRoot()
        let N = 1 << n
        var m = [[Complex<Double>]](repeating: .init(repeating: .zero, count: N), count: N)
        // H 只作用于第 q 位：i、j 在其余位相同（differ ∈ {0, 2^q}）时非零，
        // 元素 = ±h，负号当且仅当 i、j 的第 q 位均为 1。
        for i in 0..<N {
            for j in 0..<N {
                let differ = i ^ j
                if differ == 0 || differ == (1 << q) {
                    let sign = (((i >> q) & 1) == 1 && ((j >> q) & 1) == 1) ? -1.0 : 1.0
                    m[i][j] = Complex(h * sign)
                }
            }
        }
        return m
    }

    /// 输入基矢 |x⟩ 经 QFT 的输出振幅数组 ψ_out[y] = U[y][x]。
    static func outputAmplitudes(n: Int, inputIndex x: Int) -> [Complex<Double>] {
        let u = qftMatrix(n: n)
        return u.map { $0[x] }
    }

    /// 输出概率 |ψ_out[y]|²。
    static func outputProbabilities(n: Int, inputIndex x: Int) -> [Double] {
        outputAmplitudes(n: n, inputIndex: x).map { $0.real * $0.real + $0.imaginary * $0.imaginary }
    }
}

// MARK: - 模块

/// 笔记 34《量子计算与量子算法》：3 比特量子傅里叶变换——8×8 酉矩阵构造、
/// QFT² = 取反置换的谱性质、|1⟩ 输出的等幅分布（复数通道，swift-numerics）。
struct QFTModule: SimModule {

    let meta = ModuleMeta(
        id: "量子计算_QFT", title: "量子傅里叶变换 · QFT",
        subtitle: "8×8 酉矩阵 U[y][x] = e^{2πi·xy/N}/√8——QFT² = 比特取反",
        category: .quantumInfo, noteNumber: 34, tier: .seconds, difficulty: .advanced,
        keywords: ["量子傅里叶变换", "QFT", "酉矩阵", "受控相位门", "Hadamard",
                   "quantum Fourier transform", "unitary", "controlled-phase"])

    var params: [ParamSpec] {
        [
            .discrete(DiscreteSpec(
                key: "n", title: "量子比特数",
                options: [
                    .init(id: "2", title: "2 比特", subtitle: "N = 4"),
                    .init(id: "3", title: "3 比特", subtitle: "N = 8（脚本基线）"),
                    .init(id: "4", title: "4 比特", subtitle: "N = 16"),
                ],
                defaultOptionID: "3")),
            .slider(SliderSpec(key: "x", title: "输入基矢下标", symbol: "x", unit: "",
                               range: 0...15, defaultValue: 1,
                               step: 1, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .bars(BarSpec(
                title: "QFT 输出概率分布（输入 |x⟩）",
                xAxis: .init(label: "输出基矢 |y⟩"),
                yAxis: .init(label: "概率 |ψ(y)|²"))),
            .lineSeries(LineSeriesSpec(
                title: "QFT 输出振幅的 Re/Im 分量（输入 |x⟩）",
                xAxis: .init(label: "输出基矢下标 y"),
                yAxis: .init(label: "振幅 ψ(y)"),
                seriesNames: ["Re ψ(y)", "Im ψ(y)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let n = Int(input.discrete("n")) ?? 3
        let N = 1 << n
        let x = Int(input.slider("x")) % N

        let psi = QFTMath.outputAmplitudes(n: n, inputIndex: x)
        let probs = psi.map { $0.real * $0.real + $0.imaginary * $0.imaginary }
        let errUnit = QFTMath.unitaryError(QFTMath.qftMatrix(n: n))
        let errSq = QFTMath.squareError(n: n)
        let errCol = QFTMath.columnNormError(n: n)
        let errR2 = QFTMath.unitaryError(QFTMath.controlledPhase(k: 2))
        let errR3 = QFTMath.unitaryError(QFTMath.controlledPhase(k: 3))
        let errEmb = QFTMath.unitaryError(
            QFTMath.embed2q(QFTMath.controlledPhase(k: 2), n: n, control: 0, target: 1))

        let labels = (0..<N).map { "|\($0)⟩" }
        let ys = (0..<N).map(Double.init)
        let reSeries = SeriesPoints(name: "Re ψ(y)", points: zip(ys, psi.map(\.real)).map { Point(x: $0, y: $1) })
        let imSeries = SeriesPoints(name: "Im ψ(y)", points: zip(ys, psi.map(\.imaginary)).map { Point(x: $0, y: $1) })

        return SimResult(
            charts: [
                .bars(BarData(
                    spec: charts.requireBar(0),
                    bars: zip(labels, probs).map { BarItem(label: $0, value: $1) })),
                .lineSeries(LineSeriesData(
                    spec: charts.requireLineSeries(1),
                    series: [reSeries, imSeries],
                    referenceLines: [
                        .init(label: "振幅包络 +1/√N", axis: .y,
                              value: 1.0 / Double(N).squareRoot(), style: .subtle),
                        .init(label: "振幅包络 −1/√N", axis: .y,
                              value: -1.0 / Double(N).squareRoot(), style: .subtle),
                    ])),
            ],
            summary: [
                .init(id: "N", title: "矩阵维度", value: "\(N)×\(N)",
                      note: "n = \(n) 比特，2ⁿ = \(N)"),
                .init(id: "unitary", title: "酉性偏差 max|U†U − I|",
                      value: String(format: "%.2e", errUnit),
                      note: "浮点精度内严格酉"),
                .init(id: "square", title: "QFT² = 取反置换偏差",
                      value: String(format: "%.2e", errSq),
                      note: "QFT 作用两次 = 计算基取反（谱性质）"),
                .init(id: "col", title: "元素等幅偏差 max||U[y][x]|² − 1/N|",
                      value: String(format: "%.2e", errCol),
                      note: "每个元素模方恒 1/N——傅里叶变换的标志"),
                .init(id: "gates", title: "受控相位门 R₂/R₃ 酉性",
                      value: String(format: "%.2e / %.2e", errR2, errR3),
                      note: "嵌入 3 比特空间偏差 \(String(format: "%.2e", errEmb))"),
            ],
            theory: TheoryCard(
                title: "量子傅里叶变换（笔记 34）",
                formulas: [
                    "定义式：U[y][x] = e^{2πi·xy/N}/√N，N = 2ⁿ",
                    "谱性质：QFT² = 取反置换（x → −x mod N）——与基矢约定无关的判据",
                    "元素等幅：|U[y][x]|² = 1/N（每列模长 1/√N）",
                    "受控相位门 R_k = diag(1,1,1,e^{2πi/2^k})，酉性精确成立",
                ],
                reading: "柱图：输入基矢 |x⟩ 经 QFT 后各输出态等概率 1/N（|x⟩=0 时全在 |0⟩）。"
                    + "Re/Im 曲线给出复振幅实虚部——正是 e^{2πi·xy/N} 相位因子的实虚投影。"
                    + "摘要卡的三项偏差均在 1e-16 量级：QFT 是严格酉的。"))
    }
}