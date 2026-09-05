import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/量子纠缠_纠缠度量与可分性.py）

/// 量子纠缠度量与可分性判据的纯函数核。
/// 全部 4×4 复矩阵经 `AlgebraCore.eighHermitian`（zheevr，W9 基建）对角化。
enum EntanglementMath {

    /// Y⊗Y（实矩阵：σ_y 的虚部乘积给出 ±1）。
    /// 列语义：(Y⊗Y)|00⟩ = −|11⟩，(Y⊗Y)|01⟩ = |10⟩，(Y⊗Y)|10⟩ = |01⟩，(Y⊗Y)|11⟩ = −|00⟩。
    static let yy: [Double] = [
        0, 0, 0, -1,
        0, 0, 1, 0,
        0, 1, 0, 0,
        -1, 0, 0, 0,
    ]

    static let zeros: [Double] = [Double](repeating: 0, count: 16)

    // MARK: 复矩阵小工具（4×4，行主序）

    /// 复矩阵乘 C = A·B（实矩阵传零虚部即可）。
    static func matMul(_ aRe: [Double], _ aIm: [Double], _ bRe: [Double], _ bIm: [Double])
        -> (re: [Double], im: [Double]) {
        var re = [Double](repeating: 0, count: 16)
        var im = [Double](repeating: 0, count: 16)
        for i in 0..<4 {
            for j in 0..<4 {
                var sr = 0.0, si = 0.0
                for k in 0..<4 {
                    let ar = aRe[i * 4 + k], ai = aIm[i * 4 + k]
                    let br = bRe[k * 4 + j], bi = bIm[k * 4 + j]
                    sr += ar * br - ai * bi
                    si += ar * bi + ai * br
                }
                re[i * 4 + j] = sr
                im[i * 4 + j] = si
            }
        }
        return (re, im)
    }

    /// 转置。
    static func transpose(_ a: [Double]) -> [Double] {
        var t = [Double](repeating: 0, count: 16)
        for i in 0..<4 { for j in 0..<4 { t[j * 4 + i] = a[i * 4 + j] } }
        return t
    }

    /// 共轭转置 V†。
    static func conjTranspose(_ re: [Double], _ im: [Double]) -> (re: [Double], im: [Double]) {
        (transpose(re), transpose(im).map(-))
    }

    /// Hermitian 谱（只回本征值，升序）。ρ、ρ^T_B、Wootters M 均走此处。
    static func eigvals(_ re: [Double], _ im: [Double]) -> [Double] {
        AlgebraCore.eighHermitian(re, im, n: 4, wantsVectors: false).w
    }

    /// 对双比特密度矩阵在子系统 B 上取偏置转置（脚本 partial_transpose_B）：
    /// 指标 (a, b; a', b') → 转置 b ↔ b' → (a, b'; a', b)。
    static func partialTransposeB(_ re: [Double], _ im: [Double])
        -> (re: [Double], im: [Double]) {
        var outRe = [Double](repeating: 0, count: 16)
        var outIm = [Double](repeating: 0, count: 16)
        for i in 0..<4 {
            let a = i / 2, b = i % 2
            for j in 0..<4 {
                let ap = j / 2, bp = j % 2
                let dst = (2 * a + bp) * 4 + (2 * ap + b)
                outRe[dst] = re[i * 4 + j]
                outIm[dst] = im[i * 4 + j]
            }
        }
        return (outRe, outIm)
    }

    /// 负度 N = (‖ρ^T_B‖₁ − 1)/2 = (Σ|eig(ρ^T_B)| − 1)/2（脚本 negativity）。
    static func negativity(ptRe: [Double], ptIm: [Double]) -> Double {
        let w = eigvals(ptRe, ptIm)
        return (w.reduce(0) { $0 + abs($1) } - 1.0) / 2.0
    }

    /// 偏置转置最小本征值（< 0 即 NPT；脚本 min_ppt_eig）。
    static func minPptEig(ptRe: [Double], ptIm: [Double]) -> Double {
        eigvals(ptRe, ptIm)[0]
    }

    /// 单态密度矩阵 |Ψ⁻⟩⟨Ψ⁻| = (|01⟩−|10⟩)(⟨01|−⟨10|)/2（脚本 singlet_dm）。
    static func singletDM() -> (re: [Double], im: [Double]) {
        let s = 1.0 / sqrt(2.0)   // ψ = (0, 1/√2, −1/√2, 0)
        let psi = [0.0, s, -s, 0.0]
        var re = [Double](repeating: 0, count: 16)
        for i in 0..<4 { for j in 0..<4 { re[i * 4 + j] = psi[i] * psi[j] } }
        return (re, zeros)
    }

    /// Werner（退极化）态 ρ_W(p) = p|Ψ⁻⟩⟨Ψ⁻| + (1−p)I/4（脚本 werner_state 逐句）。
    static func wernerState(_ p: Double) -> (re: [Double], im: [Double]) {
        let singlet = singletDM()
        var re = [Double](repeating: 0, count: 16)
        for i in 0..<4 { re[i * 4 + i] = (1.0 - p) / 4.0 }   // (1−p)/4 · I（只填对角）
        for i in 0..<16 { re[i] += p * singlet.re[i] }
        return (re, zeros)
    }

    /// 两比特并发度（Wootters）。脚本用 R = ρ·Ỹ（一般非 Hermitian，eigvals 一般谱）；
    /// Swift 侧改用与其**谱恒等**的 Hermitian 形式：
    ///   M = ρ^{1/2} · Ỹ · ρ^{1/2}，Ỹ = Y⊗Y·ρ*·Y⊗Y（Hermitian 半正定）。
    /// 谱恒等性：ρỸ = (ρ^{1/2})(ρ^{1/2}Ỹ) 与 M = (ρ^{1/2}Ỹ)(ρ^{1/2}) 是 AB/BA 对，
    /// 同谱（det(xI−XY)=det(xI−YX)），且 M 半正定——全谱四值不丢，纯态亦稳健。
    /// C = max(0, λ₁ − λ₂ − λ₃ − λ₄)，λᵢ = √eig(M) 降序（笔记式 6.4）。
    static func concurrence(re: [Double], im: [Double]) -> Double {
        // ρ^{1/2} = V·diag(√w)·V†
        let d = AlgebraCore.eighHermitian(re, im, n: 4)
        let sqrtW = d.w.map { sqrt(max($0, 0)) }
        // A = V·diag(√w)：第 c 列乘 √w[c]
        var aRe = [Double](repeating: 0, count: 16)
        var aIm = [Double](repeating: 0, count: 16)
        for c in 0..<4 {
            for r in 0..<4 {
                aRe[r * 4 + c] = d.vRe[r * 4 + c] * sqrtW[c]
                aIm[r * 4 + c] = d.vIm[r * 4 + c] * sqrtW[c]
            }
        }
        // S = ρ^{1/2} = A·V†
        let vt = conjTranspose(d.vRe, d.vIm)
        let s = matMul(aRe, aIm, vt.re, vt.im)

        // Ỹ = Y⊗Y·ρ*·Y⊗Y（YY 实矩阵 → 传零虚部）
        let yRhoStar = matMul(yy, zeros, re, im.map(-))
        let yTilde = matMul(yRhoStar.re, yRhoStar.im, yy, zeros)

        // M = S·Ỹ·S（Hermitian 半正定；eig(M) = eig(ρỸ) ≥ 0）
        let m1 = matMul(s.re, s.im, yTilde.re, yTilde.im)
        let m = matMul(m1.re, m1.im, s.re, s.im)

        let w = eigvals(m.re, m.im)          // 升序
        // 纯态/低秩态的零本征值带 ±1e-32 级噪声：先钳 0 再开方（负数开方 = NaN）
        let lam0 = sqrt(max(w[0], 0)), lam1 = sqrt(max(w[1], 0))
        let lam2 = sqrt(max(w[2], 0)), lam3 = sqrt(max(w[3], 0))
        let c = lam3 - lam2 - lam1 - lam0
        return max(0.0, c.isNaN ? 0.0 : c)
    }

    // MARK: 纯态族纠缠熵（脚本 ent_entropy，笔记式 6.6）

    /// S(θ) = −c²log₂c² − s²log₂s²，|ψ(θ)⟩ = cosθ|00⟩ + sinθ|11⟩。
    static func entEntropy(_ theta: Double) -> Double {
        let c2 = cos(theta) * cos(theta)
        let s2 = sin(theta) * sin(theta)
        if c2 <= 0 || s2 <= 0 { return 0.0 }
        return -(c2 * log2(c2) + s2 * log2(s2))
    }

    // MARK: Ginibre 随机系综（脚本 §3）

    /// Ginibre 随机密度矩阵体：ρ = GG†/tr(GG†)，G = (randn + i·randn)/√2（4×4）。
    /// （单拆出来供测试做 PSD/trace=1 不变量断言。）
    static func ginibreRho(_ nextNormal: () -> Double)
        -> (re: [Double], im: [Double]) {
        // G 复元：实/虚部各一个标准正态，除以 √2
        var gRe = [Double](repeating: 0, count: 16)
        var gIm = [Double](repeating: 0, count: 16)
        let invSqrt2 = 1.0 / sqrt(2.0)
        for i in 0..<16 {
            gRe[i] = nextNormal() * invSqrt2
            gIm[i] = nextNormal() * invSqrt2
        }
        // ρ = G·G†（Hermitian）
        var re = [Double](repeating: 0, count: 16)
        var im = [Double](repeating: 0, count: 16)
        for i in 0..<4 {
            for j in 0..<4 {
                var sr = 0.0, si = 0.0
                for k in 0..<4 {
                    sr += gRe[i * 4 + k] * gRe[j * 4 + k] + gIm[i * 4 + k] * gIm[j * 4 + k]
                    si += gRe[i * 4 + k] * gIm[j * 4 + k] - gIm[i * 4 + k] * gRe[j * 4 + k]
                }
                re[i * 4 + j] = sr
                im[i * 4 + j] = si
            }
        }
        // 归一化 tr(ρ)
        var tr = 0.0
        for i in 0..<4 { tr += re[i * 4 + i] }
        for i in 0..<16 { re[i] /= tr; im[i] /= tr }
        return (re, im)
    }

    /// Ginibre 样本三度量：负度 / 偏置转置最小本征 / 并发度。
    static func ginibreSample(_ nextNormal: () -> Double)
        -> (neg: Double, minPpt: Double, c: Double) {
        let rho = ginibreRho(nextNormal)
        let pt = partialTransposeB(rho.re, rho.im)
        return (negativity(ptRe: pt.re, ptIm: pt.im),
                minPptEig(ptRe: pt.re, ptIm: pt.im),
                concurrence(re: rho.re, im: rho.im))
    }

    /// Box-Muller 正态对生成器（缓存正态对的第二个样本；与 EPRBellMath 同款实现）。
    ///
    /// 种子口径（docs/RNG_SEED_POLICY.md）：种子 = Python 脚本原值 20260901；
    /// PCG64 序列不复刻、分布等价（统计断言对拍）。**逐样本重播种**
    /// （seed ⊕ 金色比例散列的样本序号）——采样体在 @Sendable 闭包内串行执行，
    /// 逐样本独立流既保证两次 compute 逐位可复现，也规避跨并发域的可变捕获。
    static func makeNormalGenerator(seed: UInt64, sampleIndex: Int) -> () -> Double {
        var rng = SplitMix64(seed: seed &+ UInt64(bitPattern: Int64(sampleIndex))
            &* 0x9E3779B97F4A7C15)
        var spare: Double? = nil
        return {
            if let z = spare { spare = nil; return z }
            let u1 = max(rng.nextDouble(), Double.leastNonzeroMagnitude)
            let u2 = rng.nextDouble()
            let r = sqrt(-2 * log(u1))
            let ang = 2 * .pi * u2
            spare = r * sin(ang)
            return r * cos(ang)
        }
    }
}

// MARK: - 模块

/// 笔记 29《量子纠缠》：纠缠度量与可分性——Werner 态并发度/负度/PPT 判据、
/// 纯态族纠缠熵、Ginibre 随机两比特系综 C–N 散点（2×2 下 PPT ⇔ 可分性数值验证）。
/// 采样双模式（W9 基建首用）：交互 300（滑杆 200–500）/ 完整后台 4000。秒级档。
struct EntanglementModule: SimModule {

    /// Python 脚本口径：4000 样本、种子 20260901。
    static let fullSamples = 4000
    static let rngSeed: UInt64 = 20260901

    let meta = ModuleMeta(
        id: "量子纠缠_纠缠度量与可分性",
        title: "量子纠缠 · 纠缠度量与可分性",
        subtitle: "并发度、负度与 PPT 判据：2×2 下 PPT ⇔ 可分",
        category: .quantumInfo, noteNumber: 29, tier: .seconds, difficulty: .advanced,
        keywords: ["纠缠度量", "并发度", "负度", "PPT", "可分性", "Ginibre",
                   "纠缠熵", "密度矩阵", "部分转置", "Wootters"])

    var params: [ParamSpec] {
        [
            .discrete(DiscreteSpec(
                key: "mode", title: "采样模式",
                options: [
                    .init(id: "interactive", title: "交互采样",
                          subtitle: "滑杆指定小样本量（200–500），拖动即时反馈"),
                    .init(id: "full", title: "完整后台",
                          subtitle: "脚本口径 4000 样本，分批 + 进度 + 可取消"),
                ],
                defaultOptionID: "interactive")),
            .slider(SliderSpec(key: "nInteractive", title: "交互样本数", symbol: "N", unit: "个",
                               range: 200...500, defaultValue: 300, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "Werner 态并发度：C = max(0, (3p−1)/2)",
                xAxis: .init(label: "Werner 参数 p"),
                yAxis: .init(label: "并发度 C"),
                seriesNames: ["并发度 C（数值）"])),
            .lineSeries(LineSeriesSpec(
                title: "Werner 态负度与 PPT 判据",
                xAxis: .init(label: "Werner 参数 p"),
                yAxis: .init(label: "N 与 min eig(ρ^T_B)"),
                seriesNames: ["负度 N", "min eig(ρ^T_B)"])),
            .lineSeries(LineSeriesSpec(
                title: "纯态族 |ψ(θ)⟩ 纠缠熵（式 6.6）",
                xAxis: .init(label: "θ"),
                yAxis: .init(label: "纠缠熵 S(θ) / ebit"),
                seriesNames: ["S(θ)"])),
            .frameStack(FrameStackSpec(
                title: "Werner 态密度矩阵随 p 演化（Re ρ，4×4 帧栈）",
                xAxis: .init(label: "列基矢"),
                yAxis: .init(label: "行基矢"),
                valueLabel: "Re ρ", diverging: true)),
            .scatter(ScatterSpec(
                title: "Ginibre 随机两比特态：并发度 vs 负度（PPT ⇔ 可分性）",
                xAxis: .init(label: "负度 N"),
                yAxis: .init(label: "并发度 C"),
                seriesNames: ["PPT（可分）", "NPT（纠缠）"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        // --- 采样模式：交互滑杆小样本 vs 脚本口径全量（W9 采样双模式） ---
        let mode = input.discrete("mode")
        let n = mode == "full" ? Self.fullSamples : Int(input.slider("nInteractive"))

        // --- [1] Werner 态 201 点：C / N / min eig（脚本口径 ps = linspace(0,1,201)） ---
        let ps = Num.linspace(0, 1, count: 201)
        var cNumeric: [Double] = []
        var negWerner: [Double] = []
        var minEigWerner: [Double] = []
        cNumeric.reserveCapacity(ps.count)
        negWerner.reserveCapacity(ps.count)
        minEigWerner.reserveCapacity(ps.count)
        for p in ps {
            let rho = EntanglementMath.wernerState(p)
            let pt = EntanglementMath.partialTransposeB(rho.re, rho.im)
            cNumeric.append(EntanglementMath.concurrence(re: rho.re, im: rho.im))
            negWerner.append(EntanglementMath.negativity(ptRe: pt.re, ptIm: pt.im))
            minEigWerner.append(EntanglementMath.minPptEig(ptRe: pt.re, ptIm: pt.im))
        }

        // --- [2] 纯态族纠缠熵 201 点 ---
        let thetas = Num.linspace(0, .pi / 2, count: 201)
        let sTheta = thetas.map(EntanglementMath.entEntropy)

        // --- [3] Ginibre 系综：SamplingChannel 分批 + 进度 + 取消（W9 基建） ---
        let progress = ComputeProgress()
        let accumulator = GinibreAccumulator()
        try await SamplingChannel.run(count: n, batchSize: 200, progress: progress,
                                      phase: "Ginibre 采样") { index in
            let sample = EntanglementMath.ginibreSample(
                EntanglementMath.makeNormalGenerator(seed: Self.rngSeed, sampleIndex: index))
            accumulator.record(index, sample)
        }

        // PPT / NPT 分组（脚本 is_ppt 口径：min eig ≥ −1e-10）
        var pptPoints: [Point] = []
        var nptPoints: [Point] = []
        var nPpt = 0
        var maxCPpt = 0.0
        var minCNpt = Double.infinity
        var sumNptN = 0.0
        var sumNptC = 0.0
        for s in accumulator.all() {
            if s.minPpt >= -1e-10 {
                nPpt += 1
                maxCPpt = max(maxCPpt, s.c)
                pptPoints.append(Point(x: s.neg, y: s.c))
            } else {
                nptPoints.append(Point(x: s.neg, y: s.c))
                minCNpt = min(minCNpt, s.c)
                sumNptN += s.neg
                sumNptC += s.c
            }
        }
        let meanNptN = sumNptN / Double(max(nptPoints.count, 1))
        let meanNptC = sumNptC / Double(max(nptPoints.count, 1))

        // 散点显示抽稀（点数上限 2048）
        let pptDisplay = Num.decimate(pptPoints, maxCount: 512)
        let nptDisplay = Num.decimate(nptPoints, maxCount: 1024)

        // --- 图 1：Werner 并发度（参考线：可分阈值 p = 1/3） ---
        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "Werner 参数 p"),
                        yAxis: .init(label: "并发度 C"),
                        seriesNames: ["并发度 C（数值）"]),
            series: [.init(name: "并发度 C（数值）",
                           points: Num.strided(ps, cNumeric, stride: 2))],
            referenceLines: [
                ReferenceLine(label: "可分阈值 p = 1/3", axis: .x, value: 1.0 / 3.0, style: .subtle),
            ])

        // --- 图 2：负度 + PPT 最小本征（参考线：y = 0、p = 1/3） ---
        let chart2 = LineSeriesData(
            spec: .init(xAxis: .init(label: "Werner 参数 p"),
                        yAxis: .init(label: "N 与 min eig(ρ^T_B)"),
                        seriesNames: ["负度 N", "min eig(ρ^T_B)"]),
            series: [
                .init(name: "负度 N", points: Num.strided(ps, negWerner, stride: 2)),
                .init(name: "min eig(ρ^T_B)", points: Num.strided(ps, minEigWerner, stride: 2)),
            ],
            referenceLines: [
                ReferenceLine(label: "PPT 边界", axis: .y, value: 0),
                ReferenceLine(label: "p = 1/3", axis: .x, value: 1.0 / 3.0, style: .subtle),
            ])

        // --- 图 3：纠缠熵（参考线：θ = π/4 最大纠缠） ---
        let chart3 = LineSeriesData(
            spec: .init(xAxis: .init(label: "θ"),
                        yAxis: .init(label: "纠缠熵 S(θ) / ebit"),
                        seriesNames: ["S(θ)"]),
            series: [.init(name: "S(θ)", points: Num.strided(thetas, sTheta, stride: 2))],
            referenceLines: [
                ReferenceLine(label: "θ = π/4（最大）", axis: .x, value: .pi / 4, style: .subtle),
            ])

        // --- 图 4：密度矩阵帧栈动画（Werner ρ 实部，p = 0…1 共 11 帧；W9 热图增强首用） ---
        var frames: [FrameStackFrame] = []
        for f in 0...10 {
            let p = Double(f) / 10.0
            let rho = EntanglementMath.wernerState(p)
            var mat: [[Double]] = []
            for i in 0..<4 { mat.append(Array(rho.re[i * 4..<(i * 4 + 4)])) }
            frames.append(FrameStackFrame(label: String(format: "p = %.1f", p),
                                          value: p, values: mat))
        }
        let frameStack = FrameStackData(
            spec: FrameStackSpec(
                title: "Werner 态密度矩阵随 p 演化（Re ρ，4×4 帧栈）",
                xAxis: .init(label: "列基矢"), yAxis: .init(label: "行基矢"),
                valueLabel: "Re ρ", diverging: true),
            xTicks: ["|00⟩", "|01⟩", "|10⟩", "|11⟩"],
            yTicks: ["|00⟩", "|01⟩", "|10⟩", "|11⟩"],
            frames: frames)

        // --- 图 5：C–N 散点（按 PPT/NPT 分组着色） ---
        let chart5 = ScatterData(
            spec: .init(xAxis: .init(label: "负度 N"),
                        yAxis: .init(label: "并发度 C"),
                        seriesNames: ["PPT（可分）", "NPT（纠缠）"]),
            series: [
                .init(name: "PPT（可分）", points: pptDisplay),
                .init(name: "NPT（纠缠）", points: nptDisplay),
            ])

        // p = 1/3 精确点锚（脚本 [Werner] p=1/3 打印）
        let thrRho = EntanglementMath.wernerState(1.0 / 3.0)
        let thrPt = EntanglementMath.partialTransposeB(thrRho.re, thrRho.im)
        let thrMinEig = EntanglementMath.minPptEig(ptRe: thrPt.re, ptIm: thrPt.im)

        return SimResult(
            charts: [.lineSeries(chart1), .lineSeries(chart2), .lineSeries(chart3),
                     .frameStack(frameStack), .scatter(chart5)],
            summary: [
                .init(id: "wernerThr", title: "Werner 可分阈值",
                      value: "p* = 1/3",
                      note: String(format: "p = 1/3 精确点 min eig(ρ^T_B) = %.1e（PPT 边界）",
                                   thrMinEig)),
                .init(id: "samples", title: "Ginibre 系综",
                      value: String(format: "N = %d（%s）", n,
                                    mode == "full" ? "完整后台" : "交互"),
                      note: String(format: "PPT 可分 %d · NPT 纠缠 %d（种子 %d 固定）",
                                   nPpt, n - nPpt, Self.rngSeed)),
                .init(id: "pptLaw", title: "2×2 PPT ⇔ 可分性",
                      value: String(format: "PPT max C = %.1e", maxCPpt),
                      note: String(format: "NPT min C = %.4f > 0：全部样本无一反例（式 2.6）",
                                   minCNpt.isFinite ? minCNpt : 0)),
                .init(id: "nptMean", title: "NPT 组均值",
                      value: String(format: "⟨N⟩ = %.4f, ⟨C⟩ = %.4f", meanNptN, meanNptC),
                      note: "并发度与负度单调相关（同一纠缠序）"),
                .init(id: "entropy", title: "纠缠熵峰值",
                      value: String(format: "S(π/4) = %.4f ebit",
                                    EntanglementMath.entEntropy(.pi / 4)),
                      note: "= 1 ebit：最大纠缠单参族峰值（式 6.6）"),
            ],
            theory: TheoryCard(
                title: "纠缠度量与 Peres-Horodecki 判据",
                formulas: [
                    "Werner 态 ρ_W(p) = p|Ψ⁻⟩⟨Ψ⁻| + (1−p)I/4：C = max(0, (3p−1)/2)",
                    "负度 N = (‖ρ^T_B‖₁ − 1)/2，min eig(ρ^T_B) < 0 ⇔ NPT（纠缠）",
                    "并发度 C = max(0, λ₁−λ₂−λ₃−λ₄)，λᵢ = √eig(ρ^{1/2}Ỹρ^{1/2})（Wootters）",
                    "2×2 与 2×3：PPT ⇔ 可分（Horodecki 定理）；更高维仅必要条件",
                ],
                reading: "Werner 族在 p = 1/3 处穿过可分边界：并发度与负度同时从 0 起步，"
                    + "偏置转置谱恰好触零——PPT 判据在混态域依然锋利。Ginibre 系综"
                    + "（固定种子，交互/完整双模式重采样）给出全部样本无一反例的数值证据："
                    + "PPT 点全部落在 C = 0 线上，NPT 点全部 C > 0。帧栈动画展示密度矩阵"
                    + "随 p 从白噪声（I/4）向单态对角的演化。"))
    }
}

// MARK: - 采样累积器（SamplingChannel body 为 @Sendable；串行执行 + 锁保护可变捕获）

private final class GinibreAccumulator: @unchecked Sendable {
    struct Sample {
        let index: Int
        let neg: Double
        let minPpt: Double
        let c: Double
    }
    private let lock = NSLock()
    private var samples: [Sample] = []
    func record(_ index: Int, _ s: (neg: Double, minPpt: Double, c: Double)) {
        lock.lock()
        samples.append(Sample(index: index, neg: s.neg, minPpt: s.minPpt, c: s.c))
        lock.unlock()
    }
    /// 按样本序号排序（保持与脚本一致的采样序；串行执行下本已有序，防御性排序）。
    func all() -> [Sample] {
        lock.lock(); defer { lock.unlock() }
        return samples.sorted { $0.index < $1.index }
    }
}

// MARK: - Num 扩展（散点抽稀）

extension Num {
    /// Point 数组等距抽稀到 ≤ maxCount 点（保首末点）。
    static func decimate(_ pts: [Point], maxCount: Int) -> [Point] {
        guard pts.count > maxCount else { return pts }
        let stride = (pts.count + maxCount - 1) / maxCount
        var out = pts.enumerated().filter { $0.offset % stride == 0 }.map(\.element)
        if out.last != pts.last { out.append(pts[pts.count - 1]) }
        return out
    }
}
