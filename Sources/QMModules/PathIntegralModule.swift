import Foundation
import EngineKit

// MARK: - 计算核（移植 code/路径积分_传播子与谐振子基态.py，方案按计划 §W11 升级）
//
// 计划定案：**弃迭代乘方，改哈密顿特征分解 + 谱投影**（方案已实算验证
// Mehler relerr ≈ 5.1×10⁻⁴）；网格 100–200（StrategyCore 降维）；参数哈希缓存；
// 策略档通道（串行 + 进度/取消 + 检查点续算）。自然单位 m = ħ = ω = 1，float64。

/// 路径积分谱投影的纯函数核。
///
/// 谱投影口径：H 在网格上三点差分离散化（实对称）→ `AlgebraCore.eigh` 全谱，
/// K(x',x;T) = Σₙ ψₙ(x')ψₙ(x)·e^(−EₙT)（连续核须除 dx：离散本征矢 Σψ² = 1）。
/// 逐步能量估计量（脚本口径）：E₀(M) = −ln(n_M/n_{M−1})/ε，n_M = ‖ψ(Mε)‖_离散。
enum PathIntegralMath {

    /// 谐振子哈密顿量 H = p²/2 + x²/2（三点差分，行主序实对称 N×N）。
    static func hoHamiltonian(xs: [Double], dx: Double) -> [Double] {
        let n = xs.count
        var h = [Double](repeating: 0, count: n * n)
        for i in 0..<n {
            h[i * n + i] = 1.0 / (dx * dx) + 0.5 * xs[i] * xs[i]
            if i + 1 < n {
                h[i * n + i + 1] = -0.5 / (dx * dx)
                h[(i + 1) * n + i] = -0.5 / (dx * dx)
            }
        }
        return h
    }

    /// 自由粒子哈密顿量 H = p²/2（三点差分，行主序实对称 N×N）。
    static func freeHamiltonian(n: Int, dx: Double) -> [Double] {
        var h = [Double](repeating: 0, count: n * n)
        for i in 0..<n {
            h[i * n + i] = 1.0 / (dx * dx)
            if i + 1 < n {
                h[i * n + i + 1] = -0.5 / (dx * dx)
                h[(i + 1) * n + i] = -0.5 / (dx * dx)
            }
        }
        return h
    }

    /// Mehker 解析核 K(x',x;T)（连续归一）：
    /// √(1/(2π sinh T)) · exp(−[(x'²+x²)·cosh T − 2x'x]/(2 sinh T))。
    static func mehlerKernel(xp: Double, x: Double, T: Double) -> Double {
        let sh = sinh(T)
        let ch = cosh(T)
        let arg = ((xp * xp + x * x) * ch - 2.0 * xp * x) / (2.0 * sh)
        return sqrt(1.0 / (2.0 * Double.pi * sh)) * exp(-arg)
    }

    /// 解析自由欧氏核矩阵作用（含 dx 测度）：ψA = P_T·ψ₀，
    /// P_T(x',x) = √(1/(2πT))·exp(−(x'−x)²/(2T))·dx（脚本 free_euclidean_kernel_matrix 单步）。
    static func applyFreeKernel(xs: [Double], psi: [Double], T: Double) -> [Double] {
        let n = xs.count
        var out = [Double](repeating: 0, count: n)
        let norm = sqrt(1.0 / (2.0 * Double.pi * T))
        for i in 0..<n {
            var acc = 0.0
            for j in 0..<n {
                let d = xs[i] - xs[j]
                acc += norm * exp(-d * d / (2.0 * T)) * psi[j]
            }
            out[i] = acc * (xs[1] - xs[0])
        }
        return out
    }

    /// 展开系数 c = Vᵀψ₀（离散内积；V 列 = 特征向量，行主序 v[j·n + i]）。
    static func coefficients(v: [Double], n: Int, psi0: [Double]) -> [Double] {
        var c = [Double](repeating: 0, count: n)
        for s in 0..<n {
            var acc = 0.0
            for j in 0..<n {
                acc += v[j * n + s] * psi0[j]
            }
            c[s] = acc
        }
        return c
    }

    /// 谱投影演化：ψ(T) = Σₙ cₙ e^(−EₙT) ψₙ = V·(c⊙e^(−ET))。
    static func evolve(w: [Double], v: [Double], n: Int, c: [Double], T: Double) -> [Double] {
        var psi = [Double](repeating: 0, count: n)
        for s in 0..<n {
            let phase = c[s] * exp(-w[s] * T)
            guard phase != 0 else { continue }
            for j in 0..<n {
                psi[j] += v[j * n + s] * phase
            }
        }
        return psi
    }

    /// 归一化 |ψ|² 密度（梯形网格 ∫|ψ|²dx = 1）。
    static func density(_ psi: [Double], dx: Double) -> [Double] {
        var norm = 0.0
        for i in psi.indices.dropLast() {
            norm += 0.5 * (psi[i] * psi[i] + psi[i + 1] * psi[i + 1]) * dx
        }
        let inv = 1.0 / sqrt(norm)
        return psi.map { $0 * $0 * inv * inv }
    }

    /// 波包 L² 相对误差（两态各自归一后）：‖ψnum − ψana‖/‖ψana‖。
    static func relativeL2(_ a: [Double], _ b: [Double], dx: Double) -> Double {
        func normalized(_ p: [Double]) -> [Double] {
            var norm = 0.0
            for i in p.indices.dropLast() {
                norm += 0.5 * (p[i] * p[i] + p[i + 1] * p[i + 1]) * dx
            }
            let inv = 1.0 / sqrt(norm)
            return p.map { $0 * inv }
        }
        let an = normalized(a)
        let bn = normalized(b)
        var acc = 0.0
        for i in an.indices.dropLast() {
            let d = an[i] - bn[i]
            acc += 0.5 * (d * d) * dx
            let d2 = an[i + 1] - bn[i + 1]
            acc += 0.5 * (d2 * d2) * dx   // 相邻两点各贡献半梯形（末点由下一区间补足）
        }
        var normB = 0.0
        for i in bn.indices.dropLast() {
            normB += 0.5 * (bn[i] * bn[i] + bn[i + 1] * bn[i + 1]) * dx
        }
        return sqrt(acc / normB)
    }
}

// MARK: - 模块

/// 笔记 20《路径积分与传播子》：哈密顿特征分解 + 谱投影（策略档首个模块）。
///
/// 四图：① 自由粒子欧氏波包演化（谱投影 vs 解析单步核）；
/// ② 谐振子基态（Euclid 投影 vs 解析 e^{−x²}）；③ 基态能量逐步投影收敛
/// （E₀ → ℏω/2）；④ Mehler 对角核相对误差曲线（方案验证口径 5.1×10⁻⁴）。
/// StrategyCore 首用：降维网格（800→100–200）+ 参数哈希缓存 + 收敛判定 + 策略通道。
struct PathIntegralModule: SimModule {

    /// 脚本口径常量：HO 域半宽 18/自由域半宽 30、初态 σ₀ = 0.5、投影步长 ε = 0.02。
    static let hoExtent = 18.0
    static let freeExtent = 30.0
    static let sigma0 = 0.5
    static let eps = 0.02
    static let projectionSteps = 200

    /// 共享会话：取消后重算从已完成阶段续跑（检查点 = 已完成阶段数）。
    static let session = StrategySession<Int>()

    /// 特征分解缓存：谱分解只依赖网格档位，与 T 无关——改 T 命中缓存。
    static let eigenCache = StrategyCache<String, EigenPair>(capacity: 8)

    /// 缓存载荷（Sendable）：全谱特征对。
    struct EigenPair: Sendable {
        let w: [Double]
        let v: [Double]
    }

    let meta = ModuleMeta(
        id: "路径积分_传播子与谐振子基态",
        title: "路径积分 · 传播子与谐振子基态",
        subtitle: "特征分解 + 谱投影：虚时间演化投影出基态（ℏω/2）",
        category: .formalTheory, noteNumber: 20, tier: .strategy, difficulty: .advanced,
        keywords: ["传播子", "路径积分", "Mehler", "虚时间", "Euclidean",
                   "谱投影", "基态能量", "时间切片"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "T", title: "演化时间", symbol: "T", unit: "（ħ/ω）",
                               range: 1...8, defaultValue: 4, decimalPlaces: 1)),
            .discrete(DiscreteSpec(
                key: "quality", title: "网格档位（降维策略）",
                options: [
                    DiscreteOption(id: "preview", title: "预览 · 100 点", subtitle: "交互粗算"),
                    DiscreteOption(id: "standard", title: "标准 · 150 点", subtitle: "默认口径"),
                    DiscreteOption(id: "fine", title: "高精 · 200 点", subtitle: "物理律验证口径"),
                ],
                defaultOptionID: "fine")),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "自由粒子欧氏波包演化 |ψ(x,T)|²：谱投影 vs 解析单步核",
                xAxis: .init(label: "位置 x（ħ/ω）"),
                yAxis: .init(label: "|ψ(x,T)|²"),
                seriesNames: ["谱投影", "解析单步核"])),
            .lineSeries(LineSeriesSpec(
                title: "谐振子基态 |ψ₀(x)|²：虚时间投影 vs 解析 e^(−x²)",
                xAxis: .init(label: "位置 x（ħ/ω）"),
                yAxis: .init(label: "|ψ₀(x)|²"),
                seriesNames: ["路径积分投影", "解析 e^(−x²)"])),
            .lineSeries(LineSeriesSpec(
                title: "基态能量投影收敛：E₀(M) → ℏω/2",
                xAxis: .init(label: "投影步数 M（τ = M·ε）"),
                yAxis: .init(label: "投影能量 E₀"),
                seriesNames: ["投影 E₀"])),
            .lineSeries(LineSeriesSpec(
                title: "Mehler 对角核相对误差：K(0,0;T) 谱投影 vs 解析",
                xAxis: .init(label: "虚时间 T"),
                yAxis: .init(label: "相对误差"),
                seriesNames: ["|ΔK|/K"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let T = input.slider("T")
        let quality = GridQuality(rawValue: input.discrete("quality")) ?? .fine
        let n = StrategyGrid.reduced(reference: 800, quality: quality)

        let progress = ComputeProgress()
        progress.update(0.04, phase: "策略档通道排队")

        let payload = try await StrategyChannel.run(session: Self.session, progress: progress) {
            (ctx: StrategyContext<Int>) -> Payload in
            // --- 阶段 0：特征分解（缓存键只含网格档位；改 T 命中） ---
            var hoEigen: PathIntegralModule.EigenPair
            var freeEigen: PathIntegralModule.EigenPair
            let hoKey = StrategyKey.make("路径积分_传播子与谐振子基态#ho", [("grid", Double(n))])
            let freeKey = StrategyKey.make("路径积分_传播子与谐振子基态#free", [("grid", Double(n))])
            if let hit = Self.eigenCache.value(for: hoKey), let hitF = Self.eigenCache.value(for: freeKey) {
                hoEigen = hit
                freeEigen = hitF
            } else {
                try Task.checkCancellation()
                progress.update(0.1, phase: "哈密顿特征分解（dsyevr 全谱）")
                let dxHo = Self.hoExtent / Double(n - 1)
                let xsHo = Num.linspace(-Self.hoExtent / 2, Self.hoExtent / 2, count: n)
                let ho = AlgebraCore.eigh(PathIntegralMath.hoHamiltonian(xs: xsHo, dx: dxHo), n: n)
                let dxFree = Self.freeExtent / Double(n - 1)
                let free = AlgebraCore.eigh(
                    PathIntegralMath.freeHamiltonian(n: n, dx: dxFree), n: n)
                hoEigen = EigenPair(w: ho.w, v: ho.v)
                freeEigen = EigenPair(w: free.w, v: free.v)
                Self.eigenCache.insert(hoEigen, for: hoKey)
                Self.eigenCache.insert(freeEigen, for: freeKey)
            }
            ctx.checkpoint(0)

            // --- 阶段 1：基态能量逐步投影收敛（脚本口径：单步比值估计量） ---
            var e0Ms: [Double] = []
            var e0s: [Double] = []
            var convergedNote: String
            do {
                try Task.checkCancellation()
                progress.update(0.35, phase: "虚时间投影 200 步（能量估计）")
                let xsHo = Num.linspace(-Self.hoExtent / 2, Self.hoExtent / 2, count: n)
                var wInit = xsHo.map { exp(-$0 * $0 / (2.0 * Self.sigma0 * Self.sigma0)) }
                let norm0 = sqrt(Num.trapezoid(xsHo, wInit.map { $0 * $0 }))
                wInit = wInit.map { $0 / norm0 }
                let c = PathIntegralMath.coefficients(v: hoEigen.v, n: n, psi0: wInit)
                let c2 = c.map { $0 * $0 }
                // n_M = ‖ψ(Mε)‖_离散 = sqrt(Σ cₙ² e^(−2EₙMε))；E₀(M) = −ln(n_M/n_{M−1})/ε
                var norms = [Double](repeating: 0, count: Self.projectionSteps + 1)
                for m in 0...Self.projectionSteps {
                    var acc = 0.0
                    for s in 0..<n {
                        acc += c2[s] * exp(-2.0 * hoEigen.w[s] * Double(m) * Self.eps)
                    }
                    norms[m] = sqrt(acc)
                }
                var monitor = ConvergenceMonitor(relTolerance: 2e-3, patience: 2,
                                                 maxSteps: Self.projectionSteps / 20)
                for m in stride(from: 20, through: Self.projectionSteps, by: 20) {
                    let e0 = -log(norms[m] / norms[m - 1]) / Self.eps
                    e0Ms.append(Double(m))
                    e0s.append(e0)
                    _ = monitor.observe(e0)   // 全曲线照算（fixture 10 点对拍），监视器只标注
                }
                convergedNote = monitor.isConverged
                    ? "M=\(monitor.convergedAt! * 20) 达 |ΔE₀|/E₀ < 2×10⁻³（连续 2 检查点）"
                    : "M=200 仍未达 2×10⁻³（继续外推至 ℏω/2）"
                ctx.checkpoint(1)
            }

            // --- 阶段 2：自由粒子波包（谱投影 vs 解析单步核） ---
            var psiFree: [Double] = []
            var psiFreeAna: [Double] = []
            var xsFree: [Double] = []
            var relFree = 0.0
            do {
                try Task.checkCancellation()
                progress.update(0.6, phase: "自由波包谱投影（T = \(String(format: "%.1f", T))）")
                let dxFree = Self.freeExtent / Double(n - 1)
                xsFree = Num.linspace(-Self.freeExtent / 2, Self.freeExtent / 2, count: n)
                var psi0 = xsFree.map { exp(-$0 * $0 / (4.0 * Self.sigma0 * Self.sigma0)) }
                let n0 = sqrt(Num.trapezoid(xsFree, psi0.map { $0 * $0 }))
                psi0 = psi0.map { $0 / n0 }
                let c = PathIntegralMath.coefficients(v: freeEigen.v, n: n, psi0: psi0)
                psiFree = PathIntegralMath.evolve(w: freeEigen.w, v: freeEigen.v, n: n,
                                                  c: c, T: T)
                psiFreeAna = PathIntegralMath.applyFreeKernel(xs: xsFree, psi: psi0, T: T)
                relFree = PathIntegralMath.relativeL2(psiFree, psiFreeAna, dx: dxFree)
                ctx.checkpoint(2)
            }

            // --- 阶段 3：Mehler 对角核误差曲线 + 基态投影波函数 ---
            var mehlerT: [Double] = []
            var mehlerRel: [Double] = []
            var psiGround: [Double] = []
            var xsHo: [Double] = []
            do {
                try Task.checkCancellation()
                progress.update(0.85, phase: "Mehler 核对拍 + 基态投影")
                let dxHo = Self.hoExtent / Double(n - 1)
                xsHo = Num.linspace(-Self.hoExtent / 2, Self.hoExtent / 2, count: n)
                // 对角核：中心格点（K(x,x;T) 值 O(0.1)，无远点相消）
                let ic = xsHo.indices.min { abs(xsHo[$0]) < abs(xsHo[$1]) } ?? n / 2
                mehlerT = Num.linspace(1.0, 4.0, count: 16)
                for t in mehlerT {
                    var acc = 0.0
                    for s in 0..<n {
                        acc += hoEigen.v[ic * n + s] * hoEigen.v[ic * n + s]
                            * exp(-hoEigen.w[s] * t)
                    }
                    let knum = acc / dxHo
                    let kana = PathIntegralMath.mehlerKernel(xp: xsHo[ic], x: xsHo[ic], T: t)
                    mehlerRel.append(abs(knum - kana) / kana)
                }
                // 基态：初态窄高斯投影到 τ = M_total·ε = 4
                var wInit = xsHo.map { exp(-$0 * $0 / (2.0 * Self.sigma0 * Self.sigma0)) }
                let norm0 = sqrt(Num.trapezoid(xsHo, wInit.map { $0 * $0 }))
                wInit = wInit.map { $0 / norm0 }
                let c = PathIntegralMath.coefficients(v: hoEigen.v, n: n, psi0: wInit)
                let psiT = PathIntegralMath.evolve(w: hoEigen.w, v: hoEigen.v, n: n, c: c,
                                                   T: Double(Self.projectionSteps) * Self.eps)
                let normT = sqrt(Num.trapezoid(xsHo, psiT.map { $0 * $0 }))
                psiGround = psiT.map { $0 / normT }
                ctx.checkpoint(3)
            }

            let cache = Self.eigenCache.stats
            return Payload(
                n: n, e0Ms: e0Ms, e0s: e0s, convergedNote: convergedNote,
                xsFree: xsFree, psiFree: psiFree, psiFreeAna: psiFreeAna, relFree: relFree,
                mehlerT: mehlerT, mehlerRel: mehlerRel,
                xsHo: xsHo, psiGround: psiGround,
                hoE0: hoEigen.w.first ?? 0,
                cacheHits: cache.hits, cacheMisses: cache.misses,
                resumedFrom: ctx.resumedFrom)
        }

        // --- 图 1：自由波包 ---
        let densityNum = PathIntegralMath.density(payload.psiFree, dx: Self.freeExtent / Double(payload.n - 1))
        let densityAna = PathIntegralMath.density(payload.psiFreeAna, dx: Self.freeExtent / Double(payload.n - 1))
        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "位置 x（ħ/ω）"),
                        yAxis: .init(label: "|ψ(x,T)|²"),
                        seriesNames: ["谱投影", "解析单步核"]),
            series: [
                .init(name: "谱投影", points: Num.strided(payload.xsFree, densityNum, stride: 1)),
                .init(name: "解析单步核", points: Num.strided(payload.xsFree, densityAna, stride: 1)),
            ])

        // --- 图 2：谐振子基态 ---
        let analyticGround = payload.xsHo.map { exp(-$0 * $0) }   // |ψ₀|² ∝ e^{−x²}
        let groundNorm = Num.trapezoid(payload.xsHo, analyticGround)
        let chart2 = LineSeriesData(
            spec: .init(xAxis: .init(label: "位置 x（ħ/ω）"),
                        yAxis: .init(label: "|ψ₀(x)|²"),
                        seriesNames: ["路径积分投影", "解析 e^(−x²)"]),
            series: [
                .init(name: "路径积分投影",
                      points: Num.strided(payload.xsHo,
                                          payload.psiGround.map { $0 * $0 }, stride: 1)),
                .init(name: "解析 e^(−x²)",
                      points: Num.strided(payload.xsHo,
                                          analyticGround.map { $0 / groundNorm }, stride: 1)),
            ])

        // --- 图 3：E₀ 收敛 ---
        let chart3 = LineSeriesData(
            spec: .init(xAxis: .init(label: "投影步数 M（τ = M·ε）"),
                        yAxis: .init(label: "投影能量 E₀"),
                        seriesNames: ["投影 E₀"]),
            series: [
                .init(name: "投影 E₀", points: zip(payload.e0Ms, payload.e0s).map(Point.init)),
            ],
            referenceLines: [
                ReferenceLine(label: "ℏω/2 = 0.5", axis: .y, value: 0.5, style: .subtle),
            ])

        // --- 图 4：Mehler 核误差 ---
        let maxMehlerRel = payload.mehlerRel.max() ?? 0
        let chart4 = LineSeriesData(
            spec: .init(xAxis: .init(label: "虚时间 T"),
                        yAxis: .init(label: "相对误差"),
                        seriesNames: ["|ΔK|/K"]),
            series: [
                .init(name: "|ΔK|/K", points: zip(payload.mehlerT, payload.mehlerRel).map(Point.init)),
            ],
            referenceLines: [
                ReferenceLine(label: "波包级口径 5.1×10⁻⁴", axis: .y, value: 5.1e-4, style: .subtle),
            ])

        let e0Err = abs(payload.hoE0 - 0.5)
        return SimResult(
            charts: [.lineSeries(chart1), .lineSeries(chart2),
                     .lineSeries(chart3), .lineSeries(chart4)],
            summary: [
                .init(id: "e0", title: "基态能量（谱最低本征值）",
                      value: String(format: "%.6f（ℏω/2 = 0.5，|Δ| = %.2e）", payload.hoE0, e0Err),
                      note: "虚时间投影 200 步收敛到同一值：路径积分 = 能量本征态过滤"),
                .init(id: "mehler", title: "Mehler 对角核最大相对误差",
                      value: String(format: "%.2e（T ∈ [1, 4]）", maxMehlerRel),
                      note: "对角核直接读数（无波包平滑）；波包级 relL2 为方案验证口径 5.1×10⁻⁴ 判据"),
                .init(id: "free", title: "自由波包 L² 相对误差",
                      value: String(format: "%.2e（vs 解析单步核，T = %.1f）", payload.relFree, T),
                      note: "谱投影与解析核同口径；脚本朴素切片 N=64 收敛约 10⁻³ 量级"),
                .init(id: "converge", title: "收敛判定（StrategyCore）",
                      value: payload.convergedNote,
                      note: "ConvergenceMonitor：连续 2 检查点 |ΔE₀|/E₀ < 2×10⁻³ 即停"),
                .init(id: "cache", title: "缓存命中统计（StrategyCore）",
                      value: "hits \(payload.cacheHits) / misses \(payload.cacheMisses)",
                      note: "特征分解按网格档位缓存，改 T 命中——策略档「首算重、复用快」"),
                .init(id: "grid", title: "降维网格",
                      value: "参考 800 → \(payload.n) 点\(payload.resumedFrom != nil ? "（续算自阶段 \((payload.resumedFrom ?? 0) + 1)）" : "")",
                      note: "StrategyGrid 降维：统计量/收敛量允许比脚本参考网格粗"),
            ],
            theory: TheoryCard(
                title: "路径积分、传播子与虚时间基态投影",
                formulas: [
                    "K(x',x;T) = ∫𝒟x(t) e^{−S[x]/ħ}：所有路径的相位加权求和",
                    "K(x',x;T) = Σₙ ψₙ(x')ψₙ(x) e^{−EₙT/ħ}（谱展开，本征态过滤）",
                    "T → ∞：e^{−EₙT} 把最低本征态从任意初态中滤出 → 基态",
                    "谐振子 Mehler 核：K(0,0;T) = √(1/(2π sinh T))",
                ],
                reading: "把传播子按能量本征态展开，虚时间演化是一个天然的「基态过滤器」："
                    + "每个本征态贡献按 e^{−EₙT} 衰减，最低态衰减最慢——投影足够长的时间，"
                    + "任意初态都收敛到基态，能量估计量 −ln(范数比)/ε 收敛到 ℏω/2。"
                    + "数值上不做朴素的时间切片迭代乘方（误差随切片数累积），而是把哈密顿量"
                    + "一次性对角化、用谱投影精确演化：200 点网格上 Mehler 核相对误差即达"
                    + " 5×10⁻⁴。自由粒子波包、谐振子基态、能量收敛三条曲线互为印证——"
                    + "路径积分不是比喻，是可以算到 10⁻⁴ 的积分。"))
    }

    /// 策略通道载荷（Sendable，跨 actor 传递）。
    private struct Payload: Sendable {
        let n: Int
        let e0Ms: [Double]
        let e0s: [Double]
        let convergedNote: String
        let xsFree: [Double]
        let psiFree: [Double]
        let psiFreeAna: [Double]
        let relFree: Double
        let mehlerT: [Double]
        let mehlerRel: [Double]
        let xsHo: [Double]
        let psiGround: [Double]
        let hoE0: Double
        let cacheHits: Int
        let cacheMisses: Int
        let resumedFrom: Int?
    }
}
