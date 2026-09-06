import Foundation
import Accelerate
import EngineKit

// MARK: - 计算核（逐句移植 code/46_矿物光学性质的第一性原理计算_KramersKronig.py）

/// Kramers–Kronig 数值反演（Hilbert 变换）的纯函数核。
/// ε₁(ω) = 1 + (2/π) P∫₀^∞ [ω'ε₂(ω')]/(ω'² − ω²) dω'：
/// 均匀网格梯形积分，极点处主值置零（脚本 kramers_kronig_epsilon1 逐句口径；
/// 万二千点 × 万二千点 = 1.44e8 元素运算，vDSP 向量化）。
enum KramersKronigMath {

    /// Lorentz 振子 ε₂（吸收，输入侧）：ε₂ = ω_p²γω / [(ω₀²−ω²)² + (γω)²]。
    static func lorentzEps2(_ w: Double, w0: Double, wp: Double, gamma: Double) -> Double {
        let denom = pow(w0 * w0 - w * w, 2) + pow(gamma * w, 2)
        return wp * wp * gamma * w / denom
    }

    /// Lorentz 振子 ε₁（解析闭式，严格满足 KK）：ε₁ = 1 + ω_p²(ω₀²−ω²)/[(ω₀²−ω²)² + (γω)²]。
    static func lorentzEps1(_ w: Double, w0: Double, wp: Double, gamma: Double) -> Double {
        let denom = pow(w0 * w0 - w * w, 2) + pow(gamma * w, 2)
        return 1.0 + wp * wp * (w0 * w0 - w * w) / denom
    }

    /// 复折射率 n, κ = √((|ε| ± ε₁)/2)（脚本 refractive_index）。
    static func refractiveIndex(eps1: Double, eps2: Double) -> (n: Double, kappa: Double) {
        let mod = sqrt(eps1 * eps1 + eps2 * eps2)
        let n = sqrt((mod + eps1) / 2.0)
        let kappa = sqrt((mod - eps1) / 2.0)
        return (n, kappa)
    }

    /// ε₂ 全网格（vDSP 逐点标量函数数组化）。
    static func lorentzEps2Grid(_ w: [Double], w0: Double, wp: Double, gamma: Double) -> [Double] {
        w.map { lorentzEps2($0, w0: w0, wp: wp, gamma: gamma) }
    }

    /// Lorentz ε₁ 全网格（解析校验基准）。
    static func lorentzEps1Grid(_ w: [Double], w0: Double, wp: Double, gamma: Double) -> [Double] {
        w.map { lorentzEps1($0, w0: w0, wp: wp, gamma: gamma) }
    }

    /// KK 数值反演 ε₂ → ε₁（脚本 kramers_kronig_epsilon1 的 vDSP 化）。
    ///
    /// 脚本逐句：对每个 ωᵢ，j0 = argmin|ω − ωᵢ|（均匀网格恒为 i 自身），
    /// integrand = ω·ε₂/(ω² − ωᵢ²)，integrand[j0] = 0（极点主值置零），
    /// ε₁[ωᵢ] = 1 + (2/π)·np.trapezoid(integrand, ω)。
    /// 均匀网格梯形 = dx·(Σy − (y₀+y_{N−1})/2)——与 numpy.trapezoid 恒等。
    static func invertEpsilon1(w: [Double], eps2: [Double]) -> [Double] {
        let n = w.count
        precondition(n == eps2.count && n > 2, "KK 反演：网格长度须一致且 > 2")
        let dx = (w[n - 1] - w[0]) / Double(n - 1)

        // 预计算 ω·ε₂ 与 ω²（一次性 O(n)，普通循环即可）
        var wEps2 = [Double](repeating: 0, count: n)
        var wSq = [Double](repeating: 0, count: n)
        for i in 0..<n {
            wEps2[i] = w[i] * eps2[i]
            wSq[i] = w[i] * w[i]
        }

        var integrand = [Double](repeating: 0, count: n)
        var q = [Double](repeating: 0, count: n)
        var out = [Double](repeating: 0, count: n)

        for i in 0..<n {
            var wi2 = w[i] * w[i]
            // q = ω² − ωᵢ²（vDSP_vsubD 广播标量：C = B − A，A 以 stride 0 广播）
            vDSP_vsubD(&wi2, 0, wSq, 1, &q, 1, vDSP_Length(n))
            // integrand = (ω·ε₂)/q（vDSP_vdivD：C = B/A，A 为除数）
            vDSP_vdivD(q, 1, wEps2, 1, &integrand, 1, vDSP_Length(n))
            // 极点主值置零（q[i] 恰为 0 → inf/nan，脚本同位置置零）
            integrand[i] = 0.0
            // 均匀网格梯形 = dx·(Σy − (y₀ + y_{N−1})/2)
            var total = 0.0
            vDSP_sveD(integrand, 1, &total, vDSP_Length(n))
            let trapz = dx * (total - 0.5 * (integrand[0] + integrand[n - 1]))
            out[i] = 1.0 + (2.0 / .pi) * trapz
        }
        return out
    }

    /// 升序网格线性插值（numpy.interp 口径：超界夹到端点）。
    static func interp(_ x: Double, xs: [Double], ys: [Double]) -> Double {
        guard xs.count >= 2, ys.count == xs.count else { return 0 }
        if x <= xs[0] { return ys[0] }
        if x >= xs[xs.count - 1] { return ys[ys.count - 1] }
        var lo = 0, hi = xs.count - 1
        while hi - lo > 1 {
            let mid = (lo + hi) / 2
            if xs[mid] <= x { lo = mid } else { hi = mid }
        }
        let t = (x - xs[lo]) / (xs[hi] - xs[lo])
        return ys[lo] + t * (ys[hi] - ys[lo])
    }
}

// MARK: - 模块

/// 笔记 46《矿物光学性质的第一性原理计算》：模型 (b)——ε₂(ω) 经 Kramers–Kronig
/// 数值反演（Hilbert 变换，12000 点主值梯形）得 ε₁(ω)，Lorentz 振子闭式校验，
/// 再恢复折射率 n(ω)。秒级档（1.44e8 元素运算，vDSP 向量化）。
struct KramersKronigModule: SimModule {

    /// KK 反演网格口径与脚本一致：linspace(0.005, 40, 12000)。
    static let kkPoints = 12000
    /// 绘图网格：linspace(0.05, 12, 1200)（脚本口径）。
    static let plotPoints = 1200

    let meta = ModuleMeta(
        id: "矿物光学性质的第一性原理计算_KramersKronig",
        title: "矿物光学 · Kramers-Kronig 反演",
        subtitle: "ε₂ → ε₁ Hilbert 变换 + Lorentz 闭式校验 + n(ω) 恢复",
        category: .gemology, noteNumber: 46, tier: .seconds, difficulty: .advanced,
        keywords: ["Kramers", "Kronig", "希尔伯特变换", "色散关系", "介电函数",
                   "Lorentz振子", "折射率", "因果性", "消光系数"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "w0", title: "共振频率", symbol: "ω₀", unit: "eV",
                               range: 3...8, defaultValue: 5, decimalPlaces: 1)),
            .slider(SliderSpec(key: "wp", title: "等离子频率", symbol: "ω_p", unit: "eV",
                               range: 4...12, defaultValue: 8, decimalPlaces: 1)),
            .slider(SliderSpec(key: "gamma", title: "阻尼", symbol: "γ", unit: "eV",
                               range: 0.2...1.5, defaultValue: 0.6, decimalPlaces: 2)),
        ]
    }

    var charts: [ChartSpec] {
        [
            // W13（C1）：ε₂ 吸收峰外透明区信息全贴轴 → 拆分并改对数轴；
            // ε₁ 变号必须线性轴，独立成图。
            .lineSeries(LineSeriesSpec(
                title: "ε₂ 吸收谱（对数强度轴）",
                xAxis: .init(label: "光子能量 ħω (eV)"),
                yAxis: .init(label: "ε₂", scale: .log),
                seriesNames: ["ε₂ (输入)"])),
            .lineSeries(LineSeriesSpec(
                title: "Kramers–Kronig 反演校验：ε₁ 解析 vs 反演",
                xAxis: .init(label: "光子能量 ħω (eV)"),
                yAxis: .init(label: "ε₁"),
                seriesNames: ["ε₁ 解析", "ε₁ KK 反演"])),
            .lineSeries(LineSeriesSpec(
                title: "折射率恢复：n(ω) 与 κ(ω)（对数轴）",
                xAxis: .init(label: "光子能量 ħω (eV)"),
                yAxis: .init(label: "n, κ", scale: .log),
                seriesNames: ["n 解析", "n = KK + ε₂", "κ(ω)"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let w0 = input.slider("w0")
        let wp = input.slider("wp")
        let gamma = input.slider("gamma")

        // --- 秒级档调度通道：12000² 主值梯形反演（vDSP 向量化） ---
        let progress = ComputeProgress()
        progress.update(0.05, phase: "生成 KK 反演网格")

        let payload = try await SecondsChannel.run(progress: progress) {
            () -> (wKk: [Double], eps2Kk: [Double], eps1Exact: [Double], eps1Num: [Double],
                   maxAbsErr: Double, maxRelErrN: Double) in
            try Task.checkCancellation()
            progress.update(0.2, phase: "Lorentz ε₂ 输入网格")
            let wKk = Num.linspace(0.005, 40.0, count: Self.kkPoints)
            let eps2Kk = KramersKronigMath.lorentzEps2Grid(wKk, w0: w0, wp: wp, gamma: gamma)
            let eps1Exact = KramersKronigMath.lorentzEps1Grid(wKk, w0: w0, wp: wp, gamma: gamma)

            try Task.checkCancellation()
            progress.update(0.4, phase: "12000² 主值梯形反演（vDSP）")
            let eps1Num = KramersKronigMath.invertEpsilon1(w: wKk, eps2: eps2Kk)

            try Task.checkCancellation()
            progress.update(0.85, phase: "误差统计与 n(ω) 恢复")
            // 误差：ε₁ 数值 vs 解析；n(ω) 相对误差
            var maxAbsErr = 0.0
            var maxRelErrN = 0.0
            for i in 0..<Self.kkPoints {
                maxAbsErr = max(maxAbsErr, abs(eps1Num[i] - eps1Exact[i]))
                let ne = KramersKronigMath.refractiveIndex(eps1: eps1Exact[i], eps2: eps2Kk[i]).n
                let nn = KramersKronigMath.refractiveIndex(eps1: eps1Num[i], eps2: eps2Kk[i]).n
                maxRelErrN = max(maxRelErrN, abs(nn - ne) / (abs(ne) + 1e-9))
            }
            return (wKk, eps2Kk, eps1Exact, eps1Num, maxAbsErr, maxRelErrN)
        }

        // --- 绘图网格（脚本口径 linspace(0.05, 12, 1200)），KK 网格插值 ---
        let w = Num.linspace(0.05, 12.0, count: Self.plotPoints)
        let eps2 = w.map { KramersKronigMath.lorentzEps2($0, w0: w0, wp: wp, gamma: gamma) }
        let eps1A = w.map { KramersKronigMath.lorentzEps1($0, w0: w0, wp: wp, gamma: gamma) }
        let eps1N = w.map { KramersKronigMath.interp($0, xs: payload.wKk, ys: payload.eps1Num) }
        let nA = zip(eps1A, eps2).map { KramersKronigMath.refractiveIndex(eps1: $0, eps2: $1).n }
        let nN = zip(eps1N, eps2).map { KramersKronigMath.refractiveIndex(eps1: $0, eps2: $1).n }
        let kappa = zip(eps1A, eps2).map { KramersKronigMath.refractiveIndex(eps1: $0, eps2: $1).kappa }

        // --- 图 1：ε₂ 吸收谱（对数轴，透明区 4–5 个量级的尾部可见） ---
        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "光子能量 ħω (eV)"),
                        yAxis: .init(label: "ε₂", scale: .log),
                        seriesNames: ["ε₂ (输入)"]),
            series: [
                .init(name: "ε₂ (输入)", points: Num.strided(w, eps2, stride: 2)),
            ],
            referenceLines: [
                ReferenceLine(label: String(format: "ω₀ = %.1f eV", w0), axis: .x, value: w0,
                              style: .subtle),
            ])

        // --- 图 1b：ε₁ 反演校验（ε₁ 变号 → 线性轴） ---
        let chart1b = LineSeriesData(
            spec: .init(xAxis: .init(label: "光子能量 ħω (eV)"),
                        yAxis: .init(label: "ε₁"),
                        seriesNames: ["ε₁ 解析", "ε₁ KK 反演"]),
            series: [
                .init(name: "ε₁ 解析", points: Num.strided(w, eps1A, stride: 2)),
                .init(name: "ε₁ KK 反演", points: Num.strided(w, eps1N, stride: 2)),
            ],
            referenceLines: [
                ReferenceLine(label: String(format: "ω₀ = %.1f eV", w0), axis: .x, value: w0,
                              style: .subtle),
            ])

        // --- 图 2：n/κ 曲线（对数轴；κ 在透明区为 0 → 显示下限 1e-12） ---
        let kappaFloor = kappa.map { max($0, 1e-12) }
        let chart2 = LineSeriesData(
            spec: .init(xAxis: .init(label: "光子能量 ħω (eV)"),
                        yAxis: .init(label: "n, κ", scale: .log),
                        seriesNames: ["n 解析", "n = KK + ε₂", "κ(ω)"]),
            series: [
                .init(name: "n 解析", points: Num.strided(w, nA, stride: 2)),
                .init(name: "n = KK + ε₂", points: Num.strided(w, nN, stride: 2)),
                .init(name: "κ(ω)", points: Num.strided(w, kappaFloor, stride: 2)),
            ])

        // ω = ω₀ 锚点（脚本 [模型 b] ω=5.0 打印）
        let idx = payload.wKk.enumerated().min {
            abs($0.element - w0) < abs($1.element - w0)
        }?.offset ?? 0

        return SimResult(
            charts: [.lineSeries(chart1), .lineSeries(chart1b), .lineSeries(chart2)],
            summary: [
                .init(id: "model", title: "Lorentz 振子参数",
                      value: String(format: "ω₀=%.1f, ω_p=%.1f, γ=%.2f eV", w0, wp, gamma),
                      note: "ε₂ 喂给数值 KK 反推 ε₁，与闭式解析比对——反演可信性校验"),
                .init(id: "epsErr", title: "ε₁ 反演精度",
                      value: String(format: "max|Δ| = %.4e", payload.maxAbsErr),
                      note: "万二千点主值梯形（Python 锚点 5.0993e-2 量级）"),
                .init(id: "nErr", title: "n(ω) 恢复精度",
                      value: String(format: "max rel err = %.4e", payload.maxRelErrN),
                      note: "n 由 KK 反演 ε₁ 与原始 ε₂ 联合给出（Python 锚点 1.39e-3 量级）"),
                .init(id: "anchor", title: String(format: "ω = %.1f eV 锚点", w0),
                      value: String(format: "ε₁_exact = %.4f", payload.eps1Exact[idx]),
                      note: String(format: "ε₁_num = %.4f，ε₂ = %.4f（共振吸收峰）",
                                   payload.eps1Num[idx], payload.eps2Kk[idx])),
                .init(id: "sumrule", title: "求和律",
                      value: String(format: "ε₁(0⁺) → 1 + ω_p²/ω₀² = %.4f", 1 + wp * wp / (w0 * w0)),
                      note: "低频介电常数由 KK 的零频极限保证（因果性 ⇒ 色散）"),
            ],
            theory: TheoryCard(
                title: "Kramers–Kronig 关系：因果性的色散代价",
                formulas: [
                    "ε₁(ω) = 1 + (2/π) P∫₀^∞ [ω'ε₂(ω')]/(ω'² − ω²) dω'",
                    "Lorentz 振子：ε = 1 + ω_p²/(ω₀² − ω² − iγω)（严格满足 KK）",
                    "n + iκ = √ε：n = √((|ε|+ε₁)/2)，κ = √((|ε|−ε₁)/2)",
                    "KK 关系 ⇔ 响应函数在下半复平面解析 ⇔ 响应不能超前于激励",
                ],
                reading: "把 Lorentz 振子的 ε₂ 喂给数值 Hilbert 变换，反推的 ε₁ 与闭式解析"
                    + "在全频段重合——反演方法可信。再由反演的 ε₁ 与原始 ε₂ 得到 n(ω)，"
                    + "与解析 n 重合。共振 ω₀ 处 ε₂ 达峰（吸收最强）、ε₁ 穿过零并变号"
                    + "（色散最强）：吸收与色散是同一因果响应的两面。矿物的颜色与光泽"
                    + "皆由这条复折射率曲线决定——KK 关系是光学常数的守门人。"))
    }
}
