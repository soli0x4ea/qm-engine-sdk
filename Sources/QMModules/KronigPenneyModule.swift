import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/晶体中的电子_能带与KronigPenney.py）

/// 一维 δ 势梳 Kronig-Penney 的纯函数核（含段检测与态密度）。
enum KronigPenneyMath {

    /// ℏ²/2m 换算：3.81 eV·Å²（脚本常量，自由电子）。
    static let hbar2Over2m = 3.81

    /// KP 方程左端 f(E) = cos(αa) + (P/αa)·sin(αa)，α = √(E/(ℏ²/2m))。
    static func f(_ E: Double, P: Double, a: Double) -> Double {
        let xa = sqrt(E / hbar2Over2m) * a
        return cos(xa) + P / (xa + 1e-30) * sin(xa)
    }

    /// 简约区 k = arccos(clip(f, −1, 1))/a ∈ [0, π/a]。
    static func kReduced(_ E: Double, P: Double, a: Double) -> Double {
        let value = min(max(f(E, P: P, a: a), -1), 1)
        return acos(value) / a
    }

    /// 允许带段检测：E 升序网格上 |f| ≤ 1 的连续段 [(lo, hi)]。
    static func segments(allowed: [Bool], E: [Double]) -> [(lo: Double, hi: Double)] {
        var segs: [(lo: Double, hi: Double)] = []
        var i = 0
        let n = allowed.count
        while i < n {
            if allowed[i] {
                var j = i
                while j < n && allowed[j] { j += 1 }
                segs.append((E[i], E[j - 1]))
                i = j
            } else {
                i += 1
            }
        }
        return segs
    }

    /// np.gradient 等价（edge_order=1：内部中心差分、两端单侧）。
    static func gradient(_ y: [Double]) -> [Double] {
        let n = y.count
        guard n >= 2 else { return y }
        var out = [Double](repeating: 0, count: n)
        out[0] = y[1] - y[0]
        out[n - 1] = y[n - 1] - y[n - 2]
        if n > 2 {
            for i in 1..<(n - 1) { out[i] = (y[i + 1] - y[i - 1]) / 2 }
        }
        return out
    }

    /// 线性插值分位数（np.nanpercentile 线性法对应）。
    static func percentile(_ sorted: [Double], _ q: Double) -> Double {
        guard !sorted.isEmpty else { return 0 }
        let pos = q / 100 * Double(sorted.count - 1)
        let lo = Int(pos.rounded(.down))
        let hi = Int(pos.rounded(.up))
        guard lo != hi, hi < sorted.count else { return sorted[min(lo, sorted.count - 1)] }
        let frac = pos - Double(lo)
        return sorted[lo] + (sorted[hi] - sorted[lo]) * frac
    }

    /// 段内一维态密度 g(E) = (2/π)·|dE/dk|（脚本口径：g = (2/π)/|dk/dE|）。
    /// 输入段内 (E, k) 序列，输出同长度 g。
    static func dosInSegment(E: [Double], k: [Double]) -> [Double] {
        let dE = gradient(E)
        let dk = gradient(k)
        return zip(dk, dE).map { dkv, dEv in
            dEv == 0 ? .infinity : (2 / .pi) / abs(dkv / dEv)
        }
    }
}

// MARK: - 模块

/// 笔记 24《晶体中的电子与能带》：一维 δ 势梳 Kronig-Penney 能带
/// （20000 点扫描）+ 带边范霍夫奇点的一维态密度。秒级档。
struct KronigPenneyModule: SimModule {

    let meta = ModuleMeta(
        id: "晶体中的电子_能带与KronigPenney", title: "晶体中的电子 · Kronig-Penney 能带",
        subtitle: "δ 势梳能带方程扫描 + 1D 态密度范霍夫奇点",
        category: .manyBody, noteNumber: 24, tier: .seconds, difficulty: .basic,
        keywords: ["能带", "Kronig-Penney", "晶体", "布里渊区", "禁带", "允许带",
                   "态密度", "范霍夫", "van Hove", "布拉格反射"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "P", title: "δ 势强度", symbol: "P", unit: "",
                               range: 0.2...5, defaultValue: 1, decimalPlaces: 2)),
            .slider(SliderSpec(key: "Emax", title: "能量上限", symbol: "E_max",
                               unit: "eV", range: 40...400, defaultValue: 160,
                               scale: .log, decimalPlaces: 0)),
        ]
    }

    var charts: [ChartSpec] {
        [
            .lineSeries(LineSeriesSpec(
                title: "Kronig-Penney 能带结构（δ 势梳，简约区）",
                xAxis: .init(label: "简约区波矢 k (Å⁻¹)"),
                yAxis: .init(label: "能量 E (eV)"),
                seriesNames: ["带 1", "带 2", "带 3", "带 4", "带 5", "带 6", "带 7", "带 8"])),
            .lineSeries(LineSeriesSpec(
                title: "一维态密度 g(E)（脚本口径）",
                xAxis: .init(label: "能量 E (eV)"),
                yAxis: .init(label: "1D 态密度 g(E) (a.u.)"),
                seriesNames: ["g(E) = (2/π)·|dE/dk|（脚本口径）"])),
        ]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let P = input.slider("P")
        let emax = input.slider("Emax")
        let a = 1.0  // 晶格常数 Å（与脚本一致）
        let N = 20000

        // ---- 段 1：能带扫描（E = linspace(1e-4, Emax, 20000)，与 Python 一致）----
        let E = Num.linspace(1e-4, emax, count: N)
        var allowed = [Bool](repeating: false, count: N)
        var kRed = [Double](repeating: 0, count: N)
        for i in E.indices {
            let f = KronigPenneyMath.f(E[i], P: P, a: a)
            allowed[i] = abs(f) <= 1
            kRed[i] = acos(min(max(f, -1), 1)) / a
        }
        try Task.checkCancellation()
        let segs = KronigPenneyMath.segments(allowed: allowed, E: E)

        // 带结构系列：每段 ±k 双支拼成 V 形折线（E(−k) = E(k)）
        let maxBands = 8
        var bandSeries: [SeriesPoints] = []
        for (idx, seg) in segs.enumerated() where idx < maxBands {
            var pts: [Point] = []
            var es: [Double] = [], ks: [Double] = []
            for i in E.indices where E[i] >= seg.lo - 1e-9 && E[i] <= seg.hi + 1e-9 {
                es.append(E[i]); ks.append(kRed[i])
            }
            let stride = max(1, (es.count + 511) / 512)
            var i = es.count - 1
            while i >= 0 {  // −k 支：−π/a → 0（E 自高而低）
                pts.append(Point(x: -ks[i], y: es[i]))
                i -= stride
            }
            for i in 0..<es.count {  // +k 支：0 → π/a（E 自低而高）
                if i % stride == 0 { pts.append(Point(x: ks[i], y: es[i])) }
            }
            bandSeries.append(.init(name: "带 \(idx + 1)", points: pts))
        }

        // ---- 段 2：一维态密度 ----
        var eAll: [Double] = [], gAll: [Double] = []
        for seg in segs {
            var es: [Double] = [], ks: [Double] = []
            for i in E.indices where E[i] >= seg.lo - 1e-9 && E[i] <= seg.hi + 1e-9 {
                es.append(E[i]); ks.append(kRed[i])
            }
            guard es.count >= 4 else { continue }
            let g = KronigPenneyMath.dosInSegment(E: es, k: ks)
            eAll.append(contentsOf: es)
            gAll.append(contentsOf: g)
        }
        try Task.checkCancellation()

        // 有限值过滤 + 98 分位截断（clip 范霍夫发散峰以可视）
        var finite: [(e: Double, g: Double)] = []
        for (e, g) in zip(eAll, gAll) where g.isFinite { finite.append((e, g)) }
        let cap = KronigPenneyMath.percentile(finite.map(\.g).sorted(), 98)
        let dosPts = finite.map { Point(x: $0.e, y: min($0.g, cap)) }
        let dosStride = max(1, (dosPts.count + 511) / 512)
        let dosDisplay = stride(from: 0, to: dosPts.count, by: dosStride).map { dosPts[$0] }

        let chart0 = LineSeriesData(
            spec: .init(xAxis: .init(label: "简约区波矢 k (Å⁻¹)"),
                        yAxis: .init(label: "能量 E (eV)"),
                        seriesNames: bandSeries.map(\.name)),
            series: bandSeries,
            referenceLines: [
                ReferenceLine(label: "k = π/a（BZ 边界）", axis: .x, value: .pi / a),
                ReferenceLine(label: "−π/a", axis: .x, value: -.pi / a, style: .subtle),
            ])

        let chart1 = LineSeriesData(
            spec: .init(xAxis: .init(label: "能量 E (eV)"),
                        yAxis: .init(label: "1D 态密度 g(E) (a.u.)"),
                        seriesNames: ["g(E) = (2/π)·|dE/dk|（脚本口径）"]),
            series: [.init(name: "g(E) = (2/π)·|dE/dk|（脚本口径）", points: dosDisplay)])

        var gapText = "—"
        if segs.count >= 2 {
            gapText = String(format: "%.3f eV", segs[1].lo - segs[0].hi)
        }
        let firstTwo = segs.prefix(2).map {
            String(format: "[%.2f, %.2f]", $0.lo, $0.hi)
        }.joined(separator: " / ")
        return SimResult(
            charts: [.lineSeries(chart0), .lineSeries(chart1)],
            summary: [
                .init(id: "bands", title: "允许带（前两带）",
                      value: firstTwo.isEmpty ? "无允许带" : firstTwo,
                      note: "E 区间 (eV)"),
                .init(id: "gap", title: "第一禁带",
                      value: gapText, note: "带 1 顶 − 带 2 底"),
                .init(id: "cap", title: "g 的 98 分位截断",
                      value: String(format: "%.4f", cap),
                      note: "带内极大裁剪（98 分位）"),
                .init(id: "count", title: "允许带数",
                      value: "\(segs.count)",
                      note: String(format: "E < %.0f eV", emax)),
            ],
            theory: TheoryCard(
                title: "Kronig-Penney 能带与态密度",
                formulas: [
                    "cos(ka) + (P/αa)sin(αa) = cos(ka)：|左端| ≤ 1 为允许带",
                    "α = √(E/(ℏ²/2m))，P 为 δ 势强度；E(−k) = E(k)",
                    "禁带 = 布拉格反射全反射区：k = π/a 处能隙张开",
                    "本图 g(E) = (2/π)/|dk/dE|（脚本口径）：带边归零、带内取峰；"
                        + "教科书范霍夫形式为 (2/π)·|dk/dE|（带边发散），两者互为倒数",
                ],
                reading: "P 越大禁带越宽（势垒反射越强）。本图态密度按脚本口径"
                    + " g ∝ |dE/dk|：带边归零、带内取峰，禁带内 g = 0（无态可占）；"
                    + "与教科书 1D 范霍夫奇点（带边 |E−E_edge|^(−1/2) 发散）互为倒数——"
                    + "系基准脚本公式疑点，已在周报偏差说明中报告。"))
    }
}
