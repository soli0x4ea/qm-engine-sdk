import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// UI 图表输出审计（log-only，常设）：
/// ① 全量 63 模块 SimResult.charts 清单（数量/图种/标题）——找「无图形输出」模块；
/// ② 参数敏感性扫描：滑杆 ±20% 量程、离散选项逐一切换、多选组增删候选，
///    量化图表数值序列的最大绝对变化 / 基准量程（视觉变化占比）——
///    找「参数变化后曲线变化不明显」的模块；
/// ③ 量纲悬殊检测：线性轴下数值序列 max/min > 1e3 时提示宜用对数轴。
///
/// 门禁：与 PerformanceAuditTests 同理（策略档共享 StrategyCore 状态），
/// `QM_CHART_AUDIT=1 swift test --filter ChartAuditTests` 显式运行。
@Suite(.serialized, .disabled(if: ProcessInfo.processInfo.environment["QM_CHART_AUDIT"] != "1",
                             "与策略档模块的共享 StrategyCore 状态断言冲突；用 QM_CHART_AUDIT=1 显式运行"))
struct ChartAuditTests {

    // MARK: - 工具

    private func kindName(_ chart: ChartData) -> String {
        switch chart {
        case .lineSeries: return "lineSeries"
        case .levelDiagram: return "levelDiagram"
        case .bars: return "bars"
        case .scatter: return "scatter"
        case .heatmap: return "heatmap"
        case .dualAxisLineSeries: return "dualAxis"
        case .schematic: return "schematic"
        case .contour: return "contour"
        case .frameStack: return "frameStack"
        }
    }

    /// 图表的数值载荷（曲线 y 值 / 柱值 / 能级 / 矩阵元）。
    /// schematic 为示意元素（无数值序列），返回空。
    private func payload(_ chart: ChartData) -> [Double] {
        switch chart {
        case .lineSeries(let d): return d.series.flatMap { $0.points.map(\.y) }
        case .levelDiagram(let d): return d.levels.map(\.energy)
        case .bars(let d): return d.bars.map(\.value)
        case .scatter(let d): return d.series.flatMap { $0.points.map(\.y) }
        case .heatmap(let d): return d.values.flatMap { $0 }
        case .dualAxisLineSeries(let d): return (d.primary + d.secondary).flatMap { $0.points.map(\.y) }
        case .schematic: return []
        case .contour(let d): return d.values.flatMap { $0 }
        case .frameStack(let d): return d.frames.flatMap { $0.values.flatMap { $0 } }
        }
    }

    /// 线性轴 + 数值全正 + 动态范围悬殊 → 提示宜对数轴。
    private func dynamicRangeNote(_ chart: ChartData) -> String? {
        let linearY: Bool
        switch chart {
        case .lineSeries(let d): linearY = d.spec.yAxis.scale == .linear
        case .dualAxisLineSeries(let d): linearY = d.spec.primaryAxis.scale == .linear
        case .bars(let d): linearY = d.spec.yAxis.scale == .linear
        case .scatter(let d): linearY = d.spec.yAxis.scale == .linear
        default: return nil
        }
        guard linearY else { return nil }
        let vals = payload(chart).filter { $0.isFinite && $0 > 0 }
        guard let mx = vals.max(), let mn = vals.min(), mn > 0, mx / mn > 1e3 else { return nil }
        return "线性轴动态范围悬殊 max/min=\(String(format: "%.1e", mx / mn))（宜对数轴或分段）"
    }

    /// 两个图表数值序列的视觉变化占比（相对基准量程）。
    private func visualScore(base: [Double], perturbed: [Double]) -> Double? {
        guard !base.isEmpty, !perturbed.isEmpty else { return nil }
        let span = (base.max() ?? 0) - (base.min() ?? 0)
        guard span > 0, span.isFinite else { return nil }
        guard base.count == perturbed.count else { return 1.0 }  // 结构变化（点数不同）＝显著
        var maxDelta = 0.0
        for (a, b) in zip(base, perturbed) {
            let d = abs(a - b)
            if d.isFinite { maxDelta = max(maxDelta, d) }
        }
        return min(1.0, maxDelta / span)
    }

    // MARK: - 审计

    @Test("全量 63 模块图表输出审计：清单 + 参数敏感性 + 量纲悬殊（log-only）")
    func chartAudit() async throws {
        let registry = ModuleRegistry.shared
        registry.removeAll()
        #expect(ModuleLibrary.registerBuiltins() == 63)
        let constants = try ConstantsSet.load(.v2022)

        var noChartModules: [String] = []
        var insensitiveModules: [String] = []
        var flaggedDetails = 0

        for module in registry.all {
            let defaults = ParamValues.defaults(for: module.params)
            let base: SimResult
            do {
                base = try await module.compute(defaults, constants: constants)
            } catch {
                print("🖼 \(module.meta.id)｜\(module.meta.title)｜⚠️ 默认参数计算失败：\(error)")
                continue
            }

            let kinds = base.charts.map(kindName)
            let titles = base.charts.compactMap(\.chartTitle)
            print("🖼 \(module.meta.id)｜\(module.meta.title)｜图 \(base.charts.count) 张：\(kinds.joined(separator: ","))｜\(titles.joined(separator: " / "))")

            if base.charts.isEmpty {
                noChartModules.append("\(module.meta.id)（\(module.meta.title)）")
                continue  // 无图则无敏感性可言
            }

            // ---- 参数敏感性扫描 ----
            var moduleBest: (score: Double, desc: String) = (0, "")
            var dynamicFlags: Set<String> = []

            func probe(_ label: String, _ values: ParamValues) async {
                guard let p = try? await module.compute(values, constants: constants) else {
                    return  // 扰动出物理有效域——跳过（该方向不可达）
                }
                for (i, c) in p.charts.enumerated() {
                    if let note = dynamicRangeNote(c), i < base.charts.count {
                        // 只在基准图上标注（扰动不改变图种与量纲性质）
                        if dynamicRangeNote(base.charts[i]) != nil {
                            dynamicFlags.insert("\(titles.count > i ? titles[i] : kindName(c))：\(note)")
                        }
                    }
                    guard i < base.charts.count else { continue }
                    if let s = visualScore(base: payload(base.charts[i]),
                                           perturbed: payload(c)), s > moduleBest.score {
                        moduleBest = (s, "\(label) → 图「\(c.chartTitle ?? kindName(c))」变化 \(String(format: "%.1f", s * 100))% 量程")
                    }
                }
            }

            for spec in module.params {
                switch spec {
                case .slider(let s):
                    let span = s.range.upperBound - s.range.lowerBound
                    for frac in [0.2, -0.2] {
                        var v = defaults
                        v.sliders[s.key] = min(s.range.upperBound,
                                               max(s.range.lowerBound, s.defaultValue + frac * span))
                        if v.sliders[s.key] != s.defaultValue {
                            await probe("滑杆 \(s.key) \(frac > 0 ? "+" : "-")20% 量程", v)
                        }
                    }
                case .discrete(let s):
                    for opt in s.options where opt.id != s.defaultOptionID {
                        var v = defaults
                        v.discretes[s.key] = opt.id
                        await probe("离散 \(s.key)→\(opt.id)", v)
                    }
                case .multiCompare(let s):
                    if let first = s.candidates.first(where: { !s.defaultSelectionIDs.contains($0.id) }) {
                        var v = defaults
                        v.multiCompares[s.key] = Set(s.defaultSelectionIDs).union([first.id])
                        await probe("多选 \(s.key)+\(first.id)", v)
                    }
                    if let drop = s.defaultSelectionIDs.first {
                        var v = defaults
                        v.multiCompares[s.key] = Set(s.defaultSelectionIDs).subtracting([drop])
                        if !v.multiCompares[s.key]!.isEmpty {
                            await probe("多选 \(s.key)-\(drop)", v)
                        }
                    }
                case .constant:
                    break
                }
            }

            for f in dynamicFlags.sorted() {
                print("   ⚠️ 量纲：\(f)")
                flaggedDetails += 1
            }
            if moduleBest.score < 0.05 {
                insensitiveModules.append("\(module.meta.id)（\(module.meta.title)）——最大视觉变化仅 \(String(format: "%.1f", moduleBest.score * 100))% 量程〔\(moduleBest.desc)〕")
            } else if moduleBest.score < 0.12 {
                print("   ⚠️ 近不敏感：\(moduleBest.desc)")
                flaggedDetails += 1
            }
        }

        print("========= 图表审计汇总 =========")
        print("无图形输出模块（\(noChartModules.count)）：")
        for m in noChartModules { print("  · \(m)") }
        print("参数不敏感模块（全部扰动 < 5% 量程，共 \(insensitiveModules.count)）：")
        for m in insensitiveModules { print("  · \(m)") }
        _ = flaggedDetails
    }
}
