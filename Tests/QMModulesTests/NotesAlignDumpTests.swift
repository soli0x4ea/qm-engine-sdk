import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// W13d 笔记对齐审计 · 模拟侧图表清单导出（log-only，不进常规红绿）：
/// `QM_NOTES_DUMP=1 swift test --filter NotesAlignDump --disable-sandbox`
/// JSON 写入 `$QM_NOTES_DUMP_PATH`（默认 build/UI图表审计/sim_charts_dump.json）。
@Suite(.serialized, .disabled(if: ProcessInfo.processInfo.environment["QM_NOTES_DUMP"] != "1",
                             "导出工具：用 QM_NOTES_DUMP=1 显式运行"))
struct NotesAlignDumpTests {

    private func axis(_ a: AxisSpec?) -> [String: Any]? {
        guard let a else { return nil }
        return ["label": a.label, "scale": a.scale.rawValue]
    }

    private func range(_ vs: [Double]) -> [String: Any]? {
        let finite = vs.filter { $0.isFinite }
        guard let lo = finite.min(), let hi = finite.max() else { return nil }
        return ["min": lo, "max": hi, "n": vs.count]
    }

    private func seriesInfo(_ s: [SeriesPoints]) -> [[String: Any]] {
        s.map { series -> [String: Any] in
            var d: [String: Any] = [
                "name": series.name,
                "n": series.points.count,
            ]
            if let xr = range(series.points.map(\.x)) { d["x"] = xr }
            if let yr = range(series.points.map(\.y)) { d["y"] = yr }
            return d
        }
    }

    /// 单图 → JSON dict
    private func chartDict(_ chart: ChartData) -> [String: Any] {
        var d: [String: Any] = [:]
        switch chart {
        case .lineSeries(let c):
            d["kind"] = "lineSeries"
            d["title"] = c.spec.title ?? ""
            d["xAxis"] = axis(c.spec.xAxis) ?? [:]
            d["yAxis"] = axis(c.spec.yAxis) ?? [:]
            d["seriesNames"] = c.spec.seriesNames
            d["series"] = seriesInfo(c.series)
            if !c.referenceLines.isEmpty {
                d["refLines"] = c.referenceLines.map { "\($0.axis.rawValue)=\($0.value) \($0.label)" }
            }
            if !c.pointMarkers.isEmpty {
                d["markers"] = c.pointMarkers.compactMap(\.label)
            }
        case .dualAxisLineSeries(let c):
            d["kind"] = "dualAxis"
            d["title"] = c.spec.title ?? ""
            d["xAxis"] = axis(c.spec.xAxis) ?? [:]
            d["yAxis"] = axis(c.spec.primaryAxis) ?? [:]
            d["y2Axis"] = axis(c.spec.secondaryAxis) ?? [:]
            d["seriesNames"] = c.spec.primaryNames + c.spec.secondaryNames
            d["series"] = seriesInfo(c.primary + c.secondary)
            if !c.referenceLines.isEmpty {
                d["refLines"] = c.referenceLines.map { "\($0.axis.rawValue)=\($0.value) \($0.label)" }
            }
        case .bars(let c):
            d["kind"] = "bars"
            d["title"] = c.spec.title ?? ""
            d["xAxis"] = axis(c.spec.xAxis) ?? [:]
            d["yAxis"] = axis(c.spec.yAxis) ?? [:]
            d["seriesNames"] = c.spec.seriesNames
            d["bars"] = c.bars.map { ["label": $0.label, "value": $0.value, "series": $0.series ?? ""] }
        case .scatter(let c):
            d["kind"] = "scatter"
            d["title"] = c.spec.title ?? ""
            d["xAxis"] = axis(c.spec.xAxis) ?? [:]
            d["yAxis"] = axis(c.spec.yAxis) ?? [:]
            d["seriesNames"] = c.spec.seriesNames
            d["series"] = seriesInfo(c.series)
            if !c.referenceLines.isEmpty {
                d["refLines"] = c.referenceLines.map { "\($0.axis.rawValue)=\($0.value) \($0.label)" }
            }
        case .levelDiagram(let c):
            d["kind"] = "levelDiagram"
            d["title"] = c.spec.title ?? ""
            d["energyAxis"] = c.spec.energyAxis
            d["levels"] = c.levels.map { ["label": $0.label, "energy": $0.energy] }
            d["transitions"] = c.transitions.compactMap(\.label)
        case .heatmap(let c):
            d["kind"] = "heatmap"
            d["title"] = c.spec.title ?? ""
            d["xAxis"] = axis(c.spec.xAxis) ?? [:]
            d["yAxis"] = axis(c.spec.yAxis) ?? [:]
            d["valueLabel"] = c.spec.valueLabel
            d["xTicks"] = []   // ticks 在 Data 层，spec 无；下方 data 补
        case .contour(let c):
            d["kind"] = "contour"
            d["title"] = c.spec.title ?? ""
            d["xAxis"] = axis(c.spec.xAxis) ?? [:]
            d["yAxis"] = axis(c.spec.yAxis) ?? [:]
            d["valueLabel"] = c.spec.valueLabel
        case .schematic(let c):
            d["kind"] = "schematic"
            d["title"] = c.spec.title ?? ""
            d["xAxis"] = axis(c.spec.xAxis) ?? [:]
            d["yAxis"] = axis(c.spec.yAxis) ?? [:]
        case .frameStack(let c):
            d["kind"] = "frameStack"
            d["title"] = c.spec.title ?? ""
            d["xAxis"] = axis(c.spec.xAxis) ?? [:]
            d["yAxis"] = axis(c.spec.yAxis) ?? [:]
            d["valueLabel"] = c.spec.valueLabel
        case .bloch(let c):
            d["kind"] = "bloch"
            d["title"] = c.spec.title ?? ""
        }
        return d
    }

    /// 图 + 数据合并（spec 与 data 在 ChartData 内已合一，此处仅补 data 专属字段）
    private func enriched(_ chart: ChartData, _ d0: [String: Any]) -> [String: Any] {
        var d = d0
        switch chart {
        case .heatmap(let c):
            d["xTicks"] = c.xTicks
            d["yTicks"] = c.yTicks
            if let r = range(c.values.flatMap { $0 }) { d["valueRange"] = r }
        case .contour(let c):
            if let r = range(c.xGrid) { d["xGrid"] = r }
            if let r = range(c.yGrid) { d["yGrid"] = r }
            if let r = range(c.values.flatMap { $0 }) { d["valueRange"] = r }
        case .frameStack(let c):
            d["xTicks"] = c.xTicks
            d["yTicks"] = c.yTicks
            d["frames"] = c.frames.map(\.label)
        case .levelDiagram(let c):
            d["levels"] = c.levels.map { ["label": $0.label, "energy": $0.energy] }
        case .bloch(let c):
            d["stateLabel"] = c.stateLabel ?? ""
            d["north"] = c.northLabel
            d["south"] = c.southLabel
            d["trajectoryN"] = c.trajectory.count
        default:
            break
        }
        return d
    }

    @Test("导出全模块图表清单 JSON（log-only）")
    func dump() async throws {
        let registry = ModuleRegistry.shared
        registry.removeAll()
        #expect(ModuleLibrary.registerBuiltins() == 63)
        let constants = try ConstantsSet.load(.v2022)

        var modules: [[String: Any]] = []
        for module in registry.all {
            let defaults = ParamValues.defaults(for: module.params)
            let result: SimResult
            do {
                result = try await module.compute(defaults, constants: constants)
            } catch {
                modules.append(["id": module.meta.id, "title": module.meta.title,
                                "noteNumber": module.meta.noteNumber,
                                "error": "\(error)"])
                continue
            }
            let charts = result.charts.map { enriched($0, chartDict($0)) }
            modules.append([
                "id": module.meta.id,
                "title": module.meta.title,
                "noteNumber": module.meta.noteNumber,
                "summary": result.summary.map { "\($0.title)=\($0.value)" },
                "charts": charts,
            ])
        }

        let json = try JSONSerialization.data(
            withJSONObject: ["modules": modules],
            options: [.prettyPrinted, .sortedKeys, .withoutEscapingSlashes])
        let outPath = ProcessInfo.processInfo.environment["QM_NOTES_DUMP_PATH"]
            ?? "/Users/soli/Documents/trae/量子力学正式版/build/UI图表审计/sim_charts_dump.json"
        try json.write(to: URL(fileURLWithPath: outPath))
        print("📝 dumped \(modules.count) modules → \(outPath)")
    }
}
