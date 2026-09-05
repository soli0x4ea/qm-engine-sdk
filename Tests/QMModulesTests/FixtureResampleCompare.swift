import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// 通用 fixture 重采样对拍 helper（W7 引入，后续各周复用，故取周无关名）。
///
/// 多数模块的可视分辨率（实时档/秒级档点数）与 fixtures 管线产物（通常 512 点）不同，
/// 故采用「把 fixture 曲线重采样到模块 x 网格」的值级对拍，而非逐点同网格比较。
/// 容差按曲线空间选取：log-log 曲线用 log 空间容差，线性曲线用相对容差（含地板值防除零）。
///
/// 两条使用约定：
/// - **判据级对照**用 `expectPointwiseClose`：把移植后的计算核在 fixture 自身 x 网格上求值，
///   与 fixture y 逐点比较（无插值，容差可到 1e-9 量级）——这是最强的等价性证据。
/// - **结构级对照**用 `expectResampled`：把模块 compute 产出的曲线与 fixture 比对。
///   当 fixture 网格不足以分辨曲线结构（窄洛伦兹峰等）时，线性插值会在峰芯引入
///   显著误差，此时用 `skip` 跳过峰芯，只比较平滑段，并在测试名/注释中写明理由。

// MARK: - 图表数据提取（测试侧按声明序取，取不到返回 nil 交给 #require 报错）

func chartLineSeries(_ result: SimResult, _ index: Int) -> LineSeriesData? {
    guard index < result.charts.count, case .lineSeries(let c) = result.charts[index] else { return nil }
    return c
}

func chartBars(_ result: SimResult, _ index: Int) -> BarData? {
    guard index < result.charts.count, case .bars(let c) = result.charts[index] else { return nil }
    return c
}

func chartScatter(_ result: SimResult, _ index: Int) -> ScatterData? {
    guard index < result.charts.count, case .scatter(let c) = result.charts[index] else { return nil }
    return c
}

func chartDualAxis(_ result: SimResult, _ index: Int) -> DualAxisLineSeriesData? {
    guard index < result.charts.count,
          case .dualAxisLineSeries(let c) = result.charts[index] else { return nil }
    return c
}

func chartSchematic(_ result: SimResult, _ index: Int) -> SchematicData? {
    guard index < result.charts.count, case .schematic(let c) = result.charts[index] else { return nil }
    return c
}

func chartLevelDiagram(_ result: SimResult, _ index: Int) -> LevelDiagramData? {
    guard index < result.charts.count, case .levelDiagram(let c) = result.charts[index] else { return nil }
    return c
}

/// 升序 x 数组线性插值（超出范围夹到端点）。
func linterp(_ xs: [Double], _ ys: [Double], _ x: Double) -> Double? {
    guard xs.count >= 2, ys.count == xs.count else { return nil }
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

/// 把 fixture 曲线 (fx, fy) 重采样到模块 x 网格并逐点比较。
/// - logSpace: true 时在 (log x, log y) 空间插值，适合跨多数量级曲线（如 log-log）。
/// - skipOutside: true 时跳过模块 x 超出 fixture 范围的样本（只比重叠区）。
/// - skip: 可选谓词，返回 true 的模块样本点不参与比较（如跳过欠采样的峰芯）。
func expectResampled(
    moduleX: [Double], moduleY: [Double],
    fx: [Double], fy: [Double],
    tolerance: Double, _ ctx: String,
    logSpace: Bool = false, skipOutside: Bool = true,
    skip: ((Double) -> Bool)? = nil,
    sourceLocation: SourceLocation = #_sourceLocation
) {
    #expect(moduleX.count == moduleY.count, "\(ctx): 模块 x/y 点数不一致", sourceLocation: sourceLocation)
    guard let xlo = fx.first, let xhi = fx.last else { return }
    let fxw = logSpace ? fx.map { log10(max($0, 1e-300)) } : fx
    let fyw = logSpace ? fy.map { log10(max(abs($0), 1e-300)) } : fy
    for i in moduleX.indices {
        let x = moduleX[i]
        if skipOutside && (x < xlo || x > xhi) { continue }
        if let skip, skip(x) { continue }
        let xw = logSpace ? log10(max(x, 1e-300)) : x
        guard let yfw = linterp(fxw, fyw, xw) else { continue }
        let myY = moduleY[i]
        if logSpace {
            let myYw = log10(max(abs(myY), 1e-300))
            #expect(abs(myYw - yfw) < tolerance,
                    "\(ctx)@x≈\(String(format: "%.3e", x)): log y \(String(format: "%.4f", myYw)) vs \(String(format: "%.4f", yfw))",
                    sourceLocation: sourceLocation)
        } else {
            let denom = max(abs(yfw), 1e-6)
            #expect(abs(myY - yfw) / denom < tolerance,
                    "\(ctx)@x≈\(String(format: "%.3e", x)): \(String(format: "%.4e", myY)) vs \(String(format: "%.4e", yfw))",
                    sourceLocation: sourceLocation)
        }
    }
}
