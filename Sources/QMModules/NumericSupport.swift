import Foundation
import EngineKit

/// QMModules 内部数值工具——与 Python 脚本的 numpy 用法逐句对应。
enum Num {

    /// np.linspace：闭区间等距 count 点。
    /// 光学能量↔波长换算常数 hc（eV·nm，CODATA）；W13k 审核上收——原 4 个模块各抄一份。
    static let hcEVnm: Double = 1239.841984

    /// 均匀网格梯形积分（W13k 审核上收——原 3 个模块各抄一份，且其中一份误为矩形和）。
    static func trapz(_ y: [Double], dx: Double) -> Double {
        guard y.count >= 2 else { return 0 }
        var s = 0.5 * (y[0] + y[y.count - 1])
        for v in y.dropFirst().dropLast() { s += v }
        return s * dx
    }

    static func linspace(_ from: Double, _ to: Double, count: Int) -> [Double] {
        precondition(count >= 2, "linspace 需要 count >= 2")
        return (0..<count).map { from + (to - from) * Double($0) / Double(count - 1) }
    }

    /// np.logspace：十进制对数均匀 count 点（跨多个量级的物理量网格）。
    static func logspace(_ from: Double, _ to: Double, count: Int) -> [Double] {
        precondition(count >= 2, "logspace 需要 count >= 2")
        let logs = linspace(log10(from), log10(to), count: count)
        return logs.map { pow(10, $0) }
    }

    /// 梯形积分（非均匀网格通用，对应量子谐振子脚本的 trapz）。
    static func trapezoid(_ xs: [Double], _ ys: [Double]) -> Double {
        guard xs.count == ys.count, xs.count > 1 else { return 0 }
        var s = 0.0
        for i in 1..<xs.count {
            s += 0.5 * (ys[i - 1] + ys[i]) * (xs[i] - xs[i - 1])
        }
        return s
    }

    /// ln(k!)（对数域防大阶溢出）。
    static func lnFactorial(_ k: Int) -> Double {
        guard k > 1 else { return 0 }
        var s = 0.0
        for i in 2...k { s += log(Double(i)) }
        return s
    }

    /// 曲线抽稀（仅图表显示用；fixtures 对照在测试中按参考 x 网格直接闭式求值）。
    static func strided(_ xs: [Double], _ ys: [Double], stride: Int) -> [Point] {
        precondition(xs.count == ys.count)
        let step = max(1, stride)
        var out: [Point] = []
        out.reserveCapacity((xs.count + step - 1) / step)
        var i = 0
        while i < xs.count {
            out.append(Point(x: xs[i], y: ys[i]))
            i += step
        }
        return out
    }
}

extension ChartSpec {
    var lineSeriesSpec: LineSeriesSpec? {
        if case .lineSeries(let s) = self { return s }
        return nil
    }

    var barSpec: BarSpec? {
        if case .bars(let s) = self { return s }
        return nil
    }

    var scatterSpec: ScatterSpec? {
        if case .scatter(let s) = self { return s }
        return nil
    }

    var dualAxisLineSeriesSpec: DualAxisLineSeriesSpec? {
        if case .dualAxisLineSeries(let s) = self { return s }
        return nil
    }

    var schematicSpec: SchematicSpec? {
        if case .schematic(let s) = self { return s }
        return nil
    }

    var contourSpec: ContourSpec? {
        if case .contour(let s) = self { return s }
        return nil
    }

    var levelDiagramSpec: LevelDiagramSpec? {
        if case .levelDiagram(let s) = self { return s }
        return nil
    }

    var blochSpec: BlochSpec? {
        if case .bloch(let s) = self { return s }
        return nil
    }
}
