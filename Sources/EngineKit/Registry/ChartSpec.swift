import Foundation

// MARK: - 轴与参考线

/// 坐标轴声明。
public struct AxisSpec: Sendable, Equatable {
    /// 轴标签（含单位，如 "λ (nm)"、"E (E_ℏω)"）
    public let label: String
    public let scale: ParamScale

    public init(label: String, scale: ParamScale = .linear) {
        self.label = label
        self.scale = scale
    }
}

/// 计算产出的参考线（动态值，区别于静态声明的轴）。
public enum RefAxis: String, Sendable { case x, y }

/// 参考线视觉档（W5 追加，默认值向后兼容）：
/// threshold——阈值/峰值标注（虚线醒目，如 "λmax = 501 nm"）；
/// subtle——辅助基线（点线弱化，如 BB84 的 R = 0 零线）。
public enum ReferenceLineStyle: String, Sendable, Equatable {
    case threshold
    case subtle
}

public struct ReferenceLine: Sendable, Equatable, Identifiable {
    public let id: String
    /// 标注文案（如 "λmax = 501 nm"）
    public let label: String
    public let axis: RefAxis
    public let value: Double
    public let style: ReferenceLineStyle

    public init(id: String? = nil, label: String, axis: RefAxis, value: Double,
                style: ReferenceLineStyle = .threshold) {
        self.id = id ?? label
        self.label = label
        self.axis = axis
        self.value = value
        self.style = style
    }
}

// MARK: - 图表声明（ChartSpec）

/// ① 曲线族（Swift Charts LineMark）。
public struct LineSeriesSpec: Sendable, Equatable {
    public let title: String?
    public let xAxis: AxisSpec
    public let yAxis: AxisSpec
    /// 图例系列名（数据按此顺序着色）
    public let seriesNames: [String]

    public init(title: String? = nil, xAxis: AxisSpec, yAxis: AxisSpec,
                seriesNames: [String]) {
        self.title = title
        self.xAxis = xAxis
        self.yAxis = yAxis
        self.seriesNames = seriesNames
    }
}

/// ② 能级图（Canvas：能级横线 + 跃迁箭头）。
public struct LevelDiagramSpec: Sendable, Equatable {
    public let title: String?
    /// 能量轴标签（如 "E (eV)"、"E (E_ℏω)"）
    public let energyAxis: String
    /// 是否绘制跃迁箭头（如巴尔末系的发射线）
    public let showTransitions: Bool

    public init(title: String? = nil, energyAxis: String,
                showTransitions: Bool = true) {
        self.title = title
        self.energyAxis = energyAxis
        self.showTransitions = showTransitions
    }
}

/// ③ 柱状 / 谱线图（Swift Charts BarMark）。
public struct BarSpec: Sendable, Equatable {
    public let title: String?
    public let xAxis: AxisSpec
    public let yAxis: AxisSpec
    /// 分组系列名（可选，用于并排柱）
    public let seriesNames: [String]

    public init(title: String? = nil, xAxis: AxisSpec, yAxis: AxisSpec,
                seriesNames: [String] = []) {
        self.title = title
        self.xAxis = xAxis
        self.yAxis = yAxis
        self.seriesNames = seriesNames
    }
}

/// ④ 散点图（Swift Charts PointMark）。
public struct ScatterSpec: Sendable, Equatable {
    public let title: String?
    public let xAxis: AxisSpec
    public let yAxis: AxisSpec
    public let seriesNames: [String]

    public init(title: String? = nil, xAxis: AxisSpec, yAxis: AxisSpec,
                seriesNames: [String]) {
        self.title = title
        self.xAxis = xAxis
        self.yAxis = yAxis
        self.seriesNames = seriesNames
    }
}

/// ⑤ 热图 / 等高线底图（Canvas，W4 追加：2×2 密度矩阵起步，后续势阱/能带复用）。
public struct HeatmapSpec: Sendable, Equatable {
    public let title: String?
    public let xAxis: AxisSpec
    public let yAxis: AxisSpec
    /// 色标含义（如 "|ρ|"、"P(x,t)"），渲染成色标条标注
    public let valueLabel: String
    /// 色标是否对称发散（正负值以 0 为中心双色发散，如 Re ρ 元）
    public let diverging: Bool

    public init(title: String? = nil, xAxis: AxisSpec, yAxis: AxisSpec,
                valueLabel: String, diverging: Bool = false) {
        self.title = title
        self.xAxis = xAxis
        self.yAxis = yAxis
        self.valueLabel = valueLabel
        self.diverging = diverging
    }
}

/// ⑥ 双 Y 轴曲线族（W5 追加：量子霍尔 ρ_xy 平台 + ρ_xx dip 起步，matplotlib twinx 对应）。
/// 主轴（leading）与副轴（trailing）各携带一组系列；X 轴共享（compute 侧保证同网格）。
public struct DualAxisLineSeriesSpec: Sendable, Equatable {
    public let title: String?
    public let xAxis: AxisSpec
    public let primaryAxis: AxisSpec
    public let secondaryAxis: AxisSpec
    /// 主轴系列名（leading 侧，如 "ρ_xy/R_K"）
    public let primaryNames: [String]
    /// 副轴系列名（trailing 侧，如 "ρ_xx"）
    public let secondaryNames: [String]

    public init(title: String? = nil, xAxis: AxisSpec,
                primaryAxis: AxisSpec, secondaryAxis: AxisSpec,
                primaryNames: [String], secondaryNames: [String]) {
        self.title = title
        self.xAxis = xAxis
        self.primaryAxis = primaryAxis
        self.secondaryAxis = secondaryAxis
        self.primaryNames = primaryNames
        self.secondaryNames = secondaryNames
    }
}

/// ⑦ 示意图（W6 追加：Canvas 绘制能带填充/能级横线/标记点/箭头标注——
/// 对产生阈值 + 狄拉克海起步，Jablonski 能级复用；坐标为真实数值标度）。
public struct SchematicSpec: Sendable, Equatable {
    public let title: String?
    public let xAxis: AxisSpec
    /// 能量轴（y）标签，如 "E (MeV)"
    public let yAxis: AxisSpec

    public init(title: String? = nil, xAxis: AxisSpec, yAxis: AxisSpec) {
        self.title = title
        self.xAxis = xAxis
        self.yAxis = yAxis
    }
}

/// ⑧ 等高线图（W6 追加：从热图扩展——标量场等值线 + 分层填充色带，
/// 压缩真空 Wigner 320×320 起步；matplotlib contourf/contour 对应）。
public struct ContourSpec: Sendable, Equatable {
    public let title: String?
    public let xAxis: AxisSpec
    public let yAxis: AxisSpec
    /// 色标含义（如 "W(x, p)"）
    public let valueLabel: String
    /// 等值线级数（分层带数，如 12）
    public let levels: Int
    /// 色标是否对称发散（正负值以 0 为中心双色发散）
    public let diverging: Bool

    public init(title: String? = nil, xAxis: AxisSpec, yAxis: AxisSpec,
                valueLabel: String, levels: Int = 12, diverging: Bool = false) {
        self.title = title
        self.xAxis = xAxis
        self.yAxis = yAxis
        self.valueLabel = valueLabel
        self.levels = levels
        self.diverging = diverging
    }
}

/// ⑨ 帧栈动画（W9 追加：热图扩展——同规格小矩阵的多帧序列，
/// 密度矩阵演化/退相干轨迹起步；matplotlib animation 对应）。
public struct FrameStackSpec: Sendable, Equatable {
    public let title: String?
    public let xAxis: AxisSpec
    public let yAxis: AxisSpec
    /// 色标含义（如 "Re ρ"）
    public let valueLabel: String
    /// 色标是否对称发散
    public let diverging: Bool
    /// 自动翻帧间隔（秒，默认 0.4；渲染层 TimelineView 周期）
    public let frameInterval: Double

    public init(title: String? = nil, xAxis: AxisSpec, yAxis: AxisSpec,
                valueLabel: String, diverging: Bool = false,
                frameInterval: Double = 0.4) {
        self.title = title
        self.xAxis = xAxis
        self.yAxis = yAxis
        self.valueLabel = valueLabel
        self.diverging = diverging
        self.frameInterval = frameInterval
    }
}

/// 图表声明四类封装（W2 冻结；W4+ 追加热图/等高线、示意图、Bloch 球、帧栈动画）。
/// 声明只含静态元数据；数据由 SimResult.charts 携带，按下标一一对应。
public enum ChartSpec: Sendable, Equatable {
    case lineSeries(LineSeriesSpec)
    case levelDiagram(LevelDiagramSpec)
    case bars(BarSpec)
    case scatter(ScatterSpec)
    case heatmap(HeatmapSpec)   // W4 追加（契约仅追加）
    case dualAxisLineSeries(DualAxisLineSeriesSpec)   // W5 追加（契约仅追加）
    case schematic(SchematicSpec)   // W6 追加（契约仅追加）
    case contour(ContourSpec)   // W6 追加（契约仅追加）
    case frameStack(FrameStackSpec)   // W9 追加（契约仅追加）
}
