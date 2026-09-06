import Foundation

// MARK: - 基础几何/数据点

/// 二维数据点（显示单位）。
public struct Point: Sendable, Equatable {
    public let x: Double
    public let y: Double

    public init(x: Double, y: Double) {
        self.x = x
        self.y = y
    }
}

/// 一条命名曲线/散点系列。colorIndex 为 nil 时按系列序取色板。
public struct SeriesPoints: Sendable, Equatable, Identifiable {
    public let id: String
    public let name: String
    public let points: [Point]
    public let colorIndex: Int?

    public init(id: String? = nil, name: String, points: [Point],
                colorIndex: Int? = nil) {
        self.id = id ?? name
        self.name = name
        self.points = points
        self.colorIndex = colorIndex
    }
}

// MARK: - 图表数据（与 ChartSpec 一一对应）

/// 曲线下方区域填充（W5 追加，matplotlib fill_between 对应）：
/// points 为上边界（compute 侧已按 where 条件裁剪，如 BB84 安全区 Q ≤ 11%），
/// baseline 为下边界常量（通常 0）。
public struct AreaFill: Sendable, Equatable, Identifiable {
    public let id: String
    public let name: String
    public let points: [Point]
    public let baseline: Double
    public let colorIndex: Int?

    public init(id: String? = nil, name: String, points: [Point],
                baseline: Double = 0, colorIndex: Int? = nil) {
        self.id = id ?? name
        self.name = name
        self.points = points
        self.baseline = baseline
        self.colorIndex = colorIndex
    }
}

// MARK: - 工作点标记（W13 追加，B2/B3 复用）

/// 曲线上的当前工作点标记（参数对应的态/元素/阈值交点）：
/// 实心圆点 + 可选标注胶囊，LineSeriesData.pointMarkers 渲染。
public struct PointMarker: Sendable, Equatable, Identifiable {
    public let id: String
    public let x: Double
    public let y: Double
    public let label: String?
    /// 色板序（nil 用强调色）
    public let colorIndex: Int?

    public init(id: String? = nil, x: Double, y: Double,
                label: String? = nil, colorIndex: Int? = nil) {
        self.id = id ?? "marker(\(x),\(y))"
        self.x = x
        self.y = y
        self.label = label
        self.colorIndex = colorIndex
    }
}

public struct LineSeriesData: Sendable, Equatable {
    public let spec: LineSeriesSpec
    /// 与 spec.seriesNames 对齐
    public let series: [SeriesPoints]
    public let referenceLines: [ReferenceLine]
    /// 区域填充（W5 追加，默认空——既有模块零改动）
    public let areaFills: [AreaFill]
    /// 工作点标记（W13 追加，默认空——普适曲线族上的当前态/交点，B2/B3）
    public let pointMarkers: [PointMarker]

    public init(spec: LineSeriesSpec, series: [SeriesPoints],
                referenceLines: [ReferenceLine] = [],
                areaFills: [AreaFill] = [],
                pointMarkers: [PointMarker] = []) {
        self.spec = spec
        self.series = series
        self.referenceLines = referenceLines
        self.areaFills = areaFills
        self.pointMarkers = pointMarkers
    }
}

/// 单个能级（横线）。
public struct Level: Sendable, Equatable, Identifiable {
    public let id: String
    /// 能级标签（如 "n=0"、"$E_1$"）
    public let label: String
    public let energy: Double

    public init(id: String? = nil, label: String, energy: Double) {
        self.id = id ?? label
        self.label = label
        self.energy = energy
    }
}

/// 能级间跃迁箭头。
public struct Transition: Sendable, Equatable, Identifiable {
    public let id: String
    /// levels 数组下标
    public let fromIndex: Int
    public let toIndex: Int
    /// 箭头标注（如 "656 nm"）
    public let label: String?

    public init(id: String? = nil, fromIndex: Int, toIndex: Int,
                label: String? = nil) {
        self.id = id ?? "\(fromIndex)-\(toIndex)"
        self.fromIndex = fromIndex
        self.toIndex = toIndex
        self.label = label
    }
}

public struct LevelDiagramData: Sendable, Equatable {
    public let spec: LevelDiagramSpec
    /// 能量升序排列（渲染层按 spec 决定是否倒置，如束缚态负能级）
    public let levels: [Level]
    public let transitions: [Transition]

    public init(spec: LevelDiagramSpec, levels: [Level],
                transitions: [Transition] = []) {
        self.spec = spec
        self.levels = levels
        self.transitions = transitions
    }
}

/// 单根柱（谱线柱/统计柱）。
public struct BarItem: Sendable, Equatable, Identifiable {
    public let id: String
    /// x 轴类别标签（如 "n=0"、"656 nm"）
    public let label: String
    public let value: Double
    /// 分组系列名（并排柱时非空）
    public let series: String?

    public init(id: String? = nil, label: String, value: Double,
                series: String? = nil) {
        self.id = id ?? (series.map { "\($0)/\(label)" } ?? label)
        self.label = label
        self.value = value
        self.series = series
    }
}

public struct BarData: Sendable, Equatable {
    public let spec: BarSpec
    public let bars: [BarItem]

    public init(spec: BarSpec, bars: [BarItem]) {
        self.spec = spec
        self.bars = bars
    }
}

public struct ScatterData: Sendable, Equatable {
    public let spec: ScatterSpec
    public let series: [SeriesPoints]
    /// 参考线（W5 追加，默认空——既有模块零改动；EPR 实验时间线的 LHV/Tsirelson 界用）
    public let referenceLines: [ReferenceLine]

    public init(spec: ScatterSpec, series: [SeriesPoints],
                referenceLines: [ReferenceLine] = []) {
        self.spec = spec
        self.series = series
        self.referenceLines = referenceLines
    }
}

/// 热图数据（W4 追加）：行主序 values[y][x]，轴刻度与网格对齐。
public struct HeatmapData: Sendable, Equatable {
    public let spec: HeatmapSpec
    /// x 轴刻度标签（列；如 ["↑↑","↑↓","↓↑","↓↓"]）
    public let xTicks: [String]
    /// y 轴刻度标签（行）
    public let yTicks: [String]
    /// 行主序矩阵：values.count == yTicks.count，每行 .count == xTicks.count
    public let values: [[Double]]

    public init(spec: HeatmapSpec, xTicks: [String], yTicks: [String],
                values: [[Double]]) {
        self.spec = spec
        self.xTicks = xTicks
        self.yTicks = yTicks
        self.values = values
    }
}

/// 双 Y 轴曲线族数据（W5 追加）：主/副轴系列各一组，X 共享。
public struct DualAxisLineSeriesData: Sendable, Equatable {
    public let spec: DualAxisLineSeriesSpec
    /// 与 spec.primaryNames 对齐（leading 轴）
    public let primary: [SeriesPoints]
    /// 与 spec.secondaryNames 对齐（trailing 轴）
    public let secondary: [SeriesPoints]
    public let referenceLines: [ReferenceLine]

    public init(spec: DualAxisLineSeriesSpec, primary: [SeriesPoints],
                secondary: [SeriesPoints], referenceLines: [ReferenceLine] = []) {
        self.spec = spec
        self.primary = primary
        self.secondary = secondary
        self.referenceLines = referenceLines
    }
}

// MARK: - 示意图数据（W6 追加）

/// 示意元素：矩形区域（能带 / 禁区 / 填充海），数值坐标下的真实矩形。
public struct SchematicBand: Sendable, Equatable, Identifiable {
    public let id: String
    public let label: String?
    public let xRange: ClosedRange<Double>
    public let yRange: ClosedRange<Double>
    /// false = 仅描边（如禁区框），true = 填充（如负能海）
    public let filled: Bool
    /// 色板序（nil 按出现序取色）
    public let colorIndex: Int?

    public init(id: String? = nil, label: String? = nil,
                xRange: ClosedRange<Double>, yRange: ClosedRange<Double>,
                filled: Bool = true, colorIndex: Int? = nil) {
        self.id = id ?? label ?? "band(\(xRange.lowerBound)…\(yRange.lowerBound))"
        self.label = label
        self.xRange = xRange
        self.yRange = yRange
        self.filled = filled
        self.colorIndex = colorIndex
    }
}

/// 示意元素：标记点（电子态 / 空穴 / 粒子）。
public struct SchematicMarker: Sendable, Equatable, Identifiable {
    public enum Shape: String, Sendable, Equatable { case circle, square }
    public let id: String
    public let label: String?
    public let x: Double
    public let y: Double
    public let shape: Shape
    /// false = 空心（如狄拉克海的洞）
    public let filled: Bool
    public let colorIndex: Int?

    public init(id: String? = nil, label: String? = nil,
                x: Double, y: Double, shape: Shape = .circle,
                filled: Bool = true, colorIndex: Int? = nil) {
        self.id = id ?? label ?? "marker(\(x),\(y))"
        self.label = label
        self.x = x
        self.y = y
        self.shape = shape
        self.filled = filled
        self.colorIndex = colorIndex
    }
}

/// 示意元素：箭头标注（如 "γ + Z → e⁺ + e⁻ + Z"，anchor 处文字 + 可选箭头）。
public struct SchematicCallout: Sendable, Equatable, Identifiable {
    public let id: String
    public let text: String
    /// 文字锚点
    public let anchor: Point
    /// 箭头终点（nil = 仅文字）
    public let arrowEnd: Point?
    public let colorIndex: Int?

    public init(id: String? = nil, text: String, anchor: Point,
                arrowEnd: Point? = nil, colorIndex: Int? = nil) {
        self.id = id ?? text
        self.text = text
        self.anchor = anchor
        self.arrowEnd = arrowEnd
        self.colorIndex = colorIndex
    }
}

/// 示意图数据（W6 追加）：bands / markers / callouts / 能级横线（复用 ReferenceLine）。
public struct SchematicData: Sendable, Equatable {
    public let spec: SchematicSpec
    public let bands: [SchematicBand]
    public let markers: [SchematicMarker]
    public let callouts: [SchematicCallout]
    /// 水平 / 垂直参考线（费米能级、阈值线等）
    public let referenceLines: [ReferenceLine]

    public init(spec: SchematicSpec, bands: [SchematicBand] = [],
                markers: [SchematicMarker] = [], callouts: [SchematicCallout] = [],
                referenceLines: [ReferenceLine] = []) {
        self.spec = spec
        self.bands = bands
        self.markers = markers
        self.callouts = callouts
        self.referenceLines = referenceLines
    }
}

// MARK: - 等高线数据（W6 追加）

/// 等高线图数据：规则网格标量场，行主序 values[y][x]。
/// 网格点为采样值（非格心），等值线在网格点之间插值（matplotlib contour 语义）。
public struct ContourData: Sendable, Equatable {
    public let spec: ContourSpec
    /// x 轴网格（升序，列，如 320 点）
    public let xGrid: [Double]
    /// y 轴网格（升序，行，如 320 点）
    public let yGrid: [Double]
    /// 行主序矩阵：values.count == yGrid.count，每行 .count == xGrid.count
    public let values: [[Double]]
    /// 数据侧高亮等值线（强调色加粗，如压缩真空 Wigner 的 1σ 圈；默认空）
    public let highlightLevels: [Double]

    public init(spec: ContourSpec, xGrid: [Double], yGrid: [Double],
                values: [[Double]], highlightLevels: [Double] = []) {
        self.spec = spec
        self.xGrid = xGrid
        self.yGrid = yGrid
        self.values = values
        self.highlightLevels = highlightLevels
    }
}

// MARK: - 帧栈动画数据（W9 追加）

/// 单帧：同规格小矩阵 + 帧标签（如 "p = 0.85"）与帧参数值（Slider 绑定用）。
public struct FrameStackFrame: Sendable, Equatable, Identifiable {
    public let id: String
    /// 帧标签（如 "p = 0.85"、"t = 1.2"）
    public let label: String
    /// 帧参数数值（如 Werner p、时间 t）
    public let value: Double
    /// 行主序矩阵（与 spec 的 xTicks/yTicks 对齐）
    public let values: [[Double]]

    public init(id: String? = nil, label: String, value: Double, values: [[Double]]) {
        self.id = id ?? label
        self.label = label
        self.value = value
        self.values = values
    }
}

/// 帧栈动画数据：共享轴刻度的同规格矩阵帧序列。
public struct FrameStackData: Sendable, Equatable {
    public let spec: FrameStackSpec
    public let xTicks: [String]
    public let yTicks: [String]
    /// 帧序列（渲染层按序循环播放）
    public let frames: [FrameStackFrame]

    public init(spec: FrameStackSpec, xTicks: [String], yTicks: [String],
                frames: [FrameStackFrame]) {
        self.spec = spec
        self.xTicks = xTicks
        self.yTicks = yTicks
        self.frames = frames
    }
}

// MARK: - Bloch 球数据（W13 追加）

/// ⑩ Bloch 球（W13：Canvas 伪 3D 线框投影——经纬网格 + 态矢量/对跖点/轨迹；
/// matplotlib 3D 版对应；SceneKit 真 3D 属后续升级，契约不变）。
public struct BlochSpec: Sendable, Equatable {
    public let title: String?

    public init(title: String? = nil) {
        self.title = title
    }
}

/// 三维点（Bloch 矢量 / 轨迹采样）。
public struct Point3D: Sendable, Equatable {
    public let x: Double
    public let y: Double
    public let z: Double

    public init(x: Double, y: Double, z: Double) {
        self.x = x
        self.y = y
        self.z = z
    }
}

/// Bloch 球数据：当前态矢（单位球 |r|=1，混态允许 |r|<1）、对跖点、可选轨迹。
public struct BlochData: Sendable, Equatable {
    public let spec: BlochSpec
    /// 当前态矢量（Bloch 球坐标，r = Tr(ρσ)）
    public let state: Point3D
    /// 态矢标注（如 "|ψ⟩"）
    public let stateLabel: String?
    /// 对跖点（正交补态，可选）
    public let antipode: Point3D?
    /// 对跖点标注（如 "|ψ⊥⟩"）
    public let antipodeLabel: String?
    /// 态矢轨迹（可选，如演化路径）
    public let trajectory: [Point3D]
    /// 北/南极基矢标注（如 "|0⟩"/"|1⟩"）
    public let northLabel: String
    public let southLabel: String

    public init(spec: BlochSpec, state: Point3D, stateLabel: String? = nil,
                antipode: Point3D? = nil, antipodeLabel: String? = nil,
                trajectory: [Point3D] = [],
                northLabel: String = "|0⟩", southLabel: String = "|1⟩") {
        self.spec = spec
        self.state = state
        self.stateLabel = stateLabel
        self.antipode = antipode
        self.antipodeLabel = antipodeLabel
        self.trajectory = trajectory
        self.northLabel = northLabel
        self.southLabel = southLabel
    }
}

/// 图表数据四类封装——与模块 charts 声明数组按序对应。
public enum ChartData: Sendable, Equatable {
    case lineSeries(LineSeriesData)
    case levelDiagram(LevelDiagramData)
    case bars(BarData)
    case scatter(ScatterData)
    case heatmap(HeatmapData)   // W4 追加（契约仅追加）
    case dualAxisLineSeries(DualAxisLineSeriesData)   // W5 追加（契约仅追加）
    case schematic(SchematicData)   // W6 追加（契约仅追加）
    case contour(ContourData)   // W6 追加（契约仅追加）
    case frameStack(FrameStackData)   // W9 追加（契约仅追加）
    case bloch(BlochData)   // W13 追加（契约仅追加）

    public var chartTitle: String? {
        switch self {
        case .lineSeries(let d): return d.spec.title
        case .levelDiagram(let d): return d.spec.title
        case .bars(let d): return d.spec.title
        case .scatter(let d): return d.spec.title
        case .heatmap(let d): return d.spec.title
        case .dualAxisLineSeries(let d): return d.spec.title
        case .schematic(let d): return d.spec.title
        case .contour(let d): return d.spec.title
        case .frameStack(let d): return d.spec.title
        case .bloch(let d): return d.spec.title
        }
    }
}

// MARK: - 摘要数值卡

/// 计算结果摘要（模块页顶部数值卡，如 "λmax = 501 nm"）。
/// value 为已格式化字符串——格式化责任在 compute 侧，UI 零处理。
public struct SummaryItem: Sendable, Equatable, Identifiable {
    public let id: String
    public let title: String
    /// 格式化后的值（含单位）
    public let value: String
    /// 附加说明（如 "维恩位移 b/T"）
    public let note: String?

    public init(id: String? = nil, title: String, value: String,
                note: String? = nil) {
        self.id = id ?? title
        self.title = title
        self.value = value
        self.note = note
    }
}

// MARK: - 理论卡

/// 模块理论卡：关键公式 + 一句话导读。
/// W3 试点起随 SimResult 携带；W7 升级为 MathView LaTeX 渲染（当前等宽文本展示）。
public struct TheoryCard: Sendable, Equatable {
    public let title: String
    /// 公式行（当前为 Unicode 近似排版，W7 换 LaTeX 源串）
    public let formulas: [String]
    /// 一句话导读（结果如何读）
    public let reading: String?

    public init(title: String, formulas: [String], reading: String? = nil) {
        self.title = title
        self.formulas = formulas
        self.reading = reading
    }
}

// MARK: - 计算结果

/// compute 的唯一输出。纯值类型、Sendable，跨 actor 传递安全。
public struct SimResult: Sendable, Equatable {
    /// 与模块 charts 声明按序对应
    public let charts: [ChartData]
    /// 摘要数值卡
    public let summary: [SummaryItem]
    /// 理论卡（公式 + 导读；模块可选提供）
    public let theory: TheoryCard?
    /// 本次 compute 耗时（毫秒）——模块页"拖动即重算 · 上次 x ms"徽章用
    public let computeMillis: Double?

    public init(charts: [ChartData], summary: [SummaryItem] = [],
                theory: TheoryCard? = nil, computeMillis: Double? = nil) {
        self.charts = charts
        self.summary = summary
        self.theory = theory
        self.computeMillis = computeMillis
    }
}

// MARK: - 计算错误

/// 引擎层统一错误（模块页据此展示用户可读信息）。
public enum ComputeError: Error, Sendable, Equatable, LocalizedError {
    /// 参数超出物理有效域（如对数轴负值）
    case invalidParameter(String)
    /// 数值发散 / NaN
    case numericalDivergence(String)
    /// 策略档被取消
    case cancelled

    public var errorDescription: String? {
        switch self {
        case .invalidParameter(let why): return "参数无效：\(why)"
        case .numericalDivergence(let where_): return "计算发散：\(where_)"
        case .cancelled: return "已取消"
        }
    }
}
