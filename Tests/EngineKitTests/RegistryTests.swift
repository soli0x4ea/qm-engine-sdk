import Testing
import Foundation
@testable import EngineKit

/// W2 契约层单测：四类控件规格、ParamValues 默认值/换算、注册表分组/搜索/统计。
@Suite("Registry · SimModule 契约")
struct RegistryTests {

    // 测试用模块：覆盖四类控件 + 全部 4 类图表声明
    private struct DemoModule: SimModule {
        let meta = ModuleMeta(
            id: "test_demo", title: "测试模块 · 全控件",
            subtitle: "契约单测用", category: .oldQuantum, noteNumber: 3,
            tier: .realtime, difficulty: .basic,
            keywords: ["普朗克", "黑体"])

        var params: [ParamSpec] {
            [
                .slider(SliderSpec(key: "temperature", title: "温度", symbol: "T",
                                   unit: "K", range: 100...10000, defaultValue: 5772,
                                   scale: .log, decimalPlaces: 0)),
                .slider(SliderSpec(key: "wavelength", title: "波长", symbol: "λ",
                                   unit: "nm", range: 100...3000, defaultValue: 500)),
                .discrete(DiscreteSpec(key: "metal", title: "金属",
                                       options: [.init(id: "na", title: "钠", subtitle: "W = 2.28 eV"),
                                                 .init(id: "cs", title: "铯", subtitle: "W = 2.14 eV")],
                                       defaultOptionID: "na")),
                .multiCompare(MultiCompareSpec(key: "sources", title: "光源对比",
                                               candidates: [.init(id: "sun", label: "太阳", value: 5772),
                                                            .init(id: "bulb", label: "白炽灯", value: 2800),
                                                            .init(id: "cmb", label: "CMB", value: 2.725)],
                                               defaultSelectionIDs: ["sun", "bulb"])),
                .constant(ConstantSpec(key: "h", title: "h", note: "普朗克常量")),
            ]
        }

        var charts: [ChartSpec] {
            [
                .lineSeries(LineSeriesSpec(xAxis: .init(label: "λ (nm)"),
                                           yAxis: .init(label: "B_λ"),
                                           seriesNames: ["Planck", "Wien"])),
                .levelDiagram(LevelDiagramSpec(energyAxis: "E (E_ℏω)")),
                .bars(BarSpec(xAxis: .init(label: "n"), yAxis: .init(label: "P(n)"))),
                .scatter(ScatterSpec(xAxis: .init(label: "x"), yAxis: .init(label: "y"),
                                     seriesNames: ["a", "b"])),
            ]
        }

        func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
            SimResult(charts: [
                .lineSeries(LineSeriesData(spec: .init(xAxis: .init(label: "λ (nm)"),
                                                       yAxis: .init(label: "B_λ"),
                                                       seriesNames: ["Planck"]),
                                           series: [.init(name: "Planck",
                                                          points: [.init(x: 0, y: 0), .init(x: 1, y: 1)])])),
                .levelDiagram(LevelDiagramData(spec: .init(energyAxis: "E"),
                                               levels: [.init(label: "n=0", energy: 0.5)],
                                               transitions: [.init(fromIndex: 1, toIndex: 0, label: "501 nm")])),
                .bars(BarData(spec: .init(xAxis: .init(label: "n"), yAxis: .init(label: "P")),
                              bars: [.init(label: "n=0", value: 0.36)])),
                .scatter(ScatterData(spec: .init(xAxis: .init(label: "x"), yAxis: .init(label: "y"),
                                                 seriesNames: ["a"]),
                                     series: [.init(name: "a", points: [.init(x: 0, y: 0)])])),
            ], summary: [SummaryItem(title: "λmax", value: "501 nm", note: "维恩位移")],
               computeMillis: 0.3)
        }
    }

    @Test("ParamValues.defaults 覆盖四类控件默认值")
    func paramDefaults() {
        let module = DemoModule()
        let v = ParamValues.defaults(for: module.params)
        #expect(v.slider("temperature") == 5772)
        #expect(v.slider("wavelength") == 500)
        #expect(v.discrete("metal") == "na")
        #expect(v.multiCompare("sources") == Set(["sun", "bulb"]))
    }

    @Test("sliderComputeValue：computeUnit 为 nil 时原值，非空时经 UnitKit 换算")
    func computeUnitConversion() throws {
        // nm → μm：500 nm = 0.5 μm
        let spec = SliderSpec(key: "wl", title: "波长", unit: "nm",
                              range: 100...3000, defaultValue: 500, computeUnit: "μm")
        let v = ParamValues(sliders: ["wl": 500])
        let set = try ConstantsSet.load(.v2022)
        let converted = try v.sliderComputeValue(spec, constants: set)
        #expect(abs(converted - 0.5) < 1e-12)

        // nil computeUnit：原值
        let plain = SliderSpec(key: "t", title: "温度", unit: "K",
                               range: 1...10000, defaultValue: 300)
        let vPlain = ParamValues(sliders: ["t": 300])
        #expect(try vPlain.sliderComputeValue(plain, constants: set) == 300)
    }

    @Test("注册表：注册 / id 去重 / 分组保持分类声明序 / counts 统计")
    func registryBasics() {
        let reg = ModuleRegistry()
        reg.register(DemoModule())
        #expect(reg.all.count == 1)
        #expect(reg.module(id: "test_demo") != nil)
        #expect(reg.module(id: "nope") == nil)

        // 重复 id 忽略并返回 false
        #expect(reg.register(DemoModule()) == false)
        #expect(reg.all.count == 1)

        // counts
        let c = reg.counts
        #expect(c.total == 1 && c.realtime == 1 && c.basic == 1 && c.advanced == 0)

        // 分组：单模块 → 单分类组
        let grouped = reg.grouped
        #expect(grouped.count == 1)
        #expect(grouped[0].category == .oldQuantum)
        #expect(grouped[0].modules.count == 1)
    }

    @Test("搜索：标题 / 关键词 / 笔记编号 / 分类名命中，空查询全通过")
    func searchBehavior() {
        let reg = ModuleRegistry()
        reg.register(DemoModule())
        #expect(reg.search("黑体").count == 1)          // 关键词
        #expect(reg.search("测试模块").count == 1)       // 标题
        #expect(reg.search("03").count == 1)            // 笔记编号补零
        #expect(reg.search("3").count == 1)             // 笔记编号整数
        #expect(reg.search("薛定谔").count == 0)         // 不命中
        #expect(reg.search("").count == 1)              // 空查询全通过
        #expect(reg.search("旧量子论").count == 1)       // 分类名
    }

    @Test("compute 端到端：默认参数 → SimResult 四类图表数据 + 摘要")
    func computeEndToEnd() async throws {
        let module = DemoModule()
        let v = ParamValues.defaults(for: module.params)
        let set = try ConstantsSet.load(.v2022)
        let result = try await module.compute(v, constants: set)

        #expect(result.charts.count == 4)
        #expect(result.summary.first?.value == "501 nm")
        #expect(result.computeMillis != nil)

        guard case .lineSeries(let line) = result.charts[0] else {
            Issue.record("charts[0] 应为 lineSeries"); return
        }
        #expect(line.series.first?.points.count == 2)

        guard case .levelDiagram(let lv) = result.charts[1] else {
            Issue.record("charts[1] 应为 levelDiagram"); return
        }
        #expect(lv.levels.first?.energy == 0.5)
        #expect(lv.transitions.count == 1)
    }

    @Test("SliderSpec 对数刻度契约由 precondition 保证（此处只验证合法构造）")
    func sliderValidation() {
        let ok = SliderSpec(key: "wl", title: "λ", symbol: "λ", unit: "nm",
                            range: 100...3000, defaultValue: 500, scale: .log)
        #expect(ok.range.contains(ok.defaultValue))
    }

    // MARK: - W5 契约追加：双 Y 轴 / 区域填充 / 参考线样式（只追加不破坏）

    @Test("W5 追加：DualAxisLineSeries 声明与数据一一对应，Equatable 稳定")
    func dualAxisContract() {
        let spec = DualAxisLineSeriesSpec(
            title: "IQHE", xAxis: .init(label: "B (T)"),
            primaryAxis: .init(label: "ρ_xy/R_K"), secondaryAxis: .init(label: "ρ_xx (a.u.)"),
            primaryNames: ["ρ_xy"], secondaryNames: ["ρ_xx"])
        let data = DualAxisLineSeriesData(
            spec: spec,
            primary: [.init(name: "ρ_xy", points: [.init(x: 0, y: 0), .init(x: 1, y: 1)])],
            secondary: [.init(name: "ρ_xx", points: [.init(x: 0, y: 1), .init(x: 1, y: 0)])],
            referenceLines: [.init(label: "平台", axis: .y, value: 0.25)])

        let chart = ChartData.dualAxisLineSeries(data)
        #expect(chart.chartTitle == "IQHE")
        #expect(ChartData.dualAxisLineSeries(data) == ChartData.dualAxisLineSeries(data))
        #expect(ChartData.dualAxisLineSeries(data) != ChartData.lineSeries(
            LineSeriesData(spec: .init(xAxis: .init(label: "x"), yAxis: .init(label: "y"),
                                       seriesNames: []), series: [])))
        guard case .dualAxisLineSeries(let d) = chart else {
            Issue.record("应为 dualAxisLineSeries"); return
        }
        #expect(d.primary.count == 1 && d.secondary.count == 1)
        #expect(d.primary[0].points.count == 2)
        #expect(d.referenceLines.first?.value == 0.25)
        // ChartSpec 声明侧
        #expect(ChartSpec.dualAxisLineSeries(spec) == ChartSpec.dualAxisLineSeries(spec))
    }

    @Test("W5 追加：AreaFill 区域填充 + LineSeriesData 默认值向后兼容")
    func areaFillContract() {
        let legacy = LineSeriesData(
            spec: .init(xAxis: .init(label: "Q (%)"), yAxis: .init(label: "R"),
                       seriesNames: ["R(Q)"]),
            series: [.init(name: "R(Q)", points: [.init(x: 0, y: 1)])])
        #expect(legacy.areaFills.isEmpty, "既有构造零改动——areaFills 默认空")

        let fill = AreaFill(name: "安全区",
                            points: [.init(x: 0, y: 0), .init(x: 11, y: 0.001)],
                            baseline: 0, colorIndex: 0)
        let withFill = LineSeriesData(
            spec: legacy.spec, series: legacy.series, areaFills: [fill])
        #expect(withFill.areaFills.count == 1)
        #expect(withFill.areaFills[0].baseline == 0)
        #expect(withFill != legacy, "填充是值的一部分")
    }

    @Test("W5 追加：ReferenceLine 默认 threshold，subtle 档可区分")
    func referenceLineStyle() {
        let plain = ReferenceLine(label: "λmax = 501 nm", axis: .x, value: 501)
        #expect(plain.style == .threshold, "既有调用不传 style——默认 threshold")
        let zero = ReferenceLine(label: "R = 0", axis: .y, value: 0, style: .subtle)
        #expect(zero.style == .subtle)
        #expect(plain != zero)
    }
}
