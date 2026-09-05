import Foundation
import EngineKit

// MARK: - 计算核（逐句移植 code/黑体辐射_三律对比.py）

/// 黑体辐射三律的纯函数核。单位与 Python 一致：
/// u_λ 为谱能量密度（J/m⁴），图表按 ×1e9 转为 J·m⁻⁴·nm⁻¹。
enum BlackbodyMath {

    /// 普朗克谱能量密度（频率表示），J·s/m³。λ/ν > 0 假定由调用方保证。
    static func uNu(_ nu: Double, _ T: Double, h: Double, kB: Double, c: Double) -> Double {
        let x = h * nu / (kB * T)
        return 8 * .pi * h * pow(nu, 3) / pow(c, 3) / (exp(x) - 1)
    }

    /// 普朗克谱能量密度（波长表示），J/m⁴；λ ≤ 0 返回 0（与 Python u_lambda 一致）。
    static func uLambda(_ lam: Double, _ T: Double, h: Double, kB: Double, c: Double) -> Double {
        guard lam > 0 else { return 0 }
        let x = h * c / (lam * kB * T)
        return 8 * .pi * h * c / pow(lam, 5) / (exp(x) - 1)
    }

    /// 维恩近似（高频极限：略去玻色因子 −1），J/m⁴。
    static func uWien(_ lam: Double, _ T: Double, h: Double, kB: Double, c: Double) -> Double {
        let x = h * c / (lam * kB * T)
        return 8 * .pi * h * c / pow(lam, 5) * exp(-x)
    }

    /// 瑞利-金斯（低频极限），J/m⁴——短波端 ∝ λ⁻⁴ 发散（紫外灾难）。
    static func uRayleighJeans(_ lam: Double, _ T: Double, kB: Double, c: Double) -> Double {
        8 * .pi * kB * T / pow(lam, 4)
    }

    /// 维恩位移常数 b = hc/(x_peak·kB)；x_peak = 4.965114231 为 x·eˣ/(eˣ−1) = 5 的根。
    static func wienB(h: Double, kB: Double, c: Double) -> Double {
        h * c / (4.965114231 * kB)
    }

    /// 斯特藩-玻尔兹曼常数（由 SI 精确常量导出：σ = 2π⁵kB⁴/15h³c²）。
    static func sigmaDerived(h: Double, kB: Double, c: Double) -> Double {
        2 * pow(.pi, 5) * pow(kB, 4) / (15 * pow(h, 3) * c * c)
    }
}

// MARK: - 模块

/// 笔记 03《黑体辐射与能量量子化》：普朗克谱 vs 维恩近似 vs 瑞利-金斯。
/// 实时档：4000 点网格三条曲线，拖动温度即重算。
struct BlackbodyModule: SimModule {

    let meta = ModuleMeta(
        id: "黑体辐射_三律对比", title: "黑体辐射 · 三律对比",
        subtitle: "普朗克 / 维恩近似 / 瑞利-金斯谱与紫外灾难",
        category: .oldQuantum, noteNumber: 3, tier: .realtime, difficulty: .basic,
        keywords: ["普朗克", "Planck", "维恩", "Wien", "瑞利", "金斯", "Rayleigh",
                   "Jeans", "紫外灾难", "黑体", "斯特藩", "玻尔兹曼", "位移定律"])

    var params: [ParamSpec] {
        [
            .slider(SliderSpec(key: "T", title: "温度", symbol: "T", unit: "K",
                               range: 100...15000, defaultValue: 5772,
                               scale: .log, decimalPlaces: 0)),
            .multiCompare(MultiCompareSpec(
                key: "laws", title: "对比定律",
                candidates: [.init(id: "planck", label: "普朗克", value: 0),
                             .init(id: "wien", label: "维恩近似", value: 1),
                             .init(id: "rj", label: "瑞利-金斯", value: 2)],
                defaultSelectionIDs: ["planck", "wien", "rj"], maxSelection: 3)),
            .constant(ConstantSpec(key: "h", title: "h", note: "普朗克常量")),
            .constant(ConstantSpec(key: "kB", title: "kB", note: "玻尔兹曼常量")),
            .constant(ConstantSpec(key: "c", title: "c", note: "真空光速")),
        ]
    }

    var charts: [ChartSpec] {
        [.lineSeries(LineSeriesSpec(
            title: "黑体谱 u_λ(λ, T)：三律对比",
            xAxis: .init(label: "λ (nm)"),
            yAxis: .init(label: "u_λ (J·m⁻⁴·nm⁻¹)"),
            seriesNames: ["普朗克", "维恩近似", "瑞利-金斯"]))]
    }

    func compute(_ input: ParamValues, constants: ConstantsSet) async throws -> SimResult {
        let T = input.slider("T")
        let laws = input.multiCompare("laws")
        let h = try constants.value("h")
        let kB = try constants.value("kB")
        let c = try constants.value("c")

        // 与 Python lam_plot = linspace(5e-8, 3e-6, 4000) 一致（nm 表示）
        let grid = Num.linspace(50, 3000, count: 4000)
        let planck = grid.map { BlackbodyMath.uLambda($0 * 1e-9, T, h: h, kB: kB, c: c) * 1e9 }
        let planckPeak = planck.max() ?? 0

        var series: [SeriesPoints] = []
        if laws.contains("planck") {
            series.append(.init(name: "普朗克", points: Num.strided(grid, planck, stride: 8)))
        }
        if laws.contains("wien") {
            let wien = grid.map { BlackbodyMath.uWien($0 * 1e-9, T, h: h, kB: kB, c: c) * 1e9 }
            series.append(.init(name: "维恩近似", points: Num.strided(grid, wien, stride: 8)))
        }
        if laws.contains("rj") {
            // 与 Python 同规则：仅绘制 RJ < 4×普朗克峰值 的区段（截断紫外发散）
            var xs: [Double] = [], ys: [Double] = []
            for x in grid {
                let rj = BlackbodyMath.uRayleighJeans(x * 1e-9, T, kB: kB, c: c) * 1e9
                if rj < 4 * planckPeak { xs.append(x); ys.append(rj) }
            }
            if !xs.isEmpty {
                let step = max(1, (xs.count + 399) / 400)
                series.append(.init(name: "瑞利-金斯", points: Num.strided(xs, ys, stride: step)))
            }
        }

        let b = BlackbodyMath.wienB(h: h, kB: kB, c: c)
        let lamMaxNm = b / T * 1e9
        let sigma = BlackbodyMath.sigmaDerived(h: h, kB: kB, c: c)
        let jStar = sigma * pow(T, 4)

        let lineSpec = LineSeriesSpec(
            xAxis: .init(label: "λ (nm)"),
            yAxis: .init(label: "u_λ (J·m⁻⁴·nm⁻¹)"),
            seriesNames: ["普朗克", "维恩近似", "瑞利-金斯"])

        return SimResult(
            charts: [.lineSeries(LineSeriesData(
                spec: lineSpec,
                series: series,
                referenceLines: [ReferenceLine(
                    id: "lmax", label: String(format: "λmax = %.0f nm", lamMaxNm),
                    axis: .x, value: lamMaxNm)]))],
            summary: [
                .init(id: "lmax", title: "λmax", value: String(format: "%.1f nm", lamMaxNm),
                      note: "维恩位移 b/T"),
                .init(id: "jstar", title: "总辐出度 j*",
                      value: String(format: "%.4e W/m²", jStar),
                      note: "斯特藩-玻尔兹曼 σT⁴"),
                .init(id: "sigma", title: "σ（导出）",
                      value: String(format: "%.4e W/m²/K⁴", sigma),
                      note: "2π⁵kB⁴/15h³c²"),
            ],
            theory: TheoryCard(
                title: "黑体辐射三律",
                formulas: [
                    "普朗克：u_λ = (8πhc/λ⁵) · 1/(e^(hc/λk_BT) − 1)",
                    "维恩近似：eˣ ≫ 1 时略去 −1（高频端偏差）",
                    "瑞利-金斯：eˣ ≈ 1+x 低频展开 → u ∝ T/λ⁴（短波发散）",
                    "斯特藩-玻尔兹曼 j* = σT⁴；维恩位移 λmax·T = b",
                ],
                reading: "升温观察：峰值左移（维恩位移）、面积按 T⁴ 膨胀；短波端瑞利-金斯发散即紫外灾难，普朗克因子将其压回。"))
    }
}
