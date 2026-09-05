import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// 收尾执行包 1b：全量 63 模块性能审计（log-only，常设账本）。
///
/// 复用既有计时口径 `expectComputeUnderBudget`（best-of-N 墙钟：warmup 后取 N 次
/// 采样的最小值为判据，中位并打印），不发明新 harness。逐模块按档位名义预算打印：
/// `⏱ <id> compute 最小 X ms（best-of-N，预算 Y ms）· 中位 Z ms`。
///
/// runs/warmup 沿用各档既有预算测试惯例：
/// - 实时档 41/5（默认口径，PhotonStatistics 等绝大多数实时模块）
/// - 秒级档 5/2（DiracSpectrum / KleinGordonDensity / 谐振子先例）
/// - 策略档 3/1（KramersKronig / Entanglement / PathIntegral / 两带光学先例）；
///   策略档额外在清缓存后测一次**冷启动首算**（`❄️` 行），与热路径（缓存命中）分开披露。
///
/// **环境变量门禁**：两个策略档模块（路径积分、两带光学）持有静态
/// `StrategyCache`/`StrategySession`，其命中率断言测试要求冷启动语义；本审计若与
/// 全量套件并行运行会污染共享计数。故默认禁用，账本采集用
/// `QM_PERF_AUDIT=1 swift test --filter PerformanceAuditTests` 显式串行运行。
///
/// `enforce: false`：审计只测量打印不断言，账本行进周报告；
/// 红绿判据仍由各模块自身的预算测试（enforce 强制口径）承担。
@Suite(.serialized, .disabled(if: ProcessInfo.processInfo.environment["QM_PERF_AUDIT"] != "1",
                             "与策略档模块的共享 StrategyCore 状态断言冲突；用 QM_PERF_AUDIT=1 显式运行"))
struct PerformanceAuditTests {

    @Test("全量 63 模块性能审计：逐模块 best/中位 vs 名义预算（log-only）")
    func fullLedgerAudit() async throws {
        let registry = ModuleRegistry.shared
        registry.removeAll()
        #expect(ModuleLibrary.registerBuiltins() == 63)
        let constants = try ConstantsSet.load(.v2022)
        let budgetByTier: [ComputeTier: Double] = [
            .realtime: 16, .seconds: 2000, .strategy: 10_000]
        let samplingByTier: [ComputeTier: (runs: Int, warmup: Int)] = [
            .realtime: (41, 5), .seconds: (5, 2), .strategy: (3, 1)]
        var audited = 0
        for module in registry.all {
            let s = samplingByTier[module.meta.tier]!
            // 策略档冷启动首算：清模块缓存后单次计时（热路径由下方 best-of-N 覆盖）
            if module.meta.tier == .strategy {
                if module is TwoBandOpticalModule { TwoBandOpticalModule.jdosCache.removeAll() }
                if module is PathIntegralModule { PathIntegralModule.eigenCache.removeAll() }
                let t0 = ContinuousClock.now
                _ = try await module.compute(
                    ParamValues.defaults(for: module.params), constants: constants)
                let d = ContinuousClock.now - t0
                let c = d.components
                let ms = Double(c.seconds) * 1000 + Double(c.attoseconds) / 1e15
                print("❄️ \(module.meta.id) 冷启动首算 \(String(format: "%.2f", ms)) ms"
                      + "（清缓存单次，名义预算 \(budgetByTier[.strategy]!) ms）")
            }
            try await expectComputeUnderBudget(
                module: module,
                values: ParamValues.defaults(for: module.params),
                constants: constants,
                budgetMillis: budgetByTier[module.meta.tier]!,
                runs: s.runs, warmup: s.warmup, enforce: false)
            audited += 1
        }
        #expect(audited == 63, "账本须覆盖全部 63 模块，实际 \(audited)")
    }
}
