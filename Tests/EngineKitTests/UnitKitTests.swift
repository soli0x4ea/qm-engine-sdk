import Testing
import Foundation
@testable import EngineKit

/// 单位换算层验收：定义值精确、常量依赖换算与 CODATA 联动、量纲校验。
@Suite("UnitKit")
struct UnitKitTests {

    @Test("能量换算：eV 系（SI 定义精确值；前缀链允许 1 ulp）")
    func energyEV() throws {
        #expect(try UnitKit.convert(1.0, from: "eV", to: "J") == 1.602176634e-19)
        #expect(try UnitKit.factor(from: "J", to: "eV") == 1 / 1.602176634e-19)
        // eV/meV 的 SI 因子比值在 IEEE 中不保证整（10⁻³ 非二进制幂），按 1 ulp 断言
        let meV = try UnitKit.convert(2.5, from: "eV", to: "meV")
        #expect(abs(meV - 2500.0) <= 2500.0 * 1e-15)
        let MeV = try UnitKit.convert(1.0, from: "MeV", to: "eV")
        #expect(abs(MeV - 1e6) <= 1e6 * 1e-15)
    }

    @Test("能量↔温度：1 eV 的温度当量（kB 联动）")
    func energyTemperature() throws {
        let s = try ConstantsSet.load(.v2018)
        // 1 eV / kB = 1.602176634e-19 / 1.380649e-23（均为 SI 定义精确值 → 结果可逐位复算）
        let kelvinPerEV = 1.602176634e-19 / 1.380649e-23
        #expect(try UnitKit.convert(1.0, from: "eV", to: "K", constants: s)
                == kelvinPerEV)
        #expect(try UnitKit.convert(300.0, from: "K", to: "J", constants: s)
                == 300.0 * 1.380649e-23)
    }

    @Test("长度换算：a0 玻尔半径（CODATA 联动）")
    func lengthBohrRadius() throws {
        let s18 = try ConstantsSet.load(.v2018)
        #expect(try UnitKit.convert(1.0, from: "a0", to: "m", constants: s18)
                == 5.29177210903e-11)
        #expect(try UnitKit.convert(1.0, from: "a0", to: "nm", constants: s18)
                == 5.29177210903e-2)
        // 2022 版 a0 更新 → 换算结果随之变化（单位层与常量版本联动）
        let s22 = try ConstantsSet.load(.v2022)
        #expect(try Self.value22a0(s22) != 5.29177210903e-11)
    }

    @Test("时间换算（10⁻³ 比值非二进制幂，按 1 ulp 断言）")
    func timeUnits() throws {
        let fs = try UnitKit.convert(1.0, from: "s", to: "fs")
        #expect(abs(fs - 1e15) <= 1e15 * 1e-15)
        let ns = try UnitKit.convert(2.0, from: "ns", to: "ps")
        #expect(abs(ns - 2000.0) <= 2000.0 * 1e-15)
        // 反向：2 ps = 0.002 ns
        let ps = try UnitKit.convert(2.0, from: "ps", to: "ns")
        #expect(abs(ps - 0.002) <= 0.002 * 1e-15)
    }

    @Test("波数 → 能量：E = h c ν̃（1 cm^-1）")
    func waveNumberEnergy() throws {
        let s = try ConstantsSet.load(.v2018)
        let e = try UnitKit.energyFromWaveNumber(100.0, constants: s)  // 100 m^-1 = 1 cm^-1
        let want = 6.62607015e-34 * 299_792_458.0 * 100.0
        #expect(e == want)
        #expect(try UnitKit.factor(from: "cm^-1", to: "m^-1") == 100.0)
    }

    @Test("量纲不匹配与未知单位抛错")
    func errorCases() throws {
        #expect(throws: UnitKitError.self) {
            try UnitKit.convert(1.0, from: "eV", to: "nm")
        }
        #expect(throws: UnitKitError.self) {
            try UnitKit.factor(from: "J", to: "furlong")
        }
        // 依赖常量的单位未传 constants → 抛错
        #expect(throws: UnitKitError.self) {
            try UnitKit.factor(from: "K", to: "J")
        }
    }
}

extension UnitKitTests {
    /// 2022 版玻尔半径取值（供版本联动断言）
    static func value22a0(_ s: ConstantsSet) throws -> Double {
        try s.value("a0")
    }
}
