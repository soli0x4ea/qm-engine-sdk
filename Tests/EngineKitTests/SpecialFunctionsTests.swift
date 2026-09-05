import Testing
import Foundation
@testable import EngineKit

/// W2 特殊函数单测：连带拉盖尔 / 氢原子径向 / 勒让德，逐点对照 numpy 参考值。
/// 参考值构造：laguerre_reference.json、legendre_reference.json（numpy 求导构造）。
@Suite("Laguerre / Legendre 特殊函数")
struct SpecialFunctionsTests {

    private static func loadJSON(_ name: String) throws -> [String: Any] {
        let base = (name as NSString).deletingPathExtension
        let ext = (name as NSString).pathExtension.isEmpty ? "json" : (name as NSString).pathExtension
        let url = try #require(
            Bundle.module.url(forResource: base, withExtension: ext, subdirectory: "Fixtures")
        )
        let data = try Data(contentsOf: url)
        return try #require(try JSONSerialization.jsonObject(with: data) as? [String: Any])
    }

    // MARK: - Laguerre

    @Test("连带拉盖尔 L^α_n：α∈{1,3,5} × n=0..6 逐点 vs numpy（rel < 1e-12）")
    func assocLaguerreAgainstNumpy() throws {
        let ref = try Self.loadJSON("laguerre_reference.json")
        let xs = (ref["x"] as! [NSNumber]).map(\.doubleValue)
        let L = ref["L"] as! [String: [NSNumber]]
        let alphas = ref["alphas"] as! [Int]
        let ns = ref["n"] as! [Int]

        var checked = 0
        for a in alphas {
            for n in ns {
                let expected = (L["a\(a)_n\(n)"] ?? []).map(\.doubleValue)
                guard !expected.isEmpty else { continue }
                for (i, x) in xs.enumerated() {
                    let got = Laguerre.assocLaguerre(n, alpha: Double(a), x)
                    let exp = expected[i]
                    let denom = max(1.0, abs(exp))
                    #expect(abs(got - exp) / denom < 1e-12,
                            "L^\(a)_\(n)(\(x)): \(got) vs \(exp)")
                }
                checked += 1
            }
        }
        #expect(checked == 21)
    }

    @Test("氢原子径向 R_n,l：n=1..3 全部 l 逐点 vs numpy（rel < 1e-12）")
    func hydrogenRadialAgainstNumpy() throws {
        let ref = try Self.loadJSON("laguerre_reference.json")
        let rs = (ref["r"] as! [NSNumber]).map(\.doubleValue)
        let R = ref["R"] as! [String: [NSNumber]]

        for n in 1...3 {
            for l in 0..<n {
                let expected = (R["n\(n)l\(l)"] ?? []).map(\.doubleValue)
                guard !expected.isEmpty else { continue }
                let got = Laguerre.hydrogenRadial(n, l, rs: rs)
                for i in rs.indices {
                    let denom = max(1e-300, abs(expected[i]))
                    #expect(abs(got[i] - expected[i]) / denom < 1e-12
                            || abs(got[i] - expected[i]) < 1e-300,
                            "R_\(n)\(l)(r=\(rs[i])): \(got[i]) vs \(expected[i])")
                }
            }
        }
    }

    @Test("氢原子径向归一化：∫R²r²dr ≈ 1（n≤2，r≤20a₀ 域内）")
    func hydrogenNormalization() throws {
        let rs = Array(stride(from: 0.0, through: 30.0, by: 0.01))
        // n≤2：30a₀ 足够覆盖尾部；n=3 需 60a₀（fixtures 域截断故只验 n≤2）
        for (n, l) in [(1, 0), (2, 0), (2, 1)] {
            let R = Laguerre.hydrogenRadial(n, l, rs: rs)
            var integral = 0.0
            for i in 1..<rs.count {
                let f0 = R[i-1]*R[i-1] * rs[i-1]*rs[i-1]
                let f1 = R[i]*R[i] * rs[i]*rs[i]
                integral += 0.5 * (f0 + f1) * (rs[i] - rs[i-1])
            }
            #expect(abs(integral - 1) < 1e-6, "n\(n)l\(l) 归一化积分 = \(integral)")
        }
    }

    @Test("拉盖尔已知值抽查：L^1_0=1，L^1_1=2−x，L^0_2=1−2x+x²/2")
    func laguerreSpotChecks() {
        #expect(Laguerre.assocLaguerre(0, alpha: 1, 3.7) == 1)
        #expect(abs(Laguerre.assocLaguerre(1, alpha: 1, 2.5) - (2 - 2.5)) < 1e-15)
        #expect(abs(Laguerre.assocLaguerre(2, alpha: 0, 2.0) - (1 - 4 + 2)) < 1e-15)
    }

    // MARK: - Legendre

    @Test("勒让德 P_l：l=0..8 逐点 vs numpy（rel < 1e-12）")
    func legendreAgainstNumpy() throws {
        let ref = try Self.loadJSON("legendre_reference.json")
        let xs = (ref["x"] as! [NSNumber]).map(\.doubleValue)
        let P = ref["P"] as! [String: [NSNumber]]

        for l in 0...8 {
            let expected = (P["l\(l)"] ?? []).map(\.doubleValue)
            guard !expected.isEmpty else { continue }
            for (i, x) in xs.enumerated() {
                let got = Legendre.value(l, x)
                let denom = max(1.0, abs(expected[i]))
                #expect(abs(got - expected[i]) / denom < 1e-12,
                        "P_\(l)(\(x)): \(got) vs \(expected[i])")
            }
        }
    }

    @Test("连带勒让德 P^m_l：9 组代表组合 vs numpy 求导构造（rel < 1e-12）")
    func assocLegendreAgainstNumpy() throws {
        let ref = try Self.loadJSON("legendre_reference.json")
        let xs = (ref["x"] as! [NSNumber]).map(\.doubleValue)
        let Pm = ref["Pm"] as! [String: [NSNumber]]
        let combos = ref["combos"] as! [String]

        for key in combos {
            let expected = (Pm[key] ?? []).map(\.doubleValue)
            guard !expected.isEmpty else { continue }
            // key 形如 "l3m2"（兼容负阶 "l2m-2"）
            let parts = key.split(separator: "m")
            let l = Int(parts[0].dropFirst())!
            let m = Int(parts[1])!
            for (i, x) in xs.enumerated() {
                let got = Legendre.assoc(l, m, x)
                let denom = max(1e-3, abs(expected[i]))
                #expect(abs(got - expected[i]) / denom < 1e-12,
                        "\(key)(\(x)): \(got) vs \(expected[i])")
            }
        }
    }

    @Test("连带勒让德解析抽查：P^1_1=−√(1−x²)，P^2_2=3(1−x²)，负阶关系")
    func assocLegendreSpotChecks() {
        #expect(abs(Legendre.assoc(1, 1, 0.6) + sqrt(1 - 0.36)) < 1e-15)
        #expect(abs(Legendre.assoc(2, 2, 0.5) - 3 * (1 - 0.25)) < 1e-14)
        #expect(abs(Legendre.assoc(2, 1, 0.0)) < 1e-15)  // −1.5x → x=0 为 0
        // 负阶：P^{−1}_1 = (−1)¹·0!/2!·P^1_1 = −0.5·P^1_1
        let pm1 = Legendre.assoc(1, 1, 0.6)
        let pmm1 = Legendre.assoc(1, -1, 0.6)
        #expect(abs(pmm1 + 0.5 * pm1) < 1e-15)
    }
}
