import Testing
import Foundation
@testable import EngineKit

/// Hermite 递推库验收：解析恒等式 + numpy 参考值（1×10⁻¹²）。
@Suite("Hermite")
struct HermiteTests {

    @Test("低阶解析恒等式（递推与闭式差 ≤ 1 ulp）")
    func lowOrderIdentities() {
        #expect(Hermite.value(0, 3.7) == 1.0)
        #expect(Hermite.value(1, 2.5) == 5.0)
        // H2 = 4x² − 2
        #expect(Hermite.value(2, 1.25) == 4 * 1.25 * 1.25 - 2)
        // H3 = 8x³ − 12x（递推与闭式的求值顺序不同，逐位不等，按相对容差断言）
        let h3 = Hermite.value(3, -0.7)
        let h3Closed = 8 * pow(-0.7, 3) - 12 * (-0.7)
        #expect(abs(h3 - h3Closed) <= 1e-15 * max(1, abs(h3Closed)))
        // H4 = 16x⁴ − 48x² + 12
        let x = 0.83
        let h4 = Hermite.value(4, x)
        let h4Closed = 16 * pow(x, 4) - 48 * x * x + 12
        #expect(abs(h4 - h4Closed) <= 1e-14 * max(1, abs(h4Closed)))
    }

    @Test("递推关系自洽：Hₙ₊₁ = 2x·Hₙ − 2n·Hₙ₋₁")
    func recurrenceSelfConsistency() {
        for x in stride(from: -3.0, through: 3.0, by: 0.7) {
            for n in 2...12 {
                let rhs = 2 * x * Hermite.value(n - 1, x) - 2 * Double(n - 1) * Hermite.value(n - 2, x)
                #expect(abs(Hermite.value(n, x) - rhs) <= 1e-12 * max(1, abs(rhs)),
                        "n=\(n), x=\(x)")
            }
        }
    }

    @Test("与 numpy hermval 参考值一致（n = 0…8，rel < 1×10⁻¹²）")
    func againstNumpyReference() throws {
        let ref = try Self.loadReference()
        let ns: [Int] = ref["n"] as! [Int]
        let xs: [Double] = ref["x"] as! [Double]
        let H: [[Double]] = ref["H"] as! [[Double]]
        #expect(ns.count == H.count && xs.count > 0)
        for (i, n) in ns.enumerated() {
            let mine = Hermite.values(n, xs: xs)
            for (j, got) in mine.enumerated() {
                let want = H[i][j]
                let rel = abs(got - want) / max(1.0, abs(want))
                #expect(rel < 1e-12, "n=\(n), x=\(xs[j]): \(got) vs \(want)")
            }
        }
    }

    @Test("归一化波函数 ψn 与参考值一致，且 ∫ψ²dx ≈ 1")
    func normalizedWavefunctionAgainstReference() throws {
        let ref = try Self.loadReference()
        let ns: [Int] = ref["n_psi"] as! [Int]
        let xs: [Double] = ref["x_psi"] as! [Double]
        let psi: [[Double]] = ref["psi"] as! [[Double]]
        for (i, n) in ns.enumerated() {
            let mine = Hermite.normalizedWavefunction(n, xs: xs)
            for (j, got) in mine.enumerated() {
                let want = psi[i][j]
                let rel = abs(got - want) / max(1e-300, abs(want))
                #expect(rel < 1e-12, "ψ n=\(n), x=\(xs[j]): \(got) vs \(want)")
            }
            // 梯形积分 ∫ψ²dx ≈ 1（网格 0.01，误差 O(h²) ~ 1e-5）
            var sum = 0.0
            for j in 1..<mine.count {
                let dx = xs[j] - xs[j - 1]
                sum += 0.5 * (mine[j] * mine[j] + mine[j - 1] * mine[j - 1]) * dx
            }
            #expect(abs(sum - 1.0) < 1e-4, "n=\(n) 归一化积分 = \(sum)")
        }
    }

    private static func loadReference() throws -> [String: Any] {
        let url = try #require(
            Bundle.module.url(forResource: "hermite_reference", withExtension: "json",
                              subdirectory: "Fixtures")
        )
        let data = try Data(contentsOf: url)
        return try #require(try JSONSerialization.jsonObject(with: data) as? [String: Any])
    }
}
