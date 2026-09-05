import Testing
import Foundation
@testable import EngineKit

/// W3 包体裁决后的数值线代验收（W1「MLXEigh × dsyevr 互验」平移）：
/// SymEigh 直调 Accelerate dsyevr，以残差 / 正交性 / 解析谱三类不变量验证，
/// 精度基线沿用 W1 的 1×10⁻¹²。
@Suite("SymEigh（Accelerate dsyevr）")
struct SymEighTests {

    // MARK: - 确定性随机对称矩阵（B + Bᵀ 构造，固定种子复现）

    private static func randomSymmetric(n: Int, seed: UInt64) -> [Double] {
        var rng = SplitMix64(seed: seed)
        var b = [Double](repeating: 0, count: n * n)
        for i in 0..<(n * n) { b[i] = rng.unit() * 2 - 1 }
        var a = [Double](repeating: 0, count: n * n)
        for i in 0..<n {
            for j in 0..<n {
                a[i * n + j] = 0.5 * (b[i * n + j] + b[j * n + i])
            }
        }
        return a
    }

    /// SplitMix64：确定性 [0,1) 均匀采样（固定种子复现）
    private struct SplitMix64 {
        var state: UInt64
        init(seed: UInt64) { state = seed }
        mutating func unit() -> Double {
            state &+= 0x9E3779B97F4A7C15
            var z = state
            z = (z ^ (z >> 30)) &* 0xBF58476D1CE4E5B9
            z = (z ^ (z >> 27)) &* 0x94D049BB133111EB
            z = z ^ (z >> 31)
            return Double(z) / Double(UInt64.max)
        }
    }

    // MARK: - 不变量检查

    /// max|A·vₛ − wₛ·vₛ|（逐元素绝对残差）。
    private static func maxResidual(_ aFlat: [Double], n: Int, w: [Double], v: [Double]) -> Double {
        var residual = 0.0
        for s in 0..<n {
            for j in 0..<n {
                var av = 0.0
                for k in 0..<n { av += aFlat[j * n + k] * v[k * n + s] }
                residual = max(residual, abs(av - w[s] * v[j * n + s]))
            }
        }
        return residual
    }

    /// max|(VᵀV)ₛ₁ₛ₂ − δₛ₁ₛ₂|（正交归一偏差）。
    private static func maxOrthError(_ n: Int, _ v: [Double]) -> Double {
        var orth = 0.0
        for s1 in 0..<n {
            for s2 in 0..<n {
                var dot = 0.0
                for j in 0..<n { dot += v[j * n + s1] * v[j * n + s2] }
                orth = max(orth, abs(dot - (s1 == s2 ? 1.0 : 0.0)))
            }
        }
        return orth
    }

    // MARK: - 用例

    @Test("解析谱：diag(3, −1, 2) → w = [−1, 2, 3]，残差逐位")
    func analyticDiagonal() {
        let n = 3
        let a: [Double] = [3, 0, 0, 0, -1, 0, 0, 0, 2]
        let (w, v) = SymEigh.eighSymmetric(a, n: n)
        #expect(abs(w[0] - (-1)) < 1e-12)
        #expect(abs(w[1] - 2) < 1e-12)
        #expect(abs(w[2] - 3) < 1e-12)
        #expect(Self.maxResidual(a, n: n, w: w, v: v) < 1e-12)
    }

    @Test("解析谱：[[2,1],[1,2]] → w = [1, 3]，特征向量为 ±(1,−1)/√2 与 ±(1,1)/√2")
    func analytic2x2() {
        let a: [Double] = [2, 1, 1, 2]
        let (w, v) = SymEigh.eighSymmetric(a, n: 2)
        #expect(abs(w[0] - 1) < 1e-12)
        #expect(abs(w[1] - 3) < 1e-12)
        #expect(Self.maxResidual(a, n: 2, w: w, v: v) < 1e-12)
        // 列 0 ∝ (1, −1)/√2（分量反号），列 1 ∝ (1, 1)/√2（分量同号），符号任意
        let invSqrt2 = 1 / 2.squareRoot()
        #expect(abs(abs(v[0]) - invSqrt2) < 1e-12 && v[0] * v[2] < 0)
        #expect(abs(abs(v[1]) - invSqrt2) < 1e-12 && v[1] * v[3] > 0)
    }

    @Test("随机对称 N=50：残差与正交性 < 1×10⁻¹²，w 升序")
    func randomN50Invariants() {
        let n = 50
        let a = Self.randomSymmetric(n: n, seed: 20260904)
        let (w, v) = SymEigh.eighSymmetric(a, n: n)
        #expect(w.count == n)
        for i in 1..<n { #expect(w[i - 1] <= w[i]) }
        let scale = max(abs(w[0]), abs(w[n - 1]))
        #expect(Self.maxResidual(a, n: n, w: w, v: v) / scale < 1e-12)
        #expect(Self.maxOrthError(n, v) < 1e-12)
    }

    @Test("互验规模用例：N=200 随机对称（秒级档实际规模量级）")
    func randomN200Invariants() {
        let n = 200
        let a = Self.randomSymmetric(n: n, seed: 42609)
        let (w, v) = SymEigh.eighSymmetric(a, n: n)
        let scale = max(abs(w[0]), abs(w[n - 1]))
        #expect(Self.maxResidual(a, n: n, w: w, v: v) / scale < 1e-12)
        #expect(Self.maxOrthError(n, v) < 1e-12)
    }

    @Test("uplo='U' 与 'L' 等价（对称输入读取另一半三角）")
    func uploEquivalence() {
        let n = 30
        let a = Self.randomSymmetric(n: n, seed: 99)
        let lo = SymEigh.eighSymmetric(a, n: n, uplo: "L").w
        let up = SymEigh.eighSymmetric(a, n: n, uplo: "U").w
        for i in 0..<n { #expect(abs(lo[i] - up[i]) < 1e-12) }
    }

    @Test("W3 试点同款：谐振子三对角 H（N=800），前 6 本征值逼近解析 n+½")
    func harmonicOscillatorPreview() {
        // 与 code/量子谐振子_有限差分与相干态.py 同构：H = −½d²/dx² + ½x²，E_n = n+½
        let n = 800
        let L = 18.0
        let dx = 2.0 * L / Double(n - 1)
        var a = [Double](repeating: 0, count: n * n)
        for i in 0..<n {
            let x = -L + Double(i) * dx
            a[i * n + i] = 1.0 / (dx * dx) + 0.5 * x * x
            if i + 1 < n {
                let off = -0.5 / (dx * dx)
                a[i * n + (i + 1)] = off
                a[(i + 1) * n + i] = off
            }
        }
        let start = Date()
        let (w, _) = SymEigh.eighSymmetric(a, n: n)
        let elapsed = Date().timeIntervalSince(start)

        let analytic = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
        for (k, e) in analytic.enumerated() {
            #expect(abs(w[k] - e) / e < 5e-3, "E_\(k): \(w[k]) vs \(e)")
        }
        print("⏱ Accelerate dsyevr N=800 float64: \(String(format: "%.3f", elapsed))s")
    }
}
