import Foundation

/// 物理学家 Hermite 多项式 Hₙ(x)——谐振子模块（笔记 15）与相干态展开的主力依赖。
///
/// 实现路线与 `numpy.polynomial.hermite` 一致：三项递推
/// H₀ = 1，H₁ = 2x，Hₙ₊₁ = 2x·Hₙ − 2n·Hₙ₋₁。
/// W1 验收：与 numpy `hermval` 参考值相对误差 < 1×10⁻¹²（fixtures/reference/hermite_reference.json）。
public enum Hermite {

    /// 单点求值 Hₙ(x)
    public static func value(_ n: Int, _ x: Double) -> Double {
        precondition(n >= 0, "Hermite order n must be >= 0")
        switch n {
        case 0: return 1
        case 1: return 2 * x
        default:
            var hm2 = 1.0          // H_{n-2}
            var hm1 = 2 * x        // H_{n-1}
            var current = 0.0
            for k in 1..<(n) {
                current = 2 * x * hm1 - 2 * Double(k) * hm2
                hm2 = hm1
                hm1 = current
            }
            return current
        }
    }

    /// 向量求值 Hₙ(xs)
    public static func values(_ n: Int, xs: [Double]) -> [Double] {
        xs.map { value(n, $0) }
    }

    /// 归一化谐振子波函数 ψₙ(x)（自然单位 ħ = m = ω = 1）：
    /// ψₙ = Hₙ(x)·e^(−x²/2) / √(2ⁿ n! √π)
    /// 归一化系数以对数域累积计算，避免 2ⁿn! 大阶溢出。
    public static func normalizedWavefunction(_ n: Int, xs: [Double]) -> [Double] {
        precondition(n >= 0, "Hermite order n must be >= 0")
        // ln(2^n n!) = n ln2 + Σ ln k；归一化系数 = exp(-ln(2^n n!)/2 - ln(√π)/2)
        var lnCoef = Double(n) * log(2.0)
        if n > 1 { for k in 2...n { lnCoef += log(Double(k)) } }
        let pre = exp(-0.5 * lnCoef - 0.25 * log(.pi))
        return xs.map { pre * value(n, $0) * exp(-$0 * $0 / 2) }
    }
}
