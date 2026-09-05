import Foundation

/// 连带拉盖尔多项式 L^α_n(x) 与氢原子径向波函数——中心力场模块（笔记 17）的主力依赖。
///
/// 递推（与 scipy.special.genlaguerre 数值一致）：
/// (n+1)L^α_{n+1} = (2n+1+α−x)·L^α_n − (n+α)·L^α_{n−1}
/// 参考值构造（W1 fixtures）：L^α_n = (−1)^α d^α/dx^α L_{n+α}（numpy 多项式求导），
/// 单测逐点互验 < 1×10⁻¹²（fixtures/reference/laguerre_reference.json）。
public enum Laguerre {

    /// 单点求值 L^α_n(x)。α ≥ 0（氢原子取 α = 2l+1）。
    public static func assocLaguerre(_ n: Int, alpha: Double, _ x: Double) -> Double {
        precondition(n >= 0, "Laguerre 阶 n must be >= 0")
        precondition(alpha >= 0, "Laguerre alpha must be >= 0")
        switch n {
        case 0: return 1
        case 1: return 1 + alpha - x
        default:
            var lm2 = 1.0                    // L^α_{n-2}
            var lm1 = 1 + alpha - x          // L^α_{n-1}
            var current = 0.0
            for k in 1..<n {
                current = ((2 * Double(k) + 1 + alpha - x) * lm1
                           - (Double(k) + alpha) * lm2) / Double(k + 1)
                lm2 = lm1
                lm1 = current
            }
            return current
        }
    }

    /// 氢原子径向波函数 R_{nl}(r)，长度单位 a₀（自然单位 ħ=m=e=1）：
    /// R_{nl} = N (2r/n)^l · L^{2l+1}_{n-l-1}(2r/n) · e^{−r/n}
    /// N = √[ (2/n)³ (n−l−1)! / (2n (n+l)!) ]，对数域计算防溢出。
    public static func hydrogenRadial(_ n: Int, _ l: Int, rs: [Double]) -> [Double] {
        precondition(n >= 1, "主量子数 n >= 1")
        precondition(l >= 0 && l < n, "角量子数 0 <= l < n")

        let lnN = 1.5 * log(2.0 / Double(n))
            + 0.5 * (lnFactorial(n - l - 1) - log(Double(2 * n)) - lnFactorial(n + l))
        let norm = exp(lnN)

        return rs.map { r in
            let rho = 2.0 * r / Double(n)
            return norm * pow(rho, Double(l))
                * assocLaguerre(n - l - 1, alpha: Double(2 * l + 1), rho)
                * exp(-r / Double(n))
        }
    }

    /// ln(k!)——k ≤ 170 用 Γ 函数精确路径，更大阶直接 lgamma。
    private static func lnFactorial(_ k: Int) -> Double {
        k <= 1 ? 0 : lgamma(Double(k + 1))
    }
}
