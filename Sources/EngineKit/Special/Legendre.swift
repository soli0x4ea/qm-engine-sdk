import Foundation

/// 勒让德多项式 P_l(x) 与连带勒让德 P^m_l(x)（Condon-Shortley 约定，与 scipy.special.lpmv 一致）。
///
/// P_l：三项递推 (l+1)P_{l+1} = (2l+1)x·P_l − l·P_{l−1}
/// P^m_l：m 级升阶 + l 向上递推（Numerical Recipes 经典稳定算法，无需求导）：
///   P^m_m  = (−1)^m (2m−1)!! (1−x²)^{m/2}
///   P^{m}_{m+1} = x(2m+1) P^m_m
///   (l−m)P^m_l = (2l−1)x P^m_{l−1} − (l+m−1) P^m_{l−2}
/// W8 球谐函数（角向部分）基于本文件。参考值：numpy.polynomial.legendre 求导构造
/// （fixtures/reference/legendre_reference.json），互验 < 1×10⁻¹²。
public enum Legendre {

    /// 单点求值 P_l(x)，x ∈ [−1, 1]。
    public static func value(_ l: Int, _ x: Double) -> Double {
        precondition(l >= 0, "Legendre 阶 l must be >= 0")
        switch l {
        case 0: return 1
        case 1: return x
        default:
            var pm2 = 1.0
            var pm1 = x
            var current = 0.0
            for k in 1..<l {
                current = ((2 * Double(k) + 1) * x * pm1 - Double(k) * pm2) / Double(k + 1)
                pm2 = pm1
                pm1 = current
            }
            return current
        }
    }

    /// 连带勒让德 P^m_l(x)，Condon-Shortley 相位。|m| ≤ l。
    public static func assoc(_ l: Int, _ m: Int, _ x: Double) -> Double {
        precondition(l >= 0, "l must be >= 0")
        precondition(abs(m) <= l, "需要 |m| <= l")

        // 负阶：P^{−|m|}_l = (−1)^{|m|} (l−|m|)!/(l+|m|)! · P^{|m|}_l
        if m < 0 {
            let am = -m
            let ratio = exp(lnFactorial(l - am) - lnFactorial(l + am))
            let sign: Double = am % 2 == 0 ? 1 : -1
            return sign * ratio * assoc(l, am, x)
        }

        // 1) 升到 P^m_m
        var pmm = 1.0
        if m > 0 {
            let somx2 = sqrt(max(0.0, (1 - x) * (1 + x)))  // (1-x²)^{1/2}，端点钳 0
            var fact = 1.0
            for _ in 0..<m {
                pmm *= -fact * somx2
                fact += 2
            }
        }
        if l == m { return pmm }

        // 2) P^{m}_{m+1}
        var pmmp1 = x * Double(2 * m + 1) * pmm
        if l == m + 1 { return pmmp1 }

        // 3) l 向上递推
        var pll = 0.0
        var ll = m + 2
        repeat {
            pll = (x * Double(2 * ll - 1) * pmmp1 - Double(ll + m - 1) * pmm) / Double(ll - m)
            pmm = pmmp1
            pmmp1 = pll
            ll += 1
        } while ll <= l
        return pll
    }

    private static func lnFactorial(_ k: Int) -> Double {
        k <= 1 ? 0 : lgamma(Double(k + 1))
    }
}
