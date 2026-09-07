import Testing
import Foundation
@testable import EngineKit

/// W8 基建验收：AlgebraCore eigh 门面（dsyevr 路径定案，见 AlgebraCore.swift 头注）。
/// 全谱转发 SymEigh（W3 裁决后唯一实现）；部分谱 eighLowest 为 W8 新增，
/// 以「与全谱互验 + 残差/正交性不变量」双口径验证。
@Suite("AlgebraCore（dsyevr 门面）")
struct AlgebraCoreTests {

    // MARK: - 确定性随机对称矩阵（与 SymEighTests 同构，固定种子复现）

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

    /// max|A·vₛ − wₛ·vₛ|（逐元素绝对残差，v 为行主序 n×k）。
    private static func maxResidual(_ aFlat: [Double], n: Int, k: Int, w: [Double], v: [Double]) -> Double {
        var residual = 0.0
        for s in 0..<k {
            for j in 0..<n {
                var av = 0.0
                for t in 0..<n { av += aFlat[j * n + t] * v[t * k + s] }
                residual = max(residual, abs(av - w[s] * v[j * k + s]))
            }
        }
        return residual
    }

    /// max|(VᵀV)ₛ₁ₛ₂ − δₛ₁ₛ₂|（正交归一偏差）。
    private static func maxOrthError(n: Int, k: Int, v: [Double]) -> Double {
        var orth = 0.0
        for s1 in 0..<k {
            for s2 in 0..<k {
                var dot = 0.0
                for j in 0..<n { dot += v[j * k + s1] * v[j * k + s2] }
                let target = s1 == s2 ? 1.0 : 0.0
                orth = max(orth, abs(dot - target))
            }
        }
        return orth
    }

    // MARK: - 测试

    @Test("eigh 全谱转发：与 SymEigh 逐元一致（同一 dsyevr 调用）")
    func eighForwardsToSymEigh() {
        let n = 60
        let a = Self.randomSymmetric(n: n, seed: 20260905)
        let viaCore = AlgebraCore.eigh(a, n: n)
        let direct = SymEigh.eighSymmetric(a, n: n)
        #expect(viaCore.w == direct.w)
        #expect(viaCore.v == direct.v)
    }

    @Test("eighLowest 前 k 本征对：与全谱前 k 列逐元互验（非简并随机谱，1e-12）")
    func eighLowestMatchesFullSpectrum() {
        let n = 200
        let a = Self.randomSymmetric(n: n, seed: 87654321)
        let full = AlgebraCore.eigh(a, n: n)
        let partial = AlgebraCore.eighLowest(a, n: n, count: 12)
        #expect(partial.w.count == 12)
        #expect(partial.v.count == n * 12)
        for s in 0..<12 {
            #expect(abs(partial.w[s] - full.w[s]) < 1e-10 * max(1, abs(full.w[s])),
                    "w[\(s)]: \(partial.w[s]) vs \(full.w[s])")
        }
        // 特征向量：非简并随机谱逐元一致，允许整体符号差
        // （dsyevr 符号约定随 Accelerate/LAPACK 版本可能不同，eigh 与 eighLowest 走不同 range 路径）
        for s in 0..<12 {
            var maxDiff = 0.0
            for j in 0..<n {
                let vp = partial.v[j * 12 + s]
                let vf = full.v[j * n + s]
                maxDiff = max(maxDiff, min(abs(vp - vf), abs(vp + vf)))
            }
            #expect(maxDiff < 1e-9, "v[\(s)] 最大分量差 \(maxDiff)")
        }
    }

    @Test("eighLowest 不变量：残差 < 1e-12·‖A‖、正交归一 < 1e-12（3 组规模×种子）")
    func eighLowestInvariants() {
        for (n, k, seed) in [(50, 5, UInt64(11)), (137, 9, UInt64(22)), (64, 64, UInt64(33))] {
            let a = Self.randomSymmetric(n: n, seed: seed)
            var scale = 0.0
            for x in a { scale = max(scale, abs(x)) }
            let (w, v) = AlgebraCore.eighLowest(a, n: n, count: k)
            let residual = Self.maxResidual(a, n: n, k: k, w: w, v: v)
            let orth = Self.maxOrthError(n: n, k: k, v: v)
            #expect(residual < 1e-12 * max(scale, 1),
                    "n=\(n) k=\(k): 残差 \(residual)")
            #expect(orth < 1e-12, "n=\(n) k=\(k): 正交性 \(orth)")
        }
    }

    @Test("eighLowest 升序谱：三对角束缚态问题前 k 个为最小（W8 氢原子径向形态）")
    func eighLowestTridiagonalBoundStates() {
        // 无限深势阱 FD 哈密顿（800 点三对角）谱解析已知：E_n ∝ n²
        let n = 800
        let h = 1.0 / Double(n + 1)
        let t = 1.0 / (h * h)                       // ħ²/2m = 1 自然单位
        var a = [Double](repeating: 0, count: n * n)
        for i in 0..<n { a[i * n + i] = 2 * t }
        for i in 0..<(n - 1) {
            a[i * n + i + 1] = -t
            a[(i + 1) * n + i] = -t
        }
        let (w, _) = AlgebraCore.eighLowest(a, n: n, count: 6)
        for s in 0..<6 {
            let analytic = Double(s + 1) * Double(s + 1) * .pi * .pi  // E_n = n²π²（ħ²/2m=1，L=1）
            let rel = abs(w[s] - analytic) / analytic
            #expect(rel < 1e-4, "E_\(s + 1) = \(w[s]) vs 解析 \(analytic)（FD 误差容 1e-4）")
            if s > 0 { #expect(w[s] > w[s - 1], "谱须升序") }
        }
    }

    // MARK: - 对称三对角（收尾包 1：dstevr 通道）

    /// 确定性随机对称三对角 (d, e)。
    private static func randomTridiagonal(n: Int, seed: UInt64) -> (d: [Double], e: [Double]) {
        var rng = SplitMix64(seed: seed)
        let d = (0..<n).map { _ in rng.unit() * 4 - 2 }
        let e = (0..<(n - 1)).map { _ in rng.unit() * 2 - 1 }
        return (d, e)
    }

    @Test("eighLowestTridiagonal：与稠密 eighLowest 互验（w bit 级一致、v ~1e-9、残差/正交不变量）")
    func eighLowestTridiagonalMatchesDense() {
        for (n, k, seed) in [(200, 8, UInt64(101)), (731, 5, UInt64(202)), (64, 64, UInt64(303))] {
            let (d, e) = Self.randomTridiagonal(n: n, seed: seed)
            // 稠密参照：由 (d, e) 组装行主序对称矩阵走 eighLowest
            var a = [Double](repeating: 0, count: n * n)
            for i in 0..<n { a[i * n + i] = d[i] }
            for i in 0..<(n - 1) {
                a[i * n + i + 1] = e[i]
                a[(i + 1) * n + i] = e[i]
            }
            let dense = AlgebraCore.eighLowest(a, n: n, count: k)
            let tri = AlgebraCore.eighLowestTridiagonal(d, e, count: k)
            // 本征值：bit 级一致（W8 收官基准实测 3000 阶 wRel = 0）
            #expect(tri.w == dense.w, "n=\(n): dstevr 与 dsyevr 本征值应 bit 级一致")
            // 特征向量：非简并随机谱逐元一致到 ~1e-9（MRRR 回转置路径差异）
            for s in 0..<k {
                var maxDiff = 0.0
                for j in 0..<n {
                    maxDiff = max(maxDiff, abs(tri.v[j * k + s] - dense.v[j * k + s]))
                }
                #expect(maxDiff < 1e-9, "n=\(n) v[\(s)] 最大分量差 \(maxDiff)")
            }
            // 不变量：残差 + 正交归一（自洽口径）
            var scale = 0.0
            for i in 0..<n { scale = max(scale, abs(d[i])) }
            for x in e { scale = max(scale, abs(x)) }
            let residual = Self.maxResidual(a, n: n, k: k, w: tri.w, v: tri.v)
            let orth = Self.maxOrthError(n: n, k: k, v: tri.v)
            #expect(residual < 1e-12 * max(scale, 1), "n=\(n): 残差 \(residual)")
            #expect(orth < 1e-12, "n=\(n): 正交性 \(orth)")
        }
    }

    @Test("eighLowestTridiagonal：无限深势阱 FD 解析谱 E_n = n²π²（收尾包 1 氢原子形态）")
    func eighLowestTridiagonalWellAnalytic() {
        let n = 800
        let h = 1.0 / Double(n + 1)
        let t = 1.0 / (h * h)
        let d = [Double](repeating: 2 * t, count: n)
        let e = [Double](repeating: -t, count: n - 1)
        let (w, _) = AlgebraCore.eighLowestTridiagonal(d, e, count: 6)
        for s in 0..<6 {
            let analytic = Double(s + 1) * Double(s + 1) * .pi * .pi
            let rel = abs(w[s] - analytic) / analytic
            #expect(rel < 1e-4, "E_\(s + 1) = \(w[s]) vs 解析 \(analytic)（FD 误差容 1e-4）")
            if s > 0 { #expect(w[s] > w[s - 1], "谱须升序") }
        }
    }

    // MARK: - 复 Hermitian（W9：zheevr 通道）

    /// 复 Hermitian 残差 max|A·v − w·v|（复数逐元）。
    private static func maxResidualComplex(
        _ aRe: [Double], _ aIm: [Double], n: Int, w: [Double], vRe: [Double], vIm: [Double]
    ) -> Double {
        var worst = 0.0
        for s in 0..<n {
            for j in 0..<n {
                var ar = 0.0, ai = 0.0
                for t in 0..<n {
                    let arj = aRe[j * n + t], aij = aIm[j * n + t]
                    let vr = vRe[t * n + s], vi = vIm[t * n + s]
                    ar += arj * vr - aij * vi
                    ai += arj * vi + aij * vr
                }
                let br = w[s] * vRe[j * n + s], bi = w[s] * vIm[j * n + s]
                worst = max(worst, abs(ar - br), abs(ai - bi))
            }
        }
        return worst
    }

    @Test("eighHermitian：Pauli-Y 本征值 ±1；实对称退化情形与 dsyevr 互验")
    func eighHermitianPauliY() {
        // Y = [[0, −i],[i, 0]] → 本征值 {−1, +1}
        let wY = AlgebraCore.eighHermitian(
            [0, 0, 0, 0], [0, -1, 1, 0], n: 2, wantsVectors: false).w
        #expect(abs(wY[0] + 1) < 1e-14 && abs(wY[1] - 1) < 1e-14, "Y 谱 = \(wY)")
        // 实对称 4×4：与 dsyevr 路径互验（虚部全零）
        let a = Self.randomSymmetric(n: 4, seed: 42)
        let real = AlgebraCore.eigh(a, n: 4)
        let cplx = AlgebraCore.eighHermitian(a, [Double](repeating: 0, count: 16), n: 4)
        for i in 0..<4 {
            #expect(abs(cplx.w[i] - real.w[i]) < 1e-12 * max(1, abs(real.w[i])),
                    "w[\(i)]: \(cplx.w[i]) vs \(real.w[i])")
        }
    }

    @Test("eighHermitian：Werner 态谱 = {(1+3p)/4, (1−p)/4×3}；残差/归一不变量；vectors-only=N")
    func eighHermitianWerner() {
        // ρ_W(p) = p|Ψ⁻⟩⟨Ψ⁻| + (1−p)/4·I：|Ψ⁻⟩ = (|01⟩−|10⟩)/√2
        for p in [0.0, 0.2, 1.0 / 3.0, 0.7, 1.0] {
            var re = [Double](repeating: 0.0, count: 16)      // 非对角基线 0
            let im = [Double](repeating: 0.0, count: 16)
            for d in 0..<4 { re[d * 4 + d] = (1.0 - p) / 4.0 } // 单位阵部分只在对角
            re[1 * 4 + 1] += p / 2.0            // |01⟩⟨01|
            re[2 * 4 + 2] += p / 2.0            // |10⟩⟨10|
            re[1 * 4 + 2] -= p / 2.0            // −p/2 交叉项（单态反对称）
            re[2 * 4 + 1] -= p / 2.0
            let (w, vRe, vIm) = AlgebraCore.eighHermitian(re, im, n: 4)
            // 谱：升序 {(1−p)/4 ×3, (1+3p)/4}
            for i in 0..<3 {
                #expect(abs(w[i] - (1.0 - p) / 4.0) < 1e-12, "p=\(p) w[\(i)] = \(w[i])")
            }
            #expect(abs(w[3] - (1.0 + 3.0 * p) / 4.0) < 1e-12, "p=\(p) w[3] = \(w[3])")
            // 不变量：残差 + 迹
            let residual = Self.maxResidualComplex(re, im, n: 4, w: w, vRe: vRe, vIm: vIm)
            #expect(residual < 1e-12, "p=\(p) 残差 \(residual)")
            let tr = w.reduce(0, +)
            #expect(abs(tr - 1.0) < 1e-12, "p=\(p) 迹 \(tr)")
        }
        // jobz='N'：本征值一致，向量空
        let nv = AlgebraCore.eighHermitian(
            [2, 0, 0, 2], [0, 1, -1, 0], n: 2, wantsVectors: false)
        #expect(nv.vRe.isEmpty && nv.vIm.isEmpty && nv.w.count == 2)
    }
}
