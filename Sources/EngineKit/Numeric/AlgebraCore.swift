import Foundation
import Accelerate

// MARK: - AlgebraCore（W8 基建：秒级/策略档线代核门面）

/// 秒级/策略档模块的稠密线代门面——W8 起各 tier-3 模块经由此处调用 eigh。
///
/// **eigh 路径三指标定案（W8 落案，计划 §W8 基建项）**
///
/// | 指标 | dsyevr 路径（本门面 + SymEigh） | MLXBridge 路径（W2 历史方案） |
/// |---|---|---|
/// | 基线通过率 | W1 起承担全部 eigh 模块，与 MLX 互验 1×10⁻¹²；W3 裁决后为唯一实现，321 项既有测试全绿 | W1 互验等价，但 W3 实测无 Metal 的 iOS 模拟器直接 abort |
/// | 性能 | Accelerate LAPACK 直调，零中间层；部分本征对（RANGE='I'）只回传前 k 对 | CPU float64 执行，但运行时强制建 Metal 设备/流 |
/// | 代码量 | SymEigh + 本文件约 250 行 Swift | mlx-swift 包体 +26 MB（W3 实测），且引入 GPU 初始化副作用 |
///
/// 结论：沿用既有 dsyevr 路径（W3 包体裁决的自然延伸），本门面只做命名收敛与
/// 「前 k 个本征对」部分谱扩展；计划文本中「MLXBridge vs dsyevr 三指标定案」为
/// 陈述性记录，实际裁决已在 W3 完成。MLX 版实现见 git 历史（W2 提交 876c59d）。
///
/// 后续 tier-3 模块（W8 线代批 1/2 共 8 模块）一律经由 `eigh` / `eighLowest` 调用；
/// `SymEigh` 保留为底层实现，既有 51 模块零改动。
public enum AlgebraCore {

    /// 实对称矩阵全谱特征分解：A = V·diag(w)·Vᵀ（转发 `SymEigh.eighSymmetric`）。
    ///
    /// - Parameters:
    ///   - aFlat: 行主序 n×n 扁平数组（只读取 uplo 指示三角，语义对齐 numpy.linalg.eigh）
    /// - Returns: (w 升序全谱, v 行主序 n×n；第 i 个特征向量 = v 的第 i 列，分量 v[j·n + i])
    public static func eigh(
        _ aFlat: [Double], n: Int, uplo: String = "L"
    ) -> (w: [Double], v: [Double]) {
        SymEigh.eighSymmetric(aFlat, n: n, uplo: uplo)
    }

    /// 实对称矩阵**前 k 个最小本征对**（升序）：dsyevr RANGE='I' 部分谱路径。
    ///
    /// 秒级档降维主通道：3000² 阶问题只需前几个束缚态时，MRRR 阶段只计算/回传
    /// 前 k 对（三对角化仍为全矩阵 O(n³)，但省去其余 n−k 个本征对的回转置与正交化），
    /// 与 Python 端 `np.linalg.eigh` 全谱后取前 k 列数值一致（同一矩阵同一谱，
    /// 非简并态特征向量逐元一致到 ~1e-9，符号约定见测试侧对齐）。
    ///
    /// - Parameters:
    ///   - aFlat: 行主序 n×n 扁平数组
    ///   - k: 需要的最小本征对个数（1 ≤ k ≤ n）
    /// - Returns: (w 升序前 k 个, v 行主序 n×k；第 s 个特征向量 = 第 s 列，分量 v[j·k + s])
    public static func eighLowest(
        _ aFlat: [Double], n: Int, count k: Int, uplo: String = "L"
    ) -> (w: [Double], v: [Double]) {
        precondition(aFlat.count == n * n, "eighLowest: expected \(n * n) elements, got \(aFlat.count)")
        precondition(n >= 1, "eighLowest: n must be >= 1")
        precondition(k >= 1 && k <= n, "eighLowest: k must be in 1...n, got \(k)")
        precondition(uplo == "L" || uplo == "U", "eighLowest: uplo must be \"L\" or \"U\"")

        // 行主序 → 列主序（转置拷贝），使 uplo 三角语义与 numpy 一致（同 SymEigh）
        var a = [Double](repeating: 0, count: n * n)
        for i in 0..<n {
            for j in 0..<n {
                a[i + j * n] = aFlat[i * n + j]
            }
        }

        var jobz: Int8 = 86   // 'V'：特征值 + 特征向量
        var range: Int8 = 73  // 'I'：按序号区间 IL..IU（前 k 小）
        var uploC: Int8 = uplo == "L" ? 76 : 85
        var nInt = Int32(n)
        var lda = Int32(n)
        var vl = 0.0, vu = 0.0
        var il: Int32 = 1
        var iu = Int32(k)
        var abstol = 0.0      // 0 = 用默认精度（N*eps·|T|，LAPACK 推荐）
        var m = Int32(0)
        var w = [Double](repeating: 0, count: k)
        var z = [Double](repeating: 0, count: n * k)
        var ldz = Int32(n)
        var isuppz = [Int32](repeating: 0, count: 2 * max(1, k))
        var work = [Double](repeating: 0, count: 1)
        var lwork: Int32 = -1
        var iwork = [Int32](repeating: 0, count: 1)
        var liwork: Int32 = -1
        var info = Int32(0)

        dsyevr_(&jobz, &range, &uploC, &nInt, &a, &lda, &vl, &vu, &il, &iu, &abstol,
                &m, &w, &z, &ldz, &isuppz, &work, &lwork, &iwork, &liwork, &info)
        precondition(info == 0, "dsyevr workspace query failed, info = \(info)")
        let optimalWork = max(1, Int(work[0]))
        let optimalIWork = max(1, Int(iwork[0]))
        work = [Double](repeating: 0, count: optimalWork)
        iwork = [Int32](repeating: 0, count: optimalIWork)
        lwork = Int32(optimalWork)
        liwork = Int32(optimalIWork)

        dsyevr_(&jobz, &range, &uploC, &nInt, &a, &lda, &vl, &vu, &il, &iu, &abstol,
                &m, &w, &z, &ldz, &isuppz, &work, &lwork, &iwork, &liwork, &info)
        precondition(info == 0, "dsyevr failed, info = \(info)")
        precondition(m == Int32(k), "dsyevr returned m = \(m), expected \(k)")

        // 列主序 z（n×k，第 s 列 = 第 s 个特征向量）→ 行主序 v（n×k）
        var v = [Double](repeating: 0, count: n * k)
        for s in 0..<k {
            for j in 0..<n {
                v[j * k + s] = z[j + s * n]
            }
        }
        return (w, v)
    }

    // MARK: 对称三对角（收尾包 1 基建：FD 三对角哈密顿直达通道）

    /// 实对称**三对角**矩阵前 k 个最小本征对（升序）：`dstevr` RANGE='I' 直调。
    ///
    /// 收尾包 1 性能专项新增：有限差分哈密顿天然三对角，此前经 `eighLowest`
    /// 须先组装/拷贝 n² 稠密矩阵并付 O(n³) 三对角化；本通道直接从 d/e 表示出发
    /// （MRRR O(n²)），3000 阶实测 1441 ms → 4.8 ms。
    ///
    /// 数值口径（W8 收官基准实测）：本征值与 dsyevr 路径 **bit 级一致**；
    /// 本征向量逐元相对差 ~1×10⁻⁹（两路径回转置差异，与 eighLowest 文档声明的
    /// 非简并一致度同档），处于各模块 fixture 容差（1e-6）之内。
    ///
    /// - Parameters:
    ///   - d: 长度 n 的对角元
    ///   - e: 长度 n−1 的次对角元（内部补零到 n，LAPACK 就地约定）
    ///   - k: 需要的最小本征对个数（1 ≤ k ≤ n）
    /// - Returns: (w 升序前 k 个, v 行主序 n×k；第 s 个特征向量 = 第 s 列，分量 v[j·k + s])
    public static func eighLowestTridiagonal(
        _ d: [Double], _ e: [Double], count k: Int
    ) -> (w: [Double], v: [Double]) {
        let n = d.count
        precondition(e.count == n - 1, "eighLowestTridiagonal: e must have n-1 elements, got \(e.count)")
        precondition(n >= 1, "eighLowestTridiagonal: n must be >= 1")
        precondition(k >= 1 && k <= n, "eighLowestTridiagonal: k must be in 1...n, got \(k)")

        var dd = d
        var ee = e + [0.0]     // dstevr 就地读取 e[0..n-2]，长度按 LAPACK 约定补到 n
        var jobz: Int8 = 86    // 'V'：特征值 + 特征向量
        var range: Int8 = 73   // 'I'：按序号区间 IL..IU（前 k 小）
        var nInt = Int32(n)
        var vl = 0.0, vu = 0.0
        var il: Int32 = 1
        var iu = Int32(k)
        var abstol = 0.0       // 0 = 默认精度（N*eps·|T|）
        var m = Int32(0)
        var w = [Double](repeating: 0, count: k)
        var z = [Double](repeating: 0, count: n * k)
        var ldz = Int32(n)
        var isuppz = [Int32](repeating: 0, count: 2 * max(1, k))
        var work = [Double](repeating: 0, count: 1)
        var lwork: Int32 = -1
        var iwork = [Int32](repeating: 0, count: 1)
        var liwork: Int32 = -1
        var info = Int32(0)

        dstevr_(&jobz, &range, &nInt, &dd, &ee, &vl, &vu, &il, &iu, &abstol,
                &m, &w, &z, &ldz, &isuppz, &work, &lwork, &iwork, &liwork, &info)
        precondition(info == 0, "dstevr workspace query failed, info = \(info)")
        let optimalWork = max(1, Int(work[0]))
        let optimalIWork = max(1, Int(iwork[0]))
        work = [Double](repeating: 0, count: optimalWork)
        iwork = [Int32](repeating: 0, count: optimalIWork)
        lwork = Int32(optimalWork)
        liwork = Int32(optimalIWork)

        dstevr_(&jobz, &range, &nInt, &dd, &ee, &vl, &vu, &il, &iu, &abstol,
                &m, &w, &z, &ldz, &isuppz, &work, &lwork, &iwork, &liwork, &info)
        precondition(info == 0, "dstevr failed, info = \(info)")
        precondition(m == Int32(k), "dstevr returned m = \(m), expected \(k)")

        // 列主序 z（n×k）→ 行主序 v（n×k），与 eighLowest 同约定
        var v = [Double](repeating: 0, count: n * k)
        for s in 0..<k {
            for j in 0..<n {
                v[j * k + s] = z[j + s * n]
            }
        }
        return (w, v)
    }

    // MARK: 复 Hermitian 矩阵（W9 基建：纠缠模块密度矩阵 eigh，zheevr 直调）

    /// 复 Hermitian 矩阵特征分解：A = V·diag(w)·V†（Accelerate LAPACK `zheevr` 直调）。
    ///
    /// W9 起纠缠度量模块（4×4 密度矩阵 / 偏置转置 / Wootters 矩阵）经由本通道；
    /// 实对称矩阵请用 `eigh` / `eighLowest`（dsyevr 更快）。
    ///
    /// - Parameters:
    ///   - aRe/aIm: 行主序 n×n 实部/虚部扁平数组（只读取 uplo 指示三角，
    ///     虚部三角取符号约定与 numpy.linalg.eigh 一致：A[i,j] = conj(A[j,i]）
    ///   - wantsVectors: false 时 jobz='N' 只回传本征值（vRe/vIm 为空）
    /// - Returns: (w 升序全谱, vRe/vIm 行主序 n×n；第 i 个特征向量 = 第 i 列，
    ///   分量 (vRe[j·n + i], vIm[j·n + i])）
    public static func eighHermitian(
        _ aRe: [Double], _ aIm: [Double], n: Int, uplo: String = "L",
        wantsVectors: Bool = true
    ) -> (w: [Double], vRe: [Double], vIm: [Double]) {
        precondition(aRe.count == n * n && aIm.count == n * n,
                     "eighHermitian: expected \(n * n) elements, got \(aRe.count)/\(aIm.count)")
        precondition(n >= 1, "eighHermitian: n must be >= 1")
        precondition(uplo == "L" || uplo == "U", "eighHermitian: uplo must be \"L\" or \"U\"")

        // 行主序 → 列主序打包（uplo 三角语义对齐 numpy：只读指示三角）
        var a = [__CLPK_doublecomplex](repeating: __CLPK_doublecomplex(r: 0, i: 0), count: n * n)
        for i in 0..<n {
            for j in 0..<n {
                a[i + j * n] = __CLPK_doublecomplex(r: aRe[i * n + j], i: aIm[i * n + j])
            }
        }

        var jobz: Int8 = wantsVectors ? 86 : 78   // 'V' / 'N'
        var range: Int8 = 65                      // 'A'：全谱
        var uploC: Int8 = uplo == "L" ? 76 : 85
        var nInt = Int32(n)
        var lda = Int32(n)
        var vl = 0.0, vu = 0.0
        var il: Int32 = 0, iu: Int32 = 0
        var abstol = 0.0
        var m = Int32(0)
        var w = [Double](repeating: 0, count: n)
        var z = [__CLPK_doublecomplex](repeating: __CLPK_doublecomplex(r: 0, i: 0), count: n * n)
        var ldz = Int32(n)
        var isuppz = [Int32](repeating: 0, count: 2 * n)
        var work = [__CLPK_doublecomplex](repeating: __CLPK_doublecomplex(r: 0, i: 0), count: 1)
        var lwork: Int32 = -1
        var rwork = [Double](repeating: 0, count: 1)
        var lrwork: Int32 = -1
        var iwork = [Int32](repeating: 0, count: 1)
        var liwork: Int32 = -1
        var info = Int32(0)

        // 工作区容量查询（lwork = -1）
        zheevr_(&jobz, &range, &uploC, &nInt, &a, &lda, &vl, &vu, &il, &iu, &abstol,
                &m, &w, &z, &ldz, &isuppz, &work, &lwork, &rwork, &lrwork, &iwork, &liwork, &info)
        precondition(info == 0, "zheevr workspace query failed, info = \(info)")
        let optimalWork = max(1, Int(work[0].r))
        let optimalRWork = max(1, Int(rwork[0]))
        let optimalIWork = max(1, Int(iwork[0]))
        work = [__CLPK_doublecomplex](repeating: __CLPK_doublecomplex(r: 0, i: 0), count: optimalWork)
        rwork = [Double](repeating: 0, count: optimalRWork)
        iwork = [Int32](repeating: 0, count: optimalIWork)
        lwork = Int32(optimalWork)
        lrwork = Int32(optimalRWork)
        liwork = Int32(optimalIWork)

        zheevr_(&jobz, &range, &uploC, &nInt, &a, &lda, &vl, &vu, &il, &iu, &abstol,
                &m, &w, &z, &ldz, &isuppz, &work, &lwork, &rwork, &lrwork, &iwork, &liwork, &info)
        precondition(info == 0, "zheevr failed, info = \(info)")
        precondition(Int(m) == n, "zheevr returned m = \(m), expected \(n)")

        guard wantsVectors else { return (w, [], []) }

        // 列主序 z（n×n，第 s 列 = 第 s 个特征向量）→ 行主序 vRe/vIm
        var vRe = [Double](repeating: 0, count: n * n)
        var vIm = [Double](repeating: 0, count: n * n)
        for s in 0..<n {
            for j in 0..<n {
                let c = z[j + s * n]
                vRe[j * n + s] = c.r
                vIm[j * n + s] = c.i
            }
        }
        return (w, vRe, vIm)
    }
}
