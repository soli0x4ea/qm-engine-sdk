import Foundation
import Accelerate

/// 实对称矩阵特征分解：Accelerate LAPACK `dsyevr` 直调（W3 包体裁决定案）。
///
/// W1 曾以 mlx-swift 0.31.6 eigh（CPU float64）实现，并与 dsyevr 互验 1×10⁻¹²。
/// W3 实测发现 MLX 运行时即使 CPU 执行，流初始化也会为 GPU 设备建流并强制
/// 创建 Metal 设备（Scheduler 构造路径），在无 Metal 的 iOS 模拟器上直接 abort；
/// 且为单一 LAPACK 调用引入 26 MB 包体。故移除 MLX、直调 dsyevr——
/// 数值行为与 W1 互验基线一致（dsyevr 即当时的参考实现）。
/// MLX 版实现与裁决数据见 git 历史（W2 提交 876c59d、W3 提交记录）。
public enum SymEigh {

    /// 实对称矩阵全谱特征分解：A = V·diag(w)·Vᵀ。
    ///
    /// 输入须为对称矩阵；仅读取 uplo 指示的三角（与 numpy.linalg.eigh 语义对齐）。
    /// - Parameters:
    ///   - aFlat: 行主序 n×n 扁平数组
    ///   - uplo: "L" 读下三角 / "U" 读上三角
    /// - Returns: (w 升序, v 行主序扁平；第 i 个特征向量 = v 的第 i 列，即分量 v[j·n + i])
    public static func eighSymmetric(
        _ aFlat: [Double],
        n: Int,
        uplo: String = "L"
    ) -> (w: [Double], v: [Double]) {
        precondition(aFlat.count == n * n, "eighSymmetric: expected \(n * n) elements, got \(aFlat.count)")
        precondition(n >= 1, "eighSymmetric: n must be >= 1")
        precondition(uplo == "L" || uplo == "U", "eighSymmetric: uplo must be \"L\" or \"U\"")

        // 行主序 → 列主序（转置拷贝），使 uplo 三角语义与 numpy 一致
        var a = [Double](repeating: 0, count: n * n)
        for i in 0..<n {
            for j in 0..<n {
                a[i + j * n] = aFlat[i * n + j]
            }
        }

        var jobz: Int8 = 86   // 'V'：特征值 + 特征向量
        var range: Int8 = 65  // 'A'：全谱
        var uploC: Int8 = uplo == "L" ? 76 : 85
        var nInt = Int32(n)
        var lda = Int32(n)
        var vl = 0.0, vu = 0.0
        var il: Int32 = 0, iu: Int32 = 0
        var abstol = 0.0
        var m = Int32(0)
        var w = [Double](repeating: 0, count: n)
        var z = [Double](repeating: 0, count: n * n)
        var ldz = Int32(n)
        var isuppz = [Int32](repeating: 0, count: 2 * n)
        var work = [Double](repeating: 0, count: 1)
        var lwork: Int32 = -1
        var iwork = [Int32](repeating: 0, count: 1)
        var liwork: Int32 = -1
        var info = Int32(0)

        // 工作区容量查询（lwork = -1）
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
        precondition(m == Int32(n), "dsyevr returned m = \(m), expected \(n)")

        // 列主序 z（第 s 列 = 第 s 个特征向量）→ 行主序 v
        var v = [Double](repeating: 0, count: n * n)
        for s in 0..<n {
            for j in 0..<n {
                v[j * n + s] = z[j + s * n]
            }
        }
        return (w, v)
    }
}
