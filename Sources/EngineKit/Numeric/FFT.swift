import Foundation
import Accelerate

/// FFT 数据层——vDSP 封装，与 NumPy `fft` 语义逐一对齐。
///
/// 背景：开发计划 W4 原定为「MLX FFT 封装」；门 1 已裁决 MLX 全面移除
/// （包体 +26 MB + 模拟器无 Metal abort），Accelerate 兜底转正，零包体成本。
///
/// 语义约定（与 numpy.fft 一致）：
/// - 正变换 `fft`：X(k) = Σₙ x(n)·e^(−2πi·kn/N)，**不归一化**。
/// - 逆变换 `ifft`：x(n) = (1/N)·Σₖ X(k)·e^(+2πi·kn/N)，**含 1/N**。
/// - `fftfreq`：k ≤ (n−1)/2 取 k/(N·d)，其余取 (k−N)/(N·d)（Nyquist 入负半轴）。
/// - `fftshift` / `ifftshift`：零频居中 / 还原（分割点 ceil(n/2) / floor(n/2)）。
///
/// 实现：vDSP DFT 支持 N = 2ᵏ 与 3·2ᵏ(k≥3)（实测）；其余 N 走 Bluestein
/// chirp-z 算法（卷积化 + 2 幂 FFT），任意 N 均 O(N log N)。
public enum FFT {

    // MARK: - 复数变换

    /// np.fft.fft：正变换，不归一化。re/im 等长，N ≥ 2。
    public static func fft(real re: [Double], imag im: [Double]) -> (re: [Double], im: [Double]) {
        transform(real: re, imag: im, inverse: false)
    }

    /// np.fft.ifft：逆变换，含 1/N 归一化。
    public static func ifft(real re: [Double], imag im: [Double]) -> (re: [Double], im: [Double]) {
        transform(real: re, imag: im, inverse: true)
    }

    /// 实数输入便捷入口（虚部置零）。
    public static func fft(_ real: [Double]) -> (re: [Double], im: [Double]) {
        fft(real: real, imag: [Double](repeating: 0, count: real.count))
    }

    // MARK: - 内核

    private static func transform(
        real re: [Double], imag im: [Double], inverse: Bool
    ) -> (re: [Double], im: [Double]) {
        let n = re.count
        precondition(n >= 2 && im.count == n, "FFT 需要等长输入且 N ≥ 2")
        let out: (re: [Double], im: [Double])
        if vdspSupported(n) {
            out = vdspTransform(real: re, imag: im, inverse: inverse)
        } else {
            out = bluestein(real: re, imag: im, inverse: inverse)
        }
        guard inverse else { return out }
        // numpy ifft 归一化 1/N（vDSP 逆变换不归一化，Bluestein 内部 ifft 已含 1/M）
        let scale = 1.0 / Double(n)
        return (vDSP.multiply(scale, out.re), vDSP.multiply(scale, out.im))
    }

    /// vDSP 实测可用长度：2ᵏ，或 3·2ᵏ（k ≥ 3）。
    static func vdspSupported(_ n: Int) -> Bool {
        guard n >= 2 else { return false }
        var m = n
        while m % 2 == 0 { m /= 2 }
        if m == 1 { return true }                 // 2ᵏ
        return m == 3 && n / 3 >= 8               // 3·2ᵏ, k ≥ 3
    }

    /// vDSP DFT（zrop）：FORWARD 核 e^(−2πikn/N)，INVERSE 核 e^(+2πikn/N)，均不归一化。
    private static func vdspTransform(
        real re: [Double], imag im: [Double], inverse: Bool
    ) -> (re: [Double], im: [Double]) {
        let n = re.count
        guard let setup = vDSP_DFT_zop_CreateSetupD(
            nil, vDSP_Length(n), inverse ? .INVERSE : .FORWARD) else {
            preconditionFailure("vDSP_DFT_zop_CreateSetupD 失败 (N=\(n))")
        }
        defer { vDSP_DFT_DestroySetupD(setup) }
        var outRe = [Double](repeating: 0, count: n)
        var outIm = [Double](repeating: 0, count: n)
        re.withUnsafeBufferPointer { pr in
            im.withUnsafeBufferPointer { pi in
                outRe.withUnsafeMutableBufferPointer { qr in
                    outIm.withUnsafeMutableBufferPointer { qi in
                        vDSP_DFT_ExecuteD(setup,
                                          pr.baseAddress!, pi.baseAddress!,
                                          qr.baseAddress!, qi.baseAddress!)
                    }
                }
            }
        }
        return (outRe, outIm)
    }

    /// Bluestein chirp-z：任意 N 的 DFT 化为长度 M（≥2N−1 的 2 幂）循环卷积。
    /// kn = (k² + n² − (k−n)²)/2 ⇒ X(k) = w(k)·Σₙ a(n)·b(k−n)，w(t)=e^(∓πi·t²/N)。
    private static func bluestein(
        real re: [Double], imag im: [Double], inverse: Bool
    ) -> (re: [Double], im: [Double]) {
        let n = re.count
        var m = 1
        while m < 2 * n - 1 { m <<= 1 }

        // 调频信号 w(t) = exp(sgn·πi·t²/N)；正变换 sgn = −1。
        // t² 以 mod 2N 规约（w 周期 2N），把三角参数压回 [−π, π]——
        // 不规约时参数可达 πN ≈ 1.3e4 rad，大参数舍入使弱信号相对误差劣化 5+ 个量级。
        let sgn = inverse ? 1.0 : -1.0
        let period = 2 * n
        var wRe = [Double](repeating: 0, count: n)
        var wIm = [Double](repeating: 0, count: n)
        for t in 0..<n {
            let ang = sgn * Double.pi * Double((t * t) % period) / Double(n)
            wRe[t] = cos(ang); wIm[t] = sin(ang)
        }

        // a(n) = x(n)·w(n)，零填充至 M
        var aRe = [Double](repeating: 0, count: m)
        var aIm = [Double](repeating: 0, count: m)
        for t in 0..<n {
            aRe[t] = re[t] * wRe[t] - im[t] * wIm[t]
            aIm[t] = re[t] * wIm[t] + im[t] * wRe[t]
        }
        // b：b(0)=1；b(t)=b(M−t)=conj(w(t))（t=1..<N），构成循环卷积核
        var bRe = [Double](repeating: 0, count: m)
        var bIm = [Double](repeating: 0, count: m)
        bRe[0] = 1
        for t in 1..<n {
            bRe[t] = wRe[t]; bIm[t] = -wIm[t]
            bRe[m - t] = wRe[t]; bIm[m - t] = -wIm[t]
        }

        // c = ifft(fft(a)·fft(b))（内部走 vDSP 2 幂路径；ifft 自带 1/M）
        let A = transform(real: aRe, imag: aIm, inverse: false)
        let B = transform(real: bRe, imag: bIm, inverse: false)
        var cRe = [Double](repeating: 0, count: m)
        var cIm = [Double](repeating: 0, count: m)
        for i in 0..<m {
            cRe[i] = A.re[i] * B.re[i] - A.im[i] * B.im[i]
            cIm[i] = A.re[i] * B.im[i] + A.im[i] * B.re[i]
        }
        let c = transform(real: cRe, imag: cIm, inverse: true)
        // transform 的逆含 1/M 归一化（此处 scale=1/M，正是卷积所需）

        // X(k) = w(k)·c(k)
        var outRe = [Double](repeating: 0, count: n)
        var outIm = [Double](repeating: 0, count: n)
        for k in 0..<n {
            outRe[k] = wRe[k] * c.re[k] - wIm[k] * c.im[k]
            outIm[k] = wRe[k] * c.im[k] + wIm[k] * c.re[k]
        }
        return (outRe, outIm)
    }

    // MARK: - 频率轴与移位

    /// np.fft.fftfreq(n, d)：返回长度 n 的采样频率轴（cycles/unit）。
    /// k ≤ (n−1)/2 为正半轴，Nyquist 归入负半轴（与 numpy 一致）。
    public static func fftfreq(_ n: Int, d: Double) -> [Double] {
        precondition(n >= 1 && d > 0)
        let nd = Double(n) * d
        let cutoff = (n - 1) / 2
        return (0..<n).map { k in
            k <= cutoff ? Double(k) / nd : Double(k - n) / nd
        }
    }

    /// np.fft.fftshift：零频分量移到中央，分割点 ceil(n/2)。
    /// 例 [0,1,2,3,4] → [3,4,0,1,2]；n 偶时分割点 n/2。
    public static func fftshift<T>(_ a: [T]) -> [T] {
        let h = (a.count + 1) / 2    // ceil(n/2)
        return Array(a[h...]) + Array(a[..<h])
    }

    /// np.fft.ifftshift：fftshift 的逆，分割点 floor(n/2)。
    /// 例 [0,1,2,3,4] → [2,3,4,0,1]；n 偶时与 fftshift 相同。
    public static func ifftshift<T>(_ a: [T]) -> [T] {
        let h = a.count / 2          // floor(n/2)
        return Array(a[h...]) + Array(a[..<h])
    }
}
