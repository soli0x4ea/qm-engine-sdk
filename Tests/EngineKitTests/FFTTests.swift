import Testing
import Foundation
@testable import EngineKit

/// W4 FFT 数据层：vDSP 封装与 numpy.fft 语义对拍。
@Suite("W4 FFT 数据层")
struct FFTTests {

    private static let tol = 1e-9   // FFT 与 O(N²) 朴素求和的次序容差（远严于 §6.1 二档 1e-6）

    /// 参考：直接按定义 O(N²) 计算 DFT（X(k)=Σₙ xₙ e^(−2πikn/N)，不归一化）。
    private func dftReference(re: [Double], im: [Double]) -> (re: [Double], im: [Double]) {
        let n = re.count
        var outRe = [Double](repeating: 0, count: n)
        var outIm = [Double](repeating: 0, count: n)
        for k in 0..<n {
            var sr = 0.0, si = 0.0
            for t in 0..<n {
                let ang = -2 * Double.pi * Double(k) * Double(t) / Double(n)
                let c = cos(ang), s = sin(ang)
                sr += re[t] * c - im[t] * s
                si += re[t] * s + im[t] * c
            }
            outRe[k] = sr; outIm[k] = si
        }
        return (outRe, outIm)
    }

    @Test("fft 与定义式 DFT 逐位一致（非 2 幂 + 2 幂）")
    func fftMatchesDefinition() {
        for n in [4, 8, 12, 64, 96, 256] {
            let re = (0..<n).map { sin(Double($0) * 0.7) + Double($0 % 3) * 0.25 }
            let im = (0..<n).map { cos(Double($0) * 0.3) * 0.5 }
            let got = FFT.fft(real: re, imag: im)
            let want = dftReference(re: re, im: im)
            for k in 0..<n {
                #expect(abs(got.re[k] - want.re[k]) < Self.tol * max(1, abs(want.re[k])),
                        "n=\(n) re[\(k)]: \(got.re[k]) vs \(want.re[k])")
                #expect(abs(got.im[k] - want.im[k]) < Self.tol * max(1, abs(want.im[k])),
                        "n=\(n) im[\(k)]: \(got.im[k]) vs \(want.im[k])")
            }
        }
    }

    @Test("ifft(fft(x)) == x（1/N 归一化位置正确）")
    func roundtrip() {
        let n = 128
        let re = (0..<n).map { sin(Double($0) * 1.1) }
        let im = (0..<n).map { cos(Double($0) * 0.9) * 0.3 }
        let f = FFT.fft(real: re, imag: im)
        let back = FFT.ifft(real: f.re, imag: f.im)
        for k in 0..<n {
            #expect(abs(back.re[k] - re[k]) < 1e-12)
            #expect(abs(back.im[k] - im[k]) < 1e-12)
        }
    }

    @Test("直流分量：全 1 序列 fft → X(0)=N 其余≈0（numpy 语义）")
    func dcComponent() {
        let n = 32
        let got = FFT.fft([Double](repeating: 1, count: n))
        #expect(abs(got.re[0] - Double(n)) < 1e-12)
        for k in 1..<n {
            #expect(abs(got.re[k]) < 1e-10 && abs(got.im[k]) < 1e-10)
        }
    }

    @Test("fftfreq 与 numpy 语义一致")
    func freqAxis() {
        // np.fft.fftfreq(8, d=0.1) = [0, 1.25, 2.5, 3.75, -5, -3.75, -2.5, -1.25]
        let f = FFT.fftfreq(8, d: 0.1)
        let want = [0.0, 1.25, 2.5, 3.75, -5.0, -3.75, -2.5, -1.25]
        for k in 0..<8 { #expect(abs(f[k] - want[k]) < 1e-15, "k=\(k)") }
        // 奇数长度：np.fft.fftfreq(5, d=1) = [0, 0.2, 0.4, -0.4, -0.2]
        let f5 = FFT.fftfreq(5, d: 1)
        let want5 = [0.0, 0.2, 0.4, -0.4, -0.2]
        for k in 0..<5 { #expect(abs(f5[k] - want5[k]) < 1e-15, "k=\(k)") }
    }

    @Test("fftshift / ifftshift 奇偶长度（numpy 参照值）")
    func shifts() {
        // np.fft.fftshift([0,1,2,3]) = [2,3,0,1]；ifftshift 同（偶数）
        #expect(FFT.fftshift([0, 1, 2, 3]) == [2, 3, 0, 1])
        #expect(FFT.ifftshift([0, 1, 2, 3]) == [2, 3, 0, 1])
        // np.fft.fftshift([0,1,2,3,4]) = [3,4,0,1,2]；ifftshift = [2,3,4,0,1]
        #expect(FFT.fftshift([0, 1, 2, 3, 4]) == [3, 4, 0, 1, 2])
        #expect(FFT.ifftshift([0, 1, 2, 3, 4]) == [2, 3, 4, 0, 1])
        // 互逆
        let a = (0..<17).map { $0 }
        #expect(FFT.ifftshift(FFT.fftshift(a)) == a)
        #expect(FFT.fftshift(FFT.ifftshift(a)) == a)
    }

    @Test("4096 点高斯 FFT 性能 < 5 ms（实时档余量验证）")
    func performance4096() {
        let n = 4096
        let x = (0..<n).map { exp(-pow((Double($0) - 2048) / 200, 2)) }
        let t0 = ContinuousClock.now
        for _ in 0..<20 { _ = FFT.fft(x) }
        let ms = Double((ContinuousClock.now - t0).components.attoseconds) / 1e15
            + Double((ContinuousClock.now - t0).components.seconds) * 1000
        #expect(ms / 20 < 5, "单次 \(ms / 20) ms")
    }
}
