import Foundation
import EngineKit

/// 内置真值模块库：63 个 Python 脚本的移植模块统一在此注册。
/// SOP：新模块移植完成后在此追加一行；App 侧零改动（目录页/搜索/路由由注册表驱动）。
public enum ModuleLibrary {

    /// 注册全部内置模块；返回实际注册数（重复 id 计 0，便于单测防漏防重）。
    @discardableResult
    public static func registerBuiltins() -> Int {
        let modules: [any SimModule] = [
            // W1-W3 试点
            BlackbodyModule(),        // 笔记 03 · 实时档
            PhotoelectricModule(),    // 笔记 04 · 实时档
            PhotonStatisticsModule(), // 笔记 35 · 实时档
            HarmonicOscillatorModule(), // 笔记 15 · 秒级档（eigh 通道首用）
            // W4 早期量子论与量子态批（12 模块，fixtures 全对拍）
            QuantumBehaviorModule(),   // 笔记 01 · 实时档 · 三实验对照
            WaveParticleModule(),      // 笔记 02 · 实时档 · 波长与累积模拟
            ComptonModule(),           // 笔记 05 · 实时档 · 位移与反冲
            BalmerModule(),            // 笔记 06 · 实时档 · 巴尔末系
            WavePacketModule(),        // 笔记 07 · 实时档 · 高斯波包演化
            DoubleSlitModule(),        // 笔记 09 · 实时档 · 双缝干涉
            MeasurementTheoryModule(), // 笔记 10 · 实时档 · 指针耦合与芝诺
            UncertaintyModule(),       // 笔记 08 · 秒级档 · 等号与势阱（FFT 通道首用）
            OperatorRepresentationModule(), // 笔记 12 · 秒级档 · 傅里叶对偶（4096 点 FFT）
            TwoSpinModule(),           // 笔记 16 · 秒级档 · 单态三重态（热图首用）
            SternGerlachModule(),      // 笔记 16 · 秒级档 · 空间量子化
            RabiModule(),              // 笔记 18 · 秒级档 · 绘景等价（RK4）
            // W5 凝聚态与量子信息批（12 模块，fixtures 全对拍）
            DecoherenceRDMModule(),        // 笔记 18 · 实时档 · 退相干约化密度矩阵
            IdenticalParticlesModule(),    // 笔记 21 · 实时档 · 费米气体与 HOM
            MultiElectronAtomsModule(),    // 笔记 22 · 实时档 · Slater 屏蔽与电离能（NIST 对照）
            MolecularStructureModule(),    // 笔记 23 · 实时档 · H₂ Morse 与苯 Hückel
            KronigPenneyModule(),          // 笔记 24 · 秒级档 · 能带与态密度（20000 点扫描）
            SemiconductorModule(),          // 笔记 25 · 实时档 · 载流子与 pn 结
            QuantumStatisticsModule(),      // 笔记 26 · 实时档 · BEC 与费米简并压
            SuperconductivityModule(),     // 笔记 27 · 实时档 · BCS 能隙与约瑟夫森
            QuantumHallModule(),            // 笔记 28 · 实时档 · Landau 能级与平台（双 Y 轴首用）
            EPRBellModule(),               // 笔记 30 · 秒级档 · CHSH 与 Tsirelson 界（定种子采样）
            DecoherenceDynamicsModule(),   // 笔记 31 · 实时档 · 退相干时间尺度（三子图）
            BB84Module(),                  // 笔记 32 · 实时档 · BB84 密钥率（安全区填充首用）
            // W6 量子光学与相对论量子批（13 模块，fixtures 全对拍）
            TeleportationFidelityModule(),  // 笔记 33 · 实时档 · Werner 传态保真度
            DenseCodingModule(),            // 笔记 33 · 实时档 · 密集编码容量（Holevo 界）
            GroverModule(),                 // 笔记 34 · 实时档 · Grover sin² 律与最优迭代
            QFTModule(),                    // 笔记 34 · 秒级档 · 8×8 酉矩阵（复数通道首用）
            JaynesCummingsModule(),         // 笔记 36 · 实时档 · 真空 Rabi 振荡与缀饰态分裂
            QubitRabiModule(),              // 笔记 36 · 实时档 · 三平台退相干包络
            KleinGordonDispersionModule(),  // 笔记 37 · 实时档 · 正负能支与质量阈值
            DiracHydrogenModule(),          // 笔记 38 · 实时档 · 精细结构与薛定谔极限
            PairProductionModule(),         // 笔记 39 · 实时档 · 2mc² 阈值与反冲修正（示意图首用）
            DiracSeaHoleModule(),           // 笔记 39 · 实时档 · 狄拉克海洞 = 正电子
            CasimirModule(),               // 笔记 40 · 实时档 · d⁻⁴ 幂律与局部指数
            ModeQuantizationModule(),       // 笔记 40 · 实时档 · 模式量子化与零点能
            SecondOrderCoherenceModule(),  // 笔记 35 · 秒级档 · g⁽²⁾ 三态 + 压缩真空 Wigner（等高线首用）
            // W7 宝石学批 & 实时档收官（10 模块，fixtures 全对拍，新增 MaterialDB 数据层）
            HawkingTemperatureModule(),     // 笔记 41 · 实时档 · 霍金温度（面积律 + 普朗克单位）
            BlackHoleEntropyModule(),       // 笔记 41 · 实时档 · 黑洞熵（面积律 + 普朗克单位）
            CrystalFieldModelModule(),     // 笔记 42 · 实时档 · 晶体场高/低自旋交叉
            CrystalFieldTanabeSuganoModule(), // 笔记 42 · 秒级档 · Tanabe-Sugano d² 能级图（600 点）
            ColorVibrationalCouplingModule(), // 笔记 43 · 实时档 · 致色振动耦合饱和曲线
            ColorMechanismModelModule(),   // 笔记 43 · 实时档 · 致色模型 λ=10⁷/Δ
            FCenterHydrogenModule(),        // 笔记 44 · 实时档 · F 心类氢模型六卤化物
            ElectronPhononAbsorptionModule(), // 笔记 44 · 秒级档 · 电子-声子吸收谱（Huang-Rhys）
            FluorescenceLifetimeModule(),   // 笔记 45 · 实时档 · 荧光与寿命（A 系数 + 双指数 + Jablonski）
            RamanSpectrumModule(),         // 笔记 45 · 秒级档 · 拉曼谱（60001 点 Stokes/anti-Stokes）
            // W8 秒级档线代批 1（4 模块，fixtures 全对拍；AlgebraCore + 秒级调度通道首用）
            SchrodingerFDWellModule(),     // 笔记 13 · 秒级档 · 有限差分势阱（800² 部分谱首用）
            BarrierWellModule(),           // 笔记 14 · 秒级档 · 方势垒与有限深势阱（二分求根）
            HydrogenRadialModule(),        // 笔记 17 · 秒级档 · 氢原子径向（3000² 部分谱 ×3）
            PerturbationVariationModule(), // 笔记 19 · 秒级档 · 微扰展开与氦变分
            // W9 秒级档线代批 2（4 模块，fixtures 全对拍；eighHermitian + 采样双模式 + 帧栈/降采样首用）
            EntanglementModule(),          // 笔记 29 · 秒级档 · 纠缠度量与可分性（zheevr + 采样双模式首用）
            KleinGordonDensityModule(),    // 笔记 37 · 秒级档 · KG 概率密度振荡（400×600 降采样首用）
            DiracSpectrumModule(),         // 笔记 38 · 秒级档 · 狄拉克能谱（400 点 dsyevr 扫描）
            KramersKronigModule(),         // 笔记 46 · 秒级档 · KK 反演（12000² vDSP）
            // W11A 策略档批（fixtures 对拍；StrategyCore 降维网格/收敛判定/参数哈希缓存 + 策略通道首用）
            PathIntegralModule(),          // 笔记 20 · 策略档 · 路径积分传播子与谐振子基态（谱投影）
            // W11B 策略档批（fixtures 对拍；复用 StrategyCore/StrategyChannel + W9 KK 反演核）
            TwoBandOpticalModule(),        // 笔记 46 · 策略档 · 两带光学模型（JDOS→ε₂→KK，缓存/检查点续算）
            // W12 收官批（3D 对子：Bloch 数据模型交付，SceneKit 渲染属剩余工作；MaterialDB v2 核验）
            BlochSphereModule(),           // 笔记 11 · 实时档 · Bloch 球两能级纯态参数化
            GemologyMisuseModelModule(),   // 笔记 47 · 实时档 · 量子误用辨析（双 Bloch 对照+红宝石能级，MaterialDB ruby 行）
        ]
        var registered = 0
        for m in modules where ModuleRegistry.shared.register(m) {
            registered += 1
        }
        return registered
    }
}
