import Testing
import Foundation
@testable import QMModules
@testable import EngineKit

/// 模块注册表：W3 试点 4 + W4 十二 + W5 十二 + W6 十三 + W7 十 + W8 四 + W9 四
/// + W11A 一 + W11B 一 + W12 两 = 63 个模块（全部计划模块收官）的注册、分组、搜索、幂等性。
/// 操作共享 ModuleRegistry，须串行执行。
@Suite(.serialized)
struct ModuleLibraryTests {

    private func registered() -> ModuleRegistry {
        ModuleRegistry.shared.removeAll()
        ModuleLibrary.registerBuiltins()
        return ModuleRegistry.shared
    }

    @Test("内置 63 模块注册且重复注册幂等（SDK 计算库收官）")
    func builtinRegistration() {
        let registry = registered()
        #expect(registry.counts.total == 63)
        #expect(registry.counts.realtime == 40,
                "试点 3 + W4 七个 + W5 十个 + W6 十一个 + W7 七个 + W12 两个实时档（W8/W9 秒级档 + W11A/W11B 策略档）")
        #expect(registry.counts.basic == 44, "含 W12 Bloch 球（实时基础档）")
        #expect(registry.counts.advanced == 19, "含 W11A/W11B 策略档与 W12 误用辨析（专题档）")
        // 重复注册全部被拒
        #expect(ModuleLibrary.registerBuiltins() == 0)
        #expect(registry.counts.total == 63)
    }

    @Test("分组：十个分类按固定序，模块按笔记号排序")
    func grouping() {
        let grouped = registered().grouped
        #expect(grouped.map(\.category) == [
            .oldQuantum, .quantumStates, .schrodinger1D,
            .angularCentral, .formalTheory, .manyBody,
            .quantumInfo, .quantumOptics, .relativisticQFT, .gemology,
        ])
        // 旧量子论：W4 五模块 + 试点黑体/光电，按笔记号 1-6
        #expect(grouped[0].modules.map(\.meta.noteNumber) == [1, 2, 3, 4, 5, 6])
        // 量子态：物质波包(7)/不确定性(8)/双缝(9)/测量理论(10)/Bloch 球(11)
        #expect(grouped[1].modules.map(\.meta.noteNumber) == [7, 8, 9, 10, 11])
        #expect(grouped[1].modules.last?.meta.id == "希尔伯特空间与狄拉克符号_Bloch球")
        // 薛定谔方程与一维问题：W8 有限差分势阱(13) + 方势垒与有限深势阱(14) + 谐振子(15)
        #expect(grouped[2].modules.map(\.meta.noteNumber) == [13, 14, 15])
        #expect(grouped[2].modules.first?.meta.id == "薛定谔方程_有限差分势阱")
        // 角动量与中心力场：两自旋 + 斯特恩-盖拉赫（同笔记 16，按标题排序）+ 氢原子径向(17)
        #expect(grouped[3].modules.map(\.meta.noteNumber) == [16, 16, 17])
        #expect(grouped[3].modules.map(\.meta.id) == [
            "角动量与自旋_两自旋单态三重态", "角动量与自旋_斯特恩-盖拉赫分裂",
            "中心力场与氢原子_径向方程",
        ])
        // 形式理论：算符表象(12) + 笔记 18 两模块（同号按标题排序）+ 微扰变分(19) + 路径积分(20)
        #expect(grouped[4].modules.map(\.meta.noteNumber) == [12, 18, 18, 19, 20])
        #expect(grouped[4].modules.map(\.meta.id) == [
            "算符表象_高斯波包傅里叶对偶",
            "绘景变换与密度矩阵_两能级Rabi", "绘景变换与密度矩阵_退相干约化密度矩阵",
            "微扰论与变分法_方势阱与氦变分",
            "路径积分_传播子与谐振子基态",
        ])
        // 多体与凝聚态：W5 笔记 21-28 连续八模块
        #expect(grouped[5].modules.map(\.meta.noteNumber) == [21, 22, 23, 24, 25, 26, 27, 28])
        // 量子信息：纠缠度量(29) + EPR(30)/退相干动力学(31)/BB84(32)/传态×2(33)/量子计算×2(34)
        #expect(grouped[6].modules.map(\.meta.noteNumber) == [29, 30, 31, 32, 33, 33, 34, 34])
        #expect(grouped[6].modules.map(\.meta.id) == [
            "量子纠缠_纠缠度量与可分性",
            "EPR佯谬与Bell不等式_CHSH",
            "量子测量与退相干_退相干动力学",
            "量子密钥分发_BB84密钥率",
            "量子隐形传态_保真度",
            "量子隐形传态_密集编码容量",
            "量子计算_Grover搜索",
            "量子计算_QFT",
        ])
        // 量子光学：光场量子化×2(35) + 腔 QED×2(36)
        #expect(grouped[7].modules.map(\.meta.noteNumber) == [35, 35, 36, 36])
        #expect(grouped[7].modules.map(\.meta.id) == [
            "光场量子化_二阶相干度",
            "光场量子化_相干态与光子统计",
            "腔量子电动力学_JaynesCummings",
            "腔量子电动力学_量子比特Rabi",
        ])
        // 相对论量子与场论：KG 色散(37)/KG 密度振荡(37)/Dirac 精细结构(38)/Dirac 能谱(38)/
        // 狄拉克海×2(39)/QFT 图像×2(40)/量子引力×2(41)
        #expect(grouped[8].modules.map(\.meta.noteNumber) == [37, 37, 38, 38, 39, 39, 40, 40, 41, 41])
        #expect(grouped[8].modules.map(\.meta.id) == [
            "克莱因-戈登方程_色散关系",
            "克莱因-戈登方程_概率密度振荡",
            "狄拉克方程_氢原子精细结构",
            "狄拉克方程_能谱",
            "反粒子与狄拉克海_对产生阈值",
            "反粒子与狄拉克海_狄拉克海洞",
            "量子场论的基本图像_Casimir",
            "量子场论的基本图像_模式量子化",
            "量子引力与全息原理_霍金温度",
            "量子引力与全息原理_黑洞熵",
        ])
        // 宝石学量子专题（W7 新增第 10 类）：晶体场(42)×2 / 致色(43)×2 / 色心(44)×2 /
        // 光谱(45)×2 / W9 KK 反演(46) + W11B 两带光学模型(46) + W12 误用辨析(47)
        #expect(grouped[9].modules.map(\.meta.noteNumber)
                == [42, 42, 43, 43, 44, 44, 45, 45, 46, 46, 47])
        #expect(grouped[9].modules.map(\.meta.id) == [
            "晶体场与配位场理论_TanabeSugano",
            "晶体场与配位场理论_模型",
            "过渡金属离子致色的量子机制_振动耦合",
            "过渡金属离子致色的量子机制_模型",
            "色心与晶格缺陷的量子描述_F心类氢模型",
            "色心与晶格缺陷的量子描述_电子声子耦合吸收谱",
            "光谱学仪器的量子基础_拉曼谱",
            "光谱学仪器的量子基础_荧光与寿命",
            "矿物光学性质的第一性原理计算_KramersKronig",
            "矿物光学性质的第一性原理计算_模型",
            "宝石学中的量子误用辨析_模型",
        ])
    }

    @Test("搜索：标题 / 关键词 / 笔记编号命中")
    func search() {
        let registry = registered()
        #expect(registry.search("黑体").count == 1)
        #expect(registry.search("紫外灾难").first?.meta.id == "黑体辐射_三律对比")
        #expect(registry.search("Hermite").first?.meta.id == "量子谐振子_有限差分与相干态")
        // Mandel 命中两处：HOM 干涉（关键词 Hong-Ou-Mandel，多体分类序在前）与光场 Mandel Q
        #expect(registry.search("Mandel").count == 2)
        #expect(registry.search("Mandel").first?.meta.id == "全同粒子_费米气体与HOM")
        #expect(registry.search("Mandel").contains { $0.meta.id == "光场量子化_相干态与光子统计" })
        // W4 模块关键词命中
        #expect(registry.search("双缝").contains { $0.meta.id == "量子态与叠加原理_双缝干涉" })
        #expect(registry.search("芝诺").first?.meta.id == "测量理论_指针耦合与芝诺")
        #expect(registry.search("斯特恩").first?.meta.id == "角动量与自旋_斯特恩-盖拉赫分裂")
        // Rabi 命中三处：两能级(18) + JC/量子比特 Rabi(36)；形式理论分类序在前
        #expect(registry.search("Rabi").first?.meta.id == "绘景变换与密度矩阵_两能级Rabi")
        #expect(registry.search("Rabi").count == 3)
        // W5 模块关键词命中
        #expect(registry.search("HOM").first?.meta.id == "全同粒子_费米气体与HOM")
        #expect(registry.search("Slater").first?.meta.id == "多电子原子_Zeff与电离能")
        #expect(registry.search("Hückel").first?.meta.id == "分子结构_H2与苯")
        #expect(registry.search("Shockley").first?.meta.id == "半导体_载流子与pn结")
        #expect(registry.search("BEC").first?.meta.id == "量子统计_BEC与费米气体")
        #expect(registry.search("约瑟夫森").first?.meta.id == "超导_BCS与约瑟夫森")
        #expect(registry.search("Landau").first?.meta.id == "量子霍尔_Landau能级与平台")
        #expect(registry.search("CHSH").first?.meta.id == "EPR佯谬与Bell不等式_CHSH")
        #expect(registry.search("BB84").first?.meta.id == "量子密钥分发_BB84密钥率")
        #expect(registry.search("Caldeira").first?.meta.id == "量子测量与退相干_退相干动力学")
        #expect(registry.search("范霍夫").first?.meta.id == "晶体中的电子_能带与KronigPenney")
        #expect(registry.search("纯度").first?.meta.id == "绘景变换与密度矩阵_退相干约化密度矩阵")
        // W6 模块关键词命中
        #expect(registry.search("Grover").first?.meta.id == "量子计算_Grover搜索")
        #expect(registry.search("Holevo").first?.meta.id == "量子隐形传态_密集编码容量")
        #expect(registry.search("Werner").first?.meta.id == "量子隐形传态_保真度")
        #expect(registry.search("Jaynes").first?.meta.id == "腔量子电动力学_JaynesCummings")
        #expect(registry.search("缀饰态").first?.meta.id == "腔量子电动力学_JaynesCummings")
        #expect(registry.search("transmon").first?.meta.id == "腔量子电动力学_量子比特Rabi")
        #expect(registry.search("Casimir").first?.meta.id == "量子场论的基本图像_Casimir")
        // 搜索结果按目录序（分类→笔记号→标题）而非相关度排序：
        // "零点能" 命中 Casimir（真空零点能）与模式量子化（零点能）——同笔记 40，标题序 Casimir 在前
        #expect(registry.search("零点能").count == 2)
        #expect(registry.search("零点能").first?.meta.id == "量子场论的基本图像_Casimir")
        #expect(registry.search("零点能").contains { $0.meta.id == "量子场论的基本图像_模式量子化" })
        #expect(registry.search("狄拉克海").first?.meta.id == "反粒子与狄拉克海_狄拉克海洞")
        #expect(registry.search("正电子").first?.meta.id == "反粒子与狄拉克海_狄拉克海洞")
        // "Wigner" 命中两处：退相干动力学（猫态 Wigner，形式理论分类序在前）与二阶相干度（Wigner 函数）
        #expect(registry.search("Wigner").count == 2)
        #expect(registry.search("Wigner").first?.meta.id == "量子测量与退相干_退相干动力学")
        #expect(registry.search("Wigner").contains { $0.meta.id == "光场量子化_二阶相干度" })
        #expect(registry.search("反聚束").first?.meta.id == "光场量子化_二阶相干度")
        #expect(registry.search("精细结构").first?.meta.id == "狄拉克方程_氢原子精细结构")
        #expect(registry.search("Lamb").first?.meta.id == "狄拉克方程_氢原子精细结构")
        // W7 宝石学批关键词命中（分类序：relativisticQFT 在 gemology 前）
        #expect(registry.search("霍金").first?.meta.id == "量子引力与全息原理_霍金温度")
        #expect(registry.search("黑洞熵").first?.meta.id == "量子引力与全息原理_黑洞熵")
        #expect(registry.search("Tanabe").first?.meta.id == "晶体场与配位场理论_TanabeSugano")
        #expect(registry.search("晶体场").first?.meta.id == "晶体场与配位场理论_TanabeSugano")
        #expect(registry.search("Laporte").first?.meta.id == "过渡金属离子致色的量子机制_振动耦合")
        #expect(registry.search("F心").first?.meta.id == "色心与晶格缺陷的量子描述_F心类氢模型")
        #expect(registry.search("Huang-Rhys").first?.meta.id == "色心与晶格缺陷的量子描述_电子声子耦合吸收谱")
        #expect(registry.search("拉曼").first?.meta.id == "光谱学仪器的量子基础_拉曼谱")
        #expect(registry.search("荧光").first?.meta.id == "光谱学仪器的量子基础_荧光与寿命")
        #expect(registry.search("宝石").count == 11, "分类名「宝石学量子专题」全量命中（W12 加一）")
        // W8 秒级档线代批关键词命中
        #expect(registry.search("有限差分").count == 2, "FD 势阱(13) + 谐振子(15) 关键词")
        #expect(registry.search("隧穿").first?.meta.id == "一维势场_方势垒与有限深势阱")
        #expect(registry.search("里德伯").contains { $0.meta.id == "中心力场与氢原子_径向方程" })
        #expect(registry.search("变分").first?.meta.id == "微扰论与变分法_方势阱与氦变分")
        #expect(registry.search("上界定理").first?.meta.id == "微扰论与变分法_方势阱与氦变分")
        // W9 秒级档线代批 2 关键词命中
        #expect(registry.search("并发度").first?.meta.id == "量子纠缠_纠缠度量与可分性")
        #expect(registry.search("PPT").first?.meta.id == "量子纠缠_纠缠度量与可分性")
        #expect(registry.search("概率密度振荡").first?.meta.id == "克莱因-戈登方程_概率密度振荡")
        #expect(registry.search("质量隙").first?.meta.id == "狄拉克方程_能谱")
        #expect(registry.search("Kramers").first?.meta.id
                == "矿物光学性质的第一性原理计算_KramersKronig")
        // "希尔伯特" 双命中：KK 反演关键词 + W12 Bloch 模块（搜索按目录序：
        // category.rawValue「宝石学量子专题」<「量子态与数学结构」→ 宝石学在前）
        #expect(registry.search("希尔伯特").count == 2)
        #expect(registry.search("希尔伯特").first?.meta.id
                == "矿物光学性质的第一性原理计算_KramersKronig")
        #expect(registry.search("希尔伯特").contains {
            $0.meta.id == "希尔伯特空间与狄拉克符号_Bloch球"
        })
        // W11A 策略档关键词命中
        #expect(registry.search("传播子").first?.meta.id == "路径积分_传播子与谐振子基态")
        #expect(registry.search("Mehler").first?.meta.id == "路径积分_传播子与谐振子基态")
        #expect(registry.search("虚时间").first?.meta.id == "路径积分_传播子与谐振子基态")
        #expect(registry.search("谱投影").first?.meta.id == "路径积分_传播子与谐振子基态")
        // W11B 策略档关键词命中（"Kramers" 双命中：KK 反演在前、两带模型在后）
        #expect(registry.search("两带").first?.meta.id == "矿物光学性质的第一性原理计算_模型")
        #expect(registry.search("联合态密度").first?.meta.id == "矿物光学性质的第一性原理计算_模型")
        #expect(registry.search("JDOS").first?.meta.id == "矿物光学性质的第一性原理计算_模型")
        #expect(registry.search("吸收边").first?.meta.id == "矿物光学性质的第一性原理计算_模型")
        #expect(registry.search("Kramers").count == 2)
        #expect(registry.search("Kramers").first?.meta.id
                == "矿物光学性质的第一性原理计算_KramersKronig")
        #expect(registry.search("Kramers").contains {
            $0.meta.id == "矿物光学性质的第一性原理计算_模型"
        })
        // W12 收官批关键词命中
        // "Bloch" 双命中：误用辨析(47) + Bloch 球(11)；宝石学分类串序在前
        #expect(registry.search("Bloch").count == 2)
        #expect(registry.search("Bloch").first?.meta.id
                == "宝石学中的量子误用辨析_模型")
        #expect(registry.search("Bloch").contains {
            $0.meta.id == "希尔伯特空间与狄拉克符号_Bloch球"
        })
        // "红宝石" 双命中：致色模型(43) + 误用辨析(47)，同分类按笔记号排序
        #expect(registry.search("红宝石").count == 2)
        #expect(registry.search("红宝石").first?.meta.id
                == "过渡金属离子致色的量子机制_模型")
        #expect(registry.search("混态").first?.meta.id == "宝石学中的量子误用辨析_模型")
        #expect(registry.search("694.3").first?.meta.id == "宝石学中的量子误用辨析_模型")
        // 笔记编号 "15" / "3" / "24" / "36" 均可命中
        #expect(registry.search("15").first?.meta.noteNumber == 15)
        #expect(registry.search("03").first?.meta.noteNumber == 3)
        #expect(registry.search("24").first?.meta.noteNumber == 24)
        #expect(registry.search("36").first?.meta.noteNumber == 36)
        // "11" 为子串会误命中其它模块（如 BB84 的 "B92" 类关键词），用 contains 断言
        #expect(registry.search("11").contains { $0.meta.noteNumber == 11 })
        #expect(registry.search("47").first?.meta.noteNumber == 47)
        // id 深链定位
        #expect(registry.module(id: "光电效应_截止电压线性")?.meta.noteNumber == 4)
        #expect(registry.module(id: "物质波包_高斯波包演化")?.meta.noteNumber == 7)
        #expect(registry.module(id: "量子密钥分发_BB84密钥率")?.meta.noteNumber == 32)
        #expect(registry.module(id: "量子计算_QFT")?.meta.noteNumber == 34)
        #expect(registry.module(id: "光场量子化_二阶相干度")?.meta.noteNumber == 35)
        #expect(registry.module(id: "量子场论的基本图像_Casimir")?.meta.noteNumber == 40)
        #expect(registry.module(id: "希尔伯特空间与狄拉克符号_Bloch球")?.meta.noteNumber == 11)
        #expect(registry.module(id: "宝石学中的量子误用辨析_模型")?.meta.noteNumber == 47)
        #expect(registry.module(id: "不存在的模块") == nil)
    }
}
