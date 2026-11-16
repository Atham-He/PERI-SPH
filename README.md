# PERI-SPH

基于 **SPHinXsys / Icebreaking_PD-SPH** 的 NOSB-PD 与弱可压缩 SPH 耦合开发快照。
本仓库保存当前独立液芯分支的完整库源码、构建入口和 240 μm 胶囊压缩算例。
固体 MCC/SS 芯材分支尚未合入。

## 当前模型与验证范围

- 壳体：原 NOSB-PD 体积离散、Hughes–Winget 客观应力更新；J2 严格理想塑性。
- 液芯：WCSPH；沿用冰射流虚拟流体粒子的压力耦合，自由滑移。
- 压板：原核排斥接触；修正接触邻域核与零接触偏移核的不一致。
- 时间推进：复用原显式积分模块，按流固稳定步长的较小值共同推进。
- 首轮关闭断裂，只研究破裂前 0–20% 压缩。

已通过严格理想塑性材料检查和初始无载接触检查。**承压密闭性、整段实验曲线及空间／时间收敛尚未完成验收**，不能把此快照当成已经完成验证的胶囊求解器。

## 构建

当前液芯算例支持 **Ubuntu 24.04 / Linux glibc**；Windows 请使用 WSL。
编译器需支持 C++17。安装依赖：

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake libsimbody-dev libeigen3-dev libtbb-dev libboost-all-dev python3

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target peri_sph_liquid_capsule_240um -j4
ctest --test-dir build --output-on-failure
```

CTest 只执行轻量的严格理想塑性实现检查，不会启动大型压缩计算。
上游历史测试目录未打包，`SPHINXSYS_BUILD_TESTS` 默认关闭；本仓库的检查由 `BUILD_TESTING` 控制。

## 运行

```bash
# 完整目标为 0–20%；首次先保留一段完整状态，用于检查初始化与计算成本。
python3 tools/run_liquid_capsule.py --output runs/puf240 --max-steps 40

# 审查诊断后，从同一完整状态继续；不更换材料参数。
python3 tools/run_liquid_capsule.py --output runs/puf240 --resume
```

运行器不会自动换参数或失败重试；默认使用 12 个线程，并只对该子进程设置 32 GiB 虚拟内存预算。
运行输出、原始实验 PDF、编译产物及检查点不进入版本库。

输出包括 `diagnostics.csv`、`accepted_curve.csv`、`status.json` 和 `checkpoint.bin`。
瞬态诊断不是已经收敛的实验对照点；只有通过原生平衡与模型检查后才写入 `accepted_curve.csv`。
载荷采用加速的准静态加载与保持，**不等同于按实验真实速度模拟 48 秒**；有效性需以惯性、平衡和一次数值复核为依据。

## 固定输入与来源

依据 Zhang 等（2021）的 [液芯胶囊压缩研究](https://doi.org/10.1016/j.ijmecsci.2021.106604)，选图3(d)的 PUF 壳／石蜡油芯胶囊：

| 输入 | 数值 | 来源性质 |
|---|---:|---|
| 外径／壳厚 | 240 μm／2 μm | 论文样品与均匀壳厚近似 |
| 壳弹性模量 | 3.50 GPa | 论文压缩反演值 |
| 壳屈服强度 | 40 MPa | 论文压缩反演值 |
| 泊松比 | 0.33 | 论文采用的假设 |
| 壳密度 | 1500 kg/m³ | 外部 UF 固化树脂资料参考；不是该批胶囊实测值 |
| 液体密度／体积模量 | 830 kg/m³／1.66 GPa | 论文给定 |
| WCSPH 声速 | 1414.21 m/s | 按 sqrt(K/rho) 换算 |
| 实验加载速度 | 1 μm/s | 实验条件；与加速数值时间区分 |

壳密度参考 [NETZSCH UF 材料表](https://analyzing-testing.netzsch.com/en/polymers-netzsch-com/thermosets/uf-urea-formaldehyde-resin)。
完整填充、理想球形及首轮忽略物理黏性均为建模假设。目标论文未提供壳密度、油黏度和可直接输入的断裂参数。
E 与屈服强度来自论文对同一类压缩数据的拟合，因此本算例定位为结果复现，不宣称独立预测。

[`parameters.json`](examples/liquid_capsule_240um/parameters.json) 是来源清单，**尚不是原生程序自动读取的配置文件**；当前算例的固定常量定义在源码中。
求解器采用一致单位换算：长度 μm、时间 ns、密度单位1000 kg/m³，对应应力 GPa、力 mN、能量 nJ；没有质量或材料刚度缩放。

## 计算规模

当前离散约120.5万壳粒子、87.4万液体粒子，另有72万压板采样点。
开发机器早期实测约18.5 GiB运行内存、约2秒／步；按既定加载和最低保持时长估算，整段耗时下限约125小时。
这些是该机器上的初期估算，额外收敛、硬件和编译设置会改变耗时，不是性能保证。

## 目录

```text
SPHINXsys/src/                  耦合库源码与保留的上游模块
modules/                        上游结构模拟模块
examples/liquid_capsule_240um/   独立算例、参数来源和CMake入口
tools/                          可移植的运行入口
docs/                           来源、修改说明和验证记录
LICENSES/                       保留的上游许可文本
```

## 来源与署名

作者上游：[Luwen-Zhang/Icebreaking_PD-SPH](https://github.com/Luwen-Zhang/Icebreaking_PD-SPH)，基线提交 `34562a6c7ac11cb2816236f0249bfd8e6a202ac7`。
该工作建立于 [SPHinXsys](https://www.sphinxsys.org/)；原作者和许可声明均保留。
上游说明文档原样保存在 [`docs/upstream`](docs/upstream)。

- Shi 等，JMPS 191 (2024), 105783：[冰射流 PD-SPH 方法与实验](https://doi.org/10.1016/j.jmps.2024.105783)。
- Zhang 等，IJMS 205 (2021), 106604：[液芯胶囊准静态压缩](https://doi.org/10.1016/j.ijmecsci.2021.106604)。

具体继承关系与本地改动见 [`NOTICE`](NOTICE)、[`CHANGELOG.md`](CHANGELOG.md) 和 [`docs/source_provenance.json`](docs/source_provenance.json)。
