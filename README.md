# Pipeline 

$$
q_i = \sum_j V_j[q_j] W(\vb*{x}_i - \vb*{x}_j, h) \\
\hat{q}_i = \sum_J V_J[q_J] W(\vb*{x}_i - \vb*{x}_J, H) \\
q_I = \sum_j V_j[q_j] W(\vb*{x}_I - \vb*{x}_j, h) \\
\hat{q}_I = \sum_J V_J[q_J] W(\vb*{x}_I - \vb*{x}_J, H) \\
$$

其中 $\hat{}$ 为了标记来自邻居虚粒子 $J$ 的效果，小写 $i, j$ 指示实粒子对应的分量，大写 $I, J$ 指示虚粒子对应的分量。



### 算法

1. 生成虚粒子

2. 对于实粒子而言，
   $$
   V_i = \frac{m_i}{\rho_i}, \quad \rho_i = \sum_jV_j[\rho_j]W_{ij} = \sum_j m_j W_{ij}.
   $$
   对于虚粒子而言，
   $$
   V_I = \sum_jV_j[V_j]W(\vb*{x}_{Ij}, h) = \sum_j V_j^2W(\vb*{x}_{Ij}, h).
   $$
   计算粒子体积。

3. 离散化压强投影方程 （11， 14），通过调整对角项来获得自由边界条件（4.2）

   对于虚粒子而言
   $$
   \mathcal{D}_I(\vb*{v}^*) = -\sum_j V_j(\vb*{v}_j^* - \bar{\vb*{v}}_I^*) \vdot \grad_I W_{Ij} + \kappa\frac{\max(\rho_I - \rho_0, 0)}{\rho_0},
   $$
   其中 $\bar{\vb*{v}}_I^*$ 为虚粒子上的速度（用来减少边界上的误差）
   $$
   \bar{\vb*{v}}_I^* = \frac{\sum_j{V_j\vb*{v}_j^*W_{Ij}}}{\sum_jV_jW_{Ij}}.
   $$
   后一项是为了减少密度漂移，其中的 $\kappa$ 选取与 $m/\rho_0$ 同量级，$\rho_I = \sum_j m_jW_{Ij}$。使用了
   $$
   \grad_I W_{Ij} = \frac{\vb*{x}_{Ij}}{r_{Ij}}W_{Ij}'.
   $$
   关于压强的 Laplace 算子为
   $$
   \hat{\mathcal{L}}_I(p) = \hat{\alpha}_Ip_I - \sum_J \hat{\alpha}_{IJ}p_J,
   $$
   其中
   $$
   \hat{\alpha}_I = \sum_J{\hat{\alpha}_{IJ}},\quad \hat{\alpha}_{IJ} = 2\frac{V_J W'_{IJ}}{r_{IJ} + \eta},
   $$
   $\hat{\alpha}_I = \max(\hat{\alpha}_I, \alpha_0)$，使用的 $\alpha_0$ 为模拟开始时周围邻居中均存在实粒子和虚粒子的情况下计算得出的。

   对于固体边界条件而言，在固体区域均匀地采样固体粒子，使用流体实粒子相同的质量。固体粒子携带着固体材料的速度，参与流体虚粒子的散度计算。

4. 共轭梯度法求解线性系统
   $$
   \hat{\mathcal{L}}_I(p) = \mathcal{D}_I(\vb*{v}^*).
   $$

5. 对于实粒子计算压强效果（16）
   $$
   \mathcal{G}_i(p) = \sum_J V_Jp_J\grad_iW_{iJ} + \sum_{j^s} V_j^s \operatorname{Proj}_{\vb*{n}_{j^s}}(\vb*{v}_i - \vb*{v}_{j^s})W_{ij^s},
   $$

   $$
   \vb*{v}_i^{n+1} = \vb*{v}_i^* - \mathcal{G}_i(p).
   $$

   

6. 对于实粒子更新速度和位置，删除虚粒子





# VCL-PhysX 物理模拟库

## 总览

VCL-PhysX 提供了一套面向计算机图形学及相关领域的物理模拟框架，实现了在图形学社区中广泛应用的多种数值方法，可作为在物理模拟方向做进一步科学研究的基础。

本代码库具有以下特点：

- 采用 C++20 标准，完全面向对象、平台无关的程序设计，算法逻辑清晰。
- 基于 OpenGL 4.5 标准封装了可编程图形渲染管线，实现了基础的文字渲染引擎。
- 实现了可泛用的网格、粒子数据结构及显式、隐式计算几何类。
- 提供了包括命令行解析器、物理常量及常用数学函数在内的辅助工具。

## 目录

- [VCL-PhysX 物理模拟库](#vcl-physx-物理模拟库)
  - [总览](#总览)
  - [目录](#目录)
  - [快速入门](#快速入门)
    - [编译：Windows 平台（推荐）](#编译windows-平台推荐)
    - [编译：基于 xmake 的多平台](#编译基于-xmake-的多平台)
    - [可视化工具简介](#可视化工具简介)
  - [代码架构](#代码架构)
  - [参考文献](#参考文献)

## 快速入门

### 编译：Windows 平台（推荐）

#### 需求

- 64 位 Windows 10 操作系统
- [Git][getting-started:git] （选择将核心功能集成到系统环境变量中）
- [Visual Studio 2019][getting-started:visual-studio] Ver. 16.8 或更新的版本（包含英文语言包）

#### 安装依赖库

在合适的目录下运行 Powershell 终端。依次执行下列命令：

```
> git clone https://github.com/microsoft/vcpkg
> cd vcpkg
> .\bootstrap-vcpkg.bat
> .\vcpkg integrate install
> .\vcpkg install eigen3:x64-windows-static fmt:x64-windows-static
> .\vcpkg install glad:x64-windows-static glfw3:x64-windows-static
> .\vcpkg install yaml-cpp:x64-windows-static
```

#### 编译运行

下载代码库最新版到本地：

```
> git clone https://gitlab.com/VCL-Physics/vcl-physx
```

使用 Visual Studio 打开 `.\vcl-physx\VCL-PhsyX.sln`，切换至 Release 模式。选择“编译解决方案”。待编译完成后，进入 `.\vcl-physx\Binary\Release` 目录下，运行任一 `.\*Test.exe`。

运行 `.\Viewer.exe` 可视化模拟结果。

### 编译：基于 xmake 的多平台

#### 需求

- 64 位操作系统
  - Windows 10、Ubuntu 20.04（经过测试）
  - Mac OS X、Arch Linux、Fedora 以及其他 Unix 内核 OS（未经测试）
- 支持 C++20 的编译器（GCC 10.1+、Clang 10+、MSVC 19.25+）
- [Git][getting-started:git]
- [Xmake][getting-started:xmake]

#### 编译运行

xmake 会在第一次构建时自动下载并安装相关依赖库。

自动探测配置：

```
> git clone https://gitlab.com/VCL-Physics/vcl-physx && cd vcl-physx
> xmake -y
> xmake build examples
> xmake run SpringMassSystemTest -o output
> xmake run vcl-viewer -o output
```

自定义配置（以 Linux 为例）：

```
> git clone https://gitlab.com/VCL-Physics/vcl-physx && cd vcl-physx
> xmake f -p linux -m [debug|release] --openmp=[true|false]
> xmake -y
> xmake build examples
> xmake run SpringMassSystemTest -o output
> xmake run vcl-viewer -o output
```

把编译输出的库文件和可执行文件放到指定目录下：

```
> xmake install -o dist
> xmake install -o dist examples
```

更多相关用法，以及如何将 VCL-PhysX 集成到你的项目中，请参见 [xmake 中文文档](https://xmake.io/#/zh-cn/)。

[getting-started:git]: https://git-scm.com/downloads
[getting-started:visual-studio]: https://visualstudio.microsoft.com/
[getting-started:xmake]: https://xmake.io/#/

### 可视化工具简介

Viewer 是一款轻量级的模拟结果离线可视化工具，同时也是利用 VCL-PhysX 图形渲染库进行开发的样例程序。

当模拟程序生成输出文件到特定目录 `directory` 时，以 `-o directory` 为命令行参数调用 Viewer 将加载并可视化模拟结果（若不指定，则 `directory` 默认为当前目录下的 `output`）。Viewer 可与模拟程序同时运行，二者互不干扰。

Viewer 使用下列键位及鼠标操作：

|     操作     |               功能               |
| :----------: | :------------------------------: |
| 鼠标左键拖动 | 旋转轨道相机（仅对三维场景启用） |
| 鼠标右键拖动 |           移动相机焦点           |
|   鼠标滚轮   |  径向移动相机（放大或缩小场景）  |
|     ↑←↓→     |         改变平行光源极角         |
|      P       |           播放 / 暂停            |
|     [ ]      |         上一帧 / 下一帧          |
|      0       |        显示 / 隐藏坐标轴         |
|    1 ~ 9     |         显示 / 隐藏物体          |
|   F1 ~ F12   |    （参见画面左上角文字信息）    |
|     Esc      |             退出程序             |

## 代码架构

- `VCL-PhysX`
  - `Binary`：二进制文件存放目录
    - `Debug`：调试模式输出目录
    - `Release`：发行模式输出目录
  - `Cores`：核心代码库
    - `Geometries`：几何表示与维护类
    - `Graphics`：图形渲染引擎
    - `Materials`：材质及材料模型
    - `Physics`：物理模拟引擎
    - `Solvers`：线性系统解算器
    - `Structures`：基本数据结构
    - `Utilities`：通用辅助工具
    - `Viewer`：模拟结果离线可视化程序
  - `Demos`：样例模拟程序库
    - `*Test`：各种数值模拟方法的测试程序
  - `Properties`：Visual Studio 属性表目录

## 参考文献

1. Doyub Kim. 2016. *Fluid Engine Development*. AK Peters/CRC Press, Boca Raton, FL, USA.  
   代码架构的主要参考。
2. Robert Bridson. 2015. *Fluid simulation for computer graphics* (2nd ed.). AK Peters/CRC Press, Boca Raton, FL, USA.  
   欧拉网格流体的基础框架，参见 [`Physics/EulerianFluid`](Cores/Physics/EulerianFluid.h) 和 [`Physics/LevelSetLiquid`](Cores/Physics/LevelSetLiquid.h)。
3. Dan Koschier, Jan Bender, Barbara Solenthaler, and Matthias Teschner. 2019. *Smoothed Particle Hydrodynamics for Physically-Based Simulation of Fluids and Solids*. Eurographics 2019 Tutorial.  
   光滑粒子的数据结构、插值函数等，参见 [`Structures/SmoothedParticles`](Cores/Structures/SmoothedParticles.h) 和 [`Structures/ParticlesNearbySeacher`](Cores/Structures/ParticlesNearbySearcher.h)。  
   SPH 液体的基础框架及预测—修正不可压条件解算方法，参见 [`Physics/SmthParticleHydrodLiquid`](Cores/Physics/SmthParticleHydrodLiquid.h) 和 [`Physics/PredCorrIncomprSphLiquid`](Cores/Physics/PredCorrIncomprSphLiquid.h)。
4. Stanley Osher and Ronald P Fedkiw. 2005. *Level set methods and dynamic implicit surfaces*. Springer, New York, NY, USA.  
   快速行进法（Fast marching method）重整化水平集，参见 [`Geometries/LevelSetReinitializer`](Cores/Geometries/LevelSetReinitializer.h)。  
   移动立方体法（Marching cubes method）的等值面提取，参见 [`Geometries/LevelSetContourer`](Cores/Geometries/LevelSetContourer.h)。
5. A. Ralston. 1962. Runge-Kutta methods with minimum error bounds. *Math. Comput.* 16, 80 (1962), 431–437.  
   隆格—库塔方法的系数选取，参见 [`Physics/EulerianAdvector`](Cores/Physics/EulerianAdvector.h)。
6. A. Selle, R. Fedkiw, B. Kim, Y. Liu, and J. Rossignac. 2008. An unconditionally stable MacCormack method. *J. Sci. Comput*. 35, 2-3 (June 2008), 350–371.  
   应用麦科马克差分格式对流欧拉网格流体，参见 [`Physics/EulerianAdvector`](Cores/Physics/EulerianAdvector.h)。
7. Yen Ting Ng, Chohong Min, and Frédéric Gibou. 2009. An efficient fluid–solid coupling algorithm for single-phase flows. *J. Comput. Phys*. 228, 23 (2009), 8807-8829.  
   欧拉网格流体中更精确的固体边界条件处理，参见 [`Physics/EulerianBoundaryHelper`](Cores/Physics/EulerianBoundaryHelper.h)。
8. Yongning Zhu and Robert Bridson. 2005. Animating sand as a fluid. *ACM Trans. Graph*. 24, 3 (July 2005), 965–972.  
   PIC/FLIP 方法对流液体速度场，参见 [`Physics/ParticleInCellLiquid`](Cores/Physics/ParticleInCellLiquid.h) 和 [`Physics/FlImplicitParticleLiquid`](Cores/Physics/FlImplicitParticleLiquid.h)。
9. Chenfanfu Jiang, Craig Schroeder, Andrew Selle, Joseph Teran, and Alexey Stomakhin. 2015. The affine particle-in-cell method. *ACM Trans. Graph*. 34, 4, Article 51 (August 2015), 10 pages.  
   APIC 方法对流液体速度场，参见 [`Physics/AffineParticleInCellLiquid`](Cores/Physics/AffineParticleInCellLiquid.h)。  
   APIC 方法对流速度场，并作为物质点法的基础，参见 [`Physics/MaterialPointSubstances`](Cores/Physics/MaterialPointSubstances.h)。
10. Yuanming Hu, Yu Fang, Ziheng Ge, Ziyin Qu, Yixin Zhu, Andre Pradhana, and Chenfanfu Jiang. 2018. A moving least squares material point method with displacement discontinuity and two-way rigid body coupling. *ACM Trans. Graph*. 37, 4, Article 150 (August 2018), 14 pages.  
    移动最小二乘物质点法模拟多物理场景，参见 [`Physics/MaterialPointSubstances`](Cores/Physics/MaterialPointSubstances.h)。
11. Andre Pradhana Tampubolon, Theodore Gast, Gergely Klár, Chuyuan Fu, Joseph Teran, Chenfanfu Jiang, and Ken Museth. 2017. Multi-species simulation of porous sand and water mixtures. *ACM Trans. Graph*. 36, 4, Article 105 (July 2017), 11 pages.  
    物质点法对液体的特化，参见 [`Materials/MaterialPointLiquid`](Cores/Materials/MaterialPointLiquid.h)。
