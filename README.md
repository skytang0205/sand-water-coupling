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
  - [代码架构](#代码架构)
  - [参考文献](#参考文献)

## 快速入门

### 编译：Windows 平台（推荐）

#### 需求

- 64 位 Windows 10 操作系统
- [Git][getting-started:git] （选择将核心功能集成到系统环境变量中）
- [Visual Studio 2019 Ver. 16.8][getting-started:visual-studio] 或更新的版本（包含英文语言包）
- [Vcpkg][getting-started:vcpkg]

#### 安装依赖库

在合适的目录下运行 Powershell 终端。依次执行下列命令：

``` powershell
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

``` powershell
> git clone https://gitlab.com/VCL-Physics/vcl-physx
```

使用 Visual Studio 打开 `.\vcl-physx\VCL-PhsyX.sln`，切换至 Release 模式。选择“编译解决方案”。待编译完成后，进入 `.\vcl-physx\Binary\Release` 目录下，运行任一 `.\*Test.exe`。

运行 `.\Viewer.exe` 可视化模拟结果。

### 编译：基于 xmake 的多平台

[getting-started:git]: https://git-scm.com/downloads
[getting-started:visual-studio]: https://visualstudio.microsoft.com/
[getting-started:vcpkg]: https://github.com/microsoft/vcpkg

## 代码架构

## 参考文献

1. Doyub Kim. 2016. *Fluid Engine Development*. AK Peters/CRC Press, Boca Raton, FL, USA.
2. Robert Bridson. 2015. *Fluid simulation for computer graphics* (2nd ed.). AK Peters/CRC Press, Boca Raton, FL, USA.
