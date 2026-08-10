# MATLAB-ISAC

一个用于 ISAC（Integrated Sensing and Communication）学习和仿真的 MATLAB 小型库。

作者本人是 ISAC 初学者，目前一边学习、一边把常用的模型和算法整理进来。代码会随着学习推进持续补充和调整，也欢迎其他感兴趣的同学使用、交流和提出建议。

## 这里有什么

核心代码位于 `src/+isac/`，按功能分为以下模块：

- `+array`：阵列模型。目前包含 `ULA`，可生成导向向量。
- `+waveform`：波形模型。目前包含 `OFDM` 参数配置和 QPSK 符号生成。
- `+scene`：场景与目标模型。目前包含 `PointTarget`，描述目标的角度、距离、速度和反射系数。
- `+channel`：信道模型。目前包含远场单基地回波信道 `FarFieldMonostatic`。
- `+simulation`：仿真接口。目前包含 `MonostaticSimulator`，用于连接阵列、波形、目标和信道并生成接收数据。
- `+estimation`：参数估计。目前包含 Periodogram 的 AoA、距离和多普勒估计。

