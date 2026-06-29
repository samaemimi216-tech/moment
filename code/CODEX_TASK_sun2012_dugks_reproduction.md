# CODEX 任务说明书：将 `3D.cpp` 修正为 Wang2015-DUGKS 迁移至 Sun2012 三维瞬态导热-辐射耦合复现代码

> 本文件供 Codex / 代码修改代理直接读取。请严格按本文档修改现有 C++ 代码，目标是让当前 `3D.cpp` 能够正确复现 Sun2012 三维瞬态导热-辐射耦合基准，并保持 Wang2015 / `dugks_cavity.c` 风格的 DUGKS 有限体积-特征线更新逻辑。

---

## 0. 输入文件与目标文件

### 0.1 当前可用文件

项目中至少包含以下文件：

```text
3D.cpp
dugks_cavity(4).c
Wang2015 - A Coupled Discrete Unified Gas Kinetic Scheme for Boussinesq Flows(13).pdf
sun2012(2).pdf
A_coupled_DUGKS_for_transient_conduction_radiation_in_participating_media(2).pdf
```

### 0.2 主要修改目标

优先修改：

```text
3D.cpp
```

必要时可拆分为多个文件，但若没有明确需求，请优先保持单文件结构，降低引入新错误的风险。

### 0.3 最终目标

将当前代码修正为：

```text
D3Q27 温度场 DUGKS
+
三维角度离散 RTE 伪时间 DUGKS
+
Sun2012 hot-wall 变量体系
+
每个物理时间步内辐射场收敛
+
可输出 Sun2012 Table 1 与 Fig.3-Fig.7 所需数据
```

---

## 1. 非谈判性方程契约

修改代码时必须保证以下四个关系同时成立。

### 1.1 温度方程

代码使用 hot-wall 归一化温度：

\[
\theta=\frac{T}{T_S}
\]

因此 Sun2012 的冷壁和初始温度为：

```cpp
TE = 0.5;
```

南壁热壁为：

```cpp
TW = 1.0;
```

迁移后的能量方程必须为：

\[
\frac{\partial \theta}{\partial t^*}
=
\frac{1}{\tau_L^2}\nabla^2\theta
+
\frac{8(1-\omega)}{N_{cr}}
\left(J-\theta^4\right)
\]

其中 \(t^*\) 是 Sun2012 的无量纲时间。

### 1.2 温度 DUGKS 热松弛时间

DUGKS 恢复的热扩散系数应为：

\[
\alpha^*=\tau_g R_T=\frac{1}{\tau_L^2}
\]

因此代码必须实现：

```cpp
inline double thermal_tau() {
    return 1.0 / (RT * TAU_L * TAU_L);
}
```

禁止保留：

```cpp
return 1.0 / RT;
```

除非该版本只用于 `TAU_L = 1.0` 的临时测试，并且不得作为默认实现。

### 1.3 辐射源项系数

温度源项必须为：

\[
S_\theta=C_r(J-\theta^4)
\]

其中：

\[
C_r=\frac{8(1-\omega)}{N_{cr}}
\]

代码必须实现：

```cpp
inline double source_coeff(double Ncr, double omega) {
    return 8.0 * (1.0 - omega) / std::max(Ncr, 1.0e-14);
}
```

必须删除或废弃以下错误逻辑：

```cpp
const double SUN2012_SOURCE_SCALE = 0.125;
```

以及任何形式的：

```cpp
0.125 * TAU_L * TAU_L * (1.0 - omega) / Ncr
```

### 1.4 辐射传输方程

hot-wall 映射后的 RTE 必须为：

\[
\frac{1}{\tau_L}\mathbf s_m\cdot\nabla I_m+I_m
=
(1-\omega)\theta^4+\omega J
\]

局部平衡强度为：

\[
I_m^{eq}=(1-\omega)\theta^4+\omega J
\]

零阶角矩为：

\[
J=\sum_m w_m I_m
\]

代码中辐射传播速度必须包含光学厚度：

```cpp
ax = directions[m].x / TAU_L;
ay = directions[m].y / TAU_L;
az = directions[m].z / TAU_L;
```

禁止在 RTE 伪时间推进中直接把 `directions[m].x/y/z` 当作传播速度使用。

---

## 2. 总体修改原则

### 2.1 不要大规模重写

优先在现有 `3D.cpp` 中局部修改：

```text
thermal_tau()
Energy_equation() 或温度源项相关函数
radiation_dugks()
rec_radiation()
灰壁边界处理
时间步长计算
输出与测试函数
```

不要无必要地重写全部数据结构。

### 2.2 保留 DUGKS 基本框架

现有代码中若已经包含：

```text
g_tilde
g_bar_plus
g_tilde_plus
I_tilde
I_bar_plus
I_tilde_plus
face reconstruction
finite-volume update
```

应尽量保留，并修正方程系数、源项处理、边界条件和时间步。

### 2.3 先正确，后加速

复现阶段必须使用：

```cpp
RAD_COUPLING_INTERVAL = 1;
```

即每个物理时间步都收敛辐射场。

只有在所有基准图表通过后，才允许尝试：

```cpp
RAD_COUPLING_INTERVAL = 2, 4, 8
```

进行加速测试。

---

## 3. 具体代码修改任务

---

# Task 1：建立方程相关工具函数

在 `3D.cpp` 中新增或替换以下函数。

## 3.1 `pow4`

如果已有可保留：

```cpp
inline double pow4(double v) {
    double s = v * v;
    return s * s;
}
```

## 3.2 `thermal_tau`

替换为：

```cpp
inline double thermal_tau() {
    return 1.0 / (RT * TAU_L * TAU_L);
}
```

## 3.3 `source_coeff`

新增：

```cpp
inline double source_coeff(double Ncr, double omega) {
    return 8.0 * (1.0 - omega) / std::max(Ncr, 1.0e-14);
}
```

## 3.4 `radiation_speed`

建议新增：

```cpp
inline void radiation_speed(int m, double& ax, double& ay, double& az) {
    ax = directions[m].x / TAU_L;
    ay = directions[m].y / TAU_L;
    az = directions[m].z / TAU_L;
}
```

后续 RTE 伪时间推进、特征线回溯和 CFL 都必须使用该函数。

---

# Task 2：修正温度源项

查找当前代码中类似函数：

```cpp
Energy_equation(...)
```

或任何直接更新 `T[k][j][i]` 的源项函数。

如果保留 operator splitting 版本，至少必须把源项系数改为：

```cpp
double Cr = source_coeff(Np, omega);
double S  = Cr * (G[k][j][i] - pow4(T[k][j][i]));
```

其中当前代码中的：

```cpp
G
```

语义必须是：

```text
J = Σwm Im
```

即 hot-wall 归一化入射辐射零阶角矩。

不得再使用：

```cpp
SUN2012_SOURCE_SCALE * TAU_L * TAU_L
```

---

# Task 3：将源项并入温度 DUGKS，替代或逐步替代外部分裂

如果时间允许，优先实现含源项 DUGKS，而不是只使用 `Energy_equation()` 做 Strang splitting。

## 3.1 温度 kinetic equation

实现方程：

\[
\frac{\partial g_i}{\partial t^*}
+
\boldsymbol{\xi}_i\cdot\nabla g_i
=
\frac{g_i^{eq}-g_i}{\tau_g}
+
w_iS_\theta
\]

其中：

```cpp
geq = wt[q] * theta;
S_theta = source_coeff(Ncr, omega) * (J - pow4(theta));
```

## 3.2 宏观温度恢复

从 `g_tilde` 恢复温度时必须解：

\[
\theta
=
\sum_i\tilde g_i
+
\frac{\Delta t}{2}C_r(J-\theta^4)
\]

新增：

```cpp
double recover_theta_from_gtilde(double gsum, double J, double dt, double Cr) {
    double A = 0.5 * dt * Cr;

    double theta = std::max(0.5, gsum);

    for (int it = 0; it < 50; ++it) {
        double t2 = theta * theta;
        double t3 = t2 * theta;
        double t4 = t2 * t2;

        double F  = theta - gsum - A * J + A * t4;
        double dF = 1.0 + 4.0 * A * t3;

        double next = theta - F / dF;

        if (!std::isfinite(next) || next < 0.5) {
            next = std::max(0.5, theta - 0.05);
        } else if (next > 1.0) {
            next = std::min(1.0, theta + 0.05);
        }

        if (std::abs(next - theta) < 1.0e-12) {
            return next;
        }

        theta = next;
    }

    return theta;
}
```

然后修改 `rec_temperature()`：

```cpp
void rec_temperature(double dt, double Ncr, double omega) {
    double Cr = source_coeff(Ncr, omega);

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                double gsum = 0.0;
                for (int q = 0; q < Q; ++q) {
                    gsum += g_tilde[k][j][i][q];
                }

                T[k][j][i] = recover_theta_from_gtilde(
                    gsum,
                    G[k][j][i],
                    dt,
                    Cr
                );
            }
        }
    }
}
```

如果现有 `rec_temperature()` 没有参数，请调整其调用点。

## 3.3 `g_bar_plus` 中加入源项

在温度推进函数中，将：

```cpp
g_bar_plus = co * g_tilde + ce * eq;
```

替换为：

```cpp
double src = wt[q] * source_theta;
g_bar_plus = co * g_tilde + ce * eq + cs * src;
```

其中：

```cpp
double h = 0.5 * dt;
double tau_g = thermal_tau();

double co = (2.0 * tau_g - h) / (2.0 * tau_g + dt);
double ce = 3.0 * h / (2.0 * tau_g + dt);
double cs = 3.0 * h * tau_g / (2.0 * tau_g + dt);
```

## 3.4 界面分布中加入源项

从 `g_bar` 恢复界面物理分布时，应使用：

```cpp
double fo = 2.0 * tau_g / (2.0 * tau_g + h);
double fe = h / (2.0 * tau_g + h);
double fs = tau_g * h / (2.0 * tau_g + h);

g_face =
    fo * g_bar_face
  + fe * wt[q] * theta_face
  + fs * wt[q] * source_face;
```

如果短期内不实现界面源项，必须在注释中标明这是过渡版，且至少源项系数和 `thermal_tau()` 必须正确。

---

# Task 4：修正 RTE 伪时间 DUGKS

## 4.1 平衡强度

确保所有辐射平衡强度均为：

```cpp
Ieq = (1.0 - omega) * pow4(T[k][j][i]) + omega * G[k][j][i];
```

其中：

```text
T 是 θ
G 是 J
```

## 4.2 恢复零阶角矩 `G/J`

如果当前代码使用 `G` 表示入射辐射零阶角矩，可暂时保留变量名，但必须确保语义为：

\[
G=J=\sum_m w_m I_m
\]

由于 RTE 辅助变量不是严格守恒型，若从 `I_tilde` 恢复 `J`，建议使用：

```cpp
double recover_J_from_Itilde(double sum_Itilde, double theta, double dchi, double omega) {
    double tau_I = 1.0;
    double a = 0.5 * dchi / tau_I;
    double b = a * (1.0 - omega);

    return (sum_Itilde + b * pow4(theta)) / (1.0 + b);
}
```

然后在 `rec_radiation()` 中：

```cpp
double sum = 0.0;

for (int m = 0; m < M; ++m) {
    sum += directions[m].weight * I_tilde[k][j][i][m];
}

G[k][j][i] = recover_J_from_Itilde(
    sum,
    T[k][j][i],
    dchi,
    omega
);
```

如果当前代码已有等价公式，请检查其是否与上述一致。

## 4.3 特征线速度

所有辐射 face reconstruction 中，回溯点必须使用：

```cpp
ax = directions[m].x / TAU_L;
ay = directions[m].y / TAU_L;
az = directions[m].z / TAU_L;
```

例如：

```cpp
double departure_x = face_x - ax * h;
double departure_y = face_y - ay * h;
double departure_z = face_z - az * h;
```

不要使用：

```cpp
directions[m].x * h
```

除非 `directions[m].x` 已经被定义成 `s_x / TAU_L`，但不建议这样做，避免变量含义混乱。

---

# Task 5：修正时间步长

## 5.1 温度物理时间步

添加函数：

```cpp
double compute_thermal_dt() {
    double max_xi = 0.0;

    for (int q = 0; q < Q; ++q) {
        double mag = std::sqrt(cx[q] * cx[q] + cy[q] * cy[q] + cz[q] * cz[q]);
        max_xi = std::max(max_xi, mag);
    }

    double ds = std::min(dx, std::min(dy, dz));
    return CFL_T * ds / std::max(max_xi, 1.0e-14);
}
```

如果为了匹配 Sun2012 输出时间，需要精确命中 `0.005, 0.015, 0.050`，主循环中应裁剪最后一步：

```cpp
double dt = compute_thermal_dt();

if (time + dt > next_output_time) {
    dt = next_output_time - time;
}
```

## 5.2 辐射伪时间步

添加：

```cpp
double compute_radiation_dchi() {
    double max_a = 0.0;

    for (int m = 0; m < M; ++m) {
        double ax, ay, az;
        radiation_speed(m, ax, ay, az);

        double mag = std::sqrt(ax * ax + ay * ay + az * az);
        max_a = std::max(max_a, mag);
    }

    double ds = std::min(dx, std::min(dy, dz));
    return CFL_R * ds / std::max(max_a, 1.0e-14);
}
```

---

# Task 6：修正边界条件

## 6.1 统一法向约定

新增或明确壁面外法向：

```cpp
struct Vec3 {
    double x, y, z;
};

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 wall_normal(int axis, int face) {
    // face = 0: lower wall, face = 1: upper wall
    if (axis == 0) return face == 0 ? Vec3{-1, 0, 0} : Vec3{1, 0, 0};
    if (axis == 1) return face == 0 ? Vec3{0, -1, 0} : Vec3{0, 1, 0};
    return face == 0 ? Vec3{0, 0, -1} : Vec3{0, 0, 1};
}
```

约定：

```text
n_wall · s > 0：介质打到壁面，入射到壁面
n_wall · s < 0：壁面射入介质，边界给定方向
```

## 6.2 温度定温壁面

保持：

```cpp
double wall_T(int axis, int face) {
    if (axis == 1 && face == 0) {
        return TW; // south wall = 1.0
    }
    return TE;     // other walls = 0.5
}
```

## 6.3 黑壁边界

当：

```cpp
eps_wall == 1.0
```

入射到计算域的方向直接给：

```cpp
I_wall_out[m] = pow4(theta_wall);
```

## 6.4 灰壁反射

灰壁反射必须使用物理强度。

请将任何形式的：

```cpp
incident += bar[m] * ...
```

改为：

```cpp
double I_phys_in = fo * I_bar_face[mp] + fe * I_eq_face[mp];
incident += directions[mp].weight * I_phys_in * cos_theta;
```

完整建议函数：

```cpp
double gray_wall_outgoing_intensity(
    const Vec3& n_wall,
    double theta_wall,
    double eps_wall,
    const double* I_bar_face,
    const double* I_eq_face,
    double fo,
    double fe
) {
    double reflected = 0.0;

    for (int mp = 0; mp < M; ++mp) {
        Vec3 s{directions[mp].x, directions[mp].y, directions[mp].z};
        double cos_theta = dot(n_wall, s);

        if (cos_theta > 0.0) {
            double I_phys_in = fo * I_bar_face[mp] + fe * I_eq_face[mp];

            reflected += directions[mp].weight
                       * I_phys_in
                       * cos_theta;
        }
    }

    double I_emit = eps_wall * pow4(theta_wall);
    double I_refl = 4.0 * (1.0 - eps_wall) * reflected;

    return I_emit + I_refl;
}
```

使用时只对满足：

```cpp
dot(n_wall, s[m]) < 0.0
```

的方向赋值为该壁面出射强度。

---

# Task 7：增加基础测试模式

不要求引入测试框架。可使用编译宏或命令行参数实现。

建议支持：

```text
--test-angles
--test-black-cavity
--test-gray-cavity
--test-conduction
--test-source-newton
--test-table1
--run-sun2012
```

如果不方便解析命令行，可用宏：

```cpp
#define RUN_TEST_ANGLES
```

或在 `main()` 顶部通过简单字符串判断：

```cpp
if (argc >= 2 && std::string(argv[1]) == "--test-angles") {
    test_angles();
    return 0;
}
```

---

## 7.1 角度测试

实现函数：

```cpp
void test_angles();
```

输出：

```text
sum_w
sum_w_sx
sum_w_sy
sum_w_sz
sum_w_sx2
sum_w_sy2
sum_w_sz2
```

验收：

```text
sum_w ≈ 1
一阶矩 ≈ 0
二阶矩 ≈ 1/3
```

## 7.2 等温黑腔测试

实现：

```cpp
void test_black_cavity();
```

设置：

```text
所有 T = 0.5
所有 wall_T = 0.5
eps = 1.0
```

运行辐射伪时间收敛。

验收输出：

```text
max_abs_J_minus_theta4
max_abs_source
```

目标：

```text
< 1e-8
```

## 7.3 等温灰腔测试

实现：

```cpp
void test_gray_cavity();
```

设置：

```text
所有 T = 0.5
eps = 0.5
```

验收同黑腔。

## 7.4 纯导热测试

实现：

```cpp
void test_conduction_scaling();
```

分别测试：

```text
TAU_L = 0.1, 1.0, 5.0
```

输出每个 case 的中心点温度变化或等效扩散速度。

至少确保程序不发散，且扩散快慢趋势为：

```text
tauL=0.1 fastest
tauL=1.0 middle
tauL=5.0 slowest
```

## 7.5 Newton 源项测试

实现：

```cpp
void test_source_newton();
```

输入若干：

```text
gsum, J, dt, Cr
```

输出：

```text
theta
residual
iterations
```

目标：

```text
无 NaN
theta 合理
residual < 1e-10 或至少 < 1e-8
```

---

# Task 8：Sun2012 复现实验矩阵

增加一个集中函数：

```cpp
void run_sun2012_cases();
```

至少支持以下 case。

## 8.1 Table 1

```text
Ncr = 0.01
tauL = 1.0
omega = 0.0
eps = 1.0
t* = 0.005
```

输出：

```text
x=0.5, z=0.5
y=0.25, 0.5, 0.75
```

输出文件建议：

```text
output/table1.csv
```

格式：

```csv
y,theta
0.25,...
0.50,...
0.75,...
```

## 8.2 Fig.3

```text
Ncr = 0.01
tauL = 1.0
omega = 0.0
eps = 1.0
t* = 0.005, 0.015, 0.050, steady
```

输出：

```text
x = 0.5 plane
theta(y,z)
```

文件：

```text
output/fig3_t0005.csv
output/fig3_t0015.csv
output/fig3_t0050.csv
output/fig3_steady.csv
```

## 8.3 Fig.4

```text
tauL = 1.0
Ncr = 0.01
eps = 1.0
omega = 0.5, 0.9
```

中心线输出：

```text
x=0.5, z=0.5, y from 0 to 1
```

## 8.4 Fig.5

```text
tauL = 1.0
omega = 0.0
eps = 1.0
Ncr = 0.1, 1.0
```

## 8.5 Fig.6

```text
tauL = 1.0
omega = 0.0
Ncr = 0.01
eps_south = 0.1, 0.5
```

其余壁面可先使用：

```text
eps_other = 1.0
```

## 8.6 Fig.7

```text
omega = 0.0
Ncr = 0.01
eps = 1.0
tauL = 0.1, 5.0
```

该 case 必须使用 CFL 动态时间步和正确的 \(\tau_L\) 缩放。

---

# Task 9：输出与目录

确保程序自动创建：

```text
output/
```

若已有目录创建函数，可复用。没有则添加：

```cpp
void ensure_dir(const std::string& path) {
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    mkdir(path.c_str(), 0755);
#endif
}
```

输出 CSV 文件时必须带表头，避免后处理混乱。

中心线文件格式：

```csv
y,theta
0.000000,0.500000
...
```

中心面文件格式：

```csv
y,z,theta
0.000000,0.000000,0.500000
...
```

残差日志建议：

```csv
step,time,rad_iter,rad_err,temp_err
```

---

# Task 10：编译要求

代码应至少能使用以下命令编译：

```bash
g++ -O2 -std=c++11 3D.cpp -o sun2012_dugks
```

如果使用 OpenMP：

```bash
g++ -O2 -std=c++11 -fopenmp 3D.cpp -o sun2012_dugks
```

不要引入非标准第三方库。

---

# Task 11：运行命令建议

实现命令行参数后，应支持：

```bash
./sun2012_dugks --test-angles
./sun2012_dugks --test-black-cavity
./sun2012_dugks --test-gray-cavity
./sun2012_dugks --test-source-newton
./sun2012_dugks --test-table1
./sun2012_dugks --run-sun2012
```

如果暂时没有命令行解析，至少在 `main()` 中提供清晰的开关区域：

```cpp
const bool RUN_ANGLES = false;
const bool RUN_TABLE1 = true;
```

---

# Task 12：验收标准

修改完成后，必须满足以下最低标准。

## 12.1 编译标准

```text
g++ -O2 -std=c++11 3D.cpp -o sun2012_dugks
```

必须通过。

## 12.2 方程标准

代码中必须能 grep 到：

```text
8.0 * (1.0 - omega)
RT * TAU_L * TAU_L
directions[m].x / TAU_L
```

代码中不得再存在：

```text
SUN2012_SOURCE_SCALE
0.125 * TAU_L * TAU_L
return 1.0 / RT;
```

注意：如果 `return 1.0 / RT;` 作为注释出现可以接受，但不能作为默认逻辑。

## 12.3 基础测试标准

```text
--test-angles 通过
--test-black-cavity 中 max|J-theta^4| < 1e-8
--test-gray-cavity 中 max|J-theta^4| < 1e-8 或至少 < 1e-6
--test-source-newton 无 NaN 且残差足够小
```

## 12.4 Table 1 标准

`output/table1.csv` 必须生成。

中心线三点：

```text
y = 0.25
y = 0.50
y = 0.75
```

必须输出。

短期目标：

```text
趋势正确，数值稳定
```

最终目标：

```text
相对误差 < 1%
```

## 12.5 图表标准

至少生成：

```text
output/fig3_t0005.csv
output/fig4_omega05.csv
output/fig5_ncr01.csv
output/fig6_eps01.csv
output/fig7_tauL01.csv
```

---

# Task 13：优先级顺序

必须严格按以下顺序修改，不要跳步：

```text
1. 修正 thermal_tau()
2. 修正 source_coeff()
3. 修正 RTE 速度 s/tauL
4. 修正 radiation dchi CFL
5. 修正温度 dt CFL
6. 修正 RTE 平衡强度 Ieq
7. 修正 J 恢复公式
8. 修正灰壁反射使用物理强度
9. 添加基础测试
10. 添加 Table 1 输出
11. 添加 Fig.3-Fig.7 输出
12. 再考虑含源项 DUGKS 完全化
13. 再考虑非均匀网格
14. 再考虑 RAD_COUPLING_INTERVAL 加速
```

如果时间有限，前 10 项必须完成。

---

# Task 14：不要做的事情

不要做以下修改：

```text
不要把 cold wall 改成 0
不要把 south wall 改成 2
不要混用 H=T/Tref 和 θ=T/Ts
不要把 tauL 只用于 RTE 而不用于 thermal_tau
不要把 0.125 源项比例继续保留
不要让灰壁反射用 I_bar 直接积分
不要默认 RAD_COUPLING_INTERVAL > 1
不要一开始就跑所有 Fig.3-Fig.7 而跳过基础测试
不要引入 Eigen、Boost、Python 依赖
```

---

# Task 15：建议注释写法

在代码靠近参数定义处加入以下注释，帮助后续维护：

```cpp
/*
Sun2012 hot-wall normalization contract:

theta = T / T_s = H / 2
cold walls and initial field: theta = 0.5
south hot wall: theta = 1.0

Energy equation in Sun2012 time t*:
dtheta/dt* = (1/tauL^2) Laplacian(theta)
             + 8(1-omega)/Ncr * (J - theta^4)

Radiation equation:
(1/tauL) s · grad(I_m) + I_m
    = (1-omega) theta^4 + omega J

J = sum_m w_m I_m

Therefore:
thermal_tau = 1 / (RT * tauL^2)
source_coeff = 8 * (1-omega) / Ncr
radiation velocity = s / tauL
*/
```

---

# Task 16：最终交付

完成修改后，请输出以下内容：

```text
1. 修改了哪些函数
2. 是否删除 SUN2012_SOURCE_SCALE
3. thermal_tau 当前公式
4. source_coeff 当前公式
5. RTE 速度是否使用 s/tauL
6. 灰壁反射是否使用 I_phys
7. 已通过哪些测试
8. 生成了哪些 output/*.csv
9. 编译命令
10. 运行命令
```

---

## 附录 A：Table 1 调试排查顺序

如果 `table1.csv` 数值偏差大，请按以下顺序排查：

```text
1. T 是否代表 θ，而非 H
2. south wall 是否为 1.0
3. 其他壁是否为 0.5
4. thermal_tau 是否随 tauL 改变
5. source_coeff 是否为 8(1-omega)/Ncr
6. G/J 是否接近合理范围，例如 0.5^4 到 1^4 附近
7. RTE 是否每个物理步收敛
8. rad_err 是否真的低于 1e-6
9. dchi 在 tauL=0.1 时是否显著变小
10. 灰壁/黑壁方向判断是否正确
11. 输出时间是否准确到 t*=0.005
12. 插值中心线位置是否正确
```

---

## 附录 B：建议的最小补丁版本

如果不能一次完成完整含源项 DUGKS，至少完成以下最小补丁：

```text
1. thermal_tau = 1/(RT*tauL^2)
2. source_coeff = 8(1-omega)/Ncr
3. Energy_equation 使用 Cr*(J-theta^4)
4. RTE velocity = s/tauL
5. dchi = CFL_R*dx/max(|s/tauL|)
6. RAD_COUPLING_INTERVAL = 1
7. 灰壁反射用 I_phys
8. 添加 --test-angles 和 --test-table1
```

该最小补丁版本仍可能有 splitting error，但可以快速判断主要系统性误差是否已经消除。

---

## 附录 C：完成后的理想结构

若允许重构，建议最终形成：

```text
src/
  main.cpp
  parameters.hpp
  grid.hpp
  angular_quadrature.hpp
  thermal_dugks.hpp
  radiation_dugks.hpp
  boundary.hpp
  output.hpp
tests/
  test_angles.cpp
  test_black_cavity.cpp
  test_gray_cavity.cpp
  test_source_newton.cpp
```

但当前阶段优先保持 `3D.cpp` 可编译和可复现，不强制拆文件。
