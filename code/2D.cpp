#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#endif

using namespace std;

const double PI = 3.14159265358979323846;
const double TW = 1.0;
const double TE = 0.5;
const int NX = 41;
const int NY = 41;
const int Q = 9;
const int NTHETA = 8;
const int MAX_PHI = 16;
const int MAX_M = NTHETA * MAX_PHI;
const int NTIME = 3;
const int NCURVE = NTIME + 1;
const int CASE_out = 4;
const int NOUT = 41;
const double CFL_T = 0.50;
const double CFL_R = 0.80;
const double STEADY_error = 5.0e-8;
const int STEADY_MAX_STEPS = 30000;
const double RAD_error = 2.0e-7;
const int RAD_MAX_STEPS = 12000;
const int RAD_TRACKING_STEPS = 1;
const int RAD_COUPLING_INTERVAL = 8;
const double RAD_SLOPE_FACTOR = 0.0;

const double dx = 1.0 / NX;
const double dy = 1.0 / NY;
// Thermal DUGKS model parameters.
// alpha*=tau_t*RT=1 is kept independent of the mesh.  This avoids the
// traditional LBM diffusive scaling tau_t~dx^2 and lets dt follow CFL.
const double tau_t = 1.0;
const double RT = 1.0;
const double speed_t = 1.7320508075688772935 * sqrt(RT);

const double save_times[NTIME] = {0.005, 0.015, 0.04};
const int phi_values[CASE_out] = {4, 8, 12, 16};

const double c_tx[Q] = {
    0.0, speed_t, 0.0, -speed_t, 0.0,
    speed_t, -speed_t, -speed_t, speed_t
};
const double c_ty[Q] = {
    0.0, 0.0, speed_t, 0.0, -speed_t,
    speed_t, speed_t, -speed_t, -speed_t
};
const double w_t[Q] = {
    4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
    1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0,
    1.0 / 36.0, 1.0 / 36.0
};
const int opposite[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

// === 场变量数组 ===
double T[NY + 2][NX + 2];                     // 温度
double pre_temp[NY + 2][NX + 2];              // 上一时间步温度
double G[NY + 2][NX + 2];                     // 入射辐射
double G_pre[NY + 2][NX + 2];                 // 上一次辐射迭代的入射辐射
double Qr[NY + 2][NX + 2];                    // 辐射热流

double g_tilde[NY + 2][NX + 2][Q];            // 热传导分布函数
double g_bar_plus[NY + 2][NX + 2][Q];         // 热传导半时间步碰撞后分布函数
double g_tilde_plus[NY + 2][NX + 2][Q];       // 热传导整时间步碰撞后分布函数
double g_slope_x[NY + 2][NX + 2][Q];          // 热传导分布函数的 x 向斜率
double g_slope_y[NY + 2][NX + 2][Q];          // 热传导分布函数的 y 向斜率
double g_face_x[NY + 2][NX + 1][Q];           // 热传导 x 向界面分布函数
double g_face_y[NY + 1][NX + 2][Q];           // 热传导 y 向界面分布函数

double mu_x[MAX_M];                            // 离散方向的 x 分量
double mu_y[MAX_M];                            // 离散方向的 y 分量
double weights[MAX_M];                         // 离散角度权重
double flux_weights[MAX_M];                    // 辐射热流积分权重
double I_tilde[NY + 2][NX + 2][MAX_M];        // 辐射强度分布函数
double I_bar_plus[NY + 2][NX + 2][MAX_M];     // 辐射半时间步碰撞后分布函数
double I_tilde_plus[NY + 2][NX + 2][MAX_M];   // 辐射整时间步碰撞后分布函数
double I_slope_x[NY + 2][NX + 2][MAX_M];      // 辐射强度的 x 向斜率
double I_slope_y[NY + 2][NX + 2][MAX_M];      // 辐射强度的 y 向斜率
double I_face_x[NY + 2][NX + 1][MAX_M];       // 辐射 x 向界面分布函数
double I_face_y[NY + 1][NX + 2][MAX_M];       // 辐射 y 向界面分布函数

// === 保存结果数组 ===
double saved_temp[NCURVE][NOUT];
double fig5_temp[CASE_out][NOUT];
double fig5_Qr[CASE_out][NOUT];
double fig5_Qc[CASE_out][NOUT];

// === 运行状态 ===
int M = 0;
int thermal_step_counter = 0;
bool radiation_needs_initial_convergence = true;

// === 函数声明 ===
double minmod(double a, double b);
double pow4(double value);
void ini_angles(int nphi);
void ini_case(int nphi);
void rec_temperture();
void calculate_temperature_slopes();
void temperture(double dt);
void ini_radiation();
void rec_radiation(double omega, double pseudo_dt);
int radiation_dugks(double omega, bool force_convergence = false);
void Energy_equation(double dt, double N, double omega);
void coupled_step(double dt, double N, double omega);
double sample_cell_field(const double field[NY + 2][NX + 2], int i,
                         double y, double bottom_value, double top_value);
void save_curve(int curve_id);
void run_temp_case(const string& fig, double N, double omega, int nphi);
int run_steady_case(double N, double omega, int nphi);
void run_fig5();
void run_sun2016_table2_case();
void write_temp_csv(const string& path);
void write_fig5_csv(const string& path, double data[CASE_out][NOUT]);

int main() {
#ifdef _WIN32
    mkdir("output");
#else
    mkdir("output", 0755);
#endif
    cout << "2D coupled DUGKS conduction-radiation cases -> output\n";
    cout << "T: D2Q9 DUGKS, radiation: full-solid-angle DUGKS\n";

    run_fig5();
    cout << "wrote 5a.csv, 5b.csv, 5c.csv\n";
    run_sun2016_table2_case();

#ifndef FIG5_ONLY
    run_temp_case("6a", 0.01, 0.0, 16);
    run_temp_case("6b", 0.01, 0.5, 16);
    run_temp_case("6c", 0.01, 0.9, 16);

    run_temp_case("7a", 0.01, 0.0, 16);
    run_temp_case("7b", 0.10, 0.0, 16);
    run_temp_case("7c", 1.00, 0.0, 16);
#endif

    return 0;
}

double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return fabs(a) < fabs(b) ? a : b;
}

double pow4(double value) {
    const double square = value * value;
    return square * square;
}

void ini_angles(int nphi) {
    const double dtheta = PI / NTHETA;
    const double dphi = 2.0 * PI / nphi;
    M = NTHETA * nphi;
    int m = 0;

    for (int it = 0; it < NTHETA; ++it) {
        const double theta_l = it * dtheta;
        const double theta_r = (it + 1) * dtheta;
        const double theta_flux =
            0.5 * (theta_r - theta_l)
          - 0.25 * (sin(2.0 * theta_r) - sin(2.0 * theta_l));

        for (int ip = 0; ip < nphi; ++ip) {
            const double phi_l = ip * dphi;
            const double phi_r = (ip + 1) * dphi;
            const double solid_angle =
                (cos(theta_l) - cos(theta_r)) * dphi;
            const double direction_x =
                theta_flux * (sin(phi_r) - sin(phi_l));
            const double direction_y =
                theta_flux * (cos(phi_l) - cos(phi_r));

            mu_x[m] = direction_x / solid_angle;
            mu_y[m] = direction_y / solid_angle;
            weights[m] = solid_angle / (4.0 * PI);
            flux_weights[m] = direction_y / (4.0 * PI);
            ++m;
        }
    }
}

void ini_case(int nphi) {
    ini_angles(nphi);
    thermal_step_counter = 0;
    radiation_needs_initial_convergence = true;

    for (int j = 1; j <= NY; ++j) {
        for (int i = 1; i <= NX; ++i) {
            T[j][i] = TE;
            for (int q = 0; q < Q; ++q) {
                g_tilde[j][i][q] = w_t[q] * TE;
            }
        }
    }
    ini_radiation();
}

void rec_temperture() {
    for (int j = 1; j <= NY; ++j) {
        for (int i = 1; i <= NX; ++i) {
            double sum = 0.0;
            for (int q = 0; q < Q; ++q) sum += g_tilde[j][i][q];
            T[j][i] = sum;
        }
    }
}

void calculate_temperature_slopes() {
    for (int j = 1; j <= NY; ++j) {
        for (int i = 1; i <= NX; ++i) {
            for (int q = 0; q < Q; ++q) {
                if (i == 1 || i == NX) {
                    g_slope_x[j][i][q] = 0.0;
                } else {
                    const double left =
                        (g_bar_plus[j][i][q] - g_bar_plus[j][i - 1][q]) / dx;
                    const double right =
                        (g_bar_plus[j][i + 1][q] - g_bar_plus[j][i][q]) / dx;
                    g_slope_x[j][i][q] = minmod(left, right);
                }

                if (j == 1 || j == NY) {
                    g_slope_y[j][i][q] = 0.0;
                } else {
                    const double bottom =
                        (g_bar_plus[j][i][q] - g_bar_plus[j - 1][i][q]) / dy;
                    const double top =
                        (g_bar_plus[j + 1][i][q] - g_bar_plus[j][i][q]) / dy;
                    g_slope_y[j][i][q] = minmod(bottom, top);
                }
            }
        }
    }
}

void temperture(double dt) {
    const double h = 0.5 * dt;
    const double coefficient_old = (2.0 * tau_t - h) / (2.0 * tau_t + dt);
    const double coefficient_eq = 3.0 * h / (2.0 * tau_t + dt);
    const double face_old = 2.0 * tau_t / (2.0 * tau_t + h);
    const double face_eq = h / (2.0 * tau_t + h);

    rec_temperture();
    for (int j = 1; j <= NY; ++j) {
        for (int i = 1; i <= NX; ++i) {
            for (int q = 0; q < Q; ++q) {
                const double eq = w_t[q] * T[j][i];
                g_bar_plus[j][i][q] =
                    coefficient_old * g_tilde[j][i][q]
                  + coefficient_eq * eq;
                g_tilde_plus[j][i][q] =
                    (4.0 / 3.0) * g_bar_plus[j][i][q]
                  - (1.0 / 3.0) * g_tilde[j][i][q];
            }
        }
    }
    calculate_temperature_slopes();

    for (int j = 1; j <= NY; ++j) {
        for (int face = 0; face <= NX; ++face) {
            double bar_value[Q];
            for (int q = 0; q < Q; ++q) {
                int cell;
                if (c_tx[q] > 0.0) cell = max(1, face);
                else if (c_tx[q] < 0.0) cell = min(NX, face + 1);
                else cell = min(NX, max(1, face));

                const double cell_x = (cell - 0.5) * dx;
                const double departure_x = face * dx - c_tx[q] * h;
                bar_value[q] = g_bar_plus[j][cell][q]
                    + g_slope_x[j][cell][q] * (departure_x - cell_x)
                    - g_slope_y[j][cell][q] * c_ty[q] * h;
            }

            if (face == 0 || face == NX) {
                const double wall_temperature = TE;
                for (int q = 0; q < Q; ++q) {
                    const bool incoming =
                        (face == 0 && c_tx[q] > 0.0)
                     || (face == NX && c_tx[q] < 0.0);
                    if (incoming) {
                        bar_value[q] =
                            -bar_value[opposite[q]]
                            + 2.0 * w_t[q] * wall_temperature;
                    } else if (c_tx[q] == 0.0) {
                        bar_value[q] = w_t[q] * wall_temperature;
                    }
                }
            }

            double face_temp = 0.0;
            for (int q = 0; q < Q; ++q) face_temp += bar_value[q];
            if (face == 0 || face == NX) face_temp = TE;

            for (int q = 0; q < Q; ++q) {
                g_face_x[j][face][q] =
                    face_old * bar_value[q]
                  + face_eq * w_t[q] * face_temp;
            }
        }
    }

    for (int i = 1; i <= NX; ++i) {
        for (int face = 0; face <= NY; ++face) {
            double bar_value[Q];
            for (int q = 0; q < Q; ++q) {
                int cell;
                if (c_ty[q] > 0.0) cell = max(1, face);
                else if (c_ty[q] < 0.0) cell = min(NY, face + 1);
                else cell = min(NY, max(1, face));

                const double cell_y = (cell - 0.5) * dy;
                const double departure_y = face * dy - c_ty[q] * h;
                bar_value[q] = g_bar_plus[cell][i][q]
                    + g_slope_y[cell][i][q] * (departure_y - cell_y)
                    - g_slope_x[cell][i][q] * c_tx[q] * h;
            }

            if (face == 0 || face == NY) {
                const double wall_temperature = face == 0 ? TW : TE;
                for (int q = 0; q < Q; ++q) {
                    const bool incoming =
                        (face == 0 && c_ty[q] > 0.0)
                     || (face == NY && c_ty[q] < 0.0);
                    if (incoming) {
                        bar_value[q] =
                            -bar_value[opposite[q]]
                            + 2.0 * w_t[q] * wall_temperature;
                    } else if (c_ty[q] == 0.0) {
                        bar_value[q] = w_t[q] * wall_temperature;
                    }
                }
            }

            double face_temp = 0.0;
            for (int q = 0; q < Q; ++q) face_temp += bar_value[q];
            if (face == 0) face_temp = TW;
            if (face == NY) face_temp = TE;

            for (int q = 0; q < Q; ++q) {
                g_face_y[face][i][q] =
                    face_old * bar_value[q]
                  + face_eq * w_t[q] * face_temp;
            }
        }
    }

    for (int j = 1; j <= NY; ++j) {
        for (int i = 1; i <= NX; ++i) {
            for (int q = 0; q < Q; ++q) {
                g_tilde[j][i][q] = g_tilde_plus[j][i][q]
                    - (dt / dx) * c_tx[q]
                    * (g_face_x[j][i][q] - g_face_x[j][i - 1][q])
                    - (dt / dy) * c_ty[q]
                    * (g_face_y[j][i][q] - g_face_y[j - 1][i][q]);
            }
        }
    }
    rec_temperture();
}

void ini_radiation() {
    for (int j = 1; j <= NY; ++j) {
        for (int i = 1; i <= NX; ++i) {
            const double eq = pow4(T[j][i]);
            G[j][i] = eq;
            Qr[j][i] = 0.0;
            for (int m = 0; m < M; ++m) {
                I_tilde[j][i][m] = eq;
            }
        }
    }
}

void rec_radiation(double omega, double pseudo_dt) {
    const double moment_coefficient = 0.5 * pseudo_dt * (1.0 - omega);
    const double physical_old = 2.0 / (2.0 + pseudo_dt);
    const double physical_eq = pseudo_dt / (2.0 + pseudo_dt);

    for (int j = 1; j <= NY; ++j) {
        for (int i = 1; i <= NX; ++i) {
            double tilde_moment = 0.0;
            for (int m = 0; m < M; ++m) {
                tilde_moment += weights[m] * I_tilde[j][i][m];
            }
            const double T4 = pow4(T[j][i]);
            G[j][i] =
                (tilde_moment + moment_coefficient * T4)
                / (1.0 + moment_coefficient);
            const double source =
                (1.0 - omega) * T4 + omega * G[j][i];

            Qr[j][i] = 0.0;
            for (int m = 0; m < M; ++m) {
                const double physical =
                    physical_old * I_tilde[j][i][m] + physical_eq * source;
                Qr[j][i] += flux_weights[m] * physical;
            }
        }
    }
}

int radiation_dugks(double omega, bool force_convergence) {
    double transport_limit = 0.0;
    for (int m = 0; m < M; ++m) {
        transport_limit =
            max(transport_limit, fabs(mu_x[m]) / dx + fabs(mu_y[m]) / dy);
    }
    const double pseudo_dt = CFL_R / transport_limit;
    const double h = 0.5 * pseudo_dt;
    const double coefficient_old = (2.0 - h) / (2.0 + pseudo_dt);
    const double coefficient_eq = 3.0 * h / (2.0 + pseudo_dt);
    const double face_old = 2.0 / (2.0 + h);
    const double face_eq = h / (2.0 + h);
    const double face_moment_coefficient = 0.5 * h * (1.0 - omega);
    const bool converge =
        force_convergence || radiation_needs_initial_convergence;
    const int maximum_steps = converge ? RAD_MAX_STEPS : RAD_TRACKING_STEPS;

    for (int step = 1; step <= maximum_steps; ++step) {
        rec_radiation(omega, pseudo_dt);
        for (int j = 1; j <= NY; ++j) {
            for (int i = 1; i <= NX; ++i) {
                G_pre[j][i] = G[j][i];
                const double T4 = pow4(T[j][i]);
                const double source =
                    (1.0 - omega) * T4 + omega * G[j][i];

                for (int m = 0; m < M; ++m) {
                    I_bar_plus[j][i][m] =
                        coefficient_old * I_tilde[j][i][m]
                      + coefficient_eq * source;
                    I_tilde_plus[j][i][m] =
                        (4.0 / 3.0) * I_bar_plus[j][i][m]
                      - (1.0 / 3.0) * I_tilde[j][i][m];
                }
            }
        }

        for (int j = 1; j <= NY; ++j) {
            for (int i = 1; i <= NX; ++i) {
                for (int m = 0; m < M; ++m) {
                    if (i == 1 || i == NX) {
                        I_slope_x[j][i][m] = 0.0;
                    } else {
                        const double left =
                            (I_bar_plus[j][i][m]
                            - I_bar_plus[j][i - 1][m]) / dx;
                        const double right =
                            (I_bar_plus[j][i + 1][m]
                            - I_bar_plus[j][i][m]) / dx;
                        I_slope_x[j][i][m] =
                            RAD_SLOPE_FACTOR * minmod(left, right);
                    }

                    if (j == 1 || j == NY) {
                        I_slope_y[j][i][m] = 0.0;
                    } else {
                        const double bottom =
                            (I_bar_plus[j][i][m]
                            - I_bar_plus[j - 1][i][m]) / dy;
                        const double top =
                            (I_bar_plus[j + 1][i][m]
                            - I_bar_plus[j][i][m]) / dy;
                        I_slope_y[j][i][m] =
                            RAD_SLOPE_FACTOR * minmod(bottom, top);
                    }
                }
            }
        }

        for (int j = 1; j <= NY; ++j) {
            for (int face = 0; face <= NX; ++face) {
                double bar_value[MAX_M];
                for (int m = 0; m < M; ++m) {
                    int cell;
                    if (mu_x[m] > 0.0) cell = max(1, face);
                    else if (mu_x[m] < 0.0) cell = min(NX, face + 1);
                    else cell = min(NX, max(1, face));

                    const double cell_x = (cell - 0.5) * dx;
                    const double departure_x = face * dx - mu_x[m] * h;
                    bar_value[m] = I_bar_plus[j][cell][m]
                        + I_slope_x[j][cell][m] * (departure_x - cell_x)
                        - I_slope_y[j][cell][m] * mu_y[m] * h;
                }

                if (face == 0 || face == NX) {
                    const int cell = face == 0 ? 1 : NX;
                    const double T4 = pow4(T[j][cell]);
                    const double source =
                        (1.0 - omega) * T4 + omega * G[j][cell];
                    for (int m = 0; m < M; ++m) {
                        const bool incoming =
                            (face == 0 && mu_x[m] > 0.0)
                         || (face == NX && mu_x[m] < 0.0);
                        if (incoming) {
                            I_face_x[j][face][m] = pow4(TE);
                        } else {
                            I_face_x[j][face][m] =
                                face_old * bar_value[m] + face_eq * source;
                        }
                    }
                } else {
                    double bar_moment = 0.0;
                    for (int m = 0; m < M; ++m) {
                        bar_moment += weights[m] * bar_value[m];
                    }
                    const double face_T =
                        0.5 * (T[j][face]
                             + T[j][face + 1]);
                    const double face_T4 = pow4(face_T);
                    const double face_G =
                        (bar_moment + face_moment_coefficient * face_T4)
                        / (1.0 + face_moment_coefficient);
                    const double source =
                        (1.0 - omega) * face_T4 + omega * face_G;

                    for (int m = 0; m < M; ++m) {
                        I_face_x[j][face][m] =
                            face_old * bar_value[m] + face_eq * source;
                    }
                }
            }
        }

        for (int i = 1; i <= NX; ++i) {
            for (int face = 0; face <= NY; ++face) {
                double bar_value[MAX_M];
                for (int m = 0; m < M; ++m) {
                    int cell;
                    if (mu_y[m] > 0.0) cell = max(1, face);
                    else if (mu_y[m] < 0.0) cell = min(NY, face + 1);
                    else cell = min(NY, max(1, face));

                    const double cell_y = (cell - 0.5) * dy;
                    const double departure_y = face * dy - mu_y[m] * h;
                    bar_value[m] = I_bar_plus[cell][i][m]
                        + I_slope_y[cell][i][m] * (departure_y - cell_y)
                        - I_slope_x[cell][i][m] * mu_x[m] * h;
                }

                if (face == 0 || face == NY) {
                    const int cell = face == 0 ? 1 : NY;
                    const double wall_temperature = face == 0 ? TW : TE;
                    const double T4 = pow4(T[cell][i]);
                    const double source =
                        (1.0 - omega) * T4 + omega * G[cell][i];
                    for (int m = 0; m < M; ++m) {
                        const bool incoming =
                            (face == 0 && mu_y[m] > 0.0)
                         || (face == NY && mu_y[m] < 0.0);
                        if (incoming) {
                            I_face_y[face][i][m] = pow4(wall_temperature);
                        } else {
                            I_face_y[face][i][m] =
                                face_old * bar_value[m] + face_eq * source;
                        }
                    }
                } else {
                    double bar_moment = 0.0;
                    for (int m = 0; m < M; ++m) {
                        bar_moment += weights[m] * bar_value[m];
                    }
                    const double face_T =
                        0.5 * (T[face][i]
                             + T[face + 1][i]);
                    const double face_T4 = pow4(face_T);
                    const double face_G =
                        (bar_moment + face_moment_coefficient * face_T4)
                        / (1.0 + face_moment_coefficient);
                    const double source =
                        (1.0 - omega) * face_T4 + omega * face_G;

                    for (int m = 0; m < M; ++m) {
                        I_face_y[face][i][m] =
                            face_old * bar_value[m] + face_eq * source;
                    }
                }
            }
        }

        for (int j = 1; j <= NY; ++j) {
            for (int i = 1; i <= NX; ++i) {
                for (int m = 0; m < M; ++m) {
                    I_tilde[j][i][m] = I_tilde_plus[j][i][m]
                        - (pseudo_dt / dx) * mu_x[m]
                        * (I_face_x[j][i][m] - I_face_x[j][i - 1][m])
                        - (pseudo_dt / dy) * mu_y[m]
                        * (I_face_y[j][i][m] - I_face_y[j - 1][i][m]);
                }
            }
        }

        rec_radiation(omega, pseudo_dt);
        double error = 0.0;
        for (int j = 1; j <= NY; ++j) {
            for (int i = 1; i <= NX; ++i) {
                error = max(error, fabs(G[j][i] - G_pre[j][i]));
            }
        }
        if (converge && step > 30 && error < RAD_error) {
            radiation_needs_initial_convergence = false;
            return step;
        }
    }

    radiation_needs_initial_convergence = false;
    return maximum_steps;
}

void Energy_equation(double dt, double N, double omega) {
    const double A = (1.0 - omega) / max(N, 1.0e-14);
    for (int j = 1; j <= NY; ++j) {
        for (int i = 1; i <= NX; ++i) {
            const double source =
                A * (G[j][i] - pow4(T[j][i]));
            const double delta = dt * source;
            T[j][i] += delta;
            for (int q = 0; q < Q; ++q) {
                g_tilde[j][i][q] += w_t[q] * delta;
            }
        }
    }
}

void coupled_step(double dt, double N, double omega) {
    const bool update_radiation =
        radiation_needs_initial_convergence
        || thermal_step_counter % RAD_COUPLING_INTERVAL == 0;

    if (update_radiation) radiation_dugks(omega, false);
    Energy_equation(0.5 * dt, N, omega);
    temperture(dt);
    if (update_radiation) radiation_dugks(omega, false);
    Energy_equation(0.5 * dt, N, omega);
    rec_temperture();
    ++thermal_step_counter;
}

double sample_cell_field(const double field[NY + 2][NX + 2], int i,
                         double y, double bottom_value, double top_value) {
    if (y <= 0.0) return bottom_value;
    if (y >= 1.0) return top_value;

    const double center_index = y / dy + 0.5;
    int bottom_cell = static_cast<int>(floor(center_index));
    if (bottom_cell < 1) {
        return bottom_value
            + (field[1][i] - bottom_value) * y / (0.5 * dy);
    }
    if (bottom_cell >= NY) {
        const double y0 = (NY - 0.5) * dy;
        return field[NY][i]
            + (top_value - field[NY][i]) * (y - y0) / (1.0 - y0);
    }

    const double y0 = (bottom_cell - 0.5) * dy;
    return field[bottom_cell][i]
        + (field[bottom_cell + 1][i] - field[bottom_cell][i])
        * (y - y0) / dy;
}

void save_curve(int curve_id) {
    const int center_i = (NX + 1) / 2;
    for (int n = 0; n < NOUT; ++n) {
        const double y = static_cast<double>(n) / (NOUT - 1);
        saved_temp[curve_id][n] =
            sample_cell_field(T, center_i, y, TW, TE);
    }
}

void run_temp_case(const string& fig, double N, double omega, int nphi) {
    ini_case(nphi);
    const double base_dt =
        CFL_T / (speed_t / dx + speed_t / dy);
    double time = 0.0;

    for (int k = 0; k < NTIME; ++k) {
        while (time + 1.0e-14 < save_times[k]) {
            const double dt = min(base_dt, save_times[k] - time);
            coupled_step(dt, N, omega);
            time += dt;
        }
        radiation_dugks(omega, true);
        save_curve(k);
    }

    int steady_step = 0;
    for (; steady_step < STEADY_MAX_STEPS; ++steady_step) {
        for (int j = 1; j <= NY; ++j) {
            for (int i = 1; i <= NX; ++i) {
                pre_temp[j][i] = T[j][i];
            }
        }
        coupled_step(base_dt, N, omega);

        if (thermal_step_counter % RAD_COUPLING_INTERVAL == 0) {
            double error = 0.0;
            for (int j = 1; j <= NY; ++j) {
                for (int i = 1; i <= NX; ++i) {
                    error = max(error,
                        fabs(T[j][i] - pre_temp[j][i]));
                }
            }
            if (steady_step > 200 && error < STEADY_error) break;
        }
    }

    radiation_dugks(omega, true);
    save_curve(NTIME);
    write_temp_csv("output/" + fig + ".csv");
    cout << "wrote " << fig << ".csv, steady steps = "
         << steady_step << "\n";
}

int run_steady_case(double N, double omega, int nphi) {
    ini_case(nphi);
    const double base_dt =
        CFL_T / (speed_t / dx + speed_t / dy);
    int step = 0;

    for (; step < STEADY_MAX_STEPS; ++step) {
        for (int j = 1; j <= NY; ++j) {
            for (int i = 1; i <= NX; ++i) {
                pre_temp[j][i] = T[j][i];
            }
        }
        coupled_step(base_dt, N, omega);

        if (thermal_step_counter % RAD_COUPLING_INTERVAL == 0) {
            double error = 0.0;
            for (int j = 1; j <= NY; ++j) {
                for (int i = 1; i <= NX; ++i) {
                    error = max(error,
                        fabs(T[j][i] - pre_temp[j][i]));
                }
            }
            if (step > 200 && error < STEADY_error) break;
        }
    }
    radiation_dugks(omega, true);
    return step;
}

void run_fig5() {
    const int center_i = (NX + 1) / 2;
    double conductive_cell[NY + 2][NX + 2];

    for (int p = 0; p < CASE_out; ++p) {
        const int steps = run_steady_case(0.01, 0.0, phi_values[p]);

        for (int j = 2; j < NY; ++j) {
            conductive_cell[j][center_i] =
                -(T[j + 1][center_i]
                 - T[j - 1][center_i]) / (2.0 * dy);
        }
        conductive_cell[1][center_i] =
            -(T[2][center_i] - TW) / (1.5 * dy);
        conductive_cell[NY][center_i] =
            -(TE - T[NY - 1][center_i]) / (1.5 * dy);

        for (int n = 0; n < NOUT; ++n) {
            const double y = static_cast<double>(n) / (NOUT - 1);
            fig5_temp[p][n] =
                sample_cell_field(T, center_i, y, TW, TE);
            fig5_Qr[p][n] =
                sample_cell_field(Qr, center_i, y,
                                Qr[1][center_i],
                                Qr[NY][center_i]) / 0.01;

            if (n == 0) {
                fig5_Qc[p][n] =
                    -(T[1][center_i] - TW) / (0.5 * dy);
            } else if (n == NOUT - 1) {
                fig5_Qc[p][n] =
                    -(TE - T[NY][center_i]) / (0.5 * dy);
            } else {
                fig5_Qc[p][n] =
                    sample_cell_field(conductive_cell, center_i, y,
                                    conductive_cell[1][center_i],
                                    conductive_cell[NY][center_i]);
            }
        }
        cout << "Fig.5 phi=" << phi_values[p]
             << ", steady steps = " << steps << "\n";
    }

    write_fig5_csv("output/5a.csv", fig5_temp);
    write_fig5_csv("output/5b.csv", fig5_Qr);
    write_fig5_csv("output/5c.csv", fig5_Qc);
}

void run_sun2016_table2_case() {
    const int center_i = (NX + 1) / 2;
    const double sample_N[3] = {1.0, 0.1, 0.01};
    const double sample_y[3] = {0.3, 0.5, 0.7};
    const double wu_ou[3][3] = {
        {0.733, 0.630, 0.560},
        {0.760, 0.663, 0.590},
        {0.791, 0.725, 0.663}
    };
    const double yuen_takara[3][3] = {
        {0.737, 0.630, 0.560},
        {0.763, 0.661, 0.589},
        {0.807, 0.726, 0.653}
    };
    const double mondal_mishra[3][3] = {
        {0.738, 0.631, 0.565},
        {0.761, 0.667, 0.598},
        {0.777, 0.722, 0.672}
    };
    const double sun_zhang[3][3] = {
        {0.734, 0.630, 0.565},
        {0.757, 0.664, 0.596},
        {0.786, 0.725, 0.668}
    };

    ofstream out("output/table2_dugks.csv");
    ofstream data_out("data/table2_dugks.csv");
    ofstream table_out("output/table2.csv");
    ofstream table_data("data/table2.csv");
    ofstream tex_out("output/table2_dugks_rows.tex");
    ofstream tex_data("data/table2_dugks_rows.tex");
    out << "N,y,DUGKS,WuOu1992,YuenTakara1994,"
        << "MondalMishra2007,SunZhang2016\n";
    data_out << "N,y,DUGKS,WuOu1992,YuenTakara1994,"
             << "MondalMishra2007,SunZhang2016\n";
    table_out << "N,y,DUGKS,WuOu1992,YuenTakara1994,"
              << "MondalMishra2007,SunZhang2016\n";
    table_data << "N,y,DUGKS,WuOu1992,YuenTakara1994,"
               << "MondalMishra2007,SunZhang2016\n";
    for (int n = 0; n < 3; ++n) {
        run_steady_case(sample_N[n], 0.0, 16);
        for (int j = 0; j < 3; ++j) {
            const double value =
                sample_cell_field(T, center_i, sample_y[j], TW, TE);
            out << fixed << setprecision(2) << sample_N[n]
                << "," << fixed << setprecision(1) << sample_y[j]
                << "," << fixed << setprecision(3) << value
                << "," << wu_ou[n][j]
                << "," << yuen_takara[n][j]
                << "," << mondal_mishra[n][j]
                << "," << sun_zhang[n][j] << "\n";
            data_out << fixed << setprecision(2) << sample_N[n]
                     << "," << fixed << setprecision(1) << sample_y[j]
                     << "," << fixed << setprecision(3) << value
                     << "," << wu_ou[n][j]
                     << "," << yuen_takara[n][j]
                     << "," << mondal_mishra[n][j]
                     << "," << sun_zhang[n][j] << "\n";
            table_out << fixed << setprecision(2) << sample_N[n]
                      << "," << fixed << setprecision(1) << sample_y[j]
                      << "," << fixed << setprecision(3) << value
                      << "," << wu_ou[n][j]
                      << "," << yuen_takara[n][j]
                      << "," << mondal_mishra[n][j]
                      << "," << sun_zhang[n][j] << "\n";
            table_data << fixed << setprecision(2) << sample_N[n]
                       << "," << fixed << setprecision(1) << sample_y[j]
                       << "," << fixed << setprecision(3) << value
                       << "," << wu_ou[n][j]
                       << "," << yuen_takara[n][j]
                       << "," << mondal_mishra[n][j]
                       << "," << sun_zhang[n][j] << "\n";
            tex_out << fixed << setprecision(2) << sample_N[n]
                    << " & " << fixed << setprecision(1) << sample_y[j]
                    << " & " << fixed << setprecision(3) << value
                    << " & " << wu_ou[n][j]
                    << " & " << yuen_takara[n][j]
                    << " & " << mondal_mishra[n][j]
                    << " & " << sun_zhang[n][j] << " \\\\\n";
            tex_data << fixed << setprecision(2) << sample_N[n]
                     << " & " << fixed << setprecision(1) << sample_y[j]
                     << " & " << fixed << setprecision(3) << value
                     << " & " << wu_ou[n][j]
                     << " & " << yuen_takara[n][j]
                     << " & " << mondal_mishra[n][j]
                     << " & " << sun_zhang[n][j] << " \\\\\n";
        }
    }
    cout << "wrote table2.csv, table2_dugks.csv and table2_dugks_rows.tex\n";
}

void write_temp_csv(const string& path) {
    ofstream out(path);
    out << "y";
    for (int k = 0; k < NTIME; ++k) {
        ostringstream label;
        label << defaultfloat << setprecision(12) << save_times[k];
        out << ",t=" << label.str();
    }
    out << ",Steady\n";

    for (int n = 0; n < NOUT; ++n) {
        out << fixed << setprecision(8)
            << static_cast<double>(n) / (NOUT - 1);
        for (int k = 0; k < NCURVE; ++k) {
            out << "," << scientific << setprecision(12)
                << saved_temp[k][n];
        }
        out << "\n";
    }
}

void write_fig5_csv(const string& path, double data[CASE_out][NOUT]) {
    ofstream out(path);
    out << "y";
    for (int p = 0; p < CASE_out; ++p) {
        out << ",phi=" << phi_values[p];
    }
    out << "\n";

    for (int n = 0; n < NOUT; ++n) {
        out << fixed << setprecision(8)
            << static_cast<double>(n) / (NOUT - 1);
        for (int p = 0; p < CASE_out; ++p) {
            out << "," << scientific << setprecision(12) << data[p][n];
        }
        out << "\n";
    }
}
