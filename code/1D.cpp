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
const int NX = 100;
const int Q = 3;
const int M = 16;
const int NTIME = 4;
const int NCURVE = NTIME + 1;
const int NOUT = NX + 1;
const double CFL_T = 0.50;
const double CFL_R = 0.85;
const double FIG4_TIME = 0.05;
const double STEADY_error = 2.0e-8;
const int STEADY_MAX_STEPS = 40000;
const int RAD_MAX_STEPS = 20000;
const double RAD_error = 2.0e-9;
const int RAD_TRACKING_STEPS = 3;
const double dx = 1.0 / NX;
// Thermal DUGKS model parameters.
// Following the reference DUGKS cavity code, RT is fixed at 1/3 so the
// discrete particle speed is c=sqrt(3RT)=1.  The diffusivity is set by
// alpha*=tau_t*RT=1, independent of mesh spacing.
const double RT = 1.0 / 3.0;
const double tau_t = 1.0 / RT;
const double speed_t = 1.7320508075688772935 * sqrt(RT);

const double save_times[NTIME] = {0.0025, 0.005, 0.015, 0.04};
const double c_t[Q] = {0.0, speed_t, -speed_t};
const double w_t[Q] = {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0};
const int opposite[Q] = {0, 2, 1};

// === 场变量数组 ===
double T[NX + 2];                              // 温度
double pre_temp[NX + 2];                       // 上一时间步温度
double G[NX + 2];                              // 入射辐射
double G_pre[NX + 2];                          // 上一次辐射迭代的入射辐射
double Qr[NX + 2];                             // 辐射热流

double g_tilde[NX + 2][Q];                     // 热传导分布函数
double g_bar_plus[NX + 2][Q];                  // 热传导半时间步碰撞后分布函数
double g_tilde_plus[NX + 2][Q];                // 热传导整时间步碰撞后分布函数
double g_face[NX + 1][Q];                      // 热传导界面分布函数

double mu[M];                                  // 离散方向余弦
double weights[M];                             // 离散角度权重
double flux_weights[M];                        // 辐射热流积分权重
double I_tilde[NX + 2][M];                     // 辐射强度分布函数
double I_bar_plus[NX + 2][M];                  // 辐射半时间步碰撞后分布函数
double I_tilde_plus[NX + 2][M];                // 辐射整时间步碰撞后分布函数
double I_face[NX + 1][M];                      // 辐射界面分布函数

// === 保存结果数组 ===
double saved_temp[NCURVE][NOUT];
double Qc_out[NOUT];
double Qr_out[NOUT];
double Qt_out[NOUT];

// === 运行状态 ===
bool radiation_needs_initial_convergence = true;

// === 函数声明 ===
double minmod(double a, double b);
double pow4(double value);
void ini_angles();
void ini_case(double TE);
void rec_temperture();
void temperture(double dt, double TE);
void ini_radiation();
void rec_radiation(double omega, double pseudo_dt);
int radiation_dugks(double omega, double TE,
                    bool force_convergence = false);
void Energy_equation(double dt, double N, double omega);
void coupled_step(double dt, double N, double omega, double TE);
double sample_cell_field(const double field[NX + 2], double x,
                         double left_value, double right_value);
void save_curve(int curve_id, double TE);
void run_temp_case(const string& fig, double N, double omega,
                          double TE);
void run_flux_case(const string& fig, double N, double omega,
                   double TE);
void run_sun2016_table1_case();
void write_temp_csv(const string& path);
void write_flux_csv(const string& path);

int main() {
#ifdef _WIN32
    mkdir("output");
#else
    mkdir("output", 0755);
#endif
    ini_angles();
    cout << "T: D1Q3 DUGKS, radiation: 16-angle DUGKS\n";

    run_temp_case("2a", 0.01, 0.0, 0.5);
    run_temp_case("2b", 0.01, 0.5, 0.5);
    run_temp_case("2c", 0.01, 0.9, 0.5);

    run_temp_case("3a", 0.01, 0.0, 0.5);
    run_temp_case("3b", 0.10, 0.0, 0.5);
    run_temp_case("3c", 1.00, 0.0, 0.5);

    run_flux_case("4a", 0.10, 0.0, 0.1);
    run_flux_case("4b", 0.10, 0.5, 0.1);
    run_flux_case("4c", 0.10, 1.0, 0.1);
    run_sun2016_table1_case();

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

void ini_angles() {
    const double dtheta = PI / M;
    for (int m = 0; m < M; ++m) {
        const double theta_l = m * dtheta;
        const double theta_r = (m + 1) * dtheta;
        const double cos_l = cos(theta_l);
        const double cos_r = cos(theta_r);

        weights[m] = 0.5 * (cos_l - cos_r);
        mu[m] = 0.5 * (cos_l + cos_r);
        flux_weights[m] =
            0.25 * (sin(theta_r) * sin(theta_r)
                  - sin(theta_l) * sin(theta_l));
    }
}

void ini_case(double TE) {
    for (int i = 1; i <= NX; ++i) {
        T[i] = TE;
        for (int q = 0; q < Q; ++q) {
            g_tilde[i][q] = w_t[q] * TE;
        }
    }
    T[0] = TW;
    T[NX + 1] = TE;
    radiation_needs_initial_convergence = true;
    ini_radiation();
}

void rec_temperture() {
    for (int i = 1; i <= NX; ++i) {
        double sum = 0.0;
        for (int q = 0; q < Q; ++q) sum += g_tilde[i][q];
        T[i] = sum;
    }
}

void temperture(double dt, double TE) {
    const double h = 0.5 * dt;
    const double coefficient_old = (2.0 * tau_t - h) / (2.0 * tau_t + dt);
    const double coefficient_eq = 3.0 * h / (2.0 * tau_t + dt);
    const double face_old = 2.0 * tau_t / (2.0 * tau_t + h);
    const double face_eq = h / (2.0 * tau_t + h);
    double slope[NX + 2][Q];

    rec_temperture();

    for (int i = 1; i <= NX; ++i) {
        for (int q = 0; q < Q; ++q) {
            const double eq = w_t[q] * T[i];
            g_bar_plus[i][q] =
                coefficient_old * g_tilde[i][q] + coefficient_eq * eq;
            g_tilde_plus[i][q] =
                (4.0 / 3.0) * g_bar_plus[i][q]
              - (1.0 / 3.0) * g_tilde[i][q];
        }
    }

    for (int q = 0; q < Q; ++q) {
        slope[1][q] = 0.0;
        slope[NX][q] = 0.0;
        for (int i = 2; i < NX; ++i) {
            const double left = (g_bar_plus[i][q] - g_bar_plus[i - 1][q]) / dx;
            const double right = (g_bar_plus[i + 1][q] - g_bar_plus[i][q]) / dx;
            slope[i][q] = minmod(left, right);
        }
    }

    for (int face = 0; face <= NX; ++face) {
        double bar_value[Q];
        for (int q = 0; q < Q; ++q) {
            if (c_t[q] > 0.0) {
                const int cell = max(1, face);
                const double cell_x = (cell - 0.5) * dx;
                const double departure = face * dx - c_t[q] * h;
                bar_value[q] =
                    g_bar_plus[cell][q] + slope[cell][q] * (departure - cell_x);
            } else if (c_t[q] < 0.0) {
                const int cell = min(NX, face + 1);
                const double cell_x = (cell - 0.5) * dx;
                const double departure = face * dx - c_t[q] * h;
                bar_value[q] =
                    g_bar_plus[cell][q] + slope[cell][q] * (departure - cell_x);
            } else if (face == 0) {
                bar_value[q] = w_t[q] * TW;
            } else if (face == NX) {
                bar_value[q] = w_t[q] * TE;
            } else {
                bar_value[q] =
                    0.5 * (g_bar_plus[face][q] + g_bar_plus[face + 1][q]);
            }
        }

        if (face == 0) {
            bar_value[0] = w_t[0] * TW;
            bar_value[1] = -bar_value[opposite[1]] + 2.0 * w_t[1] * TW;
        } else if (face == NX) {
            bar_value[0] = w_t[0] * TE;
            bar_value[2] =
                -bar_value[opposite[2]] + 2.0 * w_t[2] * TE;
        }

        double face_temp = 0.0;
        for (int q = 0; q < Q; ++q) face_temp += bar_value[q];
        if (face == 0) face_temp = TW;
        if (face == NX) face_temp = TE;

        for (int q = 0; q < Q; ++q) {
            g_face[face][q] =
                face_old * bar_value[q]
              + face_eq * w_t[q] * face_temp;
        }
    }

    for (int i = 1; i <= NX; ++i) {
        for (int q = 0; q < Q; ++q) {
            g_tilde[i][q] = g_tilde_plus[i][q]
                - (dt / dx) * c_t[q] * (g_face[i][q] - g_face[i - 1][q]);
        }
    }
    rec_temperture();
    T[0] = TW;
    T[NX + 1] = TE;
}

void ini_radiation() {
    for (int i = 1; i <= NX; ++i) {
        const double eq = pow4(T[i]);
        G[i] = eq;
        Qr[i] = 0.0;
        for (int m = 0; m < M; ++m) I_tilde[i][m] = eq;
    }
}

void rec_radiation(double omega, double pseudo_dt) {
    const double moment_coefficient = 0.5 * pseudo_dt * (1.0 - omega);
    const double physical_old = 2.0 / (2.0 + pseudo_dt);
    const double physical_eq = pseudo_dt / (2.0 + pseudo_dt);

    for (int i = 1; i <= NX; ++i) {
        double tilde_moment = 0.0;
        for (int m = 0; m < M; ++m) {
            tilde_moment += weights[m] * I_tilde[i][m];
        }
        const double T4 = pow4(T[i]);
        G[i] = (tilde_moment + moment_coefficient * T4)
                 / (1.0 + moment_coefficient);
        const double source = (1.0 - omega) * T4 + omega * G[i];

        Qr[i] = 0.0;
        for (int m = 0; m < M; ++m) {
            const double physical =
                physical_old * I_tilde[i][m] + physical_eq * source;
            Qr[i] += flux_weights[m] * physical;
        }
    }
}

int radiation_dugks(double omega, double TE,
                    bool force_convergence) {
    double maximum_mu = 0.0;
    for (int m = 0; m < M; ++m) maximum_mu = max(maximum_mu, fabs(mu[m]));
    const double pseudo_dt = CFL_R * dx / maximum_mu;
    const double h = 0.5 * pseudo_dt;
    const double coefficient_old = (2.0 - h) / (2.0 + pseudo_dt);
    const double coefficient_eq = 3.0 * h / (2.0 + pseudo_dt);
    const double face_old = 2.0 / (2.0 + h);
    const double face_eq = h / (2.0 + h);
    const double face_moment_coefficient = 0.5 * h * (1.0 - omega);
    double slope[NX + 2][M];

    const bool converge =
        force_convergence || radiation_needs_initial_convergence;
    const int maximum_steps = converge ? RAD_MAX_STEPS : RAD_TRACKING_STEPS;

    for (int step = 1; step <= maximum_steps; ++step) {
        rec_radiation(omega, pseudo_dt);
        for (int i = 1; i <= NX; ++i) G_pre[i] = G[i];

        for (int i = 1; i <= NX; ++i) {
            const double T4 = pow4(T[i]);
            const double source = (1.0 - omega) * T4 + omega * G[i];
            for (int m = 0; m < M; ++m) {
                I_bar_plus[i][m] =
                    coefficient_old * I_tilde[i][m] + coefficient_eq * source;
                I_tilde_plus[i][m] =
                    (4.0 / 3.0) * I_bar_plus[i][m]
                  - (1.0 / 3.0) * I_tilde[i][m];
            }
        }

        for (int m = 0; m < M; ++m) {
            slope[1][m] = 0.0;
            slope[NX][m] = 0.0;
            for (int i = 2; i < NX; ++i) {
                const double left =
                    (I_bar_plus[i][m] - I_bar_plus[i - 1][m]) / dx;
                const double right =
                    (I_bar_plus[i + 1][m] - I_bar_plus[i][m]) / dx;
                slope[i][m] = minmod(left, right);
            }
        }

        for (int face = 0; face <= NX; ++face) {
            double bar_value[M];
            for (int m = 0; m < M; ++m) {
                if (mu[m] > 0.0) {
                    const int cell = max(1, face);
                    const double cell_x = (cell - 0.5) * dx;
                    const double departure = face * dx - mu[m] * h;
                    bar_value[m] = I_bar_plus[cell][m]
                        + slope[cell][m] * (departure - cell_x);
                } else {
                    const int cell = min(NX, face + 1);
                    const double cell_x = (cell - 0.5) * dx;
                    const double departure = face * dx - mu[m] * h;
                    bar_value[m] = I_bar_plus[cell][m]
                        + slope[cell][m] * (departure - cell_x);
                }
            }

            if (face == 0 || face == NX) {
                const int cell = face == 0 ? 1 : NX;
                const double T4 = pow4(T[cell]);
                const double source =
                    (1.0 - omega) * T4 + omega * G[cell];
                for (int m = 0; m < M; ++m) {
                    const bool incoming =
                        (face == 0 && mu[m] > 0.0)
                     || (face == NX && mu[m] < 0.0);
                    if (incoming) {
                        I_face[face][m] =
                            face == 0 ? pow4(TW) : pow4(TE);
                    } else {
                        I_face[face][m] =
                            face_old * bar_value[m] + face_eq * source;
                    }
                }
            } else {
                double bar_moment = 0.0;
                for (int m = 0; m < M; ++m) {
                    bar_moment += weights[m] * bar_value[m];
                }
                const double face_T =
                    0.5 * (T[face] + T[face + 1]);
                const double face_T4 = pow4(face_T);
                const double face_G =
                    (bar_moment + face_moment_coefficient * face_T4)
                  / (1.0 + face_moment_coefficient);
                const double source =
                    (1.0 - omega) * face_T4 + omega * face_G;

                for (int m = 0; m < M; ++m) {
                    I_face[face][m] =
                        face_old * bar_value[m] + face_eq * source;
                }
            }
        }

        for (int i = 1; i <= NX; ++i) {
            for (int m = 0; m < M; ++m) {
                I_tilde[i][m] = I_tilde_plus[i][m]
                    - (pseudo_dt / dx) * mu[m]
                    * (I_face[i][m] - I_face[i - 1][m]);
            }
        }

        rec_radiation(omega, pseudo_dt);
        double error = 0.0;
        for (int i = 1; i <= NX; ++i) {
            error = max(error, fabs(G[i] - G_pre[i]));
        }
        if (converge && step > 20 && error < RAD_error) {
            radiation_needs_initial_convergence = false;
            return step;
        }
    }
    radiation_needs_initial_convergence = false;
    return maximum_steps;
}

void Energy_equation(double dt, double N, double omega) {
    const double A = (1.0 - omega) / max(N, 1.0e-14);

    for (int i = 1; i <= NX; ++i) {
        const double source = A * (G[i] - pow4(T[i]));
        const double delta = dt * source;
        T[i] += delta;
        for (int q = 0; q < Q; ++q) g_tilde[i][q] += w_t[q] * delta;
    }
}

void coupled_step(double dt, double N, double omega, double TE) {
    radiation_dugks(omega, TE, false);
    Energy_equation(0.5 * dt, N, omega);
    temperture(dt, TE);
    radiation_dugks(omega, TE, false);
    Energy_equation(0.5 * dt, N, omega);
    rec_temperture();
}

double sample_cell_field(const double field[NX + 2], double x,
                         double left_value, double right_value) {
    if (x <= 0.0) return left_value;
    if (x >= 1.0) return right_value;

    const double center_index = x / dx + 0.5;
    int left_cell = static_cast<int>(floor(center_index));
    if (left_cell < 1) {
        const double x1 = 0.5 * dx;
        return left_value + (field[1] - left_value) * x / x1;
    }
    if (left_cell >= NX) {
        const double x0 = (NX - 0.5) * dx;
        return field[NX]
            + (right_value - field[NX]) * (x - x0) / (1.0 - x0);
    }

    const double x0 = (left_cell - 0.5) * dx;
    return field[left_cell]
        + (field[left_cell + 1] - field[left_cell]) * (x - x0) / dx;
}

void save_curve(int curve_id, double TE) {
    for (int n = 0; n < NOUT; ++n) {
        const double x = static_cast<double>(n) / NX;
        saved_temp[curve_id][n] =
            sample_cell_field(T, x, TW, TE);
    }
}

void run_temp_case(const string& fig, double N, double omega,
                          double TE) {
    ini_case(TE);
    const double base_dt = CFL_T * dx / speed_t;
    double time = 0.0;

    for (int k = 0; k < NTIME; ++k) {
        while (time + 1.0e-14 < save_times[k]) {
            const double dt = min(base_dt, save_times[k] - time);
            coupled_step(dt, N, omega, TE);
            time += dt;
        }
        save_curve(k, TE);
    }

    int steady_step = 0;
    for (; steady_step < STEADY_MAX_STEPS; ++steady_step) {
        for (int i = 1; i <= NX; ++i) pre_temp[i] = T[i];
        coupled_step(base_dt, N, omega, TE);

        double error = 0.0;
        for (int i = 1; i <= NX; ++i) {
            error = max(error, fabs(T[i] - pre_temp[i]));
        }
        if (steady_step > 100 && error < STEADY_error) break;
    }

    save_curve(NTIME, TE);
    write_temp_csv("output/" + fig + ".csv");
    cout << "wrote " << fig << ".csv, steady steps = " << steady_step << "\n";
}

void run_flux_case(const string& fig, double N, double omega,
                   double TE) {
    ini_case(TE);
    const double base_dt = CFL_T * dx / speed_t;
    double time = 0.0;

    while (time + 1.0e-14 < FIG4_TIME) {
        const double dt = min(base_dt, FIG4_TIME - time);
        coupled_step(dt, N, omega, TE);
        time += dt;
    }
    radiation_dugks(omega, TE, true);

    double conductive_cell[NX + 2];
    conductive_cell[1] =
        -(T[2] - TW) / (1.5 * dx);
    conductive_cell[NX] =
        -(TE - T[NX - 1]) / (1.5 * dx);
    for (int i = 2; i < NX; ++i) {
        conductive_cell[i] =
            -(T[i + 1] - T[i - 1]) / (2.0 * dx);
    }

    for (int n = 0; n < NOUT; ++n) {
        const double x = static_cast<double>(n) / NX;
        if (n == 0) {
            Qc_out[n] = -(T[1] - TW) / (0.5 * dx);
        } else if (n == NX) {
            Qc_out[n] =
                -(TE - T[NX]) / (0.5 * dx);
        } else {
            Qc_out[n] =
                sample_cell_field(conductive_cell, x,
                                  conductive_cell[1], conductive_cell[NX]);
        }
        Qr_out[n] =
            sample_cell_field(Qr, x, Qr[1], Qr[NX])
            / max(N, 1.0e-14);
        Qt_out[n] = Qc_out[n] + Qr_out[n];
    }

    write_flux_csv("output/" + fig + ".csv");
    cout << "wrote " << fig << ".csv\n";
}

void run_sun2016_table1_case() {
    const double N = 0.10;
    const double omega = 0.50;
    const double TE = 0.0;
    const double sample_time = 0.05;
    const double sample_x[3] = {0.25, 0.50, 0.75};
    const double sutton[3] = {0.4888, 0.1778, 0.0591};
    const double tsai_lin[3] = {0.4889, 0.1773, 0.0588};
    const double talukdar_mishra[3] = {0.4892, 0.1768, 0.0585};
    const double mishra_roy[3] = {0.4897, 0.1771, 0.0587};
    const double mondal_mishra[3] = {0.4898, 0.1769, 0.0583};
    const double sun_zhang[3] = {0.4894, 0.1771, 0.0585};

    ini_case(TE);
    const double base_dt = CFL_T * dx / speed_t;
    double time = 0.0;

    while (time + 1.0e-14 < sample_time) {
        const double dt = min(base_dt, sample_time - time);
        coupled_step(dt, N, omega, TE);
        time += dt;
    }
    radiation_dugks(omega, TE, true);

    ofstream out("output/table1_dugks.csv");
    ofstream data_out("data/table1_dugks.csv");
    ofstream table_out("output/table1.csv");
    ofstream table_data("data/table1.csv");
    ofstream tex_out("output/table1_dugks_rows.tex");
    ofstream tex_data("data/table1_dugks_rows.tex");
    out << "x,DUGKS,Sutton1990,TsaiLin1994,TalukdarMishra2001,"
        << "MishraRoy2007,MondalMishra2007,SunZhang2016\n";
    data_out << "x,DUGKS,Sutton1990,TsaiLin1994,TalukdarMishra2001,"
             << "MishraRoy2007,MondalMishra2007,SunZhang2016\n";
    table_out << "x,DUGKS,Sutton1990,TsaiLin1994,TalukdarMishra2001,"
              << "MishraRoy2007,MondalMishra2007,SunZhang2016\n";
    table_data << "x,DUGKS,Sutton1990,TsaiLin1994,TalukdarMishra2001,"
               << "MishraRoy2007,MondalMishra2007,SunZhang2016\n";
    for (int n = 0; n < 3; ++n) {
        const double value = sample_cell_field(T, sample_x[n], TW, TE);
        out << fixed << setprecision(2) << sample_x[n]
            << "," << fixed << setprecision(4) << value
            << "," << sutton[n]
            << "," << tsai_lin[n]
            << "," << talukdar_mishra[n]
            << "," << mishra_roy[n]
            << "," << mondal_mishra[n]
            << "," << sun_zhang[n] << "\n";
        data_out << fixed << setprecision(2) << sample_x[n]
                 << "," << fixed << setprecision(4) << value
                 << "," << sutton[n]
                 << "," << tsai_lin[n]
                 << "," << talukdar_mishra[n]
                 << "," << mishra_roy[n]
                 << "," << mondal_mishra[n]
                 << "," << sun_zhang[n] << "\n";
        table_out << fixed << setprecision(2) << sample_x[n]
                  << "," << fixed << setprecision(4) << value
                  << "," << sutton[n]
                  << "," << tsai_lin[n]
                  << "," << talukdar_mishra[n]
                  << "," << mishra_roy[n]
                  << "," << mondal_mishra[n]
                  << "," << sun_zhang[n] << "\n";
        table_data << fixed << setprecision(2) << sample_x[n]
                   << "," << fixed << setprecision(4) << value
                   << "," << sutton[n]
                   << "," << tsai_lin[n]
                   << "," << talukdar_mishra[n]
                   << "," << mishra_roy[n]
                   << "," << mondal_mishra[n]
                   << "," << sun_zhang[n] << "\n";
        tex_out << fixed << setprecision(2) << sample_x[n]
                << " & " << fixed << setprecision(4) << value
                << " & " << sutton[n]
                << " & " << tsai_lin[n]
                << " & " << talukdar_mishra[n]
                << " & " << mishra_roy[n]
                << " & " << mondal_mishra[n]
                << " & " << sun_zhang[n] << " \\\\\n";
        tex_data << fixed << setprecision(2) << sample_x[n]
                 << " & " << fixed << setprecision(4) << value
                 << " & " << sutton[n]
                 << " & " << tsai_lin[n]
                 << " & " << talukdar_mishra[n]
                 << " & " << mishra_roy[n]
                 << " & " << mondal_mishra[n]
                 << " & " << sun_zhang[n] << " \\\\\n";
    }
    cout << "wrote table1.csv, table1_dugks.csv and table1_dugks_rows.tex\n";
}

void write_temp_csv(const string& path) {
    ofstream out(path);
    out << "x";
    for (int k = 0; k < NTIME; ++k) {
        ostringstream label;
        label << defaultfloat << setprecision(12) << save_times[k];
        out << ",t=" << label.str();
    }
    out << ",Steady\n";

    for (int n = 0; n < NOUT; ++n) {
        out << fixed << setprecision(8) << static_cast<double>(n) / NX;
        for (int k = 0; k < NCURVE; ++k) {
            out << "," << scientific << setprecision(12)
                << saved_temp[k][n];
        }
        out << "\n";
    }
}

void write_flux_csv(const string& path) {
    ofstream out(path);
    out << "x,Qr,Qc,Qt\n";
    for (int n = 0; n < NOUT; ++n) {
        out << fixed << setprecision(8) << static_cast<double>(n) / NX
            << "," << scientific << setprecision(12) << Qr_out[n]
            << "," << Qc_out[n]
            << "," << Qt_out[n] << "\n";
    }
}
