#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _WIN32
#include <direct.h>
#endif

using namespace std;

const double PI = 3.14159265358979323846;
const double TW = 1.0;
const double TE = 0.5;
const int NX = 21;
const int NY = 21;
const int NZ = 21;
const int Q = 27;
const int NTHETA = 8;
const int MAX_PHI = 16;
const int MAX_M = NTHETA * MAX_PHI;
const int M = MAX_M;
const int MAX_TRANSIENT_TIMES = 3;
const int MAX_CURVES = MAX_TRANSIENT_TIMES + 1;
const int NOUT = 41;
const double CFL_T = 0.50;
const double CFL_R = 0.80;
const double STEADY_error = 5.0e-7;
const double RAD_error = 2.0e-7;
const int RAD_TRACKING_STEPS = 30;
const int RAD_COUPLING_INTERVAL = 4;
const double RAD_SLOPE_FACTOR = 0.0;
#ifdef SMOKE_TEST
const int RAD_MAX_STEPS = 300;
const int STEADY_MAX_STEPS = 100;
#else
const int RAD_MAX_STEPS = 12000;
const int STEADY_MAX_STEPS = 30000;
#endif

const double dx = 1.0 / NX;
const double dy = 1.0 / NY;
const double dz = 1.0 / NZ;
const double tau_t = dx * dx;
const double RT = 1.0 / tau_t;
const double speed_t = sqrt(3.0 * RT);
struct Dir {
    double x, y, z, weight;
};

Dir directions[MAX_M];

// D3Q27 velocity set
double cx[Q], cy[Q], cz[Q], wt[Q];
int opposite[Q];

// === 场变量数组 ===
double T[NZ][NY][NX];                         // 温度
double pre_temp[NZ][NY][NX];                  // 上一时间步温度
double G[NZ][NY][NX];                         // 入射辐射
double G_pre[NZ][NY][NX];                     // 上一次辐射迭代的入射辐射
double Qr[NZ][NY][NX];                        // 辐射热流

double g_tilde[NZ][NY][NX][Q];                // 热传导分布函数
double g_bar_plus[NZ][NY][NX][Q];             // 热传导半时间步碰撞后分布函数
double g_tilde_plus[NZ][NY][NX][Q];           // 热传导整时间步碰撞后分布函数
double g_slope_x[NZ][NY][NX][Q];              // 热传导分布函数的 x 向斜率
double g_slope_y[NZ][NY][NX][Q];              // 热传导分布函数的 y 向斜率
double g_slope_z[NZ][NY][NX][Q];              // 热传导分布函数的 z 向斜率
double rate_t[NZ][NY][NX][Q];                 // 热传导通量变化率

double I_tilde[NZ][NY][NX][M];                // 辐射强度分布函数
double I_bar_plus[NZ][NY][NX][M];             // 辐射半时间步碰撞后分布函数
double I_tilde_plus[NZ][NY][NX][M];           // 辐射整时间步碰撞后分布函数
double I_slope_x[NZ][NY][NX][M];              // 辐射强度的 x 向斜率
double I_slope_y[NZ][NY][NX][M];              // 辐射强度的 y 向斜率
double I_slope_z[NZ][NY][NX][M];              // 辐射强度的 z 向斜率
double rate_r[NZ][NY][NX][M];                 // 辐射通量变化率

// === 保存结果数组 ===
double saved_temp[MAX_CURVES][NOUT];
double saved_fig3[MAX_CURVES][NZ][NY];

// === 运行状态 ===
int thermal_step_counter = 0;
bool radiation_needs_initial_convergence = true;
double EW = 1.0;
double TAU_L = 1.0;

inline double pow4(double v) {
    double s = v * v; return s * s;
}
inline double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return fabs(a) < fabs(b) ? a : b;
}

void ini_velocities() {
    const double w1[3] = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
    const int v1[3] = {-1, 0, 1};
    int q = 0;
    for (int kz = 0; kz < 3; ++kz)
        for (int jy = 0; jy < 3; ++jy)
            for (int ix = 0; ix < 3; ++ix) {
                cx[q] = speed_t * v1[ix];
                cy[q] = speed_t * v1[jy];
                cz[q] = speed_t * v1[kz];
                wt[q] = w1[ix] * w1[jy] * w1[kz];
                const int ox = 2 - ix, oy = 2 - jy, oz = 2 - kz;
                opposite[q] = (oz * 3 + oy) * 3 + ox;
                ++q;
            }
}

void ini_angles() {
    const double dtheta = PI / NTHETA;
    const double dphi = 2.0 * PI / MAX_PHI;
    int m = 0;
    for (int it = 0; it < NTHETA; ++it) {
        const double tl = it * dtheta, tr = (it + 1) * dtheta;
        const double txy = 0.5 * (tr - tl)
            - 0.25 * (sin(2.0 * tr) - sin(2.0 * tl));
        const double tz = 0.5 * (sin(tr) * sin(tr) - sin(tl) * sin(tl));
        for (int ip = 0; ip < MAX_PHI; ++ip) {
            const double pl = ip * dphi, pr = (ip + 1) * dphi;
            const double sa = (cos(tl) - cos(tr)) * dphi;
            directions[m].x = txy * (sin(pr) - sin(pl)) / sa;
            directions[m].y = txy * (cos(pl) - cos(pr)) / sa;
            directions[m].z = tz * dphi / sa;
            directions[m].weight = sa / (4.0 * PI);
            ++m;
        }
    }
}

void ini_case() {
    double initial_I = pow4(TE);

    for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i) {
                T[k][j][i] = TE;
                pre_temp[k][j][i] = TE;
                G[k][j][i] = initial_I;
                G_pre[k][j][i] = initial_I;
                Qr[k][j][i] = 0.0;

                for (int q = 0; q < Q; ++q) {
                    g_tilde[k][j][i][q] = wt[q] * TE;
                    g_bar_plus[k][j][i][q] = 0.0;
                    g_tilde_plus[k][j][i][q] = 0.0;
                    g_slope_x[k][j][i][q] = 0.0;
                    g_slope_y[k][j][i][q] = 0.0;
                    g_slope_z[k][j][i][q] = 0.0;
                    rate_t[k][j][i][q] = 0.0;
                }

                for (int m = 0; m < M; ++m) {
                    I_tilde[k][j][i][m] = initial_I;
                    I_bar_plus[k][j][i][m] = initial_I;
                    I_tilde_plus[k][j][i][m] = initial_I;
                    I_slope_x[k][j][i][m] = 0.0;
                    I_slope_y[k][j][i][m] = 0.0;
                    I_slope_z[k][j][i][m] = 0.0;
                    rate_r[k][j][i][m] = 0.0;
                }
            }
}

double wall_T(int axis, int face) {
    return (axis == 1 && face == 0) ? TW : TE;
}

void rec_temperature() {
    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i) {
                T[k][j][i] = 0.0;
                for (int q = 0; q < Q; ++q)
                    T[k][j][i] += g_tilde[k][j][i][q];
            }
}

void calc_T_slopes() {
    for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i)
                for (int q = 0; q < Q; ++q) {
                    g_slope_x[k][j][i][q] = 0.0;
                    g_slope_y[k][j][i][q] = 0.0;
                    g_slope_z[k][j][i][q] = 0.0;
                }

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 1; k < NZ - 1; ++k)
        for (int j = 1; j < NY - 1; ++j)
            for (int i = 1; i < NX - 1; ++i) {
                for (int q = 0; q < Q; ++q) {
                    g_slope_x[k][j][i][q] = minmod(
                        (g_bar_plus[k][j][i][q] - g_bar_plus[k][j][i-1][q]) / dx,
                        (g_bar_plus[k][j][i+1][q] - g_bar_plus[k][j][i][q]) / dx);
                    g_slope_y[k][j][i][q] = minmod(
                        (g_bar_plus[k][j][i][q] - g_bar_plus[k][j-1][i][q]) / dy,
                        (g_bar_plus[k][j+1][i][q] - g_bar_plus[k][j][i][q]) / dy);
                    g_slope_z[k][j][i][q] = minmod(
                        (g_bar_plus[k][j][i][q] - g_bar_plus[k-1][j][i][q]) / dz,
                        (g_bar_plus[k+1][j][i][q] - g_bar_plus[k][j][i][q]) / dz);
                }
            }
}

void thermal_face(
    int axis, int a, int b, int face, double h,
    double face_old, double face_eq, double (&rate)[NZ][NY][NX][Q]
) {
    double bar[27], phys[27];
    const int cells = (axis == 0) ? NX : ((axis == 1) ? NY : NZ);

    for (int q = 0; q < Q; ++q) {
        double vel = (axis == 0) ? cx[q] : ((axis == 1) ? cy[q] : cz[q]);
        int up = (vel > 0.0) ? max(0, face - 1) : ((vel < 0.0) ? min(cells - 1, face) : min(cells - 1, max(0, face)));
        int ii = a, jj = a, kk = b;
        if (axis == 0) { ii = up; jj = a; kk = b; }
        else if (axis == 1) { ii = a; jj = up; kk = b; }
        else { ii = a; jj = b; kk = up; }

        double ctr = (up + 0.5) * dx;
        double dep = face * dx - vel * h - ctr;
        bar[q] = g_bar_plus[kk][jj][ii][q]
            + ((axis == 0) ? g_slope_x[kk][jj][ii][q]
               : (axis == 1) ? g_slope_y[kk][jj][ii][q]
                             : g_slope_z[kk][jj][ii][q]) * dep
            - ((axis == 0) ? g_slope_y[kk][jj][ii][q] * cy[q]
                                + g_slope_z[kk][jj][ii][q] * cz[q]
               : (axis == 1) ? g_slope_x[kk][jj][ii][q] * cx[q]
                                + g_slope_z[kk][jj][ii][q] * cz[q]
                             : g_slope_x[kk][jj][ii][q] * cx[q]
                                + g_slope_y[kk][jj][ii][q] * cy[q]) * h;
    }

    double fT = 0.0;
    bool bnd = (face == 0 || face == cells);
    if (bnd) {
        fT = wall_T(axis, face);
        for (int q = 0; q < Q; ++q) {
            double vel = (axis == 0) ? cx[q] : ((axis == 1) ? cy[q] : cz[q]);
            bool inc = (face == 0 && vel > 0.0) || (face == cells && vel < 0.0);
            if (inc) bar[q] = -bar[opposite[q]] + 2.0 * wt[q] * fT;
            else if (vel == 0.0) bar[q] = wt[q] * fT;
        }
    } else {
        for (int q = 0; q < Q; ++q) fT += bar[q];
    }

    for (int q = 0; q < Q; ++q)
        phys[q] = face_old * bar[q] + face_eq * wt[q] * fT;

    for (int q = 0; q < Q; ++q) {
        double vel = (axis == 0) ? cx[q] : ((axis == 1) ? cy[q] : cz[q]);
        double flux = (1.0 / dx) * vel * phys[q];
        if (face > 0) {
            int ii, jj, kk;
            if (axis == 0) { ii = face - 1; jj = a; kk = b; }
            else if (axis == 1) { ii = a; jj = face - 1; kk = b; }
            else { ii = a; jj = b; kk = face - 1; }
            rate[kk][jj][ii][q] -= flux;
        }
        if (face < cells) {
            int ii, jj, kk;
            if (axis == 0) { ii = face; jj = a; kk = b; }
            else if (axis == 1) { ii = a; jj = face; kk = b; }
            else { ii = a; jj = b; kk = face; }
            rate[kk][jj][ii][q] += flux;
        }
    }
}

void temperture(double dt) {
    double h = 0.5 * dt;
    double co = (2.0 * tau_t - h) / (2.0 * tau_t + dt);
    double ce = 3.0 * h / (2.0 * tau_t + dt);
    double fo = 2.0 * tau_t / (2.0 * tau_t + h);
    double fe = h / (2.0 * tau_t + h);

    rec_temperature();
    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i)
                for (int q = 0; q < Q; ++q) {
                    g_bar_plus[k][j][i][q] =
                        co * g_tilde[k][j][i][q] + ce * wt[q] * T[k][j][i];
                    g_tilde_plus[k][j][i][q] =
                        (4.0 / 3.0) * g_bar_plus[k][j][i][q]
                        - (1.0 / 3.0) * g_tilde[k][j][i][q];
                }
    calc_T_slopes();

    for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i)
                for (int q = 0; q < Q; ++q)
                    rate_t[k][j][i][q] = 0.0;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
            for (int face = 0; face <= NX; ++face)
                thermal_face(0, j, k, face, h, fo, fe, rate_t);
    #pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < NZ; ++k)
        for (int i = 0; i < NX; ++i)
            for (int face = 0; face <= NY; ++face)
                thermal_face(1, i, k, face, h, fo, fe, rate_t);
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i)
            for (int face = 0; face <= NZ; ++face)
                thermal_face(2, i, j, face, h, fo, fe, rate_t);

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i)
                for (int q = 0; q < Q; ++q)
                    g_tilde[k][j][i][q] =
                        g_tilde_plus[k][j][i][q] + dt * rate_t[k][j][i][q];

    rec_temperature();
}

void calc_R_slopes() {
    for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i)
                for (int m = 0; m < M; ++m) {
                    I_slope_x[k][j][i][m] = 0.0;
                    I_slope_y[k][j][i][m] = 0.0;
                    I_slope_z[k][j][i][m] = 0.0;
                }

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 1; k < NZ - 1; ++k)
        for (int j = 1; j < NY - 1; ++j)
            for (int i = 1; i < NX - 1; ++i) {
                for (int m = 0; m < M; ++m) {
                    I_slope_x[k][j][i][m] = RAD_SLOPE_FACTOR * minmod(
                        (I_bar_plus[k][j][i][m] - I_bar_plus[k][j][i-1][m]) / dx,
                        (I_bar_plus[k][j][i+1][m] - I_bar_plus[k][j][i][m]) / dx);
                    I_slope_y[k][j][i][m] = RAD_SLOPE_FACTOR * minmod(
                        (I_bar_plus[k][j][i][m] - I_bar_plus[k][j-1][i][m]) / dy,
                        (I_bar_plus[k][j+1][i][m] - I_bar_plus[k][j][i][m]) / dy);
                    I_slope_z[k][j][i][m] = RAD_SLOPE_FACTOR * minmod(
                        (I_bar_plus[k][j][i][m] - I_bar_plus[k-1][j][i][m]) / dz,
                        (I_bar_plus[k+1][j][i][m] - I_bar_plus[k][j][i][m]) / dz);
                }
            }
}

void rec_radiation(double omega, double pdt) {
    double tr = 1.0 / TAU_L;
    double mc = 0.5 * pdt * (1.0 - omega) * TAU_L;
    double po = (2.0 * tr) / (2.0 * tr + pdt);
    double pe = pdt / (2.0 * tr + pdt);

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i) {
                double tm = 0.0;
                for (int m = 0; m < M; ++m)
                    tm += directions[m].weight * I_tilde[k][j][i][m];

                G[k][j][i] =
                    (tm + mc * pow4(T[k][j][i])) / (1.0 + mc);
                double src =
                    (1.0 - omega) * pow4(T[k][j][i]) + omega * G[k][j][i];

                Qr[k][j][i] = 0.0;
                for (int m = 0; m < M; ++m) {
                    double ph = po * I_tilde[k][j][i][m] + pe * src;
                    Qr[k][j][i] +=
                        directions[m].y * directions[m].weight * ph;
                }
            }
}

void radiation_face(
    int axis, int a, int b, int face, double omega, double h,
    double fo, double fe, double fmc, double (&rate)[NZ][NY][NX][M]
) {
    double bar[MAX_M], phys[MAX_M];
    const int cells = (axis == 0) ? NX : ((axis == 1) ? NY : NZ);

    for (int m = 0; m < M; ++m) {
        double vel = (axis == 0) ? directions[m].x : ((axis == 1) ? directions[m].y : directions[m].z);
        int up = (vel > 0.0) ? max(0, face - 1) : ((vel < 0.0) ? min(cells - 1, face) : min(cells - 1, max(0, face)));
        int ii = a, jj = a, kk = b;
        if (axis == 0) { ii = up; jj = a; kk = b; }
        else if (axis == 1) { ii = a; jj = up; kk = b; }
        else { ii = a; jj = b; kk = up; }

        double ctr = (up + 0.5) * dx;
        double dep = face * dx - vel * h - ctr;
        bar[m] = I_bar_plus[kk][jj][ii][m]
            + ((axis == 0) ? I_slope_x[kk][jj][ii][m]
               : (axis == 1) ? I_slope_y[kk][jj][ii][m]
                             : I_slope_z[kk][jj][ii][m]) * dep
            - ((axis == 0) ? I_slope_y[kk][jj][ii][m] * directions[m].y
                                + I_slope_z[kk][jj][ii][m] * directions[m].z
               : (axis == 1) ? I_slope_x[kk][jj][ii][m] * directions[m].x
                                + I_slope_z[kk][jj][ii][m] * directions[m].z
                             : I_slope_x[kk][jj][ii][m] * directions[m].x
                                + I_slope_y[kk][jj][ii][m] * directions[m].y) * h;
    }

    bool bnd = (face == 0 || face == cells);
    if (bnd) {
        double wT = wall_T(axis, face);
        double wIb = pow4(wT);
        double weps = ((axis == 1 && face == 0) ? EW : 1.0);

        double Jin = 0.0;
        if (weps < 1.0) {
            for (int m = 0; m < M; ++m) {
                double vel = (axis == 0) ? directions[m].x : ((axis == 1) ? directions[m].y : directions[m].z);
                bool tw = (face == 0 && vel < 0.0) || (face == cells && vel > 0.0);
                if (tw) Jin += bar[m] * fabs(vel) * directions[m].weight * (4.0 * PI);
            }
        }

        int adj = (face == 0) ? 0 : (cells - 1);
        int ii = a, jj = a, kk = b;
        if (axis == 0) { ii = adj; jj = a; kk = b; }
        else if (axis == 1) { ii = a; jj = adj; kk = b; }
        else { ii = a; jj = b; kk = adj; }
        double src =
            (1.0 - omega) * pow4(T[kk][jj][ii]) + omega * G[kk][jj][ii];

        for (int m = 0; m < M; ++m) {
            double vel = (axis == 0) ? directions[m].x : ((axis == 1) ? directions[m].y : directions[m].z);
            bool out = (face == 0 && vel > 0.0) || (face == cells && vel < 0.0);
            if (out)
                phys[m] = weps * wIb + (1.0 - weps) * Jin / PI;
            else
                phys[m] = fo * bar[m] + fe * src;
        }
    } else {
        double bm = 0.0;
        for (int m = 0; m < M; ++m) bm += directions[m].weight * bar[m];
        int il, jl, kl, ir, jr, kr;
        if (axis == 0) { il=face-1; jl=a; kl=b; ir=face; jr=a; kr=b; }
        else if (axis == 1) { il=a; jl=face-1; kl=b; ir=a; jr=face; kr=b; }
        else { il=a; jl=b; kl=face-1; ir=a; jr=b; kr=face; }
        double fT = 0.5 * (T[kl][jl][il] + T[kr][jr][ir]);
        double fT4 = pow4(fT);
        double fG = (bm + fmc * fT4) / (1.0 + fmc);
        double src = (1.0 - omega) * fT4 + omega * fG;
        for (int m = 0; m < M; ++m)
            phys[m] = fo * bar[m] + fe * src;
    }

    for (int m = 0; m < M; ++m) {
        double vel = (axis == 0) ? directions[m].x : ((axis == 1) ? directions[m].y : directions[m].z);
        double flux = (1.0 / dx) * vel * phys[m];
        if (face > 0) {
            int ii, jj, kk;
            if (axis == 0) { ii = face-1; jj = a; kk = b; }
            else if (axis == 1) { ii = a; jj = face-1; kk = b; }
            else { ii = a; jj = b; kk = face-1; }
            rate[kk][jj][ii][m] -= flux;
        }
        if (face < cells) {
            int ii, jj, kk;
            if (axis == 0) { ii = face; jj = a; kk = b; }
            else if (axis == 1) { ii = a; jj = face; kk = b; }
            else { ii = a; jj = b; kk = face; }
            rate[kk][jj][ii][m] += flux;
        }
    }
}

int radiation_dugks(double omega, bool force_conv = false) {
    double tlim = 0.0;
    for (int m = 0; m < M; ++m)
        tlim = max(tlim, (fabs(directions[m].x)+fabs(directions[m].y)+fabs(directions[m].z))/dx);
    double pdt = CFL_R / tlim;
    double h = 0.5 * pdt;
    double tr = 1.0 / TAU_L;
    double co = (2.0*tr - h) / (2.0*tr + pdt);
    double ce = 3.0 * h / (2.0*tr + pdt);
    double fo = 2.0 / (2.0 + h);
    double fe = h / (2.0 + h);
    double fmc = 0.5 * h * (1.0 - omega) * TAU_L;
    bool cnv = force_conv || radiation_needs_initial_convergence;
    int mx = cnv ? RAD_MAX_STEPS : RAD_TRACKING_STEPS;

    for (int step = 1; step <= mx; ++step) {
        rec_radiation(omega, pdt);
        for (int k = 0; k < NZ; ++k)
            for (int j = 0; j < NY; ++j)
                for (int i = 0; i < NX; ++i)
                    G_pre[k][j][i] = G[k][j][i];

        #pragma omp parallel for collapse(3) schedule(static)
        for (int k = 0; k < NZ; ++k)
            for (int j = 0; j < NY; ++j)
                for (int i = 0; i < NX; ++i) {
                    double src =
                        (1.0 - omega) * pow4(T[k][j][i])
                        + omega * G[k][j][i];
                    for (int m = 0; m < M; ++m) {
                        I_bar_plus[k][j][i][m] =
                            co * I_tilde[k][j][i][m] + ce * src;
                        I_tilde_plus[k][j][i][m] =
                            (4.0 / 3.0) * I_bar_plus[k][j][i][m]
                            - (1.0 / 3.0) * I_tilde[k][j][i][m];
                    }
                }

        calc_R_slopes();

        for (int k = 0; k < NZ; ++k)
            for (int j = 0; j < NY; ++j)
                for (int i = 0; i < NX; ++i)
                    for (int m = 0; m < M; ++m)
                        rate_r[k][j][i][m] = 0.0;
        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < NZ; ++k)
            for (int j = 0; j < NY; ++j)
                for (int face = 0; face <= NX; ++face)
                    radiation_face(0, j, k, face, omega, h, fo, fe, fmc, rate_r);
        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < NZ; ++k)
            for (int i = 0; i < NX; ++i)
                for (int face = 0; face <= NY; ++face)
                    radiation_face(1, i, k, face, omega, h, fo, fe, fmc, rate_r);
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i)
                for (int face = 0; face <= NZ; ++face)
                    radiation_face(2, i, j, face, omega, h, fo, fe, fmc, rate_r);

        #pragma omp parallel for collapse(3) schedule(static)
        for (int k = 0; k < NZ; ++k)
            for (int j = 0; j < NY; ++j)
                for (int i = 0; i < NX; ++i)
                    for (int m = 0; m < M; ++m)
                        I_tilde[k][j][i][m] =
                            I_tilde_plus[k][j][i][m]
                            + pdt * rate_r[k][j][i][m];

        rec_radiation(omega, pdt);
        double err = 0.0;
        for (int k = 0; k < NZ; ++k)
            for (int j = 0; j < NY; ++j)
                for (int i = 0; i < NX; ++i)
                    err = max(err, fabs(G[k][j][i] - G_pre[k][j][i]));
        if (cnv && step > 30 && err < RAD_error) {
            radiation_needs_initial_convergence = false;
            return step;
        }
    }
    radiation_needs_initial_convergence = false;
    return mx;
}

void Energy_equation(double dt, double Np, double omega) {
    double coeff = (1.0 - omega) / max(Np, 1.0e-14);
    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i) {
                double delta =
                    dt * coeff * (G[k][j][i] - pow4(T[k][j][i]));
                T[k][j][i] += delta;
                for (int q = 0; q < Q; ++q)
                    g_tilde[k][j][i][q] += wt[q] * delta;
            }
}

void coupled_step(double dt, double Np, double omega) {
    bool upd = radiation_needs_initial_convergence
        || (thermal_step_counter % RAD_COUPLING_INTERVAL == 0);
    if (upd) radiation_dugks(omega, false);
    Energy_equation(0.5 * dt, Np, omega);
    temperture(dt);
    if (upd) radiation_dugks(omega, false);
    Energy_equation(0.5 * dt, Np, omega);
    rec_temperature();
    ++thermal_step_counter;
}

// ---------- output ----------
double sample_T(double y) {
    if (y <= 0.0) return TW;
    if (y >= 1.0) return TE;
    int ic = NX / 2, kc = NZ / 2;
    double ci = y / dy + 0.5;
    int j0 = (int)floor(ci);
    if (j0 < 1) {
        double y0 = 0.5 * dy;
        return TW + (T[kc][0][ic] - TW) * y / y0;
    }
    if (j0 >= NY) {
        double y0 = (NY - 0.5) * dy;
        double t0 = T[kc][NY-1][ic];
        return t0 + (TE - t0) * (y - y0) / (1.0 - y0);
    }
    double y0 = (j0 - 0.5) * dy;
    double fr = (y - y0) / dy;
    return T[kc][j0][ic]
        + fr * (T[kc][j0+1][ic] - T[kc][j0][ic]);
}

void save_curve(int cid) {
    for (int n = 0; n < NOUT; ++n)
        saved_temp[cid][n] = sample_T((double)n / (NOUT - 1));
}

void save_fig3_plane(int cid) {
    const int ic = NX / 2;
    for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
            saved_fig3[cid][k][j] = T[k][j][ic];
}

void write_curve_csv(
    const string& path, const double* transient_times, int transient_count
) {
    ofstream o(path);
    o << "y";
    for (int k = 0; k < transient_count; ++k) {
        ostringstream lb;
        lb << defaultfloat << setprecision(12) << transient_times[k];
        o << ",t=" << lb.str();
    }
    o << ",Steady\n";
    for (int n = 0; n < NOUT; ++n) {
        o << fixed << setprecision(8) << (double)n / (NOUT - 1);
        for (int k = 0; k <= transient_count; ++k)
            o << "," << scientific << setprecision(12) << saved_temp[k][n];
        o << "\n";
    }
}

double fig3_output_temperature(int cid, int kz, int jy) {
    if (jy == 0) return TW;
    if (jy == NY + 1) return TE;
    if (kz == 0 || kz == NZ + 1) return TE;
    return saved_fig3[cid][kz - 1][jy - 1];
}

void write_fig3_dat(const string& path, int cid, const string& label) {
    ofstream o(path);
    o << "TITLE = \"Sun2012 Fig.3 isotherms at x*=0.5\"\n";
    o << "VARIABLES = \"Y\", \"Z\", \"T\"\n";
    o << "ZONE T=\"" << label << "\", I=" << NY + 2
      << ", J=" << NZ + 2 << ", F=POINT\n";
    for (int kz = 0; kz < NZ + 2; ++kz) {
        const double z = kz == 0 ? 0.0
            : (kz == NZ + 1 ? 1.0 : (kz - 0.5) * dz);
        for (int jy = 0; jy < NY + 2; ++jy) {
            const double y = jy == 0 ? 0.0
                : (jy == NY + 1 ? 1.0 : (jy - 0.5) * dy);
            o << scientific << setprecision(8)
              << y << " " << z << " "
              << fig3_output_temperature(cid, kz, jy) << "\n";
        }
    }
}

void write_fig3_files(
    const double* transient_times, int transient_count
) {
    const char labels[MAX_CURVES] = {'a', 'b', 'c', 'd'};
    for (int cid = 0; cid <= transient_count; ++cid) {
        ostringstream path, zone_label;
        path << "output/3" << labels[cid] << ".dat";
        if (cid < transient_count)
            zone_label << "t=" << defaultfloat << setprecision(12)
                       << transient_times[cid];
        else
            zone_label << "Steady";
        write_fig3_dat(path.str(), cid, zone_label.str());
        cout << "  wrote " << path.str() << "\n";
    }
}

// ---------- case runners ----------
void run_one_case(
    const string& name, double Np, double omega,
    const double* transient_times, int transient_count,
    bool do_fig3 = false
) {
    double base_dt = CFL_T * dx / speed_t;
    thermal_step_counter = 0;
    radiation_needs_initial_convergence = true;
    ini_case();

    double time = 0.0;
    for (int k = 0; k < transient_count; ++k) {
        while (time + 1.0e-14 < transient_times[k]) {
            double dt = min(base_dt, transient_times[k] - time);
            coupled_step(dt, Np, omega);
            time += dt;
        }
        radiation_dugks(omega, true);
        save_curve(k);
        if (do_fig3) save_fig3_plane(k);
        cout << "  saved t=" << transient_times[k] << "\n";
    }

#ifndef SMOKE_TEST
    int ss = 0;
    for (; ss < STEADY_MAX_STEPS; ++ss) {
        for (int k = 0; k < NZ; ++k)
            for (int j = 0; j < NY; ++j)
                for (int i = 0; i < NX; ++i)
                    pre_temp[k][j][i] = T[k][j][i];

        coupled_step(base_dt, Np, omega);
        if (ss % 20 == 0) {
            double err = 0.0;
            for (int k = 0; k < NZ; ++k)
                for (int j = 0; j < NY; ++j)
                    for (int i = 0; i < NX; ++i)
                        err = max(
                            err,
                            fabs(T[k][j][i] - pre_temp[k][j][i])
                        );
            if (ss > 200 && err < STEADY_error) break;
        }
    }
    radiation_dugks(omega, true);
    save_curve(transient_count);
    if (do_fig3) save_fig3_plane(transient_count);
    cout << "  steady after " << ss << " extra steps\n";
#else
    save_curve(transient_count);
    if (do_fig3) save_fig3_plane(transient_count);
#endif

    if (do_fig3) {
        write_fig3_files(transient_times, transient_count);
    }

    if (!do_fig3) {
        write_curve_csv(
            "output/" + name + ".csv", transient_times, transient_count
        );
        cout << "wrote output/" << name << ".csv\n";
    }
}

void run_table1() {
    cout << "\n=== Computed temperatures at selected positions ===\n";
    cout << "N=0.01 omega=0 tau_L=1 eps=1 t*=0.005\n";
    EW = 1.0; TAU_L = 1.0;
    double base_dt = CFL_T * dx / speed_t;
    thermal_step_counter = 0;
    radiation_needs_initial_convergence = true;
    ini_case();

    double time = 0.0;
    while (time + 1.0e-14 < 0.005) {
        double dt = min(base_dt, 0.005 - time);
        coupled_step(dt, 0.01, 0.0);
        time += dt;
    }
    radiation_dugks(0.0, true);

    double pos[3] = {0.25, 0.50, 0.75};

    ofstream o("output/table1.csv");
    o << "y,computed_temperature\n";
    cout << "  y      Computed temperature\n";
    for (int n = 0; n < 3; ++n) {
        double val = sample_T(pos[n]);
        o << fixed << setprecision(2) << pos[n] << ","
          << scientific << setprecision(12) << val << "\n";
        cout << "  " << fixed << setprecision(2) << pos[n] << "  "
             << scientific << setprecision(6) << val << "\n";
    }
    cout << "wrote output/table1.csv\n";
}

void run_sun2012_cases() {
    const double times_005_015[2] = {0.005, 0.015};
    const double times_005_015_050[3] = {0.005, 0.015, 0.050};
    const double times_005[1] = {0.005};

#ifdef _WIN32
    mkdir("output");
#else
    mkdir("output", 0755);
#endif
    cout << "3D DUGKS conduction-radiation | D3Q27 | "
         << NX << "x" << NY << "x" << NZ << " | "
         << NTHETA << "x" << MAX_PHI << " angles\n";
    cout << "==========================================\n";

    // Table 1 output is temporarily disabled.
    // run_table1();

    // Fig.3: t*=0.005, 0.015, 0.050 and steady center-plane isotherms
    cout << "\n=== Sun2012 Fig.3: Center-Plane Isotherms ===\n";
    EW = 1.0; TAU_L = 1.0;
    run_one_case(
        "fig3", 0.01, 0.0, times_005_015_050, 3, true
    );

    // Fig.4: scattering albedo
    cout << "\n=== Sun2012 Fig.4: Scattering Albedo ===\n";
    EW = 1.0; TAU_L = 1.0;
    cout << "--- 4a: omega=0.5 ---\n";
    run_one_case("4a", 0.01, 0.5, times_005_015, 2);
    cout << "--- 4b: omega=0.9 ---\n";
    run_one_case("4b", 0.01, 0.9, times_005_015_050, 3);

    // Fig.5: conduction-radiation parameter
    cout << "\n=== Sun2012 Fig.5: N parameter ===\n";
    EW = 1.0; TAU_L = 1.0;
    cout << "--- 5a: N=0.1 ---\n";
    run_one_case("5a", 0.1, 0.0, times_005_015_050, 3);
    cout << "--- 5b: N=1.0 ---\n";
    run_one_case("5b", 1.0, 0.0, times_005_015_050, 3);

    // Fig.6: wall emissivity
    cout << "\n=== Sun2012 Fig.6: Wall Emissivity ===\n";
    TAU_L = 1.0;
    cout << "--- 6a: eps_s=0.1 ---\n";
    EW = 0.1;
    run_one_case("6a", 0.01, 0.0, times_005_015, 2);
    cout << "--- 6b: eps_s=0.5 ---\n";
    EW = 0.5;
    run_one_case("6b", 0.01, 0.0, times_005_015, 2);

    // Fig.7: optical thickness
    cout << "\n=== Sun2012 Fig.7: Optical Thickness ===\n";
    EW = 1.0;
    cout << "--- 7a: tau_L=0.1 ---\n";
    TAU_L = 0.1;
    run_one_case("7a", 0.01, 0.0, times_005_015_050, 3);
    cout << "--- 7b: tau_L=5.0 ---\n";
    TAU_L = 5.0;
    run_one_case("7b", 0.01, 0.0, times_005, 1);

    cout << "\n==========================================\n";
    cout << "All Sun2012 cases complete.\n";
    cout << "Fig.3 .dat: output/3a.dat ... output/3d.dat\n";
    cout << "Fig.4-7:    output/[4-7][ab].csv\n";
}

int main() {
    ini_velocities();
    ini_angles();
    run_sun2012_cases();
    return 0;
}
