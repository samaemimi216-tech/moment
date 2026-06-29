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
// Sun2012 figures are reported with the hot-wall-normalized temperature.
// The heated south wall is 1 and the other walls are 0.5 in output space.
const double TW = 1.0;
const double TE = 0.5;
const double OUTPUT_TEMP_SCALE = 1.0;

#ifndef GRID_N
#define GRID_N 31
#endif
const int NX = GRID_N;
const int NY = GRID_N;
const int NZ = GRID_N;
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
const int RAD_COUPLING_INTERVAL = 1;
const double RAD_SLOPE_FACTOR = 0.0;

#ifdef RAD_MAX_STEPS_OVERRIDE
const int RAD_MAX_STEPS = RAD_MAX_STEPS_OVERRIDE;
#elif defined(SMOKE_TEST)
const int RAD_MAX_STEPS = 300;
#else
const int RAD_MAX_STEPS = 12000;
#endif

#ifdef SMOKE_TEST
const int STEADY_MAX_STEPS = 100;
#else
const int STEADY_MAX_STEPS = 30000;
#endif

const double SUN2012_DT = 1.0e-4;

const double dx = 1.0 / NX;
const double dy = 1.0 / NY;
const double dz = 1.0 / NZ;

const double speed_t = 100.0;
const double RT = speed_t * speed_t / 3.0;

/*
Sun2012 plotted-temperature normalization contract:

theta = H / 2 = T / Ts is evolved and reported by the solver.
cold walls and initial field: theta = 0.5
south hot wall: theta = 1.0

Energy equation in the plotted transient time:
dtheta/dt = Laplacian(theta)
            + tauL^2 * (1-omega)/Ncr * (J - theta^4)

Radiation equation:
(1/tauL) s . grad(I_m) + I_m
    = (1-omega) theta^4 + omega J

J = sum_m w_m I_m

Therefore:
thermal_tau = 1 / RT
source_coeff = tauL^2 * (1-omega) / Ncr
radiation velocity = s / tauL
*/

struct Dir {
    double x, y, z, weight;
};

Dir directions[MAX_M];
double direction_weight_y_sum = 0.0;

double cx[Q], cy[Q], cz[Q], wt[Q];
int opposite[Q];

double T[NZ][NY][NX];
double pre_temp[NZ][NY][NX];
double G[NZ][NY][NX];
double G_pre[NZ][NY][NX];
double Qr[NZ][NY][NX];

double g_tilde[NZ][NY][NX][Q];
double g_bar_plus[NZ][NY][NX][Q];
double g_tilde_plus[NZ][NY][NX][Q];
double g_slope_x[NZ][NY][NX][Q];
double g_slope_y[NZ][NY][NX][Q];
double g_slope_z[NZ][NY][NX][Q];
double g_face_x[NZ][NY][NX + 1][Q];
double g_face_y[NZ][NY + 1][NX][Q];
double g_face_z[NZ + 1][NY][NX][Q];

double I_tilde[NZ][NY][NX][M];
double I_bar_plus[NZ][NY][NX][M];
double I_tilde_plus[NZ][NY][NX][M];
double I_slope_x[NZ][NY][NX][M];
double I_slope_y[NZ][NY][NX][M];
double I_slope_z[NZ][NY][NX][M];
double I_face_x[NZ][NY][NX + 1][M];
double I_face_y[NZ][NY + 1][NX][M];
double I_face_z[NZ + 1][NY][NX][M];

double saved_temp[MAX_CURVES][NOUT];
double saved_fig3[MAX_CURVES][NZ][NY];

int thermal_step_counter = 0;
bool radiation_needs_initial_convergence = true;
double EW = 1.0;
double TAU_L = 1.0;
bool use_uniform_wall_temperature = false;
double uniform_wall_temperature = TE;
bool use_uniform_wall_emissivity = false;
double uniform_wall_emissivity = 1.0;

inline double pow4(double v) {
    double s = v * v;
    return s * s;
}

inline double minmod(double a, double b) {
    if (a * b <= 0.0) {
        return 0.0;
    }
    if (fabs(a) < fabs(b)) {
        return a;
    } else {
        return b;
    }
}

inline double thermal_tau() {
    return 1.0 / RT;
}

inline double source_coeff(double Ncr, double omega) {
    return TAU_L * TAU_L * (1.0 - omega) / max(Ncr, 1.0e-14);
}

inline void radiation_speed(int m, double& ax, double& ay, double& az) {
    ax = directions[m].x / TAU_L;
    ay = directions[m].y / TAU_L;
    az = directions[m].z / TAU_L;
}

struct Vec3 {
    double x, y, z;
};

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 wall_normal(int axis, int face) {
    if (axis == 0) {
        return face == 0 ? Vec3{-1.0, 0.0, 0.0} : Vec3{1.0, 0.0, 0.0};
    }
    if (axis == 1) {
        return face == 0 ? Vec3{0.0, -1.0, 0.0} : Vec3{0.0, 1.0, 0.0};
    }
    return face == 0 ? Vec3{0.0, 0.0, -1.0} : Vec3{0.0, 0.0, 1.0};
}

inline double radiation_equilibrium(double theta, double J, double omega) {
    return (1.0 - omega) * pow4(theta) + omega * J;
}

double recover_J_from_Itilde(double sum_Itilde, double theta, double dchi, double omega) {
    double tau_I = 1.0;
    double a = 0.5 * dchi / tau_I;
    double b = a * (1.0 - omega);
    return (sum_Itilde + b * pow4(theta)) / (1.0 + b);
}

double recover_theta_from_gtilde(
    double gsum, double J, double dt, double Cr, int* iterations = 0
) {
    double A = 0.5 * dt * Cr;
    double theta = max(0.5, gsum);

    for (int it = 0; it < 50; ++it) {
        double t2 = theta * theta;
        double t3 = t2 * theta;
        double t4 = t2 * t2;

        double F = theta - gsum - A * J + A * t4;
        double dF = 1.0 + 4.0 * A * t3;
        double next = theta - F / dF;

        if (!isfinite(next) || next < 0.5) {
            next = max(0.5, theta - 0.05);
        } else if (next > 1.0) {
            next = min(1.0, theta + 0.05);
        }

        if (iterations) {
            *iterations = it + 1;
        }
        if (fabs(next - theta) < 1.0e-12) {
            return next;
        }
        theta = next;
    }

    return theta;
}

double compute_thermal_dt() {
    double max_xi = 0.0;

    for (int q = 0; q < Q; ++q) {
        double mag = sqrt(cx[q] * cx[q] + cy[q] * cy[q] + cz[q] * cz[q]);
        max_xi = max(max_xi, mag);
    }

    double ds = min(dx, min(dy, dz));
    return CFL_T * ds / max(max_xi, 1.0e-14);
}

double compute_radiation_dchi() {
    double max_rate = 0.0;

    for (int m = 0; m < M; ++m) {
        double ax, ay, az;
        radiation_speed(m, ax, ay, az);

        double rate = fabs(ax) / dx + fabs(ay) / dy + fabs(az) / dz;
        max_rate = max(max_rate, rate);
    }

    return CFL_R / max(max_rate, 1.0e-14);
}

void ensure_dir(const string& path) {
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    mkdir(path.c_str(), 0755);
#endif
}

string time_label(double value) {
    ostringstream label;
    label << defaultfloat << setprecision(12) << value;
    return label.str();
}

void ini_velocities() {
    const double w1[3] = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
    const int v1[3] = {-1, 0, 1};
    int q = 0;
    for (int kz = 0; kz < 3; ++kz) {
        for (int jy = 0; jy < 3; ++jy) {
            for (int ix = 0; ix < 3; ++ix) {
                cx[q] = speed_t * v1[ix];
                cy[q] = speed_t * v1[jy];
                cz[q] = speed_t * v1[kz];
                wt[q] = w1[ix] * w1[jy] * w1[kz];
                int ox = 2 - ix;
                int oy = 2 - jy;
                int oz = 2 - kz;
                opposite[q] = (oz * 3 + oy) * 3 + ox;
                ++q;
            }
        }
    }
}

void ini_angles() {
    const double dtheta = PI / NTHETA;
    const double dphi = 2.0 * PI / MAX_PHI;
    int m = 0;
    direction_weight_y_sum = 0.0;
    for (int it = 0; it < NTHETA; ++it) {
        double tl = it * dtheta;
        double tr = (it + 1) * dtheta;
        double txy = 0.5 * (tr - tl) - 0.25 * (sin(2.0 * tr) - sin(2.0 * tl));
        double tz = 0.5 * (sin(tr) * sin(tr) - sin(tl) * sin(tl));
        for (int ip = 0; ip < MAX_PHI; ++ip) {
            double pl = ip * dphi;
            double pr = (ip + 1) * dphi;
            double sa = (cos(tl) - cos(tr)) * dphi;
            directions[m].x = txy * (sin(pr) - sin(pl)) / sa;
            directions[m].y = txy * (cos(pl) - cos(pr)) / sa;
            directions[m].z = tz * dphi / sa;
            directions[m].weight = sa / (4.0 * PI);
            direction_weight_y_sum += directions[m].weight * directions[m].y;
            ++m;
        }
    }
}

void ini_case() {
    double initial_I = pow4(TE);

    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
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
                }

                for (int m = 0; m < M; ++m) {
                    I_tilde[k][j][i][m] = initial_I;
                    I_bar_plus[k][j][i][m] = initial_I;
                    I_tilde_plus[k][j][i][m] = initial_I;
                    I_slope_x[k][j][i][m] = 0.0;
                    I_slope_y[k][j][i][m] = 0.0;
                    I_slope_z[k][j][i][m] = 0.0;
                }
            }
        }
    }
}

double wall_T(int axis, int face) {
    if (use_uniform_wall_temperature) {
        return uniform_wall_temperature;
    }
    if (axis == 1 && face == 0) {
        return TW;
    } else {
        return TE;
    }
}

double wall_eps(int axis, int face) {
    if (use_uniform_wall_emissivity) {
        return uniform_wall_emissivity;
    }
    if (axis == 1 && face == 0) {
        return EW;
    }
    return 1.0;
}

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
            reflected += directions[mp].weight * I_phys_in * cos_theta;
        }
    }

    double I_emit = eps_wall * pow4(theta_wall);
    double I_refl = 4.0 * (1.0 - eps_wall) * reflected;
    return I_emit + I_refl;
}

void rec_temperature() {
    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                T[k][j][i] = 0.0;
                for (int q = 0; q < Q; ++q) {
                    T[k][j][i] += g_tilde[k][j][i][q];
                }
            }
        }
    }
}

void calc_T_slopes() {
    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                for (int q = 0; q < Q; ++q) {
                    g_slope_x[k][j][i][q] = 0.0;
                    g_slope_y[k][j][i][q] = 0.0;
                    g_slope_z[k][j][i][q] = 0.0;
                }
            }
        }
    }

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 1; k < NZ - 1; ++k) {
        for (int j = 1; j < NY - 1; ++j) {
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
    }
}

void temperture(double dt) {
    double h = 0.5 * dt;
    double tau_g = thermal_tau();
    double co = (2.0 * tau_g - h) / (2.0 * tau_g + dt);
    double ce = 3.0 * h / (2.0 * tau_g + dt);
    double fo = 2.0 * tau_g / (2.0 * tau_g + h);
    double fe = h / (2.0 * tau_g + h);

    rec_temperature();

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                for (int q = 0; q < Q; ++q) {
                    double eq = wt[q] * T[k][j][i];
                    g_bar_plus[k][j][i][q] =
                        co * g_tilde[k][j][i][q] + ce * eq;
                    g_tilde_plus[k][j][i][q] =
                        (4.0 / 3.0) * g_bar_plus[k][j][i][q]
                        - (1.0 / 3.0) * g_tilde[k][j][i][q];
                }
            }
        }
    }
    calc_T_slopes();

    #pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            for (int face = 0; face <= NX; ++face) {
                double bar[Q];
                for (int q = 0; q < Q; ++q) {
                    int cell;
                    if (cx[q] > 0.0) {
                        cell = max(0, face - 1);
                    } else if (cx[q] < 0.0) {
                        cell = min(NX - 1, face);
                    } else {
                        cell = min(NX - 1, max(0, face));
                    }

                    double cell_x = (cell + 0.5) * dx;
                    double departure_x = face * dx - cx[q] * h;
                    bar[q] = g_bar_plus[k][j][cell][q]
                        + g_slope_x[k][j][cell][q] * (departure_x - cell_x)
                        - (g_slope_y[k][j][cell][q] * cy[q]
                        + g_slope_z[k][j][cell][q] * cz[q]) * h;
                }

                if (face == 0 || face == NX) {
                    double wall_temperature = wall_T(0, face == 0 ? 0 : 1);
                    int wall_cell = (face == 0) ? 0 : NX - 1;
                    double interior_temperature = T[k][j][wall_cell];
                    for (int q = 0; q < Q; ++q) {
                        bool incoming =
                            (face == 0 && cx[q] > 0.0)
                            || (face == NX && cx[q] < 0.0);
                        if (incoming) {
                            double non_equilibrium =
                                g_bar_plus[k][j][wall_cell][q]
                                - wt[q] * interior_temperature;
                            bar[q] = wt[q] * wall_temperature
                                + non_equilibrium;
                        } else if (cx[q] == 0.0) {
                            bar[q] = wt[q] * wall_temperature;
                        }
                    }
                }

                double face_temperature = 0.0;
                for (int q = 0; q < Q; ++q) {
                    face_temperature += bar[q];
                }
                if (face == 0 || face == NX) {
                    face_temperature = wall_T(0, face == 0 ? 0 : 1);
                }

                for (int q = 0; q < Q; ++q) {
                    g_face_x[k][j][face][q] =
                        fo * bar[q] + fe * wt[q] * face_temperature;
                }
            }
        }
    }

    #pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < NZ; ++k) {
        for (int i = 0; i < NX; ++i) {
            for (int face = 0; face <= NY; ++face) {
                double bar[Q];
                for (int q = 0; q < Q; ++q) {
                    int cell;
                    if (cy[q] > 0.0) {
                        cell = max(0, face - 1);
                    } else if (cy[q] < 0.0) {
                        cell = min(NY - 1, face);
                    } else {
                        cell = min(NY - 1, max(0, face));
                    }

                    double cell_y = (cell + 0.5) * dy;
                    double departure_y = face * dy - cy[q] * h;
                    bar[q] = g_bar_plus[k][cell][i][q]
                        + g_slope_y[k][cell][i][q] * (departure_y - cell_y)
                        - (g_slope_x[k][cell][i][q] * cx[q]
                        + g_slope_z[k][cell][i][q] * cz[q]) * h;
                }

                if (face == 0 || face == NY) {
                    double wall_temperature = wall_T(1, face == 0 ? 0 : 1);
                    int wall_cell = (face == 0) ? 0 : NY - 1;
                    double interior_temperature = T[k][wall_cell][i];
                    for (int q = 0; q < Q; ++q) {
                        bool incoming =
                            (face == 0 && cy[q] > 0.0)
                            || (face == NY && cy[q] < 0.0);
                        if (incoming) {
                            double non_equilibrium =
                                g_bar_plus[k][wall_cell][i][q]
                                - wt[q] * interior_temperature;
                            bar[q] = wt[q] * wall_temperature
                                + non_equilibrium;
                        } else if (cy[q] == 0.0) {
                            bar[q] = wt[q] * wall_temperature;
                        }
                    }
                }

                double face_temperature = 0.0;
                for (int q = 0; q < Q; ++q) {
                    face_temperature += bar[q];
                }
                if (face == 0) {
                    face_temperature = wall_T(1, 0);
                }
                if (face == NY) {
                    face_temperature = wall_T(1, 1);
                }

                for (int q = 0; q < Q; ++q) {
                    g_face_y[k][face][i][q] =
                        fo * bar[q] + fe * wt[q] * face_temperature;
                }
            }
        }
    }

    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            for (int face = 0; face <= NZ; ++face) {
                double bar[Q];
                for (int q = 0; q < Q; ++q) {
                    int cell;
                    if (cz[q] > 0.0) {
                        cell = max(0, face - 1);
                    } else if (cz[q] < 0.0) {
                        cell = min(NZ - 1, face);
                    } else {
                        cell = min(NZ - 1, max(0, face));
                    }

                    double cell_z = (cell + 0.5) * dz;
                    double departure_z = face * dz - cz[q] * h;
                    bar[q] = g_bar_plus[cell][j][i][q]
                        + g_slope_z[cell][j][i][q] * (departure_z - cell_z)
                        - (g_slope_x[cell][j][i][q] * cx[q]
                        + g_slope_y[cell][j][i][q] * cy[q]) * h;
                }

                if (face == 0 || face == NZ) {
                    double wall_temperature = wall_T(2, face == 0 ? 0 : 1);
                    int wall_cell = (face == 0) ? 0 : NZ - 1;
                    double interior_temperature = T[wall_cell][j][i];
                    for (int q = 0; q < Q; ++q) {
                        bool incoming =
                            (face == 0 && cz[q] > 0.0)
                            || (face == NZ && cz[q] < 0.0);
                        if (incoming) {
                            double non_equilibrium =
                                g_bar_plus[wall_cell][j][i][q]
                                - wt[q] * interior_temperature;
                            bar[q] = wt[q] * wall_temperature
                                + non_equilibrium;
                        } else if (cz[q] == 0.0) {
                            bar[q] = wt[q] * wall_temperature;
                        }
                    }
                }

                double face_temperature = 0.0;
                for (int q = 0; q < Q; ++q) {
                    face_temperature += bar[q];
                }
                if (face == 0 || face == NZ) {
                    face_temperature = wall_T(2, face == 0 ? 0 : 1);
                }

                for (int q = 0; q < Q; ++q) {
                    g_face_z[face][j][i][q] =
                        fo * bar[q] + fe * wt[q] * face_temperature;
                }
            }
        }
    }

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                for (int q = 0; q < Q; ++q) {
                    g_tilde[k][j][i][q] = g_tilde_plus[k][j][i][q]
                        - (dt / dx) * cx[q]
                        * (g_face_x[k][j][i + 1][q] - g_face_x[k][j][i][q])
                        - (dt / dy) * cy[q]
                        * (g_face_y[k][j + 1][i][q] - g_face_y[k][j][i][q])
                        - (dt / dz) * cz[q]
                        * (g_face_z[k + 1][j][i][q] - g_face_z[k][j][i][q]);
                }
            }
        }
    }

    rec_temperature();
}

void calc_R_slopes() {
    if (RAD_SLOPE_FACTOR == 0.0) {
        return;
    }

    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                for (int m = 0; m < M; ++m) {
                    I_slope_x[k][j][i][m] = 0.0;
                    I_slope_y[k][j][i][m] = 0.0;
                    I_slope_z[k][j][i][m] = 0.0;
                }
            }
        }
    }

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 1; k < NZ - 1; ++k) {
        for (int j = 1; j < NY - 1; ++j) {
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
    }
}

void rec_radiation(double omega, double pdt) {
    double po = 2.0 / (2.0 + pdt);
    double pe = pdt / (2.0 + pdt);

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                double tm = 0.0;
                double qsum_y = 0.0;
                for (int m = 0; m < M; ++m) {
                    double weighted_I = directions[m].weight * I_tilde[k][j][i][m];
                    tm += weighted_I;
                    qsum_y += directions[m].y * weighted_I;
                }

                G[k][j][i] = recover_J_from_Itilde(tm, T[k][j][i], pdt, omega);

                double src = radiation_equilibrium(T[k][j][i], G[k][j][i], omega);
                Qr[k][j][i] = po * qsum_y + pe * src * direction_weight_y_sum;
            }
        }
    }
}

int radiation_dugks(double omega, bool force_conv = false) {
    double pdt = compute_radiation_dchi();
    double h = 0.5 * pdt;
    double co = (2.0 - h) / (2.0 + pdt);
    double ce = 3.0 * h / (2.0 + pdt);
    double fo = 2.0 / (2.0 + h);
    double fe = h / (2.0 + h);

    bool cnv = false;
    if (force_conv || radiation_needs_initial_convergence) {
        cnv = true;
    }

    int mx = RAD_TRACKING_STEPS;
    if (cnv) {
        mx = RAD_MAX_STEPS;
    }

    for (int step = 1; step <= mx; ++step) {
        rec_radiation(omega, pdt);

        for (int k = 0; k < NZ; ++k) {
            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    G_pre[k][j][i] = G[k][j][i];
                }
            }
        }

        #pragma omp parallel for collapse(3) schedule(static)
        for (int k = 0; k < NZ; ++k) {
            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    double src = radiation_equilibrium(T[k][j][i], G[k][j][i], omega);
                    for (int m = 0; m < M; ++m) {
                        I_bar_plus[k][j][i][m] = co * I_tilde[k][j][i][m] + ce * src;
                        I_tilde_plus[k][j][i][m] = (4.0 / 3.0) * I_bar_plus[k][j][i][m] - (1.0 / 3.0) * I_tilde[k][j][i][m];
                    }
                }
            }
        }

        calc_R_slopes();

        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < NZ; ++k) {
            for (int j = 0; j < NY; ++j) {
                for (int face = 0; face <= NX; ++face) {
                    double bar[MAX_M];
                    for (int m = 0; m < M; ++m) {
                        double ax, ay, az;
                        radiation_speed(m, ax, ay, az);
                        int cell;
                        if (ax > 0.0) {
                            cell = max(0, face - 1);
                        } else if (ax < 0.0) {
                            cell = min(NX - 1, face);
                        } else {
                            cell = min(NX - 1, max(0, face));
                        }

                        double cell_x = (cell + 0.5) * dx;
                        double departure_x = face * dx - ax * h;
                        bar[m] = I_bar_plus[k][j][cell][m]
                            + I_slope_x[k][j][cell][m] * (departure_x - cell_x)
                            - (I_slope_y[k][j][cell][m] * ay
                            + I_slope_z[k][j][cell][m] * az) * h;
                    }

                    if (face == 0 || face == NX) {
                        int wall_face = (face == 0) ? 0 : 1;
                        Vec3 n_wall = wall_normal(0, wall_face);
                        double wall_temperature = wall_T(0, wall_face);
                        double eps = wall_eps(0, wall_face);
                        int wall_cell = (face == 0) ? 0 : NX - 1;
                        double I_eq_face[MAX_M];
                        double interior_eq = radiation_equilibrium(
                            T[k][j][wall_cell], G[k][j][wall_cell], omega);
                        for (int m = 0; m < M; ++m) {
                            I_eq_face[m] = interior_eq;
                        }
                        double wall_out = gray_wall_outgoing_intensity(
                            n_wall, wall_temperature, eps, bar, I_eq_face, fo, fe);
                        for (int m = 0; m < M; ++m) {
                            Vec3 s{directions[m].x, directions[m].y, directions[m].z};
                            if (dot(n_wall, s) < 0.0) {
                                I_face_x[k][j][face][m] = wall_out;
                            } else {
                                I_face_x[k][j][face][m] =
                                    fo * bar[m] + fe * I_eq_face[m];
                            }
                        }
                    } else {
                        double bm = 0.0;
                        for (int m = 0; m < M; ++m) {
                            bm += directions[m].weight * bar[m];
                        }
                        double face_T = 0.5 * (T[k][j][face - 1]
                            + T[k][j][face]);
                        double face_G = recover_J_from_Itilde(bm, face_T, h, omega);
                        double src = radiation_equilibrium(face_T, face_G, omega);
                        for (int m = 0; m < M; ++m) {
                            I_face_x[k][j][face][m] = fo * bar[m] + fe * src;
                        }
                    }
                }
            }
        }

        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < NZ; ++k) {
            for (int i = 0; i < NX; ++i) {
                for (int face = 0; face <= NY; ++face) {
                    double bar[MAX_M];
                    for (int m = 0; m < M; ++m) {
                        double ax, ay, az;
                        radiation_speed(m, ax, ay, az);
                        int cell;
                        if (ay > 0.0) {
                            cell = max(0, face - 1);
                        } else if (ay < 0.0) {
                            cell = min(NY - 1, face);
                        } else {
                            cell = min(NY - 1, max(0, face));
                        }

                        double cell_y = (cell + 0.5) * dy;
                        double departure_y = face * dy - ay * h;
                        bar[m] = I_bar_plus[k][cell][i][m]
                            + I_slope_y[k][cell][i][m] * (departure_y - cell_y)
                            - (I_slope_x[k][cell][i][m] * ax
                            + I_slope_z[k][cell][i][m] * az) * h;
                    }

                    if (face == 0 || face == NY) {
                        int wall_face = (face == 0) ? 0 : 1;
                        Vec3 n_wall = wall_normal(1, wall_face);
                        double wall_temperature = wall_T(1, wall_face);
                        double eps = wall_eps(1, wall_face);
                        int wall_cell = (face == 0) ? 0 : NY - 1;
                        double I_eq_face[MAX_M];
                        double interior_eq = radiation_equilibrium(
                            T[k][wall_cell][i], G[k][wall_cell][i], omega);
                        for (int m = 0; m < M; ++m) {
                            I_eq_face[m] = interior_eq;
                        }
                        double wall_out = gray_wall_outgoing_intensity(
                            n_wall, wall_temperature, eps, bar, I_eq_face, fo, fe);

                        for (int m = 0; m < M; ++m) {
                            Vec3 s{directions[m].x, directions[m].y, directions[m].z};
                            if (dot(n_wall, s) < 0.0) {
                                I_face_y[k][face][i][m] = wall_out;
                            } else {
                                I_face_y[k][face][i][m] =
                                    fo * bar[m] + fe * I_eq_face[m];
                            }
                        }
                    } else {
                        double bm = 0.0;
                        for (int m = 0; m < M; ++m) {
                            bm += directions[m].weight * bar[m];
                        }
                        double face_T = 0.5 * (T[k][face - 1][i]
                            + T[k][face][i]);
                        double face_G = recover_J_from_Itilde(bm, face_T, h, omega);
                        double src = radiation_equilibrium(face_T, face_G, omega);
                        for (int m = 0; m < M; ++m) {
                            I_face_y[k][face][i][m] = fo * bar[m] + fe * src;
                        }
                    }
                }
            }
        }

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                for (int face = 0; face <= NZ; ++face) {
                    double bar[MAX_M];
                    for (int m = 0; m < M; ++m) {
                        double ax, ay, az;
                        radiation_speed(m, ax, ay, az);
                        int cell;
                        if (az > 0.0) {
                            cell = max(0, face - 1);
                        } else if (az < 0.0) {
                            cell = min(NZ - 1, face);
                        } else {
                            cell = min(NZ - 1, max(0, face));
                        }

                        double cell_z = (cell + 0.5) * dz;
                        double departure_z = face * dz - az * h;
                        bar[m] = I_bar_plus[cell][j][i][m]
                            + I_slope_z[cell][j][i][m] * (departure_z - cell_z)
                            - (I_slope_x[cell][j][i][m] * ax
                            + I_slope_y[cell][j][i][m] * ay) * h;
                    }

                    if (face == 0 || face == NZ) {
                        int wall_face = (face == 0) ? 0 : 1;
                        Vec3 n_wall = wall_normal(2, wall_face);
                        double wall_temperature = wall_T(2, wall_face);
                        double eps = wall_eps(2, wall_face);
                        int wall_cell = (face == 0) ? 0 : NZ - 1;
                        double I_eq_face[MAX_M];
                        double interior_eq = radiation_equilibrium(
                            T[wall_cell][j][i], G[wall_cell][j][i], omega);
                        for (int m = 0; m < M; ++m) {
                            I_eq_face[m] = interior_eq;
                        }
                        double wall_out = gray_wall_outgoing_intensity(
                            n_wall, wall_temperature, eps, bar, I_eq_face, fo, fe);
                        for (int m = 0; m < M; ++m) {
                            Vec3 s{directions[m].x, directions[m].y, directions[m].z};
                            if (dot(n_wall, s) < 0.0) {
                                I_face_z[face][j][i][m] = wall_out;
                            } else {
                                I_face_z[face][j][i][m] =
                                    fo * bar[m] + fe * I_eq_face[m];
                            }
                        }
                    } else {
                        double bm = 0.0;
                        for (int m = 0; m < M; ++m) {
                            bm += directions[m].weight * bar[m];
                        }
                        double face_T = 0.5 * (T[face - 1][j][i]
                            + T[face][j][i]);
                        double face_G = recover_J_from_Itilde(bm, face_T, h, omega);
                        double src = radiation_equilibrium(face_T, face_G, omega);
                        for (int m = 0; m < M; ++m) {
                            I_face_z[face][j][i][m] = fo * bar[m] + fe * src;
                        }
                    }
                }
            }
        }

        #pragma omp parallel for collapse(3) schedule(static)
        for (int k = 0; k < NZ; ++k) {
            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    for (int m = 0; m < M; ++m) {
                        double ax, ay, az;
                        radiation_speed(m, ax, ay, az);
                        I_tilde[k][j][i][m] = I_tilde_plus[k][j][i][m]
                            - (pdt / dx) * ax
                            * (I_face_x[k][j][i + 1][m] - I_face_x[k][j][i][m])
                            - (pdt / dy) * ay
                            * (I_face_y[k][j + 1][i][m] - I_face_y[k][j][i][m])
                            - (pdt / dz) * az
                            * (I_face_z[k + 1][j][i][m] - I_face_z[k][j][i][m]);
                    }
                }
            }
        }

        rec_radiation(omega, pdt);

        double err = 0.0;
        #pragma omp parallel for collapse(3) reduction(max:err) schedule(static)
        for (int k = 0; k < NZ; ++k) {
            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    err = max(err, fabs(G[k][j][i] - G_pre[k][j][i]));
                }
            }
        }

        if (cnv && step > 30 && err < RAD_error) {
            radiation_needs_initial_convergence = false;
            return step;
        }
    }

    radiation_needs_initial_convergence = false;
    return mx;
}

// Transitional split source update. The coefficient follows the Sun2012
// hot-wall contract; the fully source-coupled thermal DUGKS is a later step.
void Energy_equation(double dt, double Np, double omega) {
    double coeff = source_coeff(Np, omega);

    double A = dt * coeff;

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                double T_star = T[k][j][i];
                double G_val = G[k][j][i];

                // 构造隐式方程: T + A * T^4 = T_star + A * G
                double RHS = T_star + A * G_val;

                double T_guess = T_star;
                double err = 1.0;
                int iter = 0;

                // Newton-Raphson 局部迭代
                while (err > 1.0e-12) {
                    if (iter >= 50) {
                        break;
                    }

                    double T_guess_2 = T_guess * T_guess;
                    double T_guess_3 = T_guess_2 * T_guess;
                    double T_guess_4 = T_guess_3 * T_guess;

                    double f_val = T_guess + A * T_guess_4 - RHS;
                    double df_val = 1.0 + 4.0 * A * T_guess_3;

                    double T_new = T_guess - f_val / df_val;
                    if (!isfinite(T_new) || T_new <= 0.0) {
                        T_new = 0.5 * T_guess;
                    }

                    if (T_new > T_guess) {
                        err = T_new - T_guess;
                    } else {
                        err = T_guess - T_new;
                    }

                    T_guess = T_new;
                    iter++;
                }

                double delta = T_guess - T_star;
                T[k][j][i] = T_guess;
                for (int q = 0; q < Q; ++q) {
                    g_tilde[k][j][i][q] += wt[q] * delta;
                }
            }
        }
    }
}

void coupled_step(double dt, double Np, double omega) {
    bool upd = false;
    if (radiation_needs_initial_convergence || (thermal_step_counter % RAD_COUPLING_INTERVAL == 0)) {
        upd = true;
    }

    if (upd) {
        radiation_dugks(omega, true);
    }

    Energy_equation(0.5 * dt, Np, omega);
    temperture(dt);

    if (upd) {
        radiation_dugks(omega, true);
    }

    Energy_equation(0.5 * dt, Np, omega);
    rec_temperature();

    ++thermal_step_counter;
}

double sample_T(double y) {
    double h_value = TE;
    if (y <= 0.0) {
        h_value = TW;
        return OUTPUT_TEMP_SCALE * h_value;
    }
    if (y >= 1.0) {
        h_value = TE;
        return OUTPUT_TEMP_SCALE * h_value;
    }

    int ic = NX / 2;
    int kc = NZ / 2;
    double ci = y / dy - 0.5;
    int j0 = (int)floor(ci);

    if (j0 < 0) {
        double y0 = 0.5 * dy;
        h_value = TW + (T[kc][0][ic] - TW) * y / y0;
        return OUTPUT_TEMP_SCALE * h_value;
    }
    if (j0 >= NY - 1) {
        double y0 = (NY - 0.5) * dy;
        double t0 = T[kc][NY-1][ic];
        h_value = t0 + (TE - t0) * (y - y0) / (1.0 - y0);
        return OUTPUT_TEMP_SCALE * h_value;
    }

    double y0 = (j0 + 0.5) * dy;
    double fr = (y - y0) / dy;

    h_value = T[kc][j0][ic] + fr * (T[kc][j0+1][ic] - T[kc][j0][ic]);
    return OUTPUT_TEMP_SCALE * h_value;
}

void save_curve(int cid) {
    for (int n = 0; n < NOUT; ++n) {
        saved_temp[cid][n] = sample_T((double)n / (NOUT - 1));
    }
}

void save_fig3_plane(int cid) {
    const int ic = NX / 2;
    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            saved_fig3[cid][k][j] = OUTPUT_TEMP_SCALE * T[k][j][ic];
        }
    }
}

void write_curve_csv(const string& path, const double* transient_times, int transient_count) {
    ofstream o(path);
    o << "y";
    for (int k = 0; k < transient_count; ++k) {
        o << "," << time_label(transient_times[k]);
    }
    o << ",ss\n";
    for (int n = 0; n < NOUT; ++n) {
        o << fixed << setprecision(8) << (double)n / (NOUT - 1);
        for (int k = 0; k <= transient_count; ++k) {
            o << "," << scientific << setprecision(12) << saved_temp[k][n];
        }
        o << "\n";
    }
}

void write_table1_csv(const string& path) {
    ofstream o(path);
    const double yvals[3] = {0.25, 0.50, 0.75};
    o << "y,theta\n";
    for (int n = 0; n < 3; ++n) {
        o << fixed << setprecision(8) << yvals[n]
          << "," << scientific << setprecision(12) << sample_T(yvals[n]) << "\n";
    }
}

double fig3_output_temperature(int cid, int kz, int jy) {
    double boundary_sum = 0.0;
    int boundary_count = 0;

    if (jy == 0) {
        boundary_sum += OUTPUT_TEMP_SCALE * TW;
        ++boundary_count;
    }
    if (jy == NY + 1) {
        boundary_sum += OUTPUT_TEMP_SCALE * TE;
        ++boundary_count;
    }
    if (kz == 0) {
        boundary_sum += OUTPUT_TEMP_SCALE * TE;
        ++boundary_count;
    }
    if (kz == NZ + 1) {
        boundary_sum += OUTPUT_TEMP_SCALE * TE;
        ++boundary_count;
    }
    if (boundary_count > 0) {
        return boundary_sum / boundary_count;
    }
    return saved_fig3[cid][kz - 1][jy - 1];
}

string fig3_csv_label(int cid, const double* transient_times, int transient_count) {
    if (cid >= transient_count) {
        return "steady";
    }
    double t = transient_times[cid];
    if (fabs(t - 0.005) < 1.0e-12) {
        return "t0005";
    }
    if (fabs(t - 0.015) < 1.0e-12) {
        return "t0015";
    }
    if (fabs(t - 0.050) < 1.0e-12) {
        return "t0050";
    }
    ostringstream label;
    label << "t" << fixed << setprecision(6) << t;
    string text = label.str();
    for (size_t p = 0; p < text.size(); ++p) {
        if (text[p] == '.') {
            text[p] = 'p';
        }
    }
    return text;
}

void write_fig3_csv(const string& path, int cid) {
    ofstream o(path);
    o << "y,z,theta\n";
    for (int kz = 0; kz < NZ + 2; ++kz) {
        double z = 0.0;
        if (kz == 0) {
            z = 0.0;
        } else if (kz == NZ + 1) {
            z = 1.0;
        } else {
            z = (kz - 0.5) * dz;
        }

        for (int jy = 0; jy < NY + 2; ++jy) {
            double y = 0.0;
            if (jy == 0) {
                y = 0.0;
            } else if (jy == NY + 1) {
                y = 1.0;
            } else {
                y = (jy - 0.5) * dy;
            }

            o << fixed << setprecision(8) << y << "," << z << ","
              << scientific << setprecision(12)
              << fig3_output_temperature(cid, kz, jy) << "\n";
        }
    }
}

void write_fig3_dat(const string& path, int cid, const string& label) {
    ofstream o(path);
    o << "TITLE = \"Sun2012 Fig.3 isotherms at x*=0.5\"\n";
    o << "VARIABLES = \"Y\", \"Z\", \"T\"\n";
    o << "ZONE T=\"" << label << "\", I=" << NY + 2 << ", J=" << NZ + 2 << ", F=POINT\n";
    for (int kz = 0; kz < NZ + 2; ++kz) {
        double z = 0.0;
        if (kz == 0) {
            z = 0.0;
        } else if (kz == NZ + 1) {
            z = 1.0;
        } else {
            z = (kz - 0.5) * dz;
        }

        for (int jy = 0; jy < NY + 2; ++jy) {
            double y = 0.0;
            if (jy == 0) {
                y = 0.0;
            } else if (jy == NY + 1) {
                y = 1.0;
            } else {
                y = (jy - 0.5) * dy;
            }

            o << scientific << setprecision(8) << y << " " << z << " " << fig3_output_temperature(cid, kz, jy) << "\n";
        }
    }
}

void write_fig3_files(const double* transient_times, int transient_count) {
    const char labels[MAX_CURVES] = {'a', 'b', 'c', 'd'};
    for (int cid = 0; cid <= transient_count; ++cid) {
        ostringstream output_path, data_path, csv_path, zone_label;
        output_path << "output/3" << labels[cid] << ".dat";
        data_path << "data/3" << labels[cid] << ".dat";
        csv_path << "output/fig3_"
                 << fig3_csv_label(cid, transient_times, transient_count)
                 << ".csv";
        if (cid < transient_count) {
            zone_label << time_label(transient_times[cid]);
        } else {
            zone_label << "ss";
        }
        write_fig3_dat(output_path.str(), cid, zone_label.str());
        write_fig3_dat(data_path.str(), cid, zone_label.str());
        write_fig3_csv(csv_path.str(), cid);
        cout << "  wrote " << output_path.str() << ", "
             << data_path.str() << ", and " << csv_path.str() << "\n";
    }
}

void run_one_case(
    const string& name, double Np, double omega,
    const double* transient_times, int transient_count,
    bool do_fig3 = false
) {
    ensure_dir("output");
    ensure_dir("data");
    double base_dt = min(compute_thermal_dt(), SUN2012_DT);
    thermal_step_counter = 0;
    radiation_needs_initial_convergence = true;
    ini_case();

    double time = 0.0;
    for (int k = 0; k < transient_count; ++k) {
        while (time + 1.0e-14 < transient_times[k]) {
            double dt = base_dt;
            if (transient_times[k] - time < base_dt) {
                dt = transient_times[k] - time;
            }
            coupled_step(dt, Np, omega);
            time += dt;
        }
        radiation_dugks(omega, true);
        save_curve(k);
        if (do_fig3) {
            save_fig3_plane(k);
        }
        cout << "  saved t=" << transient_times[k] << "\n";
    }

#ifndef SMOKE_TEST
    int ss = 0;
    for (; ss < STEADY_MAX_STEPS; ++ss) {
        for (int k = 0; k < NZ; ++k) {
            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    pre_temp[k][j][i] = T[k][j][i];
                }
            }
        }

        coupled_step(base_dt, Np, omega);

        if (ss % 20 == 0) {
            double err = 0.0;
            for (int k = 0; k < NZ; ++k) {
                for (int j = 0; j < NY; ++j) {
                    for (int i = 0; i < NX; ++i) {
                        if (fabs(T[k][j][i] - pre_temp[k][j][i]) > err) {
                            err = fabs(T[k][j][i] - pre_temp[k][j][i]);
                        }
                    }
                }
            }
            if (ss > 200 && err < STEADY_error) {
                break;
            }
        }
    }
    radiation_dugks(omega, true);
    save_curve(transient_count);
    if (do_fig3) {
        save_fig3_plane(transient_count);
    }
    cout << "  steady after " << ss << " extra steps\n";
#else
    save_curve(transient_count);
    if (do_fig3) {
        save_fig3_plane(transient_count);
    }
#endif

    if (do_fig3) {
        write_fig3_files(transient_times, transient_count);
    }

    if (!do_fig3) {
        write_curve_csv("output/" + name + ".csv", transient_times, transient_count);
        write_curve_csv("data/" + name + ".csv", transient_times, transient_count);
        cout << "wrote output/" << name << ".csv and data/" << name << ".csv\n";
    }
}

void reset_wall_overrides() {
    use_uniform_wall_temperature = false;
    uniform_wall_temperature = TE;
    use_uniform_wall_emissivity = false;
    uniform_wall_emissivity = 1.0;
}

void run_table1_case() {
    ensure_dir("output");
    reset_wall_overrides();
    EW = 1.0;
    TAU_L = 1.0;
    thermal_step_counter = 0;
    radiation_needs_initial_convergence = true;
    ini_case();

    const double target_time = 0.005;
    const double Ncr = 0.01;
    const double omega = 0.0;
    const double base_dt = min(compute_thermal_dt(), SUN2012_DT);
    double time = 0.0;

    while (time + 1.0e-14 < target_time) {
        double dt = min(base_dt, target_time - time);
        coupled_step(dt, Ncr, omega);
        time += dt;
    }

    radiation_dugks(omega, true);
    write_table1_csv("output/table1.csv");
    cout << "wrote output/table1.csv\n";
}

void test_angles() {
    double sum_w = 0.0;
    double sum_w_sx = 0.0;
    double sum_w_sy = 0.0;
    double sum_w_sz = 0.0;
    double sum_w_sx2 = 0.0;
    double sum_w_sy2 = 0.0;
    double sum_w_sz2 = 0.0;

    for (int m = 0; m < M; ++m) {
        double w = directions[m].weight;
        sum_w += w;
        sum_w_sx += w * directions[m].x;
        sum_w_sy += w * directions[m].y;
        sum_w_sz += w * directions[m].z;
        sum_w_sx2 += w * directions[m].x * directions[m].x;
        sum_w_sy2 += w * directions[m].y * directions[m].y;
        sum_w_sz2 += w * directions[m].z * directions[m].z;
    }

    cout << scientific << setprecision(12);
    cout << "sum_w," << sum_w << "\n";
    cout << "sum_w_sx," << sum_w_sx << "\n";
    cout << "sum_w_sy," << sum_w_sy << "\n";
    cout << "sum_w_sz," << sum_w_sz << "\n";
    cout << "sum_w_sx2," << sum_w_sx2 << "\n";
    cout << "sum_w_sy2," << sum_w_sy2 << "\n";
    cout << "sum_w_sz2," << sum_w_sz2 << "\n";
}

void test_isothermal_cavity(double eps, const string& label) {
    use_uniform_wall_temperature = true;
    uniform_wall_temperature = TE;
    use_uniform_wall_emissivity = true;
    uniform_wall_emissivity = eps;
    TAU_L = 1.0;
    EW = eps;
    thermal_step_counter = 0;
    radiation_needs_initial_convergence = true;
    ini_case();

    int rad_iter = radiation_dugks(0.0, true);

    double max_abs_J_minus_theta4 = 0.0;
    double max_abs_source = 0.0;
    double Cr = source_coeff(0.01, 0.0);
    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                double theta4 = pow4(T[k][j][i]);
                double err = fabs(G[k][j][i] - theta4);
                max_abs_J_minus_theta4 = max(max_abs_J_minus_theta4, err);
                max_abs_source = max(max_abs_source, fabs(Cr * (G[k][j][i] - theta4)));
            }
        }
    }

    cout << label << "\n";
    cout << scientific << setprecision(12);
    cout << "rad_iter," << rad_iter << "\n";
    cout << "max_abs_J_minus_theta4," << max_abs_J_minus_theta4 << "\n";
    cout << "max_abs_source," << max_abs_source << "\n";
    reset_wall_overrides();
}

void test_black_cavity() {
    test_isothermal_cavity(1.0, "--test-black-cavity");
}

void test_gray_cavity() {
    test_isothermal_cavity(0.5, "--test-gray-cavity");
}

void test_conduction_scaling() {
    const double tau_values[3] = {0.1, 1.0, 5.0};
    cout << "tauL,y005_theta,time,steps\n";

    for (int n = 0; n < 3; ++n) {
        reset_wall_overrides();
        TAU_L = tau_values[n];
        EW = 1.0;
        thermal_step_counter = 0;
        radiation_needs_initial_convergence = false;
        ini_case();

        double target_time = 2.0e-4;
        double time = 0.0;
        int steps = 0;
        while (time + 1.0e-14 < target_time && steps < 2000) {
            double dt = min(compute_thermal_dt(), target_time - time);
            dt = min(dt, max(1.0e-12, 0.25 * thermal_tau()));
            temperture(dt);
            time += dt;
            ++steps;
        }

        cout << fixed << setprecision(6) << tau_values[n] << ","
             << scientific << setprecision(12) << sample_T(0.05) << ","
             << time << "," << steps << "\n";
    }

    reset_wall_overrides();
    TAU_L = 1.0;
}

void test_source_newton() {
    struct CaseData {
        double gsum, J, dt, Cr;
    };
    const CaseData cases[4] = {
        {0.50, pow4(0.50), 1.0e-4, source_coeff(0.01, 0.0)},
        {0.55, pow4(0.80), 1.0e-4, source_coeff(0.01, 0.0)},
        {0.90, pow4(0.60), 5.0e-5, source_coeff(0.10, 0.0)},
        {0.70, pow4(0.90), 1.0e-4, source_coeff(0.01, 0.5)}
    };

    cout << "case,theta,residual,iterations\n";
    for (int n = 0; n < 4; ++n) {
        int iterations = 0;
        double theta = recover_theta_from_gtilde(
            cases[n].gsum, cases[n].J, cases[n].dt, cases[n].Cr, &iterations);
        double A = 0.5 * cases[n].dt * cases[n].Cr;
        double residual = theta - cases[n].gsum - A * (cases[n].J - pow4(theta));
        cout << n << "," << scientific << setprecision(12) << theta << ","
             << residual << "," << iterations << "\n";
    }
}

void run_sun2012_cases() {
    const double times_005_015[2] = {0.005, 0.015};
    const double times_005_015_050[3] = {0.005, 0.015, 0.050};
    const double times_005[1] = {0.005};

    ensure_dir("output");
    ensure_dir("data");
    cout << "3D DUGKS conduction-radiation | D3Q27 | "
         << NX << "x" << NY << "x" << NZ << " | "
         << NTHETA << "x" << MAX_PHI << " angles\n";
    cout << "==========================================\n";

#ifdef SMOKE_TEST
    cout << "\n=== Smoke test: Sun2012 baseline short transient ===\n";
#ifdef SMOKE_SUITE
    EW = 1.0;
    TAU_L = 1.0;
    run_one_case("smoke4a", 0.01, 0.5, times_005, 1);
    EW = 1.0;
    TAU_L = 1.0;
    run_one_case("smoke5a", 0.1, 0.0, times_005, 1);
    EW = 0.1;
    TAU_L = 1.0;
    run_one_case("smoke6a", 0.01, 0.0, times_005, 1);
    EW = 1.0;
    TAU_L = 0.1;
    run_one_case("smoke7a", 0.01, 0.0, times_005, 1);
    return;
#endif
#ifdef SMOKE_CONDUCTION_TEST
    EW = 1.0;
    TAU_L = 1.0;
    run_one_case("smoke_conduction", 1.0, 1.0, times_005, 1);
    return;
#endif
#ifdef SMOKE_FIG3_TEST
    EW = 1.0;
    TAU_L = 1.0;
    run_one_case("smoke_fig3", 0.01, 0.0, times_005_015_050, 3, true);
    return;
#endif
#ifdef SMOKE_OPTICAL_TEST
    EW = 1.0;
    TAU_L = 0.1;
    run_one_case("smoke7a", 0.01, 0.0, times_005, 1);
    return;
#endif
    EW = 1.0;
    TAU_L = 1.0;
    run_one_case("smoke3d", 0.01, 0.0, times_005, 1);
    return;
#endif

    cout << "\n=== Sun2012 Table 1 ===\n";
    run_table1_case();

    cout << "\n=== Sun2012 Fig.3: Center-Plane Isotherms ===\n";
    EW = 1.0;
    TAU_L = 1.0;
    run_one_case("fig3", 0.01, 0.0, times_005_015_050, 3, true);

    cout << "\n=== Sun2012 Fig.4: Scattering Albedo ===\n";
    EW = 1.0;
    TAU_L = 1.0;
    cout << "--- 4a: omega=0.5 ---\n";
    run_one_case("fig4_omega05", 0.01, 0.5, times_005_015, 2);
    cout << "--- 4b: omega=0.9 ---\n";
    run_one_case("fig4_omega09", 0.01, 0.9, times_005_015_050, 3);

    cout << "\n=== Sun2012 Fig.5: N parameter ===\n";
    EW = 1.0;
    TAU_L = 1.0;
    cout << "--- 5a: N=0.1 ---\n";
    run_one_case("fig5_ncr01", 0.1, 0.0, times_005_015_050, 3);
    cout << "--- 5b: N=1.0 ---\n";
    run_one_case("fig5_ncr10", 1.0, 0.0, times_005_015_050, 3);

    cout << "\n=== Sun2012 Fig.6: Wall Emissivity ===\n";
    TAU_L = 1.0;
    cout << "--- 6a: eps_s=0.1 ---\n";
    EW = 0.1;
    run_one_case("fig6_eps01", 0.01, 0.0, times_005_015, 2);
    cout << "--- 6b: eps_s=0.5 ---\n";
    EW = 0.5;
    run_one_case("fig6_eps05", 0.01, 0.0, times_005_015, 2);

    cout << "\n=== Sun2012 Fig.7: Optical Thickness ===\n";
    EW = 1.0;
    cout << "--- 7a: tau_L=0.1 ---\n";
    TAU_L = 0.1;
    run_one_case("fig7_tauL01", 0.01, 0.0, times_005_015_050, 3);
    cout << "--- 7b: tau_L=5.0 ---\n";
    TAU_L = 5.0;
    run_one_case("fig7_tauL50", 0.01, 0.0, times_005, 1);

    cout << "\n==========================================\n";
    cout << "All Sun2012 cases complete.\n";
    cout << "Table 1:    output/table1.csv\n";
    cout << "Fig.3 .csv: output/fig3_t0005.csv ... output/fig3_steady.csv\n";
    cout << "Fig.4-7:    output/fig[4-7]_*.csv\n";
}

void print_usage(const char* program) {
    cout << "Usage: " << program << " [command]\n";
    cout << "Commands:\n";
    cout << "  --test-angles\n";
    cout << "  --test-black-cavity\n";
    cout << "  --test-gray-cavity\n";
    cout << "  --test-conduction\n";
    cout << "  --test-source-newton\n";
    cout << "  --test-table1\n";
    cout << "  --run-sun2012\n";
}

int main(int argc, char* argv[]) {
    ini_velocities();
    ini_angles();

    if (argc >= 2) {
        string command = argv[1];
        if (command == "--test-angles") {
            test_angles();
            return 0;
        }
        if (command == "--test-black-cavity") {
            test_black_cavity();
            return 0;
        }
        if (command == "--test-gray-cavity") {
            test_gray_cavity();
            return 0;
        }
        if (command == "--test-conduction") {
            test_conduction_scaling();
            return 0;
        }
        if (command == "--test-source-newton") {
            test_source_newton();
            return 0;
        }
        if (command == "--test-table1") {
            run_table1_case();
            return 0;
        }
        if (command == "--run-sun2012") {
            run_sun2012_cases();
            return 0;
        }
        if (command == "--help" || command == "-h") {
            print_usage(argv[0]);
            return 0;
        }

        cerr << "Unknown command: " << command << "\n";
        print_usage(argv[0]);
        return 1;
    }

    run_sun2012_cases();
    return 0;
}
