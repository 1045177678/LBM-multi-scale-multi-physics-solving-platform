// lbm_core.cpp (clean rebuild)
// -----------------------------------------------------------------------------
// 实现 lbm_core.h 中声明的接口：
//  • TimeExpr 工厂
//  • D2Q9 的 6-矩投影/重构
//  • flow/scalar 字段（BGK/MRT 简化版）
//  • 统一 Operators（速度/无滑移/Dirichlet/体力/浮力/脉冲热源/占位桥）
//  • SimulationEngine 的分层推进 r^2 子循环 + 细↔粗 矩空间限制/延拓
//  • 单位系统 init_units_from_config（发布 L0/T0/Θ0/Tref 与偏好）
// -----------------------------------------------------------------------------
// lbm_core.cpp (clean rebuild)
// -----------------------------------------------------------------------------
// Implements interfaces declared in lbm_core.h:
// • TimeExpr factory
// • 6-moment projection/reconstruction for D2Q9
// • flow/scalar fields (BGK/MRT simplified)
// • Unified Operators (velocity/no-slip/Dirichlet/body/buoyancy/pulsed heat source/occupancy bridge)
// • SimulationEngine's hierarchical boosting r^2 sub-loop + fine↔coarse moment space constraints/extensions
// • Unit system init_units_from_config (publishes L0/T0/Θ0/Tref and preferences)
// -----------------------------------------------------------------------------
#include "lbm_core.h"

#if __has_include(<nlohmann/json.hpp>)
#include <nlohmann/json.hpp>
#elif __has_include("json.hpp")
#include "json.hpp"
#else
#error "Cannot find nlohmann/json.hpp header."
#endif
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <array>
#include <unordered_map>

// ========== TimeExpr 工厂 ==========
std::unique_ptr<ITimeExpr> make_timeexpr_from_json(const Json& j) {
    if (j.is_number_float() || j.is_number_integer())
        return std::make_unique<ConstExpr>(j.get<double>());
    if (j.contains("table")) {
        auto te = std::make_unique<TableExpr>();
        te->linear = (j.value("method", std::string("linear")) != std::string("step"));
        for (auto& p : j.at("table")) te->tbl.emplace_back(p[0].get<double>(), p[1].get<double>());
        std::sort(te->tbl.begin(), te->tbl.end(), [](auto& a, auto& b) {return a.first < b.first; });
        return te;
    }
    if (j.contains("expr")) {
        auto fe = std::make_unique<FuncExpr>();
        fe->expr = j.at("expr").get<std::string>();
        if (j.contains("const")) for (auto it = j["const"].begin(); it != j["const"].end(); ++it)
            fe->consts[it.key()] = it.value().get<double>();
        return fe;
    }
    return std::make_unique<ConstExpr>(0.0);
}

// ========== D2Q9  ==========
namespace d2q9 {
    static const int Q = 9;
    static const int ex[Q] = { 0,1,0,-1,0, 1,-1,-1, 1 };
    static const int ey[Q] = { 0,0,1, 0,-1,1, 1,-1,-1 };
    static const double w[Q] = { 4.0 / 9.0,
                                1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,
                                1.0 / 36.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0 };
    static constexpr double cs2 = 1.0 / 3.0;
}

// ========== 6-M：P 与 P^+ ==========
namespace {
    bool g_proj_init = false; double P6x9[6][9]; double Pinv9x6[9][6];
    void build_projection_mats() {
        using namespace d2q9;
        for (int i = 0; i < 9; ++i) {
            double exi = ex[i], eyi = ey[i]; double ex2 = exi * exi, ey2 = eyi * eyi;
            P6x9[0][i] = 1.0;
            P6x9[1][i] = exi;
            P6x9[2][i] = eyi;
            P6x9[3][i] = (ex2 + ey2) - 2.0 * cs2;
            P6x9[4][i] = (ex2 - ey2);
            P6x9[5][i] = (exi * eyi);
        }
        double A[6][6]{}; // P P^T
        for (int r = 0; r < 6; ++r) for (int c = 0; c < 6; ++c) 
        {
           double s = 0;
           for (int k = 0; k < 9; ++k)
           {
               s += P6x9[r][k] * P6x9[c][k]; A[r][c] = s;
           }
        }
        // 6x6 inv
        double AI[6][12]{};
        for (int r = 0; r < 6; ++r) 
        {
            for (int c = 0; c < 6; ++c) AI[r][c] = A[r][c];
            AI[r][6 + r] = 1.0; 
        }
        for (int col = 0; col < 6; ++col) 
        {
            int piv = col; double best = std::fabs(AI[piv][col]);
            for (int r = col + 1; r < 6; ++r) { double v = std::fabs(AI[r][col]); if (v > best) { best = v; piv = r; } }
            if (best < 1e-12) { AI[col][col] += 1e-12; piv = col; }
            if (piv != col) for (int c = 0; c < 12; ++c) std::swap(AI[piv][c], AI[col][c]);
            double d = AI[col][col]; if (std::fabs(d) < 1e-18) d = (d >= 0 ? 1 : -1) * 1e-18;
            for (int c = 0; c < 12; ++c) AI[col][c] /= d;
            for (int r = 0; r < 6; ++r) { if (r == col) continue; double f = AI[r][col]; if (std::fabs(f) < 1e-18) continue; for (int c = 0; c < 12; ++c) AI[r][c] -= f * AI[col][c]; }
        }
        double Ainv[6][6]{}; for (int r = 0; r < 6; ++r) for (int c = 0; c < 6; ++c) Ainv[r][c] = AI[r][6 + c];
        for (int i = 0; i < 9; ++i) for (int k = 0; k < 6; ++k) { double s = 0; for (int r = 0; r < 6; ++r) s += P6x9[r][i] * Ainv[r][k]; Pinv9x6[i][k] = s; }
        g_proj_init = true;
    }
}
void project6_from_f(const double f[9], double m[6]) {
    if (!g_proj_init) build_projection_mats();
    for (int r = 0; r < 6; ++r) { double s = 0; for (int i = 0; i < 9; ++i) s += P6x9[r][i] * f[i]; m[r] = s; }
}
void reconstruct_f_from6(const double m[6], double f[9]) {
    if (!g_proj_init) build_projection_mats();
    for (int i = 0; i < 9; ++i) { double s = 0; for (int k = 0; k < 6; ++k) s += Pinv9x6[i][k] * m[k]; f[i] = s; }
}

// ========== macro ==========
static inline void macro_flow(const std::vector<double>& f, int E,
    std::vector<double>& rho,
    std::vector<double>& ux,
    std::vector<double>& uy) 
{
    using namespace d2q9; rho.assign(E, 0); ux.assign(E, 0); uy.assign(E, 0);
    for (int e = 0; e < E; ++e) {
        double r = 0, u = 0, v = 0; int off = e * Q; for (int i = 0; i < Q; ++i) { double fi = f[off + i]; r += fi; u += fi * ex[i]; v += fi * ey[i]; }
        if (r > 1e-30) { rho[e] = r; ux[e] = u / r; uy[e] = v / r; }
        else { rho[e] = 0; ux[e] = uy[e] = 0; }
    }
}
static inline void macro_scalar(const std::vector<double>& g, int E,
    const std::vector<double>& /*ux*/, const std::vector<double>& /*uy*/,
    std::vector<double>& phi) {
    using namespace d2q9; phi.assign(E, 0); for (int e = 0; e < E; ++e) { double s = 0; int off = e * Q; for (int i = 0; i < Q; ++i) s += g[off + i]; phi[e] = s; }
}

// ========== Field base class ==========
struct FieldBase : IField {
    std::string nm, kd; bool out{ true }; int L{ 0 }; const GridLevel* GL{ nullptr }; int E{ 0 };
    std::vector<double> rho, ux, uy; // flow 
    std::vector<double> phi;         // scalar 
    const std::vector<double>* FX{ nullptr };
    const std::vector<double>* FY{ nullptr };
    const std::vector<double>* SRC{ nullptr };
    FieldBase(std::string n, std::string k, int l, const GridLevel* gl) : nm(std::move(n)), kd(std::move(k)), L(l), GL(gl) { E = gl ? gl->NxC * gl->NyC : 0; }
    std::string name() const override { return nm; }
    std::string kind() const override { return kd; }
    bool output() const override { return out; }
    int level_index() const override { return L; }
    int cells() const override { return E; }
    void bindBuffers(const ForceBuffer& fb) { FX = &fb.fx; FY = &fb.fy; SRC = &fb.src; }
};

// ========== Flow: D2Q9-BGK ==========
struct FlowD2Q9_BGK : FieldBase {
    using FieldBase::FieldBase;
    double tau{ 0.8 };
    std::vector<double> f; // 

    FlowD2Q9_BGK(const std::string& name, int L, const GridLevel* gl)
        : FieldBase(name, "flow", L, gl) {
        f.assign(E * 9, 0.0);
        rho.assign(E, 1.0); ux.assign(E, 0.0); uy.assign(E, 0.0);
        using namespace d2q9; for (int e = 0; e < E; ++e) { int off = e * 9; for (int i = 0; i < 9; ++i) f[off + i] = w[i]; }
    }
    const std::vector<double>& dist() const override { return f; }
    std::vector<double>& dist() override { return f; }

    void collide(double dt) override {
        using namespace d2q9; const double omega = 1.0 / tau;
        for (int e = 0; e < E; ++e) {
            int off = e * 9; double r = rho[e], u = ux[e], v = uy[e];
            double Fx = FX ? (*FX)[e] : 0.0, Fy = FY ? (*FY)[e] : 0.0;
            const double uu = u * u + v * v;
            for (int i = 0; i < 9; ++i) {
                const double eiu = ex[i] * u + ey[i] * v;
                const double feq = w[i] * r * (1.0 + 3.0 * eiu + 4.5 * eiu * eiu - 1.5 * uu);
                // Guo force
                const double Fdot = (ex[i] * Fx + ey[i] * Fy);
                const double Fi = w[i] * (3.0 * Fdot + 9.0 * eiu * Fdot - 3.0 * (u * Fx + v * Fy));
                f[off + i] = f[off + i] - omega * (f[off + i] - feq) + (1.0 - 0.5 * omega) * Fi * dt;
            }
        }
    }
    void stream() override {
        using namespace d2q9; std::vector<double> fn(f);
        int Nx = GL->NxC, Ny = GL->NyC;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int e = j * Nx + i; int off = e * 9;
                for (int k = 0; k < 9; ++k) {
                    int ii = i - ex[k], jj = j - ey[k];
                    if (ii >= 0 && ii < Nx && jj >= 0 && jj < Ny) { int en = jj * Nx + ii; fn[off + k] = f[en * 9 + k]; }
                }
            }
        }
        f.swap(fn);
    }
    void compute_macro() override { macro_flow(f, E, rho, ux, uy); }
};

// ========== Flow: D2Q9-MRT(To be improved) ==========
struct FlowD2Q9_MRT : FieldBase {
    using FieldBase::FieldBase;
    double nu{ 0.1 }; // Kinematic viscosity（LU）
    double dt_last{ 1.0 };
    std::vector<double> f;
    double M[9][9]{}; double Minv[9][9]{}; double S[9]{};

    FlowD2Q9_MRT(const std::string& name, int L, const GridLevel* gl)
        : FieldBase(name, "flow", L, gl) {
        f.assign(E * 9, 0.0);
        rho.assign(E, 1.0); ux.assign(E, 0.0); uy.assign(E, 0.0);
        using namespace d2q9; for (int e = 0; e < E; ++e) { int off = e * 9; for (int i = 0; i < 9; ++i) f[off + i] = w[i]; }
        // (To be improved)
        const double M_[9][9] = { {1,1,1,1,1,1,1,1,1},{-4,-1,-1,-1,-1,2,2,2,2},{4,-2,-2,-2,-2,1,1,1,1},{0,1,0,-1,0,1,-1,-1,1},{0,-2,0,2,0,1,-1,-1,1},{0,0,1,0,-1,1,1,-1,-1},{0,0,-2,0,2,1,1,-1,-1},{0,1,-1,1,-1,0,0,0,0},{0,0,0,0,0,1,-1,1,-1} };
        const double Minv_[9][9] = { {1.0 / 9,-1.0 / 9,1.0 / 9,0,0,0,0,0,0},{1.0 / 9,-1.0 / 9,1.0 / 9,1.0 / 6,-1.0 / 6,0,0,1.0 / 4,0},{1.0 / 9,-1.0 / 9,1.0 / 9,0,0,1.0 / 6,-1.0 / 6,-1.0 / 4,0},{1.0 / 9,-1.0 / 9,1.0 / 9,-1.0 / 6,1.0 / 6,0,0,1.0 / 4,0},{1.0 / 9,-1.0 / 9,1.0 / 9,0,0,-1.0 / 6,1.0 / 6,-1.0 / 4,0},{1.0 / 9,2.0 / 9,1.0 / 9,1.0 / 6,1.0 / 12,1.0 / 6,1.0 / 12,0,1.0 / 4},{1.0 / 9,2.0 / 9,1.0 / 9,-1.0 / 6,-1.0 / 12,1.0 / 6,1.0 / 12,0,-1.0 / 4},{1.0 / 9,2.0 / 9,1.0 / 9,-1.0 / 6,-1.0 / 12,-1.0 / 6,-1.0 / 12,0,1.0 / 4},{1.0 / 9,2.0 / 9,1.0 / 9,1.0 / 6,1.0 / 12,-1.0 / 6,-1.0 / 12,0,-1.0 / 4} };
        for (int i = 0; i < 9; ++i) { for (int j = 0; j < 9; ++j) { M[i][j] = M_[i][j]; Minv[i][j] = Minv_[i][j]; } }
    }
    const std::vector<double>& dist() const override { return f; }
    std::vector<double>& dist() override { return f; }
    void collide(double dt) override {
        using namespace d2q9; dt_last = dt; double cs2 = d2q9::cs2; double tau = 0.5 + nu / (cs2 * dt); double s_nu = 1.0 / tau;
        double s[9] = { 1.7,1.1,1.0, s_nu, s_nu, 1.2,1.2, 1.0,1.0 }; for (int i = 0; i < 9; ++i) S[i] = s[i];
        for (int e = 0; e < E; ++e) {
            int off = e * 9; double r = 0, u = 0, v = 0; for (int i = 0; i < 9; ++i) { double fi = f[off + i]; r += fi; u += fi * ex[i]; v += fi * ey[i]; }
            if (r > 1e-30) { u /= r; v /= r; }
            else { u = v = 0; }
            rho[e] = r; ux[e] = u; uy[e] = v;
            double m[9]{}; for (int a = 0; a < 9; ++a) for (int i = 0; i < 9; ++i) m[a] += M[a][i] * f[off + i];
            double u2 = u * u + v * v; double meq[9]{}; meq[0] = r; meq[1] = -2 * r + 3 * r * u2; meq[2] = r - 3 * r * u2; meq[3] = r * u; meq[5] = r * v; meq[4] = -r * u; meq[6] = -r * v; meq[7] = r * (u * u - v * v); meq[8] = r * u * v;
            double Fx = FX ? (*FX)[e] : 0.0, Fy = FY ? (*FY)[e] : 0.0; double Fm[9]{}; Fm[3] = Fx; Fm[5] = Fy;
            for (int a = 0; a < 9; ++a) m[a] = m[a] - S[a] * (m[a] - meq[a]) + (1.0 - 0.5 * S[a]) * Fm[a] * dt;
            for (int i = 0; i < 9; ++i) { double sum = 0; for (int a = 0; a < 9; ++a) sum += Minv[i][a] * m[a]; f[off + i] = sum; }
        }
    }
    void stream() override {
        using namespace d2q9; std::vector<double> fn(f);
        int Nx = GL->NxC, Ny = GL->NyC;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int e = j * Nx + i; int off = e * 9;
                for (int k = 0; k < 9; ++k) { int ii = i - ex[k], jj = j - ey[k]; if (ii >= 0 && ii < Nx && jj >= 0 && jj < Ny) { int en = jj * Nx + ii; fn[off + k] = f[en * 9 + k]; } }
            }
        }
        f.swap(fn);
    }
    void compute_macro() override { macro_flow(f, E, rho, ux, uy); }
};

// ========== Scalar: D2Q9-BGK ==========
struct ScalarD2Q9_BGK : FieldBase {
    using FieldBase::FieldBase;
    std::string advect_by; // 
    double tau{ 0.8 };
    std::vector<double> g;

    ScalarD2Q9_BGK(const std::string& name, int L, const GridLevel* gl, std::string adv)
        : FieldBase(name, "scalar", L, gl), advect_by(std::move(adv)) {
        g.assign(E * 9, 0.0);
        phi.assign(E, 0.0); ux.assign(E, 0.0); uy.assign(E, 0.0);
    }
    const std::vector<double>& dist() const override { return g; }
    std::vector<double>& dist() override { return g; }
    void collide(double dt) override {
        using namespace d2q9; const double omega = 1.0 / tau;
        for (int e = 0; e < E; ++e) {
            int off = e * 9; double ph = phi[e], u = ux[e], v = uy[e];
            double Sval = SRC ? (*SRC)[e] : 0.0; // source term
            for (int i = 0; i < 9; ++i) {
                double eiu = ex[i] * u + ey[i] * v;
                double geq = w[i] * ph * (1.0 + 3.0 * eiu);
                double Si = w[i] * Sval * (1.0 + 3.0 * eiu);
                g[off + i] = g[off + i] - omega * (g[off + i] - geq) + (1.0 - 0.5 * omega) * Si * dt;
            }
        }
    }
    void stream() override {
        using namespace d2q9; std::vector<double> gn(g);
        int Nx = GL->NxC, Ny = GL->NyC;
        for (int j = 0; j < Ny; ++j) 
        {
            for (int i = 0; i < Nx; ++i) 
            {
                int e = j * Nx + i; int off = e * 9;
                for (int k = 0; k < 9; ++k) 
                { 
                    int ii = i - ex[k], jj = j - ey[k];
                    if (ii >= 0 && ii < Nx && jj >= 0 && jj < Ny) 
                    { 
                        int en = jj * Nx + ii; gn[off + k] = g[en * 9 + k];
                    }
                } 
            }
        }
        g.swap(gn);
    }
    void compute_macro() override { macro_scalar(g, E, ux, uy, phi); }
};
struct ScalarD2Q9_MRT : ScalarD2Q9_BGK { using ScalarD2Q9_BGK::ScalarD2Q9_BGK; };

// ========== FieldFactory ==========
static std::unordered_map<std::string, FieldFactory::Creator>& field_reg() { static std::unordered_map<std::string, FieldFactory::Creator> R; return R; }
void FieldFactory::reg(const std::string& key, Creator c) { field_reg()[key] = std::move(c); }
FieldPtr FieldFactory::create(const Json& j, const GridHierarchy& gh) {
    std::string name = j.at("name"); std::string kind = j.at("kind"); std::string lat = j.at("lattice"); std::string meth = j.at("method"); bool out = j.value("output", true); int L = j.value("level", 0); const GridLevel* gl = &gh.levels.at(L);
    std::string key = kind + "." + lat + "." + meth; auto it = field_reg().find(key); if (it == field_reg().end()) throw std::runtime_error("Unknown field key: " + key);
    FieldPtr f = it->second(j, gh); if (auto* fb = dynamic_cast<FieldBase*>(f.get())) fb->out = out; return f;
}
void register_builtin_fields() {
    FieldFactory::reg("flow.D2Q9.BGK", [](const Json& j, const GridHierarchy& gh) 
        {
            int L = j.value("level", 0); auto p = std::make_unique<FlowD2Q9_BGK>(j.at("name"), L, &gh.levels.at(L)); 
            if (j.contains("tau")) p->tau = j["tau"];
            if (j.contains("nu") && j.contains("dt"))
            {
                double nu = j["nu"], dt = j["dt"];
                p->tau = 0.5 + nu / (d2q9::cs2 * dt);
            } return p;
        }
    );
    
    //////////////
    FieldFactory::reg("flow.D2Q9.MRT", [](const Json& j, const GridHierarchy& gh) { int L = j.value("level", 0); auto p = std::make_unique<FlowD2Q9_MRT>(j.at("name"), L, &gh.levels.at(L)); if (j.contains("nu")) p->nu = j["nu"]; return p; });
   
   /////////////
    FieldFactory::reg("scalar.D2Q9.BGK", [](const Json& j, const GridHierarchy& gh) { int L = j.value("level", 0); std::string adv = j.value("advect_by", ""); auto p = std::make_unique<ScalarD2Q9_BGK>(j.at("name"), L, &gh.levels.at(L), adv); if (j.contains("tau")) p->tau = j["tau"]; if (j.contains("D") && j.contains("dt")) { double D = j["D"], dt = j["dt"]; p->tau = 0.5 + D / (d2q9::cs2 * dt); } return p; });
   
    /////////////
    FieldFactory::reg("scalar.D2Q9.MRT", [](const Json& j, const GridHierarchy& gh) { int L = j.value("level", 0); std::string adv = j.value("advect_by", ""); return std::make_unique<ScalarD2Q9_MRT>(j.at("name"), L, &gh.levels.at(L), adv); });
}

// ========== 运行期单位（供算子读取） ==========
// ========== Runtime unit (for operator reading) ==========
namespace units_runtime {
    // 运行时单位换算所需的层级尺度与输入偏好。
    // Hierarchy scale and input preferences required for runtime unit conversion.
    struct Prefs { bool force_is_phys{ true }; bool source_is_phys{ true }; bool bc_is_phys{ true }; };
    std::vector<double> L0s, T0s, Theta0s, Trefs; // base units
    Prefs prefs; // 输入偏好（phys=按SI输入，内部自动换算到LU）// Input preference (phys=input in SI, automatically converted to LU internally)
}

// ========== Operators ==========
struct OperatorBase : IOperator { int L{ 0 }; int targetLevel() const override { return L; } };

struct VelocityBC : OperatorBase {
    std::string field; std::vector<int> cells; std::unique_ptr<ITimeExpr> ux, uy;
    VelocityBC(const Json& j) { L = j.at("level"); field = j.at("field"); for (auto& c : j.at("cells")) cells.push_back((int)c); ux = make_timeexpr_from_json(j.at("ux")); uy = make_timeexpr_from_json(j.at("uy")); }
    void apply(ApplyCtx& ctx) override {
        auto* F = ctx.fieldByName(field); auto* fb = dynamic_cast<FieldBase*>(F); if (!fb || F->kind() != "flow" || F->level_index() != L) return; double t = ctx.t;
        for (int e : cells) {
            if (e < 0 || e >= fb->E) continue; double u = ux->eval(t), v = uy->eval(t);
            if (units_runtime::prefs.bc_is_phys) { double T0 = ctx.dt, L0 = ctx.gh->levels[L].dx; u *= T0 / L0; v *= T0 / L0; }
            fb->ux[e] = u; fb->uy[e] = v; fb->rho[e] = 1.0;
            if (auto* fbgk = dynamic_cast<FlowD2Q9_BGK*>(F)) { using namespace d2q9; int off = e * 9; double r = 1.0; double uu = u * u + v * v; for (int i = 0; i < 9; ++i) { double eiu = ex[i] * u + ey[i] * v; fbgk->f[off + i] = w[i] * r * (1 + 3 * eiu + 4.5 * eiu * eiu - 1.5 * uu); } }
            if (auto* fmrt = dynamic_cast<FlowD2Q9_MRT*>(F)) { using namespace d2q9; int off = e * 9; double r = 1.0; double uu = u * u + v * v; for (int i = 0; i < 9; ++i) { double eiu = ex[i] * u + ey[i] * v; fmrt->f[off + i] = w[i] * r * (1 + 3 * eiu + 4.5 * eiu * eiu - 1.5 * uu); } }
        }
    }
};
struct NoSlipBC : OperatorBase 
{ 
    std::string field; std::vector<int> cells; NoSlipBC(const Json& j) 
    {
        L = j.at("level"); field = j.at("field");
        for (auto& c : j.at("cells")) cells.push_back((int)c); 
    } 
    void apply(ApplyCtx& ctx) override 
    { 
        auto* F = ctx.fieldByName(field);
        auto* fb = dynamic_cast<FieldBase*>(F); 
        if (!fb || F->kind() != "flow" || F->level_index() != L) return;
        using namespace d2q9; 
        for (int e : cells) 
        {
            if (e < 0 || e >= fb->E) continue;
            fb->ux[e] = fb->uy[e] = 0;
            fb->rho[e] = 1; 
            if (auto* fbgk = dynamic_cast<FlowD2Q9_BGK*>(F)) 
            {
                int off = e * 9; 
                for (int i = 0; i < 9; ++i) fbgk->f[off + i] = w[i];
            } 
            if (auto* fmrt = dynamic_cast<FlowD2Q9_MRT*>(F))
            { 
                int off = e * 9;
                for (int i = 0; i < 9; ++i) fmrt->f[off + i] = w[i];
            } 
        }
    }
};
struct DirichletScalarBC : OperatorBase 
{ 
    std::string field; std::vector<int> cells; std::unique_ptr<ITimeExpr> value; 
    DirichletScalarBC(const Json& j) 
    {
        L = j.at("level"); field = j.at("field");
        for (auto& c : j.at("cells")) cells.push_back((int)c);
        value = make_timeexpr_from_json(j.at("value"));
    } 
    void apply(ApplyCtx& ctx) override 
    {
     auto* F = ctx.fieldByName(field); 
    auto* fb = dynamic_cast<FieldBase*>(F); 
    if (!fb || F->kind() != "scalar" || F->level_index() != L) return; 
    double val = value->eval(ctx.t); 
    if (units_runtime::prefs.bc_is_phys) 
    { 
        double Theta0 = (L < (int)units_runtime::Theta0s.size() ? units_runtime::Theta0s[L] : 1.0); 
        double Tref = (L < (int)units_runtime::Trefs.size() ? units_runtime::Trefs[L] : 0.0); 
        val = (val - Tref) / Theta0;
    } 
    for (int e : cells) 
    {
        if (e < 0 || e >= fb->E) continue; fb->phi[e] = val;
        if (auto* gbgk = dynamic_cast<ScalarD2Q9_BGK*>(F))
        {
            using namespace d2q9; int off = e * 9; double u = fb->ux[e], v = fb->uy[e];
            for (int i = 0; i < 9; ++i) 
            { double eiu = ex[i] * u + ey[i] * v;
            gbgk->g[off + i] = w[i] * val * (1 + 3 * eiu);
            }
        }
    } 
    } 
};
struct BodyForce : OperatorBase 
{
    std::string field; std::unique_ptr<ITimeExpr> fx, fy, src;
    BodyForce(const Json& j) 
    {
        L = j.at("level"); field = j.at("field"); if (j.contains("fx")) fx = make_timeexpr_from_json(j.at("fx"));
        if (j.contains("fy")) fy = make_timeexpr_from_json(j.at("fy"));
        if (j.contains("src")) src = make_timeexpr_from_json(j.at("src"));
    } 
void apply(ApplyCtx& ctx) override 
{ 
auto& FB = ctx.forceOfLevel(L); 
double t = ctx.t; int E = (int)FB.fx.size();
if (fx || fy) 
{ 
    double T0 = ctx.dt, L0 = ctx.gh->levels[L].dx;
    double facF = units_runtime::prefs.force_is_phys ? (T0 * T0 / L0) : 1.0;
    for (int e = 0; e < E; ++e)
    {
        if (fx) FB.fx[e] += fx->eval(t) * facF;
        if (fy) FB.fy[e] += fy->eval(t) * facF; 
    }
} 
if (src) 
{ 
    double T0 = ctx.dt; 
    double Theta0 = (L < (int)units_runtime::Theta0s.size() ? units_runtime::Theta0s[L] : 1.0);
double facS = units_runtime::prefs.source_is_phys ? (T0 / Theta0) : 1.0; for (int e = 0; e < E; ++e) FB.src[e] += src->eval(t) * facS;
} 
}
};
struct Boussinesq : OperatorBase 
{ std::string target, reads; double beta{ 0 };
Vec2 g{ 0,-9.81 }; double Tref{ 300 }; 
Boussinesq(const Json& j)
{
    L = j.at("level"); target = j.at("target"); reads = j.at("reads")[0]; beta = j.at("beta"); Tref = j.at("T_ref"); auto arr = j.at("g"); g.x = arr[0]; g.y = arr[1];
} 
void apply(ApplyCtx& ctx) override 
{ auto* Tf = ctx.fieldByName(reads); 
auto* Fflow = ctx.fieldByName(target);
auto* Tfb = dynamic_cast<FieldBase*>(Tf); 
auto* Ffb = dynamic_cast<FieldBase*>(Fflow);
if (!Tfb || !Ffb || Tfb->level_index() != L || Ffb->level_index() != L) return; 
auto& FB = ctx.forceOfLevel(L);
int E = Tfb->E; double T0 = ctx.dt, L0 = ctx.gh->levels[L].dx; 
double Theta0 = (L < (int)units_runtime::Theta0s.size() ? units_runtime::Theta0s[L] : 1.0);
double gx = g.x * (units_runtime::prefs.force_is_phys ? (T0 * T0 / L0) : 1.0); 
double gy = g.y * (units_runtime::prefs.force_is_phys ? (T0 * T0 / L0) : 1.0);
for (int e = 0; e < E; ++e) 
{
    double dT_si = (Tfb->phi.empty() ? 0.0 : (Tfb->phi[e] * Theta0)); 
    FB.fx[e] += beta * dT_si * gx; FB.fy[e] += beta * dT_si * gy; 
}
} 
};

struct HeatPulse : OperatorBase 
{
    std::string field; std::unique_ptr<ITimeExpr> power; double spotR{ 1.0 };
    HeatPulse(const Json& j) 
    { 
        L = j.at("level"); field = j.at("field"); power = make_timeexpr_from_json(j.at("power"));
        if (j.contains("spot_radius")) spotR = j.at("spot_radius"); 
    }
    void apply(ApplyCtx& ctx) override 
    {
        auto* F = ctx.fieldByName(field); auto* fb = dynamic_cast<FieldBase*>(F);
    if (!fb || F->kind() != "scalar" || F->level_index() != L) return;
    auto& FB = ctx.forceOfLevel(L); double P = power->eval(ctx.t);
    if (units_runtime::prefs.source_is_phys) 
    { 
        double T0 = ctx.dt; double Theta0 = (L < (int)units_runtime::Theta0s.size() ? units_runtime::Theta0s[L] : 1.0); P *= T0 / Theta0; 
    } 
    double cx = fb->GL->x0c + (fb->GL->NxC - 1) * 0.5 * fb->GL->dx, cy = fb->GL->y0c + (fb->GL->NyC - 1) * 0.5 * fb->GL->dy;
    for (int e = 0; e < fb->E; ++e) 
    { 
        Vec2 c = fb->GL->cellCenter(e);
        double r2 = (c.x - cx) * (c.x - cx) + (c.y - cy) * (c.y - cy);
        if (r2 <= spotR * spotR) FB.src[e] += P; 
    }
    }
};
struct MethodBridge : OperatorBase { void apply(ApplyCtx&) override {/* 细↔粗桥接在引擎 advanceLevel 的对齐处完成   */ } };

static std::unordered_map<std::string, OperatorFactory::Creator>& op_reg() { static std::unordered_map<std::string, OperatorFactory::Creator> R; return R; }
void OperatorFactory::reg(const std::string& tp, Creator c) { op_reg()[tp] = std::move(c); }
std::unique_ptr<IOperator> OperatorFactory::create(const Json& j) { std::string tp = j.at("type"); auto it = op_reg().find(tp); if (it == op_reg().end()) throw std::runtime_error("Unknown operator: " + tp); return it->second(j); }
void register_builtin_operators() {
    OperatorFactory::reg("velocity_bc", [](const Json& j) { return std::make_unique<VelocityBC>(j); });
    OperatorFactory::reg("noslip_bc", [](const Json& j) { return std::make_unique<NoSlipBC>(j); });
    OperatorFactory::reg("dirichlet_bc", [](const Json& j) { return std::make_unique<DirichletScalarBC>(j); });
    OperatorFactory::reg("body_force", [](const Json& j) { return std::make_unique<BodyForce>(j); });
    OperatorFactory::reg("boussinesq", [](const Json& j) { return std::make_unique<Boussinesq>(j); });
    OperatorFactory::reg("heat_pulse", [](const Json& j) { return std::make_unique<HeatPulse>(j); });
    OperatorFactory::reg("method_bridge", [](const Json& j) { (void)j; return std::make_unique<MethodBridge>(); });
}

// ========== 引擎实现 ==========
// ========== Engine Implementation ==========
IField* SimulationEngine::fieldByName(const std::string& nm) const { for (auto& f : fields) if (f->name() == nm) return f.get(); return nullptr; }
void SimulationEngine::buildForceBuffers() 
{ 
    forces.resize(gh.levels.size());
    for (size_t L = 0; L < gh.levels.size(); ++L) forces[L].resize(gh.levels[L].NxC * gh.levels[L].NyC);
    for (auto& f : fields) if (auto* fb = dynamic_cast<FieldBase*>(f.get())) fb->bindBuffers(forces[fb->level_index()]);
}

static inline int level_ratio(const GridHierarchy& gh, int Lc) 
{
    if (Lc + 1 >= (int)gh.levels.size()) return 0; 
    const auto& C = gh.levels[Lc]; const auto& F = gh.levels[Lc + 1]; double rx = C.dx / F.dx, ry = C.dy / F.dy; int r = (int)std::lround(rx); 
    if (r <= 0 || std::fabs(rx - r) > 1e-9 || std::fabs(ry - r) > 1e-9) return 0; return r; 
}
static inline std::vector<double>& dist_of(IField* f) 
{
    if (auto* p = dynamic_cast<FlowD2Q9_BGK*>(f)) return p->f; 
    if (auto* p = dynamic_cast<FlowD2Q9_MRT*>(f)) return p->f;
    if (auto* p = dynamic_cast<ScalarD2Q9_BGK*>(f)) return p->g;
    if (auto* p = dynamic_cast<ScalarD2Q9_MRT*>(f)) return p->g;
    throw std::runtime_error("dist_of: unknown field type");
}

void SimulationEngine::advanceLevel(int L) 
{
    if ((int)forces.size() != (int)gh.levels.size()) buildForceBuffers();
    forces[L].clear();
    ApplyCtx ctx;
    ctx.t = (L < (int)level_time_phys.size() ? level_time_phys[L] : 0.0);
    ctx.dt = (L < (int)level_Ct.size() ? level_Ct[L] : 1.0);
    ctx.level = L; ctx.gh = &gh; 
    ctx.fieldByName = [&](const std::string& n) {return this->fieldByName(n); };
    ctx.forceOfLevel = [&](int LL)->ForceBuffer& { return this->forces.at(LL); };
    for (auto& op : operators) if (op->targetLevel() == L) op->apply(ctx);
    for (auto& f : fields) if (f->kind() == "flow" && f->level_index() == L) f->compute_macro();
    for (auto& f : fields) 
    {
        if (f->kind() != "scalar") continue; 
        auto* s = dynamic_cast<ScalarD2Q9_BGK*>(f.get());
        if (!s || s->level_index() != L) continue; if (s->advect_by.empty()) continue;
        auto* ff = fieldByName(s->advect_by);
        auto* fb = dynamic_cast<FieldBase*>(ff);
        if (fb && fb->level_index() == L) { s->ux = fb->ux; s->uy = fb->uy;
        }
    }
    for (auto& f : fields) if (f->level_index() == L) f->collide(1.0);
    for (auto& f : fields) if (f->level_index() == L) f->stream();
    int r = level_ratio(gh, L); 
    if (r > 0) 
    {
        int sub = r * r;
        for (int s = 0; s < sub; ++s) advanceLevel(L + 1); 
        compress_fine_to_coarse(L + 1, L); expand_coarse_to_fine(L, L + 1);
    }
    for (auto& f : fields) if (f->level_index() == L) f->compute_macro();
    if (L >= (int)level_time_phys.size()) level_time_phys.resize(L + 1, 0.0); level_time_phys[L] += ctx.dt;
}

void SimulationEngine::compress_fine_to_coarse(int Lf, int Lc) {
    for (const auto& ov : gh.overlaps) 
    {
        if (ov.coarseLevel != Lc || ov.fineLevel != Lf) continue;
        for (auto& fco : fields) 
        { 
            auto* Fco = dynamic_cast<FieldBase*>(fco.get()); 
            if (!Fco || Fco->level_index() != Lc) continue; 
            IField* ff = nullptr; for (auto& fx : fields)
            {
                if (fx->name() == fco->name() && fx->level_index() == Lf) 
                { 
                    ff = fx.get(); break;
                } 
            } 
            if (!ff) continue;
            auto* Ffi = dynamic_cast<FieldBase*>(ff);
            if (!Ffi) continue; 
            const int Ec = Fco->E; std::vector<std::array<double, 6>> ms(Ec); 
            std::vector<int> cnt(Ec, 0); auto& distF = dist_of(ff); 
            for (size_t k = 0; k < ov.fineIds.size(); ++k) 
            {
                int fe = ov.fineIds[k]; int ce = ov.coarseIds[k];
                if (ce < 0 || ce >= Ec || fe < 0) continue; 
                double f9[9]; int off = fe * 9;
                for (int i = 0; i < 9; ++i) f9[i] = distF[off + i];
                double m6[6]; project6_from_f(f9, m6);
                for (int j = 0; j < 6; ++j) ms[ce][j] += m6[j];
                cnt[ce]++;
            }
            auto& distC = dist_of(fco.get());
            for (int ce = 0; ce < Ec; ++ce) 
            {
                if (cnt[ce] == 0) continue;
                double m6[6]; 
                for (int j = 0; j < 6; ++j) m6[j] = ms[ce][j] / cnt[ce];
                double f9[9]; reconstruct_f_from6(m6, f9); int off = ce * 9;
                for (int i = 0; i < 9; ++i) distC[off + i] = f9[i];
            }
        } 
    }
}

void SimulationEngine::expand_coarse_to_fine(int Lc, int Lf) {
    const auto& C = gh.levels[Lc]; const auto& F = gh.levels[Lf];
    for (const auto& ov : gh.overlaps) 
    {
        if (ov.coarseLevel != Lc || ov.fineLevel != Lf) continue; 
        for (auto& fco : fields) 
        {
            auto* Fco = dynamic_cast<FieldBase*>(fco.get()); 
            if (!Fco || Fco->level_index() != Lc) continue;
            IField* ff = nullptr; 
            for (auto& fx : fields) 
            { 
                if (fx->name() == fco->name() && fx->level_index() == Lf) 
                { ff = fx.get(); break; 
                }
            }
            if (!ff) continue; 
            auto* Ffi = dynamic_cast<FieldBase*>(ff);
            if (!Ffi) continue;
            const int Ec = Fco->E;
            std::vector<std::array<double, 6>> mc(Ec);
            {
                auto& distC = dist_of(fco.get());
                for (int ce = 0; ce < Ec; ++ce)
                { 
                    double f9[9]; int off = ce * 9; 
                    for (int i = 0; i < 9; ++i) f9[i] = distC[off + i]; 
                    double m6[6]; project6_from_f(f9, m6);
                    for (int j = 0; j < 6; ++j) mc[ce][j] = m6[j]; 
                } 
            }
            auto& distF = dist_of(ff); int NxC = C.NxC, NyC = C.NyC;
            for (size_t k = 0; k < ov.fineIds.size(); ++k) 
            {
                int fe = ov.fineIds[k];
                if (fe < 0) continue;
                Vec2 cf = F.cellCenter(fe);
                double xi = (cf.x - C.x0c) / C.dx, yi = (cf.y - C.y0c) / C.dy;
                int i0 = (int)std::floor(xi), j0 = (int)std::floor(yi); 
                double tx = xi - i0, ty = yi - j0;
                if (i0 < 0) { i0 = 0; tx = 0; } if (j0 < 0) { j0 = 0; ty = 0; }
                if (i0 >= NxC - 1) { i0 = NxC - 2; tx = 1; } if (j0 >= NyC - 1) { j0 = NyC - 2; ty = 1; }
                int e00 = j0 * NxC + i0, e10 = e00 + 1, e01 = (j0 + 1) * NxC + i0, e11 = e01 + 1; double m6[6]; 
                for (int j = 0; j < 6; ++j) 
                {
                    double m00 = mc[e00][j], m10 = mc[e10][j], m01 = mc[e01][j], m11 = mc[e11][j];
                    m6[j] = (1 - tx) * (1 - ty) * m00 + tx * (1 - ty) * m10 + (1 - tx) * ty * m01 + tx * ty * m11; 
                }
                double f9[9];  reconstruct_f_from6(m6, f9); 
                int off = fe * 9;
                for (int i = 0; i < 9; ++i) distF[off + i] = f9[i]; 
            } 
        }
    }
}

// ========== 单位系统初始化 ==========
// ========== Unit system initialization ==========
void SimulationEngine::init_units_from_config(const Json& cfg) {
    level_Ct.assign(gh.levels.size(), 1.0);
    if (level_time_phys.size() != gh.levels.size()) level_time_phys.assign(gh.levels.size(), 0.0);
    if (!cfg.contains("units")) return; const auto& U = cfg.at("units");
    struct Ulev { double L0 = 1.0, T0 = 0.0, Theta0 = 1.0, T_ref = 0.0; }; std::vector<Ulev> units(gh.levels.size());
    if (U.contains("levels")) {
        const auto& LV = U.at("levels"); 
        for (size_t L = 0; L < gh.levels.size() && L < LV.size(); ++L) 
        {
            const auto& a = LV[L];
            if (a.contains("CL")) units[L].L0 = a["CL"].get<double>(); 
            if (a.contains("T0")) units[L].T0 = a["T0"].get<double>();
            if (a.contains("Theta0")) units[L].Theta0 = a["Theta0"].get<double>();
            if (a.contains("T_ref")) units[L].T_ref = a["T_ref"].get<double>(); 
        }
    }
    bool have_flow = U.contains("flow") && U["flow"].contains("tau") && U["flow"].contains("nu_phys");
    double tau_flow = have_flow ? U["flow"]["tau"].get<double>() : 0.8;
    double nu_phys = have_flow ? U["flow"]["nu_phys"].get<double>() : 1.0e-6;
    // 
    if (have_flow) { for (size_t L = 0; L < gh.levels.size(); ++L) { double L0 = units[L].L0; double T0 = (L0 * L0) * (tau_flow - 0.5) * d2q9::cs2 / std::max(1e-30, nu_phys); units[L].T0 = T0; level_Ct[L] = T0; } }
    else { // 
        double T0_0 = (units[0].T0 > 0 ? units[0].T0 : 1.0);
        for (size_t L = 0; L < gh.levels.size(); ++L) 
        {
            if (L == 0) { units[L].T0 = T0_0; level_Ct[L] = T0_0;
            }
            else 
            {
                double r = gh.levels[L - 1].dx / gh.levels[L].dx; int ri = (int)std::lround(r);
                double T0 = T0_0 / std::max(1, ri * ri); units[L].T0 = T0; level_Ct[L] = T0; 
            }
        }
    }
   //
    using namespace units_runtime; L0s.resize(gh.levels.size()); T0s.resize(gh.levels.size()); Theta0s.resize(gh.levels.size()); Trefs.resize(gh.levels.size());
    for (size_t L = 0; L < gh.levels.size(); ++L) { L0s[L] = units[L].L0; T0s[L] = units[L].T0; Theta0s[L] = units[L].Theta0; Trefs[L] = units[L].T_ref; }
    // 
    if (U.contains("inputs")) {
        std::string f = U["inputs"].value("force", "phys"); std::string s = U["inputs"].value("source", "phys"); std::string b = U["inputs"].value("bc", "phys");
        units_runtime::prefs.force_is_phys = (f == "phys");
        units_runtime::prefs.source_is_phys = (s == "phys");
        units_runtime::prefs.bc_is_phys = (b == "phys");
    }
}
