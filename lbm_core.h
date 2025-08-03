// lbm_core.h (v1.0)
// -----------------------------------------------------------------------------
// A lightweight 2D D2Q9 LBM mini-framework unifying multiple fields (flow/scalar),
// multi-scale support, and time-varying operators.
// - Header file: interface declarations + minimal implementations
//   (TimeExpr inline implementation contained here).
// - Companion files: lbm_core.cpp / main.cpp
//
// Design goals:
//   1) Zero-modification core for adding new fields: registration via FieldFactory / OperatorFactory.
//   2) Operator parameters can vary over time: TimeExpr (constant / table / expression).
//   3) Unified 2D D2Q9 (both flow and scalar support BGK/MRT).
//   4) Output in Tecplot .dat format (variables enumerated automatically by active fields).
//      Complex numerical kernels reside in lbm_core.cpp.
// -----------------------------------------------------------------------------
// lbm_core.h (v1.0)
// -----------------------------------------------------------------------------
// 统一多场（flow / scalar）+ 多尺度 + 可时变算子的 LBM 小框架（2D, D2Q9）
// - 头文件：接口声明 + 轻量实现（TimeExpr 内联实现在本头内）
// - 配套文件：lbm_core.cpp / main.cpp
//
// 设计目标：
//   1) 新增场“零改动”核心：通过 FieldFactory / OperatorFactory 注册。
//   2) 任意算子参数可随时间：TimeExpr（const / table / expr）。
//   3) 2D 统一 D2Q9（flow / scalar 均可 BGK/MRT）。
//   4) 输出保持 Tecplot .dat（按激活场自动选择输出变量）。
//      复杂数值核在 lbm_core.cpp。
// -----------------------------------------------------------------------------


#pragma once
#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <unordered_map>
#include <cassert>
#include <cmath>


#if __has_include(<nlohmann/json_fwd.hpp>)
#include <nlohmann/json_fwd.hpp>
#elif __has_include(<nlohmann/json.hpp>)
#include <nlohmann/json.hpp>
#else
#error "Cannot find nlohmann/json headers. Install via NuGet/vcpkg or use single_include."
#endif
using Json = nlohmann::json;

// ========================= 基础类型与工具 =========================
// ========================= Basic types and tools ==========================
struct Vec2 { double x{ 0 }, y{ 0 }; };

struct GridLevel {
    // 来自 mesh_meta_cell.json 的几何定义（cell-center 坐标系即可推算节点）
    // Geometry definition from mesh_meta_cell.json (nodes can be inferred from the cell-center coordinate system)
    int NxC{ 0 }, NyC{ 0 };        // 单元数（cells）
    int NxN{ 0 }, NyN{ 0 };        // 节点数（nodes）= cells + 1
    double dx{ 1.0 }, dy{ 1.0 };   // 格子宽度   lattice spacing
    double x0c{ 0.5 }, y0c{ 0.5 }; // cell-center 原点(第一个 cell 中心坐标)/cell-center origin (the coordinates of the first cell center)
    double x0n{ 0.0 }, y0n{ 0.0 }; // node origin

    // Convenience functions
    inline int cells() const { return NxC * NyC; }
    inline int nodes() const { return NxN * NyN; }

    // 计算第 e 个单元中心坐标（0-based, 行优先）
    //Calculate the coordinates of the center of the e-th cell (0-based, row-first)
    inline Vec2 cellCenter(int e) const {
        int i = e % NxC; int j = e / NxC;
        return { x0c + i * dx, y0c + j * dy };
    }
};

struct Overlap { 
    // 粗细层映射（coarse ↔ fine）
    //Coarse - fine layer mapping(coarse ↔ fine)
    int coarseLevel{ 0 }, fineLevel{ 1 }, ratio{ 2 };
    std::vector<int> coarseIds;
    std::vector<int> fineIds;
};

struct GridHierarchy {
    std::vector<GridLevel> levels;
    std::vector<Overlap> overlaps; // 多个粗细对 // Multiple thickness pairs
};

// ========================= TimeExpr：时变表达式 =========================
// 统一时间函数接口：任意参数都可以是 ITimeExpr（常数/表格/表达式）
// ======================== = TimeExpr: Time - varying expression ==========================
// Unified time function interface: any parameter can be an ITimeExpr (constant/table/expression)
struct ITimeExpr {
    virtual double eval(double t) const = 0;
    virtual ~ITimeExpr() = default;
};

// const
struct ConstExpr : ITimeExpr {
    double v{ 0 };
    explicit ConstExpr(double val) : v(val) {}
    double eval(double) const override { return v; }
};

// 表格插值（线性/阶梯）//Table interpolation (linear/step)
struct TableExpr : ITimeExpr {
    std::vector<std::pair<double, double>> tbl; // (t, value)
    bool linear{ true }; // true=linear, false=step
    double eval(double t) const override {
        if (tbl.empty()) return 0.0;
        if (t <= tbl.front().first) return tbl.front().second;
        if (t >= tbl.back().first)  return tbl.back().second;
        // 寻找 t 所在分段
        // Find the segment where t is located
        int lo = 0, hi = (int)tbl.size() - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) >> 1;
            if (tbl[mid].first <= t) lo = mid; else hi = mid;
        }
        if (!linear) return tbl[lo].second; // 阶梯//step
        // 线性插值
        // Linear interpolation
        double t0 = tbl[lo].first, v0 = tbl[lo].second;
        double t1 = tbl[hi].first, v1 = tbl[hi].second;
        double a = (t - t0) / (t1 - t0);
        return v0 * (1 - a) + v1 * a;
    }
};

// 简易函数表达式：为保持独立性与可替换性，这里仅实现极轻量版解析：
//  - 支持 + - * / ^ ()
//  - 支持内置函数：sin cos exp log abs sqrt min max
//  - 支持常量字典 consts{"A":...,"t0":...}
//  - 支持符号 t
// 若需更强表达式，可在后续替换为 MuParser / TinyExpr，只需保持接口不变。
// Simple function expressions: to maintain modularity and replaceability, this is a minimal parser:
//  - supports operators +, -, *, /, ^, and parentheses ( )
//  - supports built-in functions: sin, cos, exp, log, abs, sqrt, min, max
//  - supports a constants dictionary consts{"A":..., "t0":...}
//  - supports the variable t
// For more advanced parsing, you can later swap in MuParser or TinyExpr, provided the interface remains the same.

struct FuncExpr : ITimeExpr {
    std::string expr;                     // 原始表达式字符串// Original expression string
    std::unordered_map<std::string, double> consts; 

    // 递归下降解析（极简、足够本项目时序用途；错误输入不会抛异常，仅尽力解析）
    // Recursive descent parsing (minimal, sufficient for this project's timing needs; does not throw on invalid input, simply attempts best-effort parsing)
    mutable const char* p{ nullptr };

    double eval(double t) const override {
        p = expr.c_str();
        return parseExpr(t);
    }

    static bool isAlpha(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_'; }
    static bool isNum(char c) { return (c >= '0' && c <= '9') || c == '.'; }
    void eatWS() const { while (*p == ' ' || *p == '\t') ++p; }

    double parseExpr(double t) const { //    + -
        double v = parseTerm(t); eatWS();
        while (*p == '+' || *p == '-') { char op = *p++; double r = parseTerm(t); v = (op == '+') ? v + r : v - r; eatWS(); }
        return v;
    }
    double parseTerm(double t) const { //    * /
        double v = parsePow(t); eatWS();
        while (*p == '*' || *p == '/') { char op = *p++; double r = parsePow(t); v = (op == '*') ? v * r : v / r; eatWS(); }
        return v;
    }
    double parsePow(double t) const { //    ^
        double v = parseFactor(t); eatWS();
        if (*p == '^') { ++p; double r = parsePow(t); v = std::pow(v, r); }
        return v;
    }
    double parseIdentOrNumber(double t) const {
        eatWS();
        // number
        if (isNum(*p) || ((*p == '+' || *p == '-') && isNum(*(p + 1)))) {
            char* end = nullptr; double v = std::strtod(p, &end); p = end; return v;
        }
        // ident / func
        if (isAlpha(*p)) {
            const char* s = p; while (isAlpha(*p) || isNum(*p)) ++p; std::string id(s, p);
            eatWS();
            if (id == "t") return t;
            auto it = consts.find(id); if (it != consts.end()) return it->second;
            // func(
            if (*p == '(') {
                ++p; double a = parseExpr(t); eatWS(); double b = 0; if (*p == ',') { ++p; b = parseExpr(t); }
                eatWS(); if (*p == ')') ++p;
                if (id == "sin") return std::sin(a);
                if (id == "cos") return std::cos(a);
                if (id == "exp") return std::exp(a);
                if (id == "log") return std::log(std::max(1e-300, a));
                if (id == "abs") return std::fabs(a);
                if (id == "sqrt")return std::sqrt(std::max(0.0, a));
                if (id == "min") return std::min(a, b);
                if (id == "max") return std::max(a, b);
            }
            // 未知标识 -> 0
            //unkown
            return 0.0;
        }
        // 其他 -> 0
        //other
        return 0.0;
    }
    double parseFactor(double t) const {
        eatWS();
        if (*p == '(') { ++p; double v = parseExpr(t); eatWS(); if (*p == ')') ++p; return v; }
        if (*p == '+') { ++p; return parseFactor(t); }
        if (*p == '-') { ++p; return -parseFactor(t); }
        return parseIdentOrNumber(t);
    }
};

// 工厂方法：从 JSON 构造 ITimeExpr
// Factory method: construct an ITimeExpr from JSON
std::unique_ptr<ITimeExpr> make_timeexpr_from_json(const Json& j);

// ========================= Field 抽象与注册 =========================
// ========================= Field Abstraction & Registration =========================
struct IField {
    virtual ~IField() = default;
    virtual std::string name() const = 0;   // fieldname "flow","T"）
    virtual std::string kind() const = 0;   // "flow" | "scalar"
    virtual bool output() const = 0;        // output door

    // 数值步骤（在 LU 中推进；引擎会把 TimeExpr 的物理量换到 LU）
    // Numerical stepping (performed in lattice units; the engine converts TimeExpr’s physical quantities into lattice units)
    virtual void collide(double dt) = 0;    // （Guo/MRT）
    virtual void stream() = 0;              // pull streaming
    virtual void compute_macro() = 0;     

    // 访问底层缓冲（用于接口 E/C 的矩空间投影/重构）
    // 统一假设 D2Q9：每单元 9 个分布，dist().size() == cells()*9
    // Access the underlying buffer (used for moment-space projection/reconstruction in the E/C interfaces)
    // Assumes D2Q9 universally: 9 distributions per cell, so dist().size() == cells()*9
    virtual std::vector<double>& dist() = 0;
    virtual const std::vector<double>& dist() const = 0;


    virtual int level_index() const = 0;
    virtual int cells() const = 0;
    virtual int q() const { return 9; }
};
using FieldPtr = std::unique_ptr<IField>;

struct FieldFactory {
    using Creator = std::function<FieldPtr(const Json&, const GridHierarchy&)>;
    static void reg(const std::string& key, Creator c);
    static FieldPtr create(const Json& j, const GridHierarchy& gh);
};

// ========================= buffer =========================
struct ForceBuffer { // buffer
    std::vector<double> fx, fy;   // force
    std::vector<double> src;      // source
    void resize(int cells) { fx.assign(cells, 0.0); fy.assign(cells, 0.0); src.assign(cells, 0.0); }
    void clear() { std::fill(fx.begin(), fx.end(), 0.0); std::fill(fy.begin(), fy.end(), 0.0); std::fill(src.begin(), src.end(), 0.0); }
};

// ========================= // Operator Abstraction & Registration =========================
struct ApplyCtx {
    double t{ 0 }, dt{ 1 };
    int level{ 0 };
    const GridHierarchy* gh{ nullptr };
    // 访问场与缓冲的回调由调度器注入（避免对具体实现的依赖）
    // Callbacks for accessing fields and buffers are injected by the scheduler (to avoid dependencies on concrete implementations)
    std::function<IField* (const std::string&)> fieldByName; // 读其他场（如 boussinesq 读 T）//read other field
    std::function<ForceBuffer& (int)> forceOfLevel;           //  force/src 
    
};

struct IOperator {
    virtual ~IOperator() = default;
    virtual int targetLevel() const = 0;
    virtual void apply(ApplyCtx& ctx) = 0; // 
};

struct OperatorFactory {
    using Creator = std::function<std::unique_ptr<IOperator>(const Json&)>;
    static void reg(const std::string& type, Creator c);
    static std::unique_ptr<IOperator> create(const Json& j);
};

// ========================= 调度与输出 =========================
// Scheduling & Output
struct SimulationIO {
    int output_every{ 100 };
    std::string out_basename{ "result" };
};

struct SimulationEngine {
    // ----------------- 网格 / 场 / 算子 / I-O -----------------
    //mesh/field/operator/IO
    GridHierarchy gh;                       // 网格层级（从 mesh_meta 读取）//// Grid hierarchy (loaded from mesh_meta)
    std::vector<FieldPtr> fields;           // all field
    std::vector<std::unique_ptr<IOperator>> operators; //
    std::vector<ForceBuffer> forces;       
    SimulationIO io;

    // ----------------- 单位系统（SI↔LU，按需激活） -----------------
    // ----------------- Unit system (SI ↔ LU, activated as needed) -----------------

    struct BaseUnits { //
        double L0{ 1.0 }; // m
        double M0{ 1.0 }; // kg（
        double T0{ 1.0 }; // s
        double I0{ 1.0 }; // A
        double Theta0{ 1.0 }; // K
        double N0{ 1.0 }; // mol
        double J0{ 1.0 }; // cd）
        double T_ref{ 0.0 }; // 
        bool   use_density_phys{ false }; // 
    };
    std::vector<BaseUnits> units;           // 

    // 输入单位偏好：TimeExpr 默认按何种单位提供（可在 config 里覆盖）
    // Input unit preference: default units for TimeExpr (overridable in config)

    struct InputUnitPrefs { bool force_is_phys{ true }; bool source_is_phys{ true }; bool bc_is_phys{ true }; } input_units;

 
    std::vector<double> level_time_phys;    // 
    std::vector<double> level_Ct;           // 


    IField* fieldByName(const std::string& nm) const;


    void buildForceBuffers();

    // —— 单位系统：初始化与常用换算 ——
    // 从 config 中读取 units 段并结合网格 spacing 计算每层 L0/T0：
    // 支持：units.flow.{tau,nu_phys} 用于定标 T0；units.levels[].CL 写入 L0；
    //       units.levels[].Theta0/T_ref；inputs{force,source,bc} 决定算子默认输入单位。
    // —— Unit system: initialization and common conversions ——
// Read the “units” section from the config and, combined with grid spacing, compute L0 (length scale) and T0 (time scale) for each level:
// Supports units.flow.{tau, nu_phys} to calibrate T0; units.levels[].CL to set L0;
//          units.levels[].Theta0/T_ref; inputs.{force, source, bc} determine the default units for operators’ inputs
    void init_units_from_config(const Json& cfg);

    // PU→LU ：
    inline double vel_phys2lu(int L, double u_si)   const { return u_si * (units[L].T0 / units[L].L0); }
    inline double acc_phys2lu(int L, double a_si)   const { return a_si * (units[L].T0 * units[L].T0 / units[L].L0); }
    inline double scalar_val_phys2lu(int L, double T_si) const { return (T_si - units[L].T_ref) / units[L].Theta0; }
    inline double scalar_src_phys2lu(int L, double dTdt_si) const { return (dTdt_si / units[L].Theta0) * units[L].T0; }

    // 旧版：全局同步一步（仍保留以兼容旧调用），内部调用 advanceLevel(0)
   // Compatibility reserved
    void stepOnce(double t_phys, double /*dt_unused*/) { (void)t_phys; advanceLevel(0); }

   //new
    void advanceLevel(int L);

    // 桥接：矩空间 压缩/扩展（只在时间对齐点调用）
//Bridge: Moment Space Compression / Expansion(only called at time alignment points)
    void compress_fine_to_coarse(int Lf, int Lc); // fine→c
    void expand_coarse_to_fine(int Lc, int Lf);   // c→f

    // Tecplot .dat
    void writeTecplot(int step, double t_phys) const;
};

// ========================= 6 矩投影/重构：接口声明 =========================
// 统一 D2Q9 的 6 个“水动力矩”集合：
//   m = [phi_or_rho, jx, jy, e, pxx, pxy]^T
// 其中对 flow：phi_or_rho = rho；对 scalar：phi_or_rho = phi。
// 注意：具体的 P 与 P^+ 常量矩阵在 lbm_core.cpp 中定义。
// ========================= 6-Moment Projection/Reconstruction: Interface Declaration ========================
// Unify the six "hydrodynamic moments" set of D2Q9:
// m = [phi_or_rho, jx, jy, e, pxx, pxy]^T
// Where for flow: phi_or_rho = rho; for scalar: phi_or_rho = phi.
// Note: The specific P and P^+ constant matrices are defined in lbm_core.cpp.
void project6_from_f(const double f[9], double m[6]);
void reconstruct_f_from6(const double m[6], double f[9]);

// ========================= 工具：注册 Key 约定 =========================
// // ========================= Tool: Registration Key Convention =========================
// FieldFactory key: kind+"."+lattice+"."+method，如："flow.D2Q9.BGK", "scalar.D2Q9.MRT"
// OperatorFactory key: type，如："velocity_bc","noslip_bc","dirichlet_bc","body_force",
//                         "boussinesq","heat_pulse","method_bridge"

// ========================= 内部：注册函数（在 lbm_core.cpp 定义） =========================
// ========================= Internal: registration function (defined in lbm_core.cpp) ==========================
void register_builtin_fields();
void register_builtin_operators();

