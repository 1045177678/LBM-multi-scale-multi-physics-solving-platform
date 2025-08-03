// main.cpp (v1.0)
// -----------------------------------------------------------------------------
// 作用：
//   - 读取 mesh_meta_cell.json（若不存在则构造一个单层默认网格）
//   - 读取 config.json（新 schema v1.0：fields + operators + io）
//   - 构建 SimulationEngine，注册内置 Field/Operator，运行步进并输出 Tecplot
// 依赖：
//   - 单头 JSON 库：nlohmann/json（https://github.com/nlohmann/json）
//     将 json.hpp 放到与本项目同级的 include 目录，或直接拷到工程中。
//   - C++17
// -----------------------------------------------------------------------------
// main.cpp (v1.0)
// -----------------------------------------------------------------------------
// Purpose:
//   - Load mesh_meta_cell.json (or construct a single‐level default mesh if it doesn’t exist)
//   - Load config.json (new schema v1.0: fields + operators + io)
//   - Build the SimulationEngine, register built-in Fields/Operators, execute timestepping, and export Tecplot output
// Dependencies:
//   - Single-header JSON library: nlohmann/json (https://github.com/nlohmann/json)
//     Place json.hpp in your project’s include directory, or copy it directly into the project.
//   - Requires C++17
// -----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include "lbm_core.h"
#include "lbm_tecout.h"

#if __has_include(<nlohmann/json.hpp>)
#include <nlohmann/json.hpp>
#elif __has_include("json.hpp")
#include "json.hpp"
#else
#error "Cannot find nlohmann/json.hpp header."
#endif
using json = nlohmann::json;

// -------------------------- 工具：读取文件到 json --------------------------
// -------------------------- Tool: read file to json --------------------------
static json load_json(const std::string& fn) { std::ifstream ifs(fn); if (!ifs) throw std::runtime_error("Cannot open " + fn); json j; ifs >> j; return j; }

// -------------------------- 从 mesh_meta_cell.json 解析层级 -----------------
// -------------------------- Parse hierarchy from mesh_meta_cell.json -----------------
static GridHierarchy load_mesh_meta_or_default(const std::string& fn) {
    GridHierarchy gh; std::ifstream ifs(fn);
    if (!ifs) {
        // 默认：单层 11x11 cells，dx=1，原点 (0,0)
        // Default: single layer 11x11 cells, dx=1, origin (0,0)
        GridLevel L; L.NxC = 10; L.NyC = 10; L.NxN = 11; L.NyN = 11; L.dx = 1; L.dy = 1; L.x0n = 0; L.y0n = 0; L.x0c = 0.5; L.y0c = 0.5; gh.levels.push_back(L);
        std::cerr << "[warn] mesh_meta_cell.json not found, using default single level 10x10.\n";
        return gh;
    }
    json j; ifs >> j; auto arr = j.at("levels");
    for (size_t i = 0; i < arr.size(); ++i) {
        auto Lj = arr[i]; GridLevel L;
        L.NxN = Lj.at("dims_node")[0]; L.NyN = Lj.at("dims_node")[1];
        L.NxC = Lj.at("dims_cell")[0]; L.NyC = Lj.at("dims_cell")[1];
        L.x0n = Lj.at("origin_node")[0]; L.y0n = Lj.at("origin_node")[1];
        L.x0c = Lj.at("origin_cell")[0]; L.y0c = Lj.at("origin_cell")[1];
        L.dx = Lj.at("spacing")[0];    L.dy = Lj.at("spacing")[1];
        gh.levels.push_back(L);
    }
    // overlaps（option）
    if (j.contains("overlaps")) {
        for (auto& oj : j["overlaps"]) {
            Overlap ov; ov.coarseLevel = oj.at("coarse"); ov.fineLevel = oj.at("fine"); ov.ratio = oj.at("ratio");
            if (oj.contains("mapping")) {
                for (auto& pair : oj["mapping"]) { ov.coarseIds.push_back(pair[0]); ov.fineIds.push_back(pair[1]); }
            }
            gh.overlaps.push_back(std::move(ov));
        }
    }
    return gh;
}

// --------------------------  fields / operators / io -------------------
static void build_engine_from_config(SimulationEngine& eng, const json& cfg) {
    //  register
    register_builtin_fields();
    register_builtin_operators();

    // fields
    for (auto& fj : cfg.at("fields")) {
        auto f = FieldFactory::create(fj, eng.gh);
        eng.fields.push_back(std::move(f));
    }
    // operators
    for (auto& oj : cfg.at("operators")) {
        auto op = OperatorFactory::create(oj);
        eng.operators.push_back(std::move(op));
    }
    // io
    if (cfg.contains("io")) {
        eng.io.output_every = cfg["io"].value("output_every", 100);
        eng.io.out_basename = cfg["io"].value("out_basename", std::string("result"));
    }

    eng.buildForceBuffers();
}


int main(int argc, char** argv) {
    try {
        std::string mesh_meta = "mesh_meta_cell.json"; // 由网格生成器导出
        std::string config_fn = (argc > 1 ? argv[1] : "config.json");

        //
        SimulationEngine eng; eng.gh = load_mesh_meta_or_default(mesh_meta);

        //
        json cfg = load_json(config_fn);
        build_engine_from_config(eng, cfg);
        // ini units:   units.flow.{tau,nu_phys} and levels[].CL to  Ct(T0)
        eng.init_units_from_config(cfg);

        // time step
        int steps = cfg.value("steps", 1000);
        int out_every = eng.io.output_every;
        double t_phys = 0.0;
        for (int n = 0; n < steps; ++n) {
            eng.stepOnce(0.0, 0.0); // 
            t_phys = (eng.level_time_phys.empty() ? 0.0 : eng.level_time_phys[0]);
            if (n % out_every == 0) WriteTecplot(eng, n, t_phys);
        }

        std::cout << "Done. Steps=" << steps
            << " t_phys=" << t_phys
            << ". Output every " << out_every << " steps.\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "[FATAL] " << e.what() << "\n"; return 1;
    }
}
