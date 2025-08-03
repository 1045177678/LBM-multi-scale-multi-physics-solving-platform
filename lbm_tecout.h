// lbm_tecout.h
// 独立的 Tecplot 导出（避免受 lbm_core.cpp 行数限制影响）
// 用法：在 main.cpp 中 #include "lbm_tecout.h"，然后调用 WriteTecplot(eng, step, t_phys)
// lbm_tecout.h
// Independent Tecplot export (to avoid line-count limitations in lbm_core.cpp)
// Usage: #include "lbm_tecout.h" in main.cpp, then call WriteTecplot(eng, step, t_phys)
#pragma once
#include <string>

struct SimulationEngine; // 前置声明（由 lbm_core.h 定义）// Forward declarations (defined in lbm_core.h)


// 写出 Tecplot FEQUADRILATERAL（BLOCK）文件，变量：X,Y,(rho,ux,uy), 所有标量场
// 文件名：eng.io.out_basename + "_" + step + ".dat"
// Write Tecplot FEQUADRILATERAL (BLOCK) file with variables: X, Y, (rho, ux, uy), and all scalar fields
// Filename: eng.io.out_basename + "_" + step + ".dat"
void WriteTecplot(const SimulationEngine& eng, int step, double t_phys);
