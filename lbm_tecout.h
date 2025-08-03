// lbm_tecout.h
// ������ Tecplot ������������ lbm_core.cpp ��������Ӱ�죩
// �÷����� main.cpp �� #include "lbm_tecout.h"��Ȼ����� WriteTecplot(eng, step, t_phys)
// lbm_tecout.h
// Independent Tecplot export (to avoid line-count limitations in lbm_core.cpp)
// Usage: #include "lbm_tecout.h" in main.cpp, then call WriteTecplot(eng, step, t_phys)
#pragma once
#include <string>

struct SimulationEngine; // ǰ���������� lbm_core.h ���壩// Forward declarations (defined in lbm_core.h)


// д�� Tecplot FEQUADRILATERAL��BLOCK���ļ���������X,Y,(rho,ux,uy), ���б�����
// �ļ�����eng.io.out_basename + "_" + step + ".dat"
// Write Tecplot FEQUADRILATERAL (BLOCK) file with variables: X, Y, (rho, ux, uy), and all scalar fields
// Filename: eng.io.out_basename + "_" + step + ".dat"
void WriteTecplot(const SimulationEngine& eng, int step, double t_phys);
