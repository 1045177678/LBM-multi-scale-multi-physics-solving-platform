
#include "lbm_core.h"
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>

// D2Q9 
namespace d2q9_local {
    static const int    ex[9] = {0,1,0,-1,0, 1,-1,-1, 1};
    static const int    ey[9] = {0,0,1, 0,-1,1, 1,-1,-1};
    static const double w [9] = {4.0/9.0,
                                 1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,
                                 1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
}

static inline void flow_macro_from_f9(const double f9[9], double& rho, double& ux, double& uy){
    using namespace d2q9_local;
    double s0=0, sx=0, sy=0;
    for(int i=0;i<9;++i){ s0 += f9[i]; sx += f9[i]*ex[i]; sy += f9[i]*ey[i]; }
    rho = s0;
    if (rho>1e-30){ ux = sx / rho; uy = sy / rho; }
    else { ux=0; uy=0; }
}

static inline double scalar_from_g9(const double g9[9]){
    double s=0; for(int i=0;i<9;++i) s+=g9[i]; return s;
}

void WriteTecplot(const SimulationEngine& eng, int step, double /*t_phys*/)
{
    //
    bool anyFlow=false; std::vector<std::string> scalarNames;
    for (const auto& fp : eng.fields) {
        if (!fp->output()) continue;
        if (fp->kind()=="flow") anyFlow=true;
        else if (fp->kind()=="scalar") scalarNames.push_back(fp->name());
    }
    std::sort(scalarNames.begin(), scalarNames.end());
    scalarNames.erase(std::unique(scalarNames.begin(), scalarNames.end()), scalarNames.end());

    // 
    std::string fn = eng.io.out_basename + "_" + std::to_string(step) + ".dat";
    std::ofstream ofs(fn);
    if(!ofs){ std::cerr<<"[io] cannot open "<<fn<<""; return; }

    ofs << "TITLE = \"Cell-Center Mesh\"";
    ofs << "VARIABLES = \"X\",\"Y\"";
    if(anyFlow) ofs << ",\"rho\",\"ux\",\"uy\"";
    for (const auto& s : scalarNames) ofs << ",\"" << s << "\"";
    ofs << "";

    //by  level
    for (size_t L=0; L<eng.gh.levels.size(); ++L) {
        const auto& GL = eng.gh.levels[L];
        int NxN=GL.NxN, NyN=GL.NyN, NxC=GL.NxC, NyC=GL.NyC;
        int N = NxN*NyN, E = NxC*NyC;

        ofs << "ZONE T=\"Level" << L << "\", N="<<N<<", E="<<E
            << ", DATAPACKING=BLOCK";
        int cellVarCount = (anyFlow?3:0) + (int)scalarNames.size();
        if (cellVarCount>0) ofs << ", VARLOCATION=([3-" << (2+cellVarCount) << "]=CELLCENTERED)";
        ofs << ", ZONETYPE=FEQUADRILATERAL";
        // X block
        for(int j=0;j<NyN;++j)
            for(int i=0;i<NxN;++i)
                ofs << (GL.x0n + i*GL.dx) << "";
        // Y block
        for(int j=0;j<NyN;++j)
            for(int i=0;i<NxN;++i)
                ofs << (GL.y0n + j*GL.dy) << "";

        // flow
        const IField* flowF = nullptr;
        for (const auto& fp : eng.fields)
            if (fp->output() && fp->kind()=="flow" && fp->level_index()==(int)L) { flowF = fp.get(); break; }

        if (anyFlow) {
            std::vector<double> rho(E,0.0), ux(E,0.0), uy(E,0.0);
            if (flowF) {
                const auto& f = flowF->dist();
                for (int e=0; e<E; ++e){
                    double f9[9]; int off = e*9;
                    for(int i=0;i<9;++i) f9[i]=f[off+i];
                    flow_macro_from_f9(f9, rho[e], ux[e], uy[e]);
                }
            }
            for (int e=0;e<E;++e) ofs << rho[e] << "";
            for (int e=0;e<E;++e) ofs << ux[e]  << "";
            for (int e=0;e<E;++e) ofs << uy[e]  << "";
        }

        // b)  scalar 
        for (const auto& sname : scalarNames) {
            const IField* scF = nullptr;
            for (const auto& fp : eng.fields)
                if (fp->output() && fp->kind()=="scalar"
                    && fp->name()==sname && fp->level_index()==(int)L) { scF = fp.get(); break; }

            std::vector<double> phi(E,0.0);
            if (scF){
                const auto& g = scF->dist();
                for (int e=0; e<E; ++e){
                    double g9[9]; int off=e*9;
                    for(int i=0;i<9;++i) g9[i]=g[off+i];
                    phi[e] = scalar_from_g9(g9);
                }
            }
            for (int e=0;e<E;++e) ofs << phi[e] << "";
        }

        // FEQUADRILATERAL connection£¨1-based£©
        for (int j=0;j<NyC;++j){
            for (int i=0;i<NxC;++i){
                int n0 = j*NxN + i;
                int n1 = n0 + 1;
                int n3 = (j+1)*NxN + i;
                int n2 = n3 + 1;
                ofs << (n0+1) << " " << (n1+1) << " " << (n2+1) << " " << (n3+1) << "";
            }
        }
    }
}
