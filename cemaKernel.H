#ifndef cemaKernel_H
#define cemaKernel_H

// CEMA の固有値解析 (Eigen を使う。fvCFD.H とは別の翻訳単位にする)

#include <vector>

struct cemaResult
{
    double lambdaRe;   // CEM の固有値の実部 [1/s]
    double lambdaIm;   // CEM の固有値の虚部 [1/s]
    double phiOmega;   // b_e . omega
    double phiS;       // b_e . s
    double alpha;      // phiS/phiOmega (propagation index)
};

// J: 列優先 n x n のヤコビアン (J[i + j*n] = d f_i / d y_j)
// omega, s: 長さ n の chem / non-chem 項
// nConserved: 保存モードとして除外する |lambda| 最小の固有値の数
cemaResult cemaAnalyse
(
    const double* J,
    const double* omega,
    const double* s,
    const int n,
    const int nConserved,
    std::vector<double>* eigRe = nullptr,
    std::vector<double>* eigIm = nullptr
);

#endif
