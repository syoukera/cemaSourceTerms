#include "cemaKernel.H"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

cemaResult cemaAnalyse
(
    const double* J,
    const double* omega,
    const double* s,
    const int n,
    const int nConserved,
    std::vector<double>* eigRe,
    std::vector<double>* eigIm
)
{
    // J^T の右固有ベクトル = J の左固有ベクトル (b^T J = lambda b^T)
    const Eigen::Map<const Eigen::MatrixXd> Jmat(J, n, n);
    const Eigen::EigenSolver<Eigen::MatrixXd> es(Jmat.transpose(), true);

    const Eigen::VectorXcd lambda = es.eigenvalues();

    if (eigRe && eigIm)
    {
        eigRe->resize(n);
        eigIm->resize(n);
        for (int i = 0; i < n; i++)
        {
            (*eigRe)[i] = lambda(i).real();
            (*eigIm)[i] = lambda(i).imag();
        }
    }

    // |lambda| の小さい順に nConserved 個を保存モードとして除外
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort
    (
        order.begin(),
        order.end(),
        [&](const int a, const int b)
        {
            return std::abs(lambda(a)) < std::abs(lambda(b));
        }
    );

    // 残りの中で実部最大のものを CEM とする
    int e = -1;
    for (int k = std::min(std::max(nConserved, 0), n - 1); k < n; k++)
    {
        const int i = order[k];
        if (e < 0 || lambda(i).real() > lambda(e).real())
        {
            e = i;
        }
    }

    // 複素共役対の場合は実部を使う (alpha は比なので正規化は不要)
    const Eigen::VectorXd be = es.eigenvectors().col(e).real();

    const Eigen::Map<const Eigen::VectorXd> omegaVec(omega, n);
    const Eigen::Map<const Eigen::VectorXd> sVec(s, n);

    cemaResult r;
    r.lambdaRe = lambda(e).real();
    r.lambdaIm = lambda(e).imag();
    r.phiOmega = be.dot(omegaVec);
    r.phiS = be.dot(sVec);

    const double vSmall = std::numeric_limits<double>::min();
    if (std::abs(r.phiOmega) > vSmall)
    {
        r.alpha = r.phiS/r.phiOmega;
    }
    else
    {
        // phiOmega ~ 0: 符号付きの大きな値
        const double sgn = (r.phiS >= 0 ? 1.0 : -1.0)*(r.phiOmega >= 0 ? 1.0 : -1.0);
        r.alpha = sgn*1e30;
    }

    return r;
}
