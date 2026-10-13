// 拡散項はソルバー (Y-hEqn_DD.H / Y-hEqn.H) の左辺から
// ddt 項と対流項を除いたものと同じ離散式で評価する

const volScalarField mut(turbulence->mut());

// Soret 項を持つ種 (H, H2) の判定と熱拡散比
auto hasSoret = [&](const label i)
{
    return thermalDiffusion && (Y[i].name() == "H" || Y[i].name() == "H2");
};
auto TDRatio = [&](const label i) -> const volScalarField&
{
    return Y[i].name() == "H" ? TDRatio_H : TDRatio_H2;
};

// --------------------------------------------------
// 補正速度フラックス phiUc = (Σ ρ Vk) · Sf [kg/s]
// Y-hEqn_DD.H と同じく全化学種 (不活性種を含む) で和をとる
// differentialDiffusion off のときはゼロ
// --------------------------------------------------
DiffError *= 0;
if (differentialDiffusion)
{
    forAll(Y, i)
    {
        const volScalarField& Yi = Y[i];
        const volScalarField& Diffi = Diff[i];

        DiffError += rho*Diffi*fvc::grad(Yi) + Yi*rho*Diffi*fvc::grad(MW)/MW;
        if (hasSoret(i))
        {
            DiffError +=
                -rho*Diffi*TDRatio(i)*fvc::grad(T)*Yi/(T*(X[i] + SMALL));
        }
    }
}
const surfaceScalarField phiUc("phiUc", fvc::interpolate(DiffError) & mesh.Sf());

const dimensionedScalar zeroYTerm("zero", dimMass/dimVolume/dimTime, 0.0);

forAll(Y, i)
{
    if (composition.active(i))
    {
        volScalarField& Yi = Y[i];
        const word& name = Yi.name();

        // --------------------------------------------------
        // 対流項: ∇·(ρu Yi) [kg/m³/s]
        // --------------------------------------------------
        volScalarField conv_Yi("conv_" + name, mvConvection->fvcDiv(phi, Yi));

        // --------------------------------------------------
        // 拡散項: ∇·ji [kg/m³/s]（左辺に置かれた項なので符号はそのまま）
        //   diffFick  : -∇·((ρDi + μt/Sct)∇Yi)   (off: -∇·((α + μt/Sct)∇Yi))
        //   diffCorr  :  ∇·(phiUc Yi)             補正速度
        //   diffMW    : -∇·(ρDi Yi/MW ∇MW)
        //   diffSoret : -∇·(-ρDi Θi Yi/(T Xi) ∇T) H, H2 のみ
        // fvm 項は (mat & Yi) = (A*Yi - b)/V で陽的評価する
        // --------------------------------------------------
        volScalarField diffFick_Yi
        (
            IOobject("diffFick_" + name, runTime.timeName(), mesh),
            mesh,
            zeroYTerm
        );
        volScalarField diffCorr_Yi
        (
            IOobject("diffCorr_" + name, runTime.timeName(), mesh),
            mesh,
            zeroYTerm
        );
        volScalarField diffMW_Yi
        (
            IOobject("diffMW_" + name, runTime.timeName(), mesh),
            mesh,
            zeroYTerm
        );
        volScalarField diffSoret_Yi
        (
            IOobject("diffSoret_" + name, runTime.timeName(), mesh),
            mesh,
            zeroYTerm
        );

        if (differentialDiffusion)
        {
            const volScalarField& Diffi = Diff[i];

            diffFick_Yi = (-fvm::laplacian(rho*Diffi + mut/Sct, Yi)) & Yi;
            diffCorr_Yi = mvConvection->fvmDiv(phiUc, Yi) & Yi;
            diffMW_Yi = -fvc::laplacian(rho*Diffi*Yi/MW, MW);

            if (hasSoret(i))
            {
                diffSoret_Yi = -fvc::laplacian
                (
                    -rho*Diffi*TDRatio(i)*Yi/(T*(X[i] + SMALL)),
                    T
                );
            }
        }
        else
        {
            diffFick_Yi = (-fvm::laplacian(mixAlpha + mut/Sct, Yi)) & Yi;
        }

        volScalarField diff_Yi
        (
            "diff_" + name,
            diffFick_Yi + diffCorr_Yi + diffMW_Yi + diffSoret_Yi
        );

        // --------------------------------------------------
        // non-chemical = -conv - diff
        // --------------------------------------------------
        volScalarField nonChem_Yi("nonChem_" + name, - conv_Yi - diff_Yi);

        conv_Yi.write();
        diff_Yi.write();
        nonChem_Yi.write();

        if (writeDiffTerms)
        {
            diffFick_Yi.write();
            diffCorr_Yi.write();
            diffMW_Yi.write();
            diffSoret_Yi.write();
        }

        Info<< name << nl
            << "  nonChem  [kg/m3/s] min/max : "
            << gMin(nonChem_Yi)  << " / " << gMax(nonChem_Yi)  << nl
            << endl;

        Info << "conv min/max " << gMin(conv_Yi) << " / " << gMax(conv_Yi) << endl;
        Info << "diff min/max " << gMin(diff_Yi) << " / " << gMax(diff_Yi) << endl;
    }
}
