volScalarField& he = thermo.he();

volScalarField conv_enthalpy("conv_enthalpy", mvConvection->fvcDiv(phi, he));

// --------------------------------------------------
// 拡散項: ∇·q [W/m³]（左辺に置かれた項なので符号はそのまま）
//   diffCond  : -∇·((α + μt/Prt)∇h)
//   diffsh    :  Σ ∇·(α hs,i ∇Yi)          differentialDiffusion on のみ
//   diffFlux  : -∇·(Σ hs,i ρ Vi Yi)         differentialDiffusion on のみ
// --------------------------------------------------
volScalarField diffCond_enthalpy
(
    "diffCond_enthalpy",
    (-fvm::laplacian(mixAlpha + mut/Prt, he)) & he
);

volScalarField diffFlux_enthalpy
(
    IOobject("diffFlux_enthalpy", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
);

vk *= 0;
diffsh *= 0;

if (differentialDiffusion)
{
    // 化学種 i の顕エンタルピー hs,i(p, T) を Hsi に設定
    auto setHsi = [&](const label i)
    {
        forAll(Hsi, celli)
        {
            Hsi[celli] = composition.Hs(i, p[celli], T[celli]);
        }

        volScalarField::Boundary& HsiBf = Hsi.boundaryFieldRef();

        forAll(HsiBf, patchi)
        {
            scalarField& HsiPatch = HsiBf[patchi];
            const scalarField& pp = p.boundaryField()[patchi];
            const scalarField& Tp = T.boundaryField()[patchi];

            forAll(HsiPatch, facei)
            {
                HsiPatch[facei] = composition.Hs(i, pp[facei], Tp[facei]);
            }
        }
    };

    surfaceScalarField sumOfOtherSpecies("sumOfOtherSpecies", 0*phiUc);
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex)
        {
            const volScalarField& Yi = Y[i];
            const volScalarField& Diffi = Diff[i];

            setHsi(i);

            // 化学種拡散フラックス (補正速度込み) [kg/s]
            volVectorField speciesFluxVol
            (
                rho*Diffi*fvc::grad(Yi)
              + Yi*rho*Diffi*fvc::grad(MW)/MW
              - Yi*DiffError
            );
            if (hasSoret(i))
            {
                speciesFluxVol +=
                    -rho*Diffi*TDRatio(i)*fvc::grad(T)*Yi/(T*(X[i] + SMALL));
            }
            const surfaceScalarField speciesFlux
            (
                fvc::interpolate(speciesFluxVol) & mesh.Sf()
            );

            vk += fvc::interpolate(Hsi)*speciesFlux;
            sumOfOtherSpecies += speciesFlux;
            diffsh += fvc::laplacian(mixAlpha*Hsi, Yi);
            Yt += Yi;
        }
    }

    // 不活性種: ソルバーと同じくフラックスは -Σ(他種), 質量分率は 1 - Yt
    setHsi(inertIndex);
    vk += fvc::interpolate(Hsi)*-sumOfOtherSpecies;
    diffsh += fvc::laplacian(mixAlpha*Hsi, 1 - Yt);

    diffFlux_enthalpy = -fvc::div(vk);
}

volScalarField diff_enthalpy
(
    "diff_enthalpy",
    diffCond_enthalpy + diffsh + diffFlux_enthalpy
);

// --------------------------------------------------
// 圧力項 [W/m³]: EEqn 左辺の -dpdt (e のときは +div(phiv,p)) を右辺へ移したもの
// pEqn.H と同じく保存済みの p, p_0 から fvc::ddt(p) で評価する
// (ソルバーではこの値が次のステップの EEqn で使われる)
// --------------------------------------------------
volScalarField dpdt_enthalpy("dpdt_enthalpy", 0*conv_enthalpy);
if (he.name() == "e")
{
    dpdt_enthalpy = -fvc::div
    (
        fvc::absolute(phi/fvc::interpolate(rho), U),
        p,
        "div(phiv,p)"
    );
}
else if (thermo.dpdt())
{
    dpdt_enthalpy = fvc::ddt(p);
}

// --------------------------------------------------
// 運動エネルギー項 [W/m³]: ddt(ρK) + ∇·(φK)（左辺に置かれた項）
// K = 0.5|U|² は解像された流れの運動エネルギー (laminar でもゼロでない)
// ρ の旧時刻値は保存されていないので、連続の式 ∂ρ/∂t = -∇·φ を使い
//   ddt(ρK) = ρ ddt(K) + K_0 ∂ρ/∂t = ρ ddt(K) - K_0 ∇·φ
// として評価する（K_0 は保存済みの U_0 から計算）
// --------------------------------------------------
K.oldTime() = 0.5*magSqr(U.oldTime());

volScalarField kin_enthalpy
(
    "kin_enthalpy",
    rho*fvc::ddt(K) - K.oldTime()*fvc::div(phi) + fvc::div(phi, K)
);

volScalarField nonChem_enthalpy
(
    "nonChem_enthalpy",
    -conv_enthalpy - diff_enthalpy + dpdt_enthalpy - kin_enthalpy
);

// ソルバー出力 mixAlpha との比較用
volScalarField mixAlpha_post("mixAlpha_post", mixAlpha);
mixAlpha_post.write();

rho.write();
conv_enthalpy.write();
diff_enthalpy.write();
dpdt_enthalpy.write();
kin_enthalpy.write();
nonChem_enthalpy.write();

if (writeDiffTerms)
{
    diffCond_enthalpy.write();
    volScalarField("diffsh_enthalpy", diffsh).write();
    diffFlux_enthalpy.write();
}

Info<< "energy terms" << nl
    << "  nonChem_enthalpy [kg/m/s^3] min/max : "
    << gMin(nonChem_enthalpy) << " / " << gMax(nonChem_enthalpy) << nl
    << endl;
