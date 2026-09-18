#include "fvCFD.H"
#include "fvOptions.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "psiReactionThermophysicalTransportModel.H"

int main(int argc,char *argv[])
{
    argList::validArgs.append("time-folder");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (args.args().size() != 2)
    {
        FatalErrorInFunction
            << "Usage: " << args.executable() << " <time-folder>"
            << exit(FatalError);
    }

    const word timeFolder(args.args()[1]);

    instantList times = runTime.times();
    label timeI = -1;
    forAll(times,i)
    {
        if (times[i].name() == timeFolder)
        {
            timeI = i;
            break;
        }
    }

    if (timeI == -1)
    {
        FatalErrorInFunction
            << "Time folder " << timeFolder << " was not found."
            << exit(FatalError);
    }

    runTime.setTime(times[timeI], timeI);

    #include "createFields.H"
    #include "createFieldRefs.H"

    Info<< "cemaSourceTerms start" << endl;

    // instantList times = runTime.times();
    // instantList times = timeSelector::select0(runTime,args);
    
    // Info << "T min/max " << gMin(T)  << " " << gMax(T) << endl;
    Info << "T min/max " << gMin(thermo.T())  << " " << gMax(thermo.T()) << endl;
    // Info << "p min/max " << gMin(p)  << " " << gMax(p) << endl;
    Info << "p min/max " << gMin(thermo.p())  << " " << gMax(thermo.p()) << endl;
    Info << "rho min/max " << gMin(rho)  << " " << gMax(rho) << endl;
    // Info << "U min/max " << gMin(U)  << " " << gMax(U) << endl;
    Info << "phi min/max " << gMin(phi)  << " " << gMax(phi) << endl;

    tmp<volScalarField> tcp = thermo.Cp();
    Info << "Cp min/max " << gMin(tcp()) << " " << gMax(tcp()) << endl;
    Info << "he min/max " << gMin(thermo.he()) << " " << gMax(thermo.he()) << endl;

    
    forAll(Y,i)
    {
        Info << Y[i].name() << " min/max " << gMin(Y[i])  << " " << gMax(Y[i]) << endl;
    }

    // ========================================================
    // CEMA: YEqn の chemical / non-chemical 項の計算・出力
    // ========================================================
    // --- 対流スキームの生成 ---
    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );

    thermo.correct();

    runTime.setDeltaT(1e-12);
    reaction->correct();

    forAll(Y, i)
    {
        if (composition.active(i))
        {
            volScalarField& Yi = Y[i];
            const word& name = Yi.name();

            // --------------------------------------------------
            // chemical項: ω̇i [kg/m³/s]
            // R(Yi) は右辺項 → -(mat & Yi)/V が陽的評価
            // --------------------------------------------------

            // NOTE:
            // この source() ベースの抽出は combustionModel = laminar / PaSR / EDC 系
            // （fvm::Sp() を使わず Su += chemistryPtr_->RR(i) のみで構成される実装）
            // でのみ厳密に正しい。
            // 将来 singleStepCombustion 系（infinitelyFastChemistry, diffusion, FSD 等）
            // や semiImplicit 処理を使うモデルに切り替える場合は、
            // fvm::Sp() による陰的項が加わるため本コードの前提が崩れる。
            // その際は (tRi() & Yi) ベースの実装、または diag()/hasDiag() を
            // チェックした上での再検証が必要。
            tmp<fvScalarMatrix> tRi = reaction->R(Yi);

            volScalarField chem_Yi
            (
                IOobject
                (
                    "chem_" + name,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            );
            const scalarField& source = tRi().source();
            const scalarField& Vcells = mesh.V();
            forAll(chem_Yi, cellI)
            {
                chem_Yi[cellI] = -source[cellI] / Vcells[cellI];
            }
            tRi.clear();

            volScalarField Qdot(reaction->Qdot());
            Info << "Qdot min/max " << gMin(Qdot) << " / " << gMax(Qdot) << endl;

            // --------------------------------------------------
            // 対流項: ∇·(ρu Yi) [kg/m³/s]
            // fvc::div で陽的に計算
            // --------------------------------------------------
            volScalarField conv_Yi
            (
                IOobject
                (
                    "conv_" + name,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            );
            conv_Yi = mvConvection->fvcDiv(phi, Yi);

            // --------------------------------------------------
            // 拡散項: ∇·ji [kg/m³/s]
            // divj は fvm::laplacian 相当の純線形行列
            // source() はゼロ → (mat & Yi) で陽的評価する
            //
            // (mat & Yi) = A*Yi - H(Yi) = 残差ベクトル [kg/s]
            // 左辺に置かれた項なので符号はそのまま: +∇·ji
            // --------------------------------------------------
            tmp<fvScalarMatrix> tdiffMat = thermophysicalTransport->divj(Yi);

            volScalarField diff_Yi
            (
                IOobject
                (
                    "diff_" + name,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            );
            {
                diff_Yi = tdiffMat() & Yi;
            }
            tdiffMat.clear();

            // --------------------------------------------------
            // non-chemical = -conv - diff
            // --------------------------------------------------
            volScalarField nonChem_Yi
            (
                IOobject
                (
                    "nonChem_" + name,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                - conv_Yi - diff_Yi
                // - conv_Yi // for check conv_Yi
            );

            chem_Yi.write();
            diff_Yi.write();
            nonChem_Yi.write();


            // Info<< name
            //     << "  chem   [kg/m3/s] min/max: "
            //     << gMin(chem_Yi)    << " / " << gMax(chem_Yi)    << nl
            //     << "  nonChem[kg/m3/s] min/max: "
            //     << gMin(nonChem_Yi) << " / " << gMax(nonChem_Yi) << nl
            //     << endl;

            // --------------------------------------------------
            // 残差: chem - nonChem
            // ddt = 0になって収束していればゼロになるはず
            // --------------------------------------------------
            volScalarField residual_Yi
            (
                IOobject
                (
                    "residual_" + name,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                chem_Yi - nonChem_Yi
            );

            // 残差の統計（ゼロに近いほど整合している）
            const scalar resMax  = gMax(mag(residual_Yi)());
            const scalar resMean = residual_Yi.weightedAverage(mesh.V()).value();
            const scalar chemMax = gMax(mag(chem_Yi)());

            Info<< name << nl
                // << "  ddt      [kg/m3/s] min/max : "
                // << gMin(ddt_Yi)      << " / " << gMax(ddt_Yi)      << nl
                << "  chem     [kg/m3/s] min/max : "
                << gMin(chem_Yi)     << " / " << gMax(chem_Yi)     << nl
                << "  nonChem  [kg/m3/s] min/max : "
                << gMin(nonChem_Yi)  << " / " << gMax(nonChem_Yi)  << nl
                << "  residual [kg/m3/s] max|mean|: "
                << resMax << " / " << resMean << nl
                << "  relative residual (resMax/chemMax): "
                << (chemMax > SMALL ? resMax/chemMax : 0.0) << nl
                << endl;

            
            Info << "conv min/max " << gMin(conv_Yi) << " / " << gMax(conv_Yi) << endl;
            Info << "diff min/max " << gMin(diff_Yi) << " / " << gMax(diff_Yi) << endl;

            // Info<< "Qdot = " << reaction->Qdot()()[0] << endl;   // 発熱速度との整合性チェック

        }
    }

    // ========================================================
    // CEMA: EEqn の chemical / non-chemical 項の計算・出力
    // ========================================================
    volScalarField& he = thermo.he();

    volScalarField chem_h
    (
        IOobject
        (
            "chem_h",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        reaction->Qdot()
    );

    volScalarField conv_h
    (
        IOobject
        (
            "conv_h",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", chem_h.dimensions(), 0.0)
    );
    conv_h = mvConvection->fvcDiv(phi, he);

    tmp<fvScalarMatrix> tdiffMatE = thermophysicalTransport->divq(he);
    volScalarField diff_h
    (
        IOobject
        (
            "diff_h",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", chem_h.dimensions(), 0.0)
    );
    diff_h = tdiffMatE() & he;
    tdiffMatE.clear();

    volScalarField nonChem_h
    (
        IOobject
        (
            "nonChem_h",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -conv_h - diff_h
    );

    volScalarField residual_h
    (
        IOobject
        (
            "residual_h",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        chem_h - nonChem_h
    );

    chem_h.write();
    conv_h.write();
    diff_h.write();
    nonChem_h.write();
    residual_h.write();
    he.write();

    Info<< "energy terms" << nl
        << "  chem_h [kg/m/s^3] min/max : "
        << gMin(chem_h) << " / " << gMax(chem_h) << nl
        << "  nonChem_h [kg/m/s^3] min/max : "
        << gMin(nonChem_h) << " / " << gMax(nonChem_h) << nl
        << "  residual_h [kg/m/s^3] max|mean| : "
        << gMax(mag(residual_h)()) << " / "
        << residual_h.weightedAverage(mesh.V()).value() << nl
        << endl;

    return 0;
}
