#include "fvCFD.H"
#include "fvOptions.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "psiReactionThermophysicalTransportModel.H"

int main(int argc,char *argv[])
{
    argList::validArgs.append("time-folder");
    
    argList::addBoolOption
    (
        "calculateChem",
        "Calculate chemical source terms"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (args.args().size() != 2)
    {
        FatalErrorInFunction
            << "Usage: " << args.executable()
            << " <time-folder> [-calculateChem]"
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

    forAll(Y, i)
    {
        if (composition.active(i))
        {
            volScalarField& Yi = Y[i];
            const word& name = Yi.name();

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
            );

            conv_Yi.write();
            diff_Yi.write();
            nonChem_Yi.write();

            Info<< name << nl
                << "  nonChem  [kg/m3/s] min/max : "
                << gMin(nonChem_Yi)  << " / " << gMax(nonChem_Yi)  << nl
                << endl;

            Info << "conv min/max " << gMin(conv_Yi) << " / " << gMax(conv_Yi) << endl;
            Info << "diff min/max " << gMin(diff_Yi) << " / " << gMax(diff_Yi) << endl;
        }
    }

    // #include <iomanip>   // setprecision 用（必要なら）

    // // --- 化学種ごとの Hf [J/kg] (Tstd 基準) と分子量 [kg/kmol] ---
    // {
    //     OFstream hfFile(runTime.path()/runTime.timeName()/"Hf_species.dat");
    //     hfFile << "# species  Hf[J/kg]  W[kg/kmol]" << nl;
    //     hfFile.precision(16);

    //     forAll(Y, i)
    //     {
    //         const scalar Hfi = composition.Hf(i);   // [J/kg]
    //         const scalar Wi  = composition.Wi(i);   // [kg/kmol]

    //         hfFile << Y[i].name() << " " << Hfi << " " << Wi << nl;
    //         Info<< Y[i].name() << " Hf = " << Hfi << " J/kg, W = " << Wi << endl;
    //     }
    // }

    // ========================================================
    // CEMA: EEqn の chemical / non-chemical 項の計算・出力
    // ========================================================
    volScalarField& he = thermo.he();

    volScalarField conv_enthalpy
    (
        IOobject
        (
            "conv_enthalpy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        // dimensionedScalar("zero", chem_enthalpy.dimensions(), 0.0)
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    );
    conv_enthalpy = mvConvection->fvcDiv(phi, he);

    tmp<fvScalarMatrix> tdiffMatE = thermophysicalTransport->divq(he);
    volScalarField diff_enthalpy
    (
        IOobject
        (
            "diff_enthalpy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        // dimensionedScalar("zero", chem_enthalpy.dimensions(), 0.0)
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    );
    diff_enthalpy = tdiffMatE() & he;
    tdiffMatE.clear();

    volScalarField nonChem_enthalpy
    (
        IOobject
        (
            "nonChem_enthalpy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -conv_enthalpy - diff_enthalpy
    );
    
    volScalarField alphaEff(thermophysicalTransport->alphaEff());
    alphaEff.write();
    
    // volScalarField divPhi("divPhi", fvc::div(phi));
    // divPhi.write();

    rho.write();
    // he.write();
    // conv_enthalpy.write();
    // diff_enthalpy.write();
    nonChem_enthalpy.write();
    
    Info<< "energy terms" << nl
        << "  nonChem_enthalpy [kg/m/s^3] min/max : "
        << gMin(nonChem_enthalpy) << " / " << gMax(nonChem_enthalpy) << nl
        << endl;

    // flag for calculation of chemical source terms
    // default: false
    // const bool flagCalculateChem = args.found("calculateChem");
    const bool flagCalculateChem = args.optionFound("calculateChem");

    if (flagCalculateChem) {

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

                chem_Yi.write();
            }
        }

        volScalarField chem_enthalpy
        (
            IOobject
            (
                "chem_enthalpy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            reaction->Qdot()
        );

        chem_enthalpy.write();
        // residual_enthalpy.write();

        Info<< "energy terms" << nl
            << "  chem_enthalpy [kg/m/s^3] min/max : "
            << gMin(chem_enthalpy) << " / " << gMax(chem_enthalpy) << nl
            << endl;
    }

    return 0;
}
