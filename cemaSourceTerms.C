#include "fvCFD.H"
#include "fvOptions.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "psiReactionThermophysicalTransportModel.H"

int main(int argc,char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // runTime.setTime(instant(0.1, "0.1"), 0);

    instantList times = runTime.times();
    label timeI = -1;
    forAll(times,i)
    {
        if (mag(times[i].value() - 0.05) < SMALL)
        {
            timeI = i;
            break;
        }
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
            {
                const scalarField& Vcells = mesh.V();
                volScalarField chemResidField = tRi() & Yi;  // tmp<volScalarField>を実体化
                const scalarField& chemResid = chemResidField.primitiveField();
                forAll(chem_Yi, cellI)
                {
                    chem_Yi[cellI] = -chemResid[cellI] / Vcells[cellI];
                }
            }
            tRi.clear();

            // // --------------------------------------------------
            // // chemical項: ω̇i [kg/m³/s]
            // // --------------------------------------------------
            // tmp<fvScalarMatrix> tRi = reaction->R(Yi);
            // const scalarField& srcRi = tRi->source();  // [kg/s] per cell
            
            // // 対角成分の寄与が無視できるか確認
            // Info << "R diag min/max: " 
            //     << gMin(tRi->diag()) << " / " << gMax(tRi->diag()) << endl;


            // volScalarField chem_Yi
            // (
            //     IOobject
            //     (
            //         "chem_" + name,
            //         runTime.timeName(),
            //         mesh,
            //         IOobject::NO_READ,
            //         IOobject::AUTO_WRITE
            //     ),
            //     mesh,
            //     dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            // );
            // // セル体積で割って [kg/m³/s] に変換
            // forAll(chem_Yi, cellI)
            // {
            //     chem_Yi[cellI] = srcRi[cellI] / mesh.V()[cellI];
            // }
            // tRi.clear();

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
                mvConvection->fvcDiv(phi, Yi)
            );

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
                const scalarField& Vcells = mesh.V();
                volScalarField diffResidField = tdiffMat() & Yi;  // tmp<volScalarField>を実体化
                const scalarField& diffResid = diffResidField.primitiveField();
                forAll(diff_Yi, cellI)
                {
                    diff_Yi[cellI] = diffResid[cellI] / Vcells[cellI];
                }
            }
            tdiffMat.clear();

            // --------------------------------------------------
            // 拡散項: ∇·ji [kg/m³/s]
            // divj が返す fvScalarMatrix から source() で取り出す
            // --------------------------------------------------
            // tmp<fvScalarMatrix> tdiffMat = thermophysicalTransport->divj(Yi);
            // const scalarField& srcDiff = tdiffMat->source();  // [kg/s] per cell

            // Info << "diff diag min/max: " << gMin(tdiffMat->diag()) << " / "
            //     << gMax(tdiffMat->diag()) << endl;


            // volScalarField diff_Yi
            // (
            //     IOobject
            //     (   
            //         "diff_" + name,
            //         runTime.timeName(),
            //         mesh,
            //         IOobject::NO_READ,
            //         IOobject::AUTO_WRITE
            //     ),
            //     mesh,
            //     dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            // );
            // forAll(diff_Yi, cellI)
            // {
            //     diff_Yi[cellI] = srcDiff[cellI] / mesh.V()[cellI];
            // }
            // tdiffMat.clear();

            // // --------------------------------------------------
            // // ddt
            // // --------------------------------------------------
            // volScalarField ddtYi
            // (
            //     IOobject
            //     (
            //         "ddt_"+name,
            //         runTime.timeName(),
            //         mesh
            //     ),
            //     fvc::ddt(rho,Yi)
            // );

            // Info
            //     << "ddt min/max "
            //     << gMin(ddtYi)
            //     << " / "
            //     << gMax(ddtYi)
            //     << endl;

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

            chem_Yi.write();
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

        }
    }

    return 0;
}
