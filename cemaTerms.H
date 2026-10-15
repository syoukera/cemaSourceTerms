// CEMA: 各セルで化学ヤコビアン J = dω/dy (y = [T, Y_0..Y_{N-2}], pyJac) の
// 固有値解析を行い、CEM の左固有ベクトル b_e へ chem 項 ω と non-chem 項 s を
// 射影して propagation index α = (b_e·s)/(b_e·ω) を求める
//
//   ω: -cemaChem pyjac    : pyJac dydt の瞬時値
//      -cemaChem openfoam : chem_<Yi>/ρ, (chem_enthalpy - Σ hs_i chem_<Yi>)/(ρ cp)
//   s: -cemaFrame lagrangian : 拡散項 (+ dpdt, kin) のみ (ρ Dy/Dt = ... の形)
//      -cemaFrame eulerian   : さらに非保存形の対流項 -(∇·(φY) - Y∇·φ) を加える
//   温度の項は he = hs = Σ Y_i hs_i(T) から ρ cp dT = d(ρhs) - ρ Σ hs_i dY_i

if (he.name() != "h")
{
    FatalErrorInFunction
        << "CEMA requires energy sensibleEnthalpy (he = h), not " << he.name()
        << exit(FatalError);
}

const label nCema = pyjac_nsp();
const label nConserved =
    args.optionLookupOrDefault<label>("cemaNConserved", pyjac_n_conserved());

Info<< "CEMA: " << nCema << " pyJac species, chem = " << cemaChem
    << ", frame = " << cemaFrame << ", nConserved = " << nConserved << endl;

// --- pyJac の種 k → OpenFOAM の Y の index ---
labelList ofIndex(nCema, -1);
for (label k = 0; k < nCema; k++)
{
    const word name(pyjac_species_name(k));
    if (!composition.species().found(name))
    {
        FatalErrorInFunction
            << "pyJac species " << name << " is not in the OpenFOAM species "
            << composition.species() << exit(FatalError);
    }
    ofIndex[k] = composition.species()[name];
}
if (nCema != Y.size())
{
    FatalErrorInFunction
        << "Number of species differs: pyJac " << nCema
        << ", OpenFOAM " << Y.size() << exit(FatalError);
}
if (ofIndex[nCema - 1] != inertIndex)
{
    WarningInFunction
        << "The last pyJac species " << pyjac_species_name(nCema - 1)
        << " (excluded from the state vector) is not the inert species "
        << inertSpecie << endl;
}

forAll(Y, i)
{
    if (!convY.set(i) || !diffY.set(i))
    {
        FatalErrorInFunction
            << "Non-chemical terms of " << Y[i].name() << " are not available"
            << " (inactive species)" << exit(FatalError);
    }
    if (cemaChem == "openfoam" && !chemY.set(i))
    {
        FatalErrorInFunction
            << "chem_" << Y[i].name() << " is not available" << exit(FatalError);
    }
}

// --- 出力フィールド ---
auto cemaField = [&](const word& name, const dimensionSet& dims)
{
    return new volScalarField
    (
        IOobject
        (
            name,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dims, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
};

PtrList<volScalarField> cemaFields(8);
cemaFields.set(0, cemaField("cema_lambda", dimless/dimTime));
cemaFields.set(1, cemaField("cema_lambdaLog", dimless));
cemaFields.set(2, cemaField("cema_phiOmega", dimless));
cemaFields.set(3, cemaField("cema_phiS", dimless));
cemaFields.set(4, cemaField("cema_alpha", dimless));
cemaFields.set(5, cemaField("cema_mode", dimless));
cemaFields.set(6, cemaField("cema_dTdtChem", dimTemperature/dimTime));
cemaFields.set(7, cemaField("cema_dTdtNonChem", dimTemperature/dimTime));

volScalarField& cemaLambda = cemaFields[0];
volScalarField& cemaLambdaLog = cemaFields[1];
volScalarField& cemaPhiOmega = cemaFields[2];
volScalarField& cemaPhiS = cemaFields[3];
volScalarField& cemaAlpha = cemaFields[4];
volScalarField& cemaMode = cemaFields[5];
volScalarField& cemaDTdtChem = cemaFields[6];
volScalarField& cemaDTdtNonChem = cemaFields[7];

// --- 非化学項 s (セル値, 質量あたりにする前) ---
const volScalarField divPhi(fvc::div(phi));
const volScalarField cpField(thermo.Cp());

// ρ D(hs)/Dt の非化学項 [W/m³]
volScalarField sH("sH", -diff_enthalpy + dpdt_enthalpy - kin_enthalpy);
if (cemaFrame == "eulerian")
{
    sH -= conv_enthalpy - he*divPhi;
}

const bool useOpenFOAMChem = (cemaChem == "openfoam");

const label nY = Y.size();
List<scalar> y(nCema), J(nCema*nCema), dy(nCema);
List<scalar> omega(nCema), s(nCema);
List<scalar> hs(nY), sYc(nY), wYc(nY);

// 固有値をログに出す (ローカルの T 最大・最小セル)
const label cellTmax = findMax(T.primitiveField());
const label cellTmin = findMin(T.primitiveField());

forAll(T, celli)
{
    const scalar rhoc = rho[celli];
    const scalar Tc = T[celli];
    const scalar pc = p[celli];
    const scalar cpc = cpField[celli];

    // 状態ベクトル y = [T, Y_0..Y_{N-2}] (pyJac の順序)
    y[0] = Tc;
    for (label k = 0; k < nCema - 1; k++)
    {
        y[k + 1] = Y[ofIndex[k]][celli];
    }

    pyjac_jacob(pc, y.begin(), J.begin());

    forAll(Y, i)
    {
        hs[i] = composition.Hs(i, pc, Tc);
    }

    // --- 非化学項 s: dY_i/dt [1/s] ---
    // 不活性種はソルバーと同じく他種の和の符号反転 (Σ dY/dt = 0)
    scalar sumOther = 0;
    forAll(Y, i)
    {
        if (i == inertIndex) continue;

        scalar ri = -diffY[i][celli];
        if (cemaFrame == "eulerian")
        {
            ri -= convY[i][celli] - Y[i][celli]*divPhi[celli];
        }
        sYc[i] = ri/rhoc;
        sumOther += sYc[i];
    }
    sYc[inertIndex] = -sumOther;

    scalar sumHsS = 0;
    forAll(Y, i)
    {
        sumHsS += hs[i]*sYc[i];
    }

    s[0] = (sH[celli] - rhoc*sumHsS)/(rhoc*cpc);
    for (label k = 0; k < nCema - 1; k++)
    {
        s[k + 1] = sYc[ofIndex[k]];
    }

    // --- 化学項 ω ---
    if (useOpenFOAMChem)
    {
        scalar sumOtherW = 0;
        forAll(Y, i)
        {
            if (i == inertIndex) continue;
            wYc[i] = chemY[i][celli]/rhoc;
            sumOtherW += wYc[i];
        }
        wYc[inertIndex] = -sumOtherW;

        scalar sumHsW = 0;
        forAll(Y, i)
        {
            sumHsW += hs[i]*wYc[i];
        }

        omega[0] =
            (chemEnthalpyPtr()[celli] - rhoc*sumHsW)/(rhoc*cpc);
        for (label k = 0; k < nCema - 1; k++)
        {
            omega[k + 1] = wYc[ofIndex[k]];
        }
    }
    else
    {
        pyjac_dydt(pc, y.begin(), dy.begin());
        omega = dy;
    }

    std::vector<double> eigRe, eigIm;
    const bool printEig = (celli == cellTmax || celli == cellTmin);

    const cemaResult r = cemaAnalyse
    (
        J.begin(),
        omega.begin(),
        s.begin(),
        nCema,
        nConserved,
        printEig ? &eigRe : nullptr,
        printEig ? &eigIm : nullptr
    );

    cemaLambda[celli] = r.lambdaRe;
    cemaLambdaLog[celli] = sign(r.lambdaRe)*Foam::log10(scalar(1) + mag(r.lambdaRe));
    cemaPhiOmega[celli] = r.phiOmega;
    cemaPhiS[celli] = r.phiS;
    cemaAlpha[celli] = r.alpha;
    cemaDTdtChem[celli] = omega[0];
    cemaDTdtNonChem[celli] = s[0];

    // 0: 非爆発 (λ_e < 0), 1: 自着火 (|α| < 1),
    // 2: 拡散支援着火 (α > 1), 3: 消炎 (α < -1)
    if (r.lambdaRe < 0)
    {
        cemaMode[celli] = 0;
    }
    else if (r.alpha > 1)
    {
        cemaMode[celli] = 2;
    }
    else if (r.alpha < -1)
    {
        cemaMode[celli] = 3;
    }
    else
    {
        cemaMode[celli] = 1;
    }

    if (printEig)
    {
        Pout<< "CEMA eigenvalues at cell " << celli << " (T = " << Tc
            << " K), lambda_e = " << r.lambdaRe << " + " << r.lambdaIm
            << "i:" << nl;
        for (label i = 0; i < nCema; i++)
        {
            Pout<< "    " << eigRe[i] << " + " << eigIm[i] << "i" << nl;
        }
        Pout<< endl;
    }
}

forAll(cemaFields, fi)
{
    cemaFields[fi].correctBoundaryConditions();
    cemaFields[fi].write();

    Info<< cemaFields[fi].name() << " min/max : "
        << gMin(cemaFields[fi]) << " / " << gMax(cemaFields[fi]) << endl;
}
Info<< endl;
