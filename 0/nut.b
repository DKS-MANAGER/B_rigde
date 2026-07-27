/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Web:      www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0.0;

boundaryField
{
    inlet
    {
        type            calculated;
        value           uniform 0.0;
    }

    outlet
    {
        type            calculated;
        value           uniform 0.0;
    }

    bed
    {
        type            nutkRoughWallFunction;
        Ks              uniform 0.00068; // Grade II sand
        Cs              uniform 0.5;
        value           uniform 0.0;
    }

    top
    {
        type            slip;
    }

    bridge
    {
        type            nutkWallFunction;
        value           uniform 0.0;
    }

    frontAndBack
    {
        type            empty;
    }
}

// ************************************************************************* //
