#!/usr/bin/env python3
import pathlib,os, acts, acts.examples, acts.examples.itk
import argparse
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    TruthSeedRanges,
    addCKFTracks,
    TrackSelectorConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addAthenaAmbiguityResolution,
    AthenaAmbiguityResolutionConfig,
    addAmbiguityResolutionML,
    AmbiguityResolutionMLConfig,
    addVertexFitting,
    VertexFinder,
)

parser = argparse.ArgumentParser(description="Full chain with the ITk detector")

parser.add_argument("--events", "-n", help="Number of events", type=int, default=100)
parser.add_argument("--geo_dir", help="Path to the ITk geometry", type=str, default="/homeijclab/chakkappai/Acts/acts-itk")

parser.add_argument(
    "--geant4", help="Use Geant4 instead of fatras", action="store_true"
)
parser.add_argument(
    "--ttbar",
    help="Use Pythia8 (ttbar, pile-up 200) instead of particle gun",
    action="store_true",
)
parser.add_argument(
    "--MLSolver",
    help="Use the Ml Ambiguity Solver instead of the classical one",
    action="store_true",
)
parser.add_argument(
    "--AthenaSolver",
    help="Use the Athena Ambiguity Solver instead of the classical one",
    action="store_true",
)
parser.add_argument(
    "--GreedySolver",
    help="Use the Greedy Ambiguity Solvera and then Athena Ambiguity solver",
    action="store_true",
)

args = vars(parser.parse_args())

ambiguity_MLSolver = args["MLSolver"]
athena_ambiguity_resolution = args["AthenaSolver"]
greedy_ambiguity_resolution = args["GreedySolver"]
geo_dir = pathlib.Path(args["geo_dir"])

ttbar_pu200 = False
u = acts.UnitConstants
outputDir = pathlib.Path.cwd() / "itk_output"
# acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls

detector, trackingGeometry, decorators = acts.examples.itk.buildITkGeometry(geo_dir)
field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(
    events=args["events"],
    numThreads=-1,
    outputDir=str(outputDir),
)

if not ttbar_pu200:
    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        EtaConfig(-4.0, 4.0, uniform=True),
        ParticleConfig(2, acts.PdgParticle.eMuon, randomizeCharge=True),
        rnd=rnd,
    )
else:
    addPythia8(
        s,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=200,
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
            mean=acts.Vector4(0, 0, 0, 0),
        ),
        rnd=rnd,
        outputDirRoot=outputDir,
    )

addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    preSelectParticles=ParticleSelectorConfig(
        rho=(0.0 * u.mm, 28.0 * u.mm),
        absZ=(0.0 * u.mm, 1.0 * u.m),
        eta=(-4.0, 4.0),
        pt=(150 * u.MeV, None),
        removeNeutral=True,
    )
    if ttbar_pu200
    else ParticleSelectorConfig(),
    outputDirRoot=outputDir,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir / "itk-hgtd/itk-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
)

addSeeding(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-4.0, 4.0), nHits=(9, None))
    if ttbar_pu200
    else TruthSeedRanges(),
    seedingAlgorithm=SeedingAlgorithm.Default,
    *acts.examples.itk.itkSeedingAlgConfig(
        acts.examples.itk.InputSpacePointsType.PixelSpacePoints
    ),
    initialSigmas=[
        1 * u.mm,
        1 * u.mm,
        1 * u.degree,
        1 * u.degree,
        0.1 / u.GeV,
        1 * u.ns,
    ],
    initialVarInflation=[1.0] * 6,
    geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json",
    outputDirRoot=outputDir,
)

addCKFTracks(
    s,
    trackingGeometry,
    field,
    trackSelectorConfig=(
        # fmt: off
        TrackSelectorConfig(absEta=(None, 2.0), pt=(0.9 * u.GeV, None), nMeasurementsMin=9, maxHoles=2, maxSharedHits=2),
        TrackSelectorConfig(absEta=(None, 2.6), pt=(0.4 * u.GeV, None), nMeasurementsMin=8, maxHoles=2, maxSharedHits=2),
        TrackSelectorConfig(absEta=(None, 4.0), pt=(0.4 * u.GeV, None), nMeasurementsMin=7, maxHoles=2, maxSharedHits=2),
        # fmt: on
    ),
    outputDirRoot=outputDir,
)


if ambiguity_MLSolver:
    addAmbiguityResolutionML(
        s,
        AmbiguityResolutionMLConfig(
            maximumSharedHits=3, maximumIterations=1000000, nMeasurementsMin=7
        ),
        outputDirRoot=outputDir,
        # outputDirCsv=outputDir,
        onnxModelFile=os.path.dirname(__file__)
        + "/MLAmbiguityResolution/duplicateClassifier.onnx",
    )

elif greedy_ambiguity_resolution:
    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(
            maximumSharedHits=3,
            maximumIterations=1000000,
            nMeasurementsMin=7,
        ),
        outputDirRoot=outputDir,
        writeCovMat=True,
        # outputDirCsv=outputDir,
    )

    addAthenaAmbiguityResolution(
        s,
        outputDirRoot=outputDir,
        writeCovMat=True,
        # outputDirCsv=outputDir,
    )

else:
    addAthenaAmbiguityResolution(
        s,
        outputDirRoot=outputDir,
        writeCovMat=True,
    )


# addVertexFitting(
#     s,
#     field,
#     vertexFinder=VertexFinder.Iterative,
#     outputDirRoot=outputDir,
# )

s.run()
