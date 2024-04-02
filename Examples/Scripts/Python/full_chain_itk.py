#!/usr/bin/env python3
import pathlib
import os
import acts
import acts.examples
import acts.examples.itk
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
    addSeedFilterML,
    SeedFilterMLDBScanConfig,
)

parser = argparse.ArgumentParser(description="Full chain with the ITk detector")

parser.add_argument("--events", "-n", help="Number of events", type=int, default=100)
parser.add_argument(
    "--geo_dir",
    help="Path to the ITk geometry",
    type=str,
    default="/homeijclab/chakkappai/Acts/acts-itk",
)
parser.add_argument(
    "--ambi_config",
    help="Path to the ambiguity resolution config",
    type=str,
    default="/ACTS_itk/ambiguity_resolution_config.json",
)
parser.add_argument(
    "--out_dir",
    help="Path to the output directory",
    type=str,
    default=pathlib.Path.cwd() / "itk_output",
)

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
    "--seedFilter_ML",
    help="Use the Ml Seed Filter instead of the classical one",
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
seedFilter_ML = args["seedFilter_ML"]
geo_dir = pathlib.Path(args["geo_dir"])
ambi_config = pathlib.Path(args["ambi_config"])
ambi_config = str(ambi_config)
ttbar_pu200 = args["ttbar"]

u = acts.UnitConstants
outputDir = pathlib.Path.cwd() / "itk_output"
# acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls

detector, trackingGeometry, decorators = acts.examples.itk.buildITkGeometry(geo_dir)
field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(
    events=args["events"],
    numThreads=1,
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
    preSelectParticles=(
        ParticleSelectorConfig(
            rho=(0.0 * u.mm, 28.0 * u.mm),
            absZ=(0.0 * u.mm, 1.0 * u.m),
            eta=(-4.0, 4.0),
            pt=(150 * u.MeV, None),
            removeNeutral=True,
        )
        if ttbar_pu200
        else ParticleSelectorConfig()
    ),
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
    (
        TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-4.0, 4.0), nHits=(9, None))
        if ttbar_pu200
        else TruthSeedRanges()
    ),
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
if seedFilter_ML:
    addSeedFilterML(
        s,
        SeedFilterMLDBScanConfig(
            epsilonDBScan=0.03, minPointsDBScan=2, minSeedScore=0.1
        ),
        onnxModelFile=os.path.dirname(__file__)
        + "/MLAmbiguityResolution/seedDuplicateClassifier.onnx",
        outputDirRoot=outputDir,
        # outputDirCsv=outputDir,
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
        AthenaAmbiguityResolutionConfig(
            minScore=100,
            minScoreSharedTracks=400,
            maxShared=2,
            maxSharedTracksPerMeasurement=15,
            pTMax=1400,
            pTMin=0.5,
            phiMax=3.14,
            phiMin=-3.14,
            etaMax=2.7,
            etaMin=-2.7,
            useAmbigFcn=False,
        ),
        outputDirRoot=outputDir,
        AmbiVolumeFile=ambi_config,
        writeCovMat=True,
        # outputDirCsv=outputDir,
    )

else:
    addAthenaAmbiguityResolution(
        s,
        AthenaAmbiguityResolutionConfig(
            minScore=0,
            minScoreSharedTracks=1,
            maxShared=2,
            maxSharedTracksPerMeasurement=2,
            pTMax=1400,
            pTMin=0.5,
            phiMax=3.14,
            phiMin=-3.14,
            etaMax=4,
            etaMin=-4,
            useAmbigFcn=True,
        ),
        outputDirRoot=outputDir,
        AmbiVolumeFile=ambi_config,
        writeCovMat=True,
        # outputDirCsv=outputDir,
    )


# addVertexFitting(
#     s,
#     field,
#     vertexFinder=VertexFinder.Iterative,
#     outputDirRoot=outputDir,
# )

s.run()
