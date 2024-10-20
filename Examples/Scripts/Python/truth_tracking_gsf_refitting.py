#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants


def runRefittingGsf(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    s: acts.examples.Sequencer = None,
):
    s = runTruthTrackingKalman(
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        outputDir=outputDir,
        s=s,
    )

    gsfOptions = {
        "betheHeitlerApprox": acts.examples.AtlasBetheHeitlerApprox.makeDefault(),
        "maxComponents": 12,
        "componentMergeMethod": acts.examples.ComponentMergeMethod.maxWeight,
        "mixtureReductionAlgorithm": acts.examples.MixtureReductionAlgorithm.KLDistance,
        "weightCutoff": 1.0e-4,
        "level": acts.logging.INFO,
    }

    s.addAlgorithm(
        acts.examples.RefittingAlgorithm(
            acts.logging.INFO,
            inputTracks="kf_tracks",
            outputTracks="gsf_refit_tracks",
            fit=acts.examples.makeGsfFitterFunction(
                trackingGeometry, field, **gsfOptions
            ),
        )
    )

    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=acts.logging.INFO,
            inputTracks="gsf_refit_tracks",
            inputParticles="truth_seeds_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="refit_track_particle_matching",
            outputParticleTrackMatching="refit_particle_track_matching",
        )
    )

    s.addWriter(
        acts.examples.RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="gsf_refit_tracks",
            inputParticles="truth_seeds_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_gsf_refit.root"),
        )
    )

    s.addWriter(
        acts.examples.RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="truth_seeds_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            filePath=str(outputDir / "tracksummary_gsf_refit.root"),
        )
    )

    s.addWriter(
        acts.examples.TrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="truth_seeds_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_gsf_refit.root"),
        )
    )

    return s


if __name__ == "__main__":
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    # ODD
    from acts.examples.odd import getOpenDataDetector

    detector, trackingGeometry, decorators = getOpenDataDetector()
    digiConfigFile = (
        srcdir / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json"
    )

    ## GenericDetector
    # detector, trackingGeometry, _ = acts.examples.GenericDetector.create()
    # digiConfigFile = (
    #     srcdir
    #     / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
    # )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runRefittingGsf(
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=digiConfigFile,
        outputDir=Path.cwd(),
    ).run()
