package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link GermlineCNVCaller}.
 */
public final class GermlineCNVCallerIntegrationTest extends CommandLineProgramTest {
    private static final String gCNVSimDataDir = toolsTestDir + "copynumber/gcnv-sim-data/";
    private static final File[] testCountFiles = IntStream.range(0, 20)
            .mapToObj(n -> new File(gCNVSimDataDir + String.format("SAMPLE_%03d_counts.tsv", n)))
            .toArray(File[]::new);
    private static final File contigPloidyCallsOutputDir = new File(gCNVSimDataDir + "contig-ploidy-calls/");
    private static final File simIntervalListSubsetFile = new File(gCNVSimDataDir + "sim_intervals_subset.interval_list");
    private final File tempOutputDir = createTempDir("test-germline-cnv");

    /**
     * Run the tool in the COHORT mode for all 20 samples on a small subset of intervals
     */
    @Test(groups = {"python"})
    public void testCohortWithoutIntervalAnnotations() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(testCountFiles).forEach(argsBuilder::addInput);
        argsBuilder.addArgument(GermlineCNVCaller.RUN_MODE_LONG_NAME, GermlineCNVCaller.RunMode.COHORT.name())
                .addArgument("L", simIntervalListSubsetFile.getAbsolutePath())
                .addArgument(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME,
                        contigPloidyCallsOutputDir.getAbsolutePath())
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempOutputDir.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-germline-cnv-cohort")
                .addArgument(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .addArgument(StandardArgumentDefinitions.VERBOSITY_NAME, "DEBUG");
        runCommandLine(argsBuilder);
    }

    /**
     * Run the tool in CASE mode for the first 5 samples using the model generated by
     * {@link #testCohortWithoutIntervalAnnotations()}
     */
    @Test(groups = {"python"}, dependsOnMethods = "testCohortWithoutIntervalAnnotations")
    public void testCase() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(testCountFiles, 0, 5).forEach(argsBuilder::addInput);
        argsBuilder.addArgument(GermlineCNVCaller.RUN_MODE_LONG_NAME, GermlineCNVCaller.RunMode.CASE.name())
                .addArgument(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME,
                        contigPloidyCallsOutputDir.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.MODEL_LONG_NAME,
                        new File(tempOutputDir, "test-germline-cnv-cohort-model").getAbsolutePath())
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempOutputDir.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-germline-cnv-case")
                .addArgument(StandardArgumentDefinitions.VERBOSITY_NAME, "DEBUG");
        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"}, expectedExceptions = IllegalArgumentException.class)
    public void testCaseWithoutModel() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(testCountFiles, 0, 5).forEach(argsBuilder::addInput);
        argsBuilder.addArgument(GermlineCNVCaller.RUN_MODE_LONG_NAME, GermlineCNVCaller.RunMode.CASE.name())
                .addArgument(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME,
                        contigPloidyCallsOutputDir.getAbsolutePath())
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempOutputDir.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-germline-cnv-case")
                .addArgument(StandardArgumentDefinitions.VERBOSITY_NAME, "DEBUG");
        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"}, enabled = false)
    public void testCohortWithInputModel() {
    }

    @Test(groups = {"python"}, enabled = false)
    public void testCohortWithAnnotatedIntervals() {
    }
}
