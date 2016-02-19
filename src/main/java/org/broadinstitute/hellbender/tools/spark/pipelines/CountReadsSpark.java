package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.PrintStream;

@CommandLineProgramProperties(summary = "Counts reads in the input SAM/BAM",
        oneLineSummary = "CountReads on Spark",
        programGroup = SparkProgramGroup.class)
public final class CountReadsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();

        final long count = reads.count();
        System.out.println(count);

        try ( final PrintStream ps = new PrintStream(BucketUtils.createFile(out, getAuthenticatedGCSOptions())) ) {
            ps.print(count);
        }
    }
}
