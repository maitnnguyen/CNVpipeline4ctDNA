#!/usr/bin/env anduril
//$OPT --threads 20
//$OPT --pipe "tee $ANDURIL_EXECUTION_DIR/_log"
//$OPT --pipe "anduril-pager --ls"
//$OPT --wrapper slurm-prefix
//$PRE export ANDURIL_NODELIST="evm01 evm04 evm05 evm06 evm07 evm08 evm09 evm10 evmfull01 evmfull02"
//$PRE export ANDURIL_SELECTNODE=plain

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.microarray._
import anduril.sequencing._

// *Note: in coverage file for each sample, chrY should be removed
// data contain input file has at least 6 columns: patient | sample | bamfiles | type | batch | TP53
// vcf file by sample is the output from GATK, after filtering 

object ichorCNA{

    val files = INPUT(path="/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/int/bamlistset7.csv")
    val fltbam = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/targetMar01/result_CNV_ctDNA/fltBamList/out.csv")
    val germlines = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/old_pipe/result_pon/varsPerNormalCSV/out.csv")
    val SNVperSample = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/old_pipe/result_variantCalling/varsBySamplePassCSV/out.csv")
    val intervals = INPUT(path="/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v4/int/ROIS7.interval")
    val interval_GCannot = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v4/int/ROIS7_annot.intervals")
    val ctDNAinfo = INPUT(path = "/mnt/storageBig8/work/joikkone/ctDNA_s8/ctDNA_allBatches_baseciInfo.csv")
    val centromere = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2021/PipelineScript/CNV/oncosysCNA/int/GRCh38_centromere_acen.txt")
    val tp53_info = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/manuscript/int/TP53_vaf_SAGE.csv")

    // step1: calculate raw coverage from bamfiles
    val Coverage = NamedMap[BashEvaluate]("Coverage")
    val CoverageArray = NamedMap[BinaryFile]("CoverageArray")
    // step2: correct GC and mappability
    val CorrectBias = NamedMap[BashEvaluate]("CorrectBias")
    // step3: making PoN
    //val makePoNfromBDNA = NamedMap[BashEvaluate]("makePoNfromBDNA")
    //val makePoNfromZeroTPL = NamedMap[BashEvaluate]("makePoNfromZeroTPL")
    // step4: CNV calling
    val CNVcall = NamedMap[BashEvaluate]("CNVcall")
    val LogRArray = NamedMap[BinaryFile]("LogRArray")
    val SegmArray = NamedMap[BinaryFile]("SegmArray")

    /// Load coverage from PureCN pipeline to this for ichorCNA
    /*for ( row <- iterCSV(files)) {
            //val patient = row("patient")
            val bamfile = INPUT(row("file"))
            //val vcffile = row("vcf")
            val sample = row("sample")

            Coverage(sample) = BashEvaluate(var1 = bamfile,
                                         var2 = intervals,
                                         script = """
                                         Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/S1_coverage_cal.R @var1@ @var2@ @out1@
                                         """
            )
            Coverage(sample)._filename("out1","out1.txt")
            Coverage(sample)._custom("memory") = "2G"
            CoverageArray(sample) = Coverage(sample).out1
            
    } */
    val CoverageList = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v4/process/result_PureCN_S7/RawCoverageList/out.csv") //Array2CSV(in = CoverageArray)
    //val CoverageList = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v3/result_PureCN_v1/CoverageList/out.csv")

    val TissueSamples = REvaluate(table1 = CoverageList,
                            script = """
                                suppressMessages({
                                    library(dplyr)
                                })
                                    table1 <- table1 %>%
                                            filter(!grepl('BDNA', Key))

                                    table.out = data.frame(table1)
                            """
                            )
    // Create PoN based on BDNA samples
    val makePoNfromBDNA = BashEvaluate(var1 = CoverageList,
                                var2 = interval_GCannot,
                                var3 = centromere,
                                script = """
                                Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v4/code/R_scripts/MakingPoN_ichor.R \
                                -f @var1@ -g @var2@ -c @var3@ -o @out1@

                                """)
    makePoNfromBDNA._filename("out1","normDB.rds")
    makePoNfromBDNA._custom("memory") = "3G"

    val makePoNfromTPL = BashEvaluate(var1 = CoverageList,
                                var2 = interval_GCannot,
                                var3 = centromere,
                                var4 = tp53_info,
                                script = """
                                Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v4/code/R_scripts/MakingPoN_ichor.R \
                                -f @var1@ -g @var2@ -c @var3@ --tp53 @var4@ --pon "TPL" -o @out1@ 

                                """)
    makePoNfromTPL._filename("out1","normDB.rds")
    makePoNfromTPL._custom("memory") = "3G"

    val CombinePon = BashEvaluate(var1 = CoverageList,
                                var2 = interval_GCannot,
                                var3 = centromere,
                                var4 = tp53_info,
                                script = """
                                Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v4/code/R_scripts/PoNCombine_ichor.R \
                                -f @var1@ -g @var2@ -c @var3@ --tp53 @var4@ --pon "TPL" -o @out1@ 

                                """)
    CombinePon._filename("out1","normDB.rds")
    CombinePon._custom("memory") = "3G"


    for ( row <- iterCSV(TissueSamples.table)) {
            //val patient = row("patient")
            val covfile = INPUT(row("File"))
            val sample = row("Key")

            CNVcall(sample) = BashEvaluate(var1 = covfile, 
                                        var2 = interval_GCannot,
                                        var3 = centromere,
                                        var4 = makePoNfromBDNA.out1, 
                                        var5 = makePoNfromTPL.out1,
                                        var6 = CombinePon.out1,
                                        script = """
                                            Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2021/PipelineScript/CNV/oncosysCNA/scripts/segmentforSample.R \
                                            -i @var1@ -g @var2@ -c @var3@ -p @var4@ -o @out1@ -O @folder1@ \
                                            --normal "seq(.01,1,by=.01)" --ploidy "seq(1,6)" --maxCN "20"

                                            Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2021/PipelineScript/CNV/oncosysCNA/scripts/segmentforSample.R \
                                            -i @var1@ -g @var2@ -c @var3@ -p @var5@ -o @out2@ -O @folder2@ \
                                            --normal "seq(.01,1,by=.01)" --ploidy "seq(1,6)" --maxCN "20"

                                            Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2021/PipelineScript/CNV/oncosysCNA/scripts/segmentforSample.R \
                                            -i @var1@ -g @var2@ -c @var3@ -p @var6@ -o @out2@ -O @folder3@ \
                                            --normal "seq(.01,1,by=.01)" --ploidy "seq(1,6)" --maxCN "20"

                                            """)

            CNVcall(sample)._filename("out1","res.rds")
            CNVcall(sample)._filename("out2","out2.rds")
            CNVcall(sample)._filename("out2","out3.rds")
            SegmArray(sample) = CNVcall(sample).out1
        }
        val Segmlist = Array2CSV(in = SegmArray)

}