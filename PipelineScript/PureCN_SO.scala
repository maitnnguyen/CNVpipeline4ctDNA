#!/usr/bin/env anduril
//$OPT --threads 20
//$OPT --pipe "tee $ANDURIL_EXECUTION_DIR/_log"
//$OPT --pipe "anduril-pager --ls"
//$OPT --wrapper slurm-prefix
//$PRE export ANDURIL_SELECTNODE=plain

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.microarray._
import anduril.sequencing._

// *Note: in coverage file for each sample, chrY should be removed
// data contain input file has at least 6 columns: patient | sample | bamfiles | type | batch | TP53
// vcf file by sample is the output from GATK, after filtering 

object purecn{

    val files = INPUT(path="/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/old_pipe/allOldBamlist.csv")
    val fltbam = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v3/result_PureCN_oldseq_v1/FilterBamList/out.csv")
    val germlines = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/old_pipe/result_pon/varsPerNormalCSV/out.csv")
    val intervals = INPUT(path="/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v4/int/ROISO.interval")
    val interval_GCannot = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v4/int/ROISO_annot.intervals")
    val ctDNAinfo = INPUT(path = "/mnt/storageBig8/work/joikkone/ctDNA_s8/ctDNA_allBatches_baseciInfo.csv")
    val tp53_info = INPUT(path = "/mnt/storageBig8/work/giovama/projects/ctDNA/Oseq/analysis/tumor_fraction/TP53_vaf_SAGE.csv")
    val target_bed = INPUT(path = "/mnt/storageBig8/work/nguyenma/projects/ctDNA/Oseq/result/ROI4OldSeq_final.bed")

    // step1: calculate raw coverage from bamfiles
    val Coverage = NamedMap[BashEvaluate]("Coverage")
    val CoverageArray = NamedMap[BinaryFile]("CoverageArray")
    val RawCoverageArray = NamedMap[BinaryFile]("RawCoverageArray")
    // step2: correct GC and mappability
    val CorrectBias = NamedMap[BashEvaluate]("CorrectBias")
    // step3: making PoN
    //val makePoNfromBDNA = NamedMap[BashEvaluate]("makePoNfromBDNA")
    //val makePoNfromZeroTPL = NamedMap[BashEvaluate]("makePoNfromZeroTPL")
    // step4: CNV calling
    val CNVcall = NamedMap[BashEvaluate]("CNVcall")
    val LogRArray = NamedMap[BinaryFile]("LogRArray")
    val SegmArray = NamedMap[BinaryFile]("SegmArray")

    // For filter short reads
    val FltCoverage = NamedMap[BashEvaluate]("FltCoverage")
    val FltRawCoverageArray = NamedMap[BinaryFile]("FltRawCoverageArray")
    val FltCorrectBias = NamedMap[BashEvaluate]("FltCorrectBias")
    val FltCoverageArray = NamedMap[BinaryFile]("FltCoverageArray")
    val FltCNVcall = NamedMap[BashEvaluate]("FltCNVcall")
    val FltLogRArray = NamedMap[BinaryFile]("FltLogRArray")
    val FltSegmArray = NamedMap[BinaryFile]("FltSegmArray")
    
    
    for ( row <- iterCSV(files)) {
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
            RawCoverageArray(sample) = Coverage(sample).out1
            
            CorrectBias(sample) = BashEvaluate(var1 = Coverage(sample).out1,
                                            var2 = interval_GCannot,
                                            script = """
                                                Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/S2_correctBias.R @var1@ @var2@ @out1@
                                                """)
            CorrectBias(sample)._filename("out1","out1.txt")
            CoverageArray(sample) = CorrectBias(sample).out1

    }
    val RawCoverageList = Array2CSV(in = RawCoverageArray)
    val CoverageList = Array2CSV(in = CoverageArray)

    val cov_vcf = REvaluate(table1 = CoverageList,
                            table2 = germlines,
                            table3 = tp53_info,
                            script = """
                            suppressMessages({
                                library(dplyr)
                                library(tidyr)
                                library(stringr)})

                            table1$patient = str_split_fixed(table1$Key,'_',n=2)[,1]
                            names(table1)=c('Key','covfile','patient')
                            table1 = table1[-grep('BDNA', table1$covfile),]

                            table2$patient = str_split_fixed(table2$Key,'_',n=2)[,1]
                            names(table2) = c('Key','vcf','patient')
                            table = left_join(table1, table2[,c('vcf','patient')], by = 'patient')

                            # get info from ctDNAinfo to map TP53 value for sample
                            table3 <- table3 %>% #$Key = table3$sample
                                dplyr::rename(Key = sample, TP53 = vaf)
                            
                            table = inner_join(table, table3[,c('Key', 'TP53')], by = 'Key')
                            #names(table)[5] = 'TP53'

                            # out
                            table.out = table[-grep('BDNA',table$vcf),]
                            """
                    )

    val makePoNfromBDNA = BashEvaluate(var1 = CoverageList,
                                script = """
                                    Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/S3_makingPoN.R @var1@ @out1@

                                    """)
    makePoNfromBDNA._filename("out1","normDB.rds")

    val makePoNfromZeroTPL = BashEvaluate(var1 = cov_vcf.table,
                                var2 = CoverageList,
                                var3 = tp53_info,
                                script = """
                                    Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v4/code/R_scripts/MakingPoN.R @var2@ @out1@ @var3@
                                    
                                    """)
    makePoNfromZeroTPL._filename("out1","normDB.rds")

    val CombinePon = BashEvaluate(var1 = cov_vcf.table,
                                var2 = CoverageList,
                                var3 = tp53_info,
                                script = """
                                    Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/PureCN/update_v4/code/R_scripts/PoNCombine.R @var2@ @out1@ @var3@
                                    
                                    """)
    CombinePon._filename("out1","normDB.rds")

    for ( row <- iterCSV(cov_vcf.table)) {
            //val patient = row("patient")
            val covfile = INPUT(row("covfile"))
            val vcffile = INPUT(row("vcf"))
            val sample = row("Key")

            CNVcall(sample) = BashEvaluate(var1 = covfile, var2 = vcffile,
                                        var3 = makePoNfromBDNA.out1, var4 = interval_GCannot,
                                        var5 = makePoNfromZeroTPL.out1,
                                        var6 = CombinePon.out1,
                                        script = """
                                            Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/CNVcall.R \
                                            @var1@ @var2@ @var3@ @var4@  @out1@ @out2@ @out3@
                                            
                                            Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/CNVcall.R \
                                            @var1@ @var2@ @var5@ @var4@ @folder1@/logR.csv @folder1@/segment.csv @folder1@/result.rds

                                            Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/CNVcall.R \
                                            @var1@ @var2@ @var6@ @var4@ @folder2@/logR.csv @folder2@/segment.csv @folder2@/result.rds
                                            """)
            CNVcall(sample)._filename("out1","logR.csv")
            CNVcall(sample)._filename("out2","segment.csv")
            CNVcall(sample)._filename("out3","result.rds")
            LogRArray(sample) = CNVcall(sample).out1
            SegmArray(sample) = CNVcall(sample).out2
        }
        val LogRlist = Array2CSV(in = LogRArray)
        val Segmlist = Array2CSV(in = SegmArray)
/*
        val logRlistbyPatient = REvaluate(table1 = LogRlist,
                                    script = """
                                        suppressMessages({
                                        library(dplyr)
                                        library(tidyr)
                                        library(stringr)})
                                        table1$patient = str_split_fixed(table1$Key, '_', n=2)[,1]
                                        table.out = table1
                                        array.out = split(table.out[,c("Key","File")],table.out$patient)
                                        """)
        val logRbyPatient = Array2CSV(logRlistbyPatient.outArray)

        for (row <- iterCSV(logRbyPatient.out)){
            val patient = row("Key")
            val logRfile = INPUT(row("File"))

            PWSSeg(patient) = BashEvaluate(var1 = logRfile,
                                        script = """
                                            Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/S5_piecewiseSegment.R @var1@ @out1@ @out2@

                                            """)
            PWSSeg(patient)._filename("out1","segment.csv")
            PWSSeg(patient)._filename("out2","interpolate.csv")
        }*/

    // CNV result for short reads only
    for ( row <- iterCSV(fltbam)) {
            //val patient = row("patient")
            val bamfile = INPUT(row("File"))
            //val vcffile = row("vcf")
            val sample = row("Key")

            FltCoverage(sample) = BashEvaluate(var1 = bamfile,
                                         var2 = intervals,
                                         script = """
                                         Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/S1_coverage_cal.R @var1@ @var2@ @out1@
                                         """
            )
            FltCoverage(sample)._filename("out1","out1.txt")
            FltCoverage(sample)._custom("memory") = "2G"
            FltRawCoverageArray(sample) = FltCoverage(sample).out1
            
            FltCorrectBias(sample) = BashEvaluate(var1 = Coverage(sample).out1,
                                            var2 = interval_GCannot,
                                            script = """
                                                Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/S2_correctBias.R @var1@ @var2@ @out1@
                                                """)
            FltCorrectBias(sample)._filename("out1","out1.txt")
            FltCoverageArray(sample) = FltCorrectBias(sample).out1

    }
    val FltRawCoverageList = Array2CSV(in = FltRawCoverageArray)
    val FltCoverageList = Array2CSV(in = FltCoverageArray)

    val fltcov_vcf = REvaluate(table1 = FltCoverageList,
                            table2 = germlines,
                            table3 = tp53_info,
                            script = """
                            suppressMessages({
                                library(dplyr)
                                library(tidyr)
                                library(stringr)})

                            table1$patient = str_split_fixed(table1$Key,'_',n=2)[,1]
                            names(table1)=c('Key','covfile','patient')
                            table1 = table1[-grep('BDNA', table1$covfile),]

                            table2$patient = str_split_fixed(table2$Key,'_',n=2)[,1]
                            names(table2) = c('Key','vcf','patient')
                            table = left_join(table1, table2[,c('vcf','patient')], by = 'patient')

                            # get info from ctDNAinfo to map TP53 value for sample
                            table3 <- table3 %>% #$Key = table3$sample
                                dplyr::rename(Key = sample, TP53 = vaf)
                            
                            table = inner_join(table, table3[,c('Key', 'TP53')], by = 'Key')
                            #names(table)[5] = 'TP53'

                            # out
                            table.out = table[-grep('BDNA',table$vcf),]
                            """
                    )

    for ( row <- iterCSV(fltcov_vcf.table)) {
            //val patient = row("patient")
            val covfile = INPUT(row("covfile"))
            val vcffile = INPUT(row("vcf"))
            val sample = row("Key")

            FltCNVcall(sample) = BashEvaluate(var1 = covfile, var2 = vcffile,
                                        var3 = makePoNfromBDNA.out1, var4 = interval_GCannot,
                                        var5 = makePoNfromZeroTPL.out1,
                                        var6 = CombinePon.out1,
                                        script = """

                                            Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/CNVcall.R \
                                            @var1@ @var2@ @var3@ @var4@  @out1@ @out2@ @out3@
                                            
                                            Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/CNVcall.R \
                                            @var1@ @var2@ @var5@ @var4@ @folder1@/logR.csv @folder1@/segment.csv @folder1@/result.rds

                                            Rscript /mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/CNVCalling/PureCN/CNVcall.R \
                                            @var1@ @var2@ @var6@ @var4@ @folder2@/logR.csv @folder2@/segment.csv @folder2@/result.rds

                                            """)
            FltCNVcall(sample)._filename("out1","logR.csv")
            FltCNVcall(sample)._filename("out2","segment.csv")
            FltCNVcall(sample)._filename("out3","result.rds")
            FltLogRArray(sample) = FltCNVcall(sample).out1
            FltSegmArray(sample) = FltCNVcall(sample).out2
        }
        val FltLogRlist = Array2CSV(in = FltLogRArray)
        val FltSegmlist = Array2CSV(in = FltSegmArray)

}