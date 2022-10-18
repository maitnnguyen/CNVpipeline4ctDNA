#!/usr/bin/env anduril
//$OPT --threads 20
//$OPT --pipe "tee $ANDURIL_EXECUTION_DIR/_log"
//$OPT --pipe "anduril-pager --ls"
//$OPT --wrapper slurm-prefix
//$PRE export ANDURIL_SELECTNODE=plain
//$PRE export ANDURIL_NODELIST="evm01 evm02 evm03 evm04 evm05 evm06 evm07 evm08 evm09 evm10 evmfull01 evmfull02"
//$PRE export JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64"
//$PRE export PATH="/mnt/csc-gc8/work/joikkone/miniconda2/bin:/mnt/csc-gc8/work/joikkone/miniconda2/envs/snakemake/bin:$PATH"

import anduril.builtin._
import anduril.tools._
import anduril.anima._
import anduril.microarray._
import anduril.sequencing._
import org.anduril.runtime._

object variantCall{

    val reference = INPUT(path = "/mnt/storageBig8/resources/references_annotations/GRCh38d1vd1/GRCh38.d1.vd1.fa")
	val dbsnp = INPUT(path = "/mnt/csc-gc8/resources/pipeline_resources/gatk/dbsnp_154.hg38.vcf.gz")
    val cosmic = INPUT(path = "/mnt/csc-gc8/resources/references_annotations/variant_annotations/COSMIC/cosmic92/CosmicCodingMuts.vcf.gz")
    val bams = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/int/bamsSet7.csv")
    val ctDNAInfo = INPUT(path = "/mnt/storageBig8/work/joikkone/ctDNA_s8/ctDNA_allBatches_baseciInfo.csv")
    val bamfiles = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/int/bamlistset7.csv")
	//val getPairScript = INPUT(path = "get_matchBams.py")
    val intervals = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/peakcall/BDNAset7/TargetEvaluation/finalTarget4ctDNA.interval_list")
    val PoN = INPUT(path = "/mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/old_pipe/result_pon/ponOVCActDNA/ponOVCActDNA.vcf.gz")
    val recalTables = NamedMap[BashEvaluate]("recalTables")
    val recalBams = NamedMap[BashEvaluate]("recalBams")    
	val recalBamsArray = NamedMap[BinaryFile]("recalBamsArray") 
    val varsBySample = NamedMap[BashEvaluate]("varsBySample")
    val varsBySampleMutect2 = NamedMap[BinaryFile]("varsBySampleMutect2")
    val varsBySamplePass = NamedMap[BashEvaluate]("varsBySamplePass")
    val varsBySamplePassArray = NamedMap[BinaryFile]("varsBySamplePassArray")
    val arrayIn = NamedMap[CSV2Array]("arrayIn")
    val bamArray = NamedMap[CSV2Array]("bamArray")
    val mergedByPatient = NamedMap[BashEvaluate]("mergedByPatient")   
	val varsHC = NamedMap[BashEvaluate]("varsHC")
  val varsHCexMA = NamedMap[BashEvaluate]("varsHCexMA")
  val varsAnno = NamedMap[Annovar]("varsAnno")
  val varsAnnoTable = NamedMap[BashEvaluate]("varsAnnoTable")
  val varsAnnoTableOut = NamedMap[CSVTransformer]("varsAnnoTableOut")
  val varsAnnoTableFildbSNP = NamedMap[CSVDplyr]("varsAnnoTableFildbSNP")
  val varsAnnoTableFilstrand = NamedMap[CSVDplyr]("varsAnnoTableFilstrand")
  val varsAnnoTableShortFormat = NamedMap[CSVDplyr]("varsAnnoTableShortFormat")
  val variantsDPexonic = NamedMap[CSVTransformer]("variantsDPexonic")
  val variantsDP = NamedMap[CSVTransformer]("variantsDP")
  val variantsDPexonicFilstrand = NamedMap[CSVTransformer]("variantsDPexonicFilstrand")
  val variantsDPFilstrand = NamedMap[CSVTransformer]("variantsDPFilstrand")  
  val CADD = NamedMap[BashEvaluate]("CADD")

   /* for (row <- iterCSV(bams)){
        val key = row("Key")
        val bam = INPUT(row("File"))
        
        recalTables(key) = BashEvaluate(var1 = bam,
                                      var2 = reference,
                                      var3 = dbsnp,
                                      var4 = intervals,
                                      script = """java -Xmx4g -jar /opt/share/gatk-3.7/GenomeAnalysisTK.jar -T BaseRecalibrator -R @var2@ -I @var1@ -knownSites @var3@ -L @var4@ -fixMisencodedQuals -o @out1@""")
        recalTables(key)._filename("out1", "recal.table")

        recalBams(key) = BashEvaluate(var1 = bam,
                                      var2 = reference,
                                      var3 = recalTables(key).out1,
                                      script = """java -Xmx4g -jar /opt/share/gatk-3.7/GenomeAnalysisTK.jar -T PrintReads -R @var2@ -I @var1@ -BQSR @var3@ -fixMisencodedQuals -o @out1@""")
        recalBams(key)._filename("out1","out1.bam")
		recalBamsArray(key) = recalBams(key).out1	             
        }
   
    val recalBamsCSV = Array2CSV(in = recalBamsArray) 
    val bamsMatched = BashEvaluate(var1 = recalBamsCSV,
								   var2 = getPairScript,
								   script = """python @var2@ @var1@ @out1@""") 
                   */
    // filtering for mutation call for read_depth: 100
        
    for (row <- iterCSV(bams)){
        val tumorKey = row("cancerKey")
        val normalKey = row("bloodKey")
        val bamTumor = INPUT(row("cancerFile"))
        val bamNormal = INPUT(row("bloodFile"))        
                
        varsBySample(tumorKey) = BashEvaluate(var1 = bamTumor,
                                         var2 = bamNormal,
                                         var3 = reference,
                                         var4 = dbsnp,
                                         var5 = cosmic,
                                         var6 = intervals,
                                         //var7 = PoN,                                        
                                         script = """java -Xmx4g -jar /opt/share/gatk-3.7/GenomeAnalysisTK.jar -T MuTect2 -R \
                                         @var3@ --cosmic @var5@ --dbsnp @var4@ --normal_panel /mnt/storageBig8/work/nguyenma/ctDNA/2020/BGI/old_pipe/result_pon/ponOVCActDNA/ponOVCActDNA.vcf.gz \
                                         -L @var6@ -I:tumor @var1@ -I:normal @var2@ -o @out1@ --enable_strand_artifact_filter --max_alt_alleles_in_normal_count 2 --max_alt_alleles_in_normal_qscore_sum 200 \
                                         --max_alt_allele_in_normal_fraction 0.07 --enable_clustered_read_position_filter""")
        varsBySample(tumorKey)._filename("out1", "out1.vcf")
        varsBySample(tumorKey)._custom("cpu") = "1"
        varsBySample(tumorKey)._custom("memory") = "2G"
        varsBySampleMutect2(tumorKey) = varsBySample(tumorKey).out1
    
        varsBySamplePass(tumorKey) = BashEvaluate(var1 = varsBySample(tumorKey).out1,
                                             script = """
                                                      grep "^#\|PASS" @var1@ > @out1@
                                                      """)
        varsBySamplePass(tumorKey)._filename("out1","out1.vcf")
        varsBySamplePassArray(tumorKey) = varsBySamplePass(tumorKey).out1
        }
    
    val varsBySamplePassCSV = Array2CSV(in = varsBySamplePassArray)
    val varsBySampleSplitByPatient = REvaluate(table1 = varsBySamplePassCSV,
						table2 = ctDNAInfo,
                                               script = """
                                                        library(tidyr); library(dplyr); library(stringr)
							t = filter(table2, table2$type %in% c('tissue','plasma','blood'))
							t$name = gsub('._DNA*','',t$sample)
							tab = filter(table1, table1$Key %in% t$name)
                                                        table.out <- tab %>% extract(Key,into="Patient",remove=F)
                                                        array.out <- split(table.out[,c("Key","File")],table.out$Patient)
                                                        """)
    val varsByPatient = Array2CSV(in = varsBySampleSplitByPatient.outArray)

    val bamsPatient = REvaluate(table1 = bamfiles,
                                  table2 = ctDNAInfo,
                                  script = """
                                        library(tidyr); library(dplyr); library(stringr)
                                        t = filter(table2, table2$type %in% c('tissue', 'plasma','blood'))
                                        t$name = gsub('._DNA*','',t$sample)
                                        tab = filter(table1, table1$sample %in% t$name)
					names(tab) = c('File','Key','Patient')
                                        table.out = tab #%>% extract(patient,into="Patient",remove=F)
                                        array.out = split(table.out[,c("Key","File")],table.out$Patient)
                                           """) 

    val bamsByPatient = Array2CSV(in = bamsPatient.outArray)

    // merge all VCFs from the same patient 
    for (row <- iterCSV(varsByPatient)){
        val patient = row("Key")
        arrayIn(patient) = CSV2Array(in = varsBySampleSplitByPatient.outArray(patient),
                                     keys = "column")
        bamArray(patient) = CSV2Array(in = bamsPatient.outArray(patient),
                                      keys = "column")
        mergedByPatient(patient) = BashEvaluate(var1 = reference,
                                                array1 = arrayIn(patient),
                                                script = "java -Xmx4g -jar /opt/share/gatk-3.7/GenomeAnalysisTK.jar -T CombineVariants -R @var1@ -o @out1@ -genotypeMergeOptions UNIQUIFY " +
                                             """ $( paste -d ' ' <(getarraykeys array1) <(getarrayfiles array1)  | sed 's,^, --variant:,' | tr -d '\\\n' ) --disable_auto_index_creation_and_locking_when_reading_rods""" )
		mergedByPatient(patient)._filename("out1","out1.vcf")


		varsHC(patient) = BashEvaluate(var1 = reference, 
						 var2 = mergedByPatient(patient).out1,
						 var3 = dbsnp,
						 array1 = bamArray(patient),
						 script = """
						 java -Xmx4g -jar /opt/share/gatk-3.7/GenomeAnalysisTK.jar -T HaplotypeCaller -R @var1@ -stand_call_conf 0 --dbsnp @var3@ -A StrandAlleleCountsBySample -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES \
             -alleles @var2@  -mmq 0 -o @out1@ $( getarrayfiles array1  | sed 's,^, -I ,' | tr -d '\\\n' ) --disable_auto_index_creation_and_locking_when_reading_rods""" )	
		varsHC(patient)._filename("out1","out1.vcf")
		// filter out triallelic site (in Mutect, this is automatically done)
		varsHCexMA(patient) = BashEvaluate(var1 = reference,
										   var2  = varsHC(patient).out1,
										   script = """
													java -Xmx4g -jar /opt/share/gatk-3.7/GenomeAnalysisTK.jar -T SelectVariants -R @var1@ -V @var2@ -o @out1@ --restrictAllelesTo BIALLELIC --disable_auto_index_creation_and_locking_when_reading_rods
												    """)
		varsHCexMA(patient)._filename("out1","out1.vcf")		
		
		// annotation, includes both snps and indels
        varsAnno(patient) = Annovar(vcfIn = varsHCexMA(patient).out1,
                                    annovarPath="/mnt/csc-gc8/opt/annovar_20191024/annovar",
                                    annovardb="/mnt/csc-gc8/resources/annovardb/",
                                    buildver="hg38",
                                    inputType="vcf",
                                    protocol="wgEncodeGencodeBasicV34,refGene,cytoBand,avsnp154,kaviar_20150923,cosmic92,cosmic92_custom,clinvar_20210102,dbnsfp40b2a_interpro_domain",
                                    operation="g,g,r,f,f,f,f,f,f",                                 
                                    threads = 6)

      CADD(patient) = BashEvaluate(var1 = varsHCexMA(patient).out1,
                script = """
                PATH="/mnt/csc-gc8/work/joikkone/miniconda2/bin:/mnt/csc-gc8/work/joikkone/miniconda2/envs/snakemake/bin:/usr/lib/jvm/java-11-openjdk-amd64/bin/:/mnt/csc-gc8/work/joikkone/anduril/bin:/mnt/csc-gc8/work/joikkone/anduril/utils:/mnt/csc-gc8/work/joikkone/anduril/bin/scala/bin:/mnt/csc-gc8/work/joikkone/anduril/bin/sbt/bin:/mnt/csc-gc8/work/joikkone/anduril:/opt/share/bin:/opt/anduril-dev/bin:/opt/anduril-dev/utils:/opt/anduril-dev/bin/scala/bin:/opt/anduril-dev/bin/sbt/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/usr/local/repos/tools/bin"
            
            /mnt/storageBig8/resources/utils/wgs_scripts/variants/cadd_annotations.sh --input @var1@ -o @out1@ -T 2

            java -jar /opt/share/snpEff/snpEff/SnpSift.jar extractFields -e "." @out1@ \
                CHROM POS ID REF ALT FILTER \
                CADD_raw CADD_phred > @out2@
        
                """
            )
      CADD(patient)._filename("out1","out1.vcf")
      CADD(patient)._filename("out2","out2.csv")
        
        varsAnnoTable(patient) = BashEvaluate(var1 = varsAnno(patient).vcfOut,
		                                      var2 = reference,
         			                          script = """java -Xmx8g -jar /opt/share/gatk-3.7/GenomeAnalysisTK.jar -T VariantsToTable -R @var2@ --allowMissingData -V @var1@ -o @out1@ -F CHROM -F POS -F REF -F ALT -F ID -F Filter -F cytoBand -F Func.refGene -F Gene.refGene -F GeneDetail.refGene -F ExonicFunc.refGene -F AAChange.refGene \
                                              -F Func.wgEncodeGencodeBasicV34 -F Gene.wgEncodeGencodeBasicV34 -F GeneDetail.wgEncodeGencodeBasicV34 -F ExonicFunc.wgEncodeGencodeBasicV34 \
                                              -F AAChange.wgEncodeGencodeBasicV34 -F avsnp154 -F cosmic92 -F cosmic92_custom  \
                                              -F clinvar_20210102 -F Kaviar_AF -F Kaviar_AC -F Kaviar_AN -F SIFT_score -F SIFT_pred -F Polyphen2_HDIV_score -F Polyphen2_HDIV_pred \
                                              -F Polyphen2_HVAR_score -F Polyphen2_HVAR_pred -F MutationTaster_score -F MutationTaster_pred -F MutationAssessor_score -F MutationAssessor_pred \
                                              -F CADD_phred -F DANN_score -GF AD -GF DP -GF SAC -F FS """)

	    varsAnnoTable(patient)._filename("out1","variants_annotated.table")        

     	// remove variants without variant allele reads in any sample
        // add VAF
		varsAnnoTableOut(patient) = CSVTransformer(csv1=varsAnnoTable(patient).out1,
    	                                           transform1="""
	                                                      csv <- as.data.frame(csv1[,grepl("\\.AD",colnames(csv1))])
                                                          VA <- apply(csv,2,function(x)as.numeric(sapply(strsplit(x,","),function(x)x[2])))
                                                          RA <- apply(csv,2,function(x)as.numeric(sapply(strsplit(x,","),function(x)x[1])))
                                                          colnames(VA) <- gsub("\\.AD","",colnames(VA))
                                                          colnames(RA) <- gsub("\\.AD","",colnames(RA))
                                                          AF <- VA / (VA + RA)
                                                          colnames(AF) <- paste0(colnames(AF),".AF")
                                                          AF <- apply(AF,2,function(x)round(x,3))
                                                          csv1 <- data.frame(csv1,AF)
                                                          csv1$Truncal <- ifelse(rowSums(AF>0)==ncol(AF),"Truncal","Heterogenous")
                                                          idx <- rowSums(VA,na.rm=T)>0
                                                          csv1[idx,]
                                                          """)			
        
     // filter variants in dbSNP except those in cosmic
      varsAnnoTableFildbSNP(patient) = CSVDplyr(csv1 = varsAnnoTableOut(patient),
                                              function1 = """filter(avsnp154 == "." |  cosmic92 != ".")""")     

      // variants should have at least 1 read in each strand for variant allele
    // filter strand bias based on fisher test FS>20
  
      varsAnnoTableFilstrand(patient) = CSVDplyr(csv1 = varsAnnoTableFildbSNP(patient),
                                               script = """source("/mnt/storageBig8/work/nguyenma/ctDNA/2020/pipeline_script/MutationCalling/strandFilter.R")""",
                                               //function1 = """filter(rowSums(select(csv1, ends_with(".AF") )>0.05)>0)""",
                                               function1 = """filterStrand()""")
  			                                         
		varsAnnoTableShortFormat(patient) = CSVDplyr(csv1= varsAnnoTableFildbSNP(patient),
                                                 script="library(tidyr)",
                                                 function1="""mutate(samples = paste(gsub("\\.AD","",colnames(csv1)[grepl("\\.AD",colnames(csv1))]),collapse=";"))""",
                                                 function2="""unite("readCounts",ends_with("AD"),sep=";")""",
                                                 function3="""select(-ends_with("DP"))""",
                                                 function4="""select(-ends_with("AD"))""",
                                                 function5="""select(-ends_with("AF"))""")

		variantsDPexonic(patient) = CSVTransformer(csv1=varsAnnoTableFildbSNP(patient),
                                       transform1="""
                                                  #csv1 <- csv1[csv1$Func.refGene %in% "exonic",]
                                                  csv <- csv1[,grepl("\\.AD",colnames(csv1))]
                                                  VA <- apply(csv,2,function(x)as.numeric(sapply(strsplit(x,","),function(x)x[2])))
                                                  RA <- apply(csv,2,function(x)as.numeric(sapply(strsplit(x,","),function(x)x[1])))
                                                  colnames(VA) <- gsub("\\.AD","",colnames(VA))
                                                  colnames(RA) <- gsub("\\.AD","",colnames(RA))
                                                  csv <- VA / (VA + RA)
                                                  csv <- na.omit(csv)
                                                  csv[rowSums(csv)>0,]
                                                  """)

	variantsDP(patient) = CSVTransformer(csv1=varsAnnoTableFildbSNP(patient),
                                       transform1="""
                                                  csv <- csv1[,grepl("\\.AD",colnames(csv1))]
                                                  VA <- apply(csv,2,function(x)as.numeric(sapply(strsplit(x,","),function(x)x[2])))
                                                  RA <- apply(csv,2,function(x)as.numeric(sapply(strsplit(x,","),function(x)x[1])))
                                                  colnames(VA) <- gsub("\\.AD","",colnames(VA))
                                                  colnames(RA) <- gsub("\\.AD","",colnames(RA))
                                                  csv <- VA / (VA + RA)
                                                  csv <- na.omit(csv)
                                                  csv[rowSums(csv)>0,]
                                                  """)

    variantsDPexonicFilstrand(patient) = CSVTransformer(csv1=varsAnnoTableFilstrand(patient),
                                       transform1="""
                                                  #csv1 <- csv1[csv1$Func.refGene %in% "exonic",]
                                                  csv <- csv1[,grepl("\\.AD",colnames(csv1))]
                                                  VA <- apply(csv,2,function(x)as.numeric(sapply(strsplit(x,","),function(x)x[2])))
                                                  RA <- apply(csv,2,function(x)as.numeric(sapply(strsplit(x,","),function(x)x[1])))
                                                  colnames(VA) <- gsub("\\.AD","",colnames(VA))
                                                  colnames(RA) <- gsub("\\.AD","",colnames(RA))
                                                  csv <- VA / (VA + RA)
                                                  csv <- na.omit(csv)
                                                  csv[rowSums(csv)>0,]
                                                  """)

    variantsDPFilstrand(patient) = CSVTransformer(csv1=varsAnnoTableFilstrand(patient),
                                       transform1="""
                                                  csv <- csv1[,grepl("\\.AD",colnames(csv1))]
                                                  VA <- apply(csv,2,function(x)as.numeric(sapply(strsplit(x,","),function(x)x[2])))
                                                  RA <- apply(csv,2,function(x)as.numeric(sapply(strsplit(x,","),function(x)x[1])))
                                                  colnames(VA) <- gsub("\\.AD","",colnames(VA))
                                                  colnames(RA) <- gsub("\\.AD","",colnames(RA))
                                                  csv <- VA / (VA + RA)
                                                  csv <- na.omit(csv)
                                                  csv[rowSums(csv)>0,]
                                                  """)
  }

	val varsHeatmaps = REvaluate(inArray= variantsDPexonic,
                           script="""
                                  library(pheatmap)
                                  library(viridis)
                                  setwd(document.dir)
                                  for( i in names(array)){
                                    freq <- array[[i]]
                                    bin <- freq
                                    bin[bin>0] <- 1
                                    rownames(bin) <- paste0("snv",c(1:nrow(bin)))
                                    rownames(freq) <- paste0("snv",c(1:nrow(freq)))
                                    bin <- bin[do.call(order,bin),]
                                    bin <- as.data.frame(apply(bin,2,rev))
                                    freq <- freq[rownames(bin),]
                                    tbin <- t(bin)
                                    tfreq <- t(freq)
                                    pheatmap(tbin,show_rownames=T,show_colnames=F,cluster_rows=F,cellheight=10,cellwidth=3,color=c("grey85","#4575B4"),
                                             border_color="white",cluster_cols=F,filename=paste0(i,"_heatmap.pdf"))
                                    pheatmap(tfreq,show_rownames=T,show_colnames=F,cluster_rows=F,cellheight=10,cellwidth=3,color=rev(magma(12)),
                                             breaks=c(0,0.03,0.06,0.1,0.15,0.20,0.25,0.30,0.35,0.45,1),cluster_cols=F,filename=paste0(i,"_vaf_heatmap.pdf"))
                                                              }
                                                    table.out <- data.frame()
                                  """)

	val varsFolder = Array2Folder(in = varsAnnoTableOut,
								  fileMode = "@key@.csv")
    val varsFildbSNPFolder = Array2Folder(in = varsAnnoTableFildbSNP,
                                  fileMode = "@key@.csv")

    val allVars = CSVListJoin(in = varsAnnoTableShortFormat)
    
 //    }
}
/*
    OUTPUT(makeArray(recalBamsArray))
}    
*/
