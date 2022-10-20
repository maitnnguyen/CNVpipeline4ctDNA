# CNVpipeline4ctDNA
A pipeline for calling CNV events in ctDNA

Details for each folder:
1. PipelineScript
- Scala script to run the CNAs workflow in Anduril
- There are 2 scripts for 2 workflows based on: PureCN and ichorCNA

2. data:
- Segmentation from CNAs calling pipelines
- Clinical data
- TP53 MAF information

3. materials:
- 2 targeted files for Regions Of Interest (ROIs)
- GC and mappability annotation for ROI (ROIs)
- Telomere coordinations of genome hg38
- Protein-coding gene list from ROIs (ROIs) for defining CNAs status (gain/loss) on gene-level

4. analyses:
- Scripts for processing data
- Notebooks for running analysis and visualization

5. figures:
- Main figures in Manuscript
