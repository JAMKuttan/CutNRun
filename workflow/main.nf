#!/usr/bin/env nextflow

/*

BICF Cut-N-Run Analysis Workflow
#### Homepage / Documentation
https://github.com/JAMKuttan/CutNRun
Licensed under MIT (https://github.com/JAMKuttan/CutNRun/LICENSE.md)


*/

// Path to an input file, or a pattern for multiple inputs
// Note - $baseDir is the location of this workflow file main.nf

// Define Input variables
params.reads = "$baseDir/../test_data/*.fastq.gz"
params.pairedEnd = false
params.designFile = "$baseDir/../test_data/design.txt"
params.genome = 'GRCm38'
params.cutoffRatio = 1.2
params.outDir= "$baseDir/output"
params.extendReadsLen = 100
params.topPeakCount = 600
params.astrocyte = false
params.skipDiff = false
params.skipMotif = false
params.skipPlotProfile = false
params.references = "$baseDir/../docs/references.md"
params.multiqc =  "$baseDir/conf/multiqc_config.yaml"
params.ci = false
params.dev = false


// Assign variables if astrocyte
if (params.astrocyte) {
  print("Running under astrocyte")
  referenceLocation = "/project/shared/bicf_workflow_ref"
  if (params.genome == 'GRCh37') {
    params.bwaIndex = "$referenceLocation/human/$params.genome"
    params.chromSizes = "$referenceLocation/human/$params.genome/genomefile.txt"
    params.fasta = "$referenceLocation/human/$params.genome/genome.fa"
    params.gtf = "$referenceLocation/human/$params.genome/gencode.v19.chr_patch_hapl_scaff.annotation.gtf"
    params.geneNames = "$referenceLocation/human/$params.genome/genenames.txt"
    params.genomeSize = 'hs'
  } else if (params.genome == 'GRCm38') {
    params.bwaIndex = "$referenceLocation/mouse/$params.genome"
    params.chromSizes = "$referenceLocation/mouse/$params.genome/genomefile.txt"
    params.fasta = "$referenceLocation/mouse/$params.genome/genome.fa"
    params.gtf = "$referenceLocation/mouse/$params.genome/gencode.vM20.annotation.gtf"
    params.geneNames = "$referenceLocation/mouse/$params.genome/genenames.txt"
    params.genomeSize = 'mm'
  } else if (params.genome == 'GRCh38') {
    params.bwaIndex = "$referenceLocation/human/$params.genome"
    params.chromSizes = "$referenceLocation/human/$params.genome/genomefile.txt"
    params.fasta = "$referenceLocation/human/$params.genome/genome.fa"
    params.gtf = "$referenceLocation/human/$params.genome/gencode.v25.chr_patch_hapl_scaff.annotation.gtf"
    params.geneNames = "$referenceLocation/human/$params.genome/genenames.txt"
    params.genomeSize = 'hs'
  }
} else {
    params.bwaIndex = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
    params.genomeSize = params.genome ? params.genomes[ params.genome ].genomesize ?: false : false
    params.chromSizes = params.genome ? params.genomes[ params.genome ].chromsizes ?: false : false
    params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
    params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
    params.geneNames = params.genome ? params.genomes[ params.genome ].geneNames ?: false : false
}


// Define List of Files
readsList = Channel
  .fromPath( params.reads )
  .flatten()
  .map { file -> [ file.getFileName().toString(), file.toString() ].join("\t")}
  .collectFile( name: 'fileList.tsv', newLine: true )

// Define regular variables
pairedEnd = params.pairedEnd
designFile = Channel.fromPath(params.designFile)
genomeSize = params.genomeSize
genome = params.genome
chromSizes = params.chromSizes
fasta = params.fasta
cutoffRatio = params.cutoffRatio
outDir = params.outDir
extendReadsLen = params.extendReadsLen
topPeakCount = params.topPeakCount
skipDiff = params.skipDiff
skipMotif = params.skipMotif
skipPlotProfile = params.skipPlotProfile
references = params.references
multiqc = params.multiqc
gtfFile = params.gtf
geneNames = params.geneNames


// Check design file for errors
process checkDesignFile {

  publishDir "$outDir/design", mode: 'copy'

  input:

  file designFile
  file readsList

  output:

  file("design.tsv") into designFilePaths

  script:

  if (pairedEnd) {
    """
    module load python/3.7.x-anaconda
    python3 $baseDir/scripts/check_design.py -d $designFile -f $readsList -p
    """
  }
  else {
    """
    module load python/3.7.x-anaconda
    python $baseDir/scripts/check_design.py -d $designFile -f $readsList
    """
  }

}


// Define channel for raw reads
if (pairedEnd) {
  rawReads = designFilePaths
    .splitCsv(sep: '\t', header: true)
    .map { row -> [ row.sample_id, [row.fastq_read1, row.fastq_read2], row.experiment_id, row.biosample, row.factor, row.treatment, row.replicate, row.control_id ] }
} else {
rawReads = designFilePaths
  .splitCsv(sep: '\t', header: true)
  .map { row -> [ row.sample_id, [row.fastq_read1], row.experiment_id, row.biosample, row.factor, row.treatment, row.replicate, row.control_id ] }
}

// Trim raw reads using trimgalore
process trimReads {

  tag "$sampleId-$replicate"
  publishDir "$outDir/${task.process}/${sampleId}", mode: 'copy'

  input:

  set sampleId, reads, experimentId, biosample, factor, treatment, replicate, controlId from rawReads

  output:

  set sampleId, file('*.fq.gz'), experimentId, biosample, factor, treatment, replicate, controlId into trimmedReads
  file('*trimming_report.txt') into trimgaloreResults
  file('version_*.txt') into trimReadsVersions

  script:

  if (pairedEnd) {
    """
    module load python/3.7.x-anaconda
    module load trimgalore/0.6.4
    python3 $baseDir/scripts/trim_reads.py -f ${reads[0]} ${reads[1]} -s $sampleId -p
    """
  }
  else {
    """
    module load python/3.7.x-anaconda
    module load trimgalore/0.6.4
    python3 $baseDir/scripts/trim_reads.py -f ${reads[0]} -s $sampleId
    """
  }

}


// Align trimmed reads using bowtie2
process alignReads {

  queue '128GB,256GB,256GBv1'
  tag "$sampleId-$replicate"
  publishDir "$outDir/${task.process}/${sampleId}", mode: 'copy'

  input:

  set sampleId, reads, experimentId, biosample, factor, treatment, replicate, controlId from trimmedReads
  file index from bwaIndex.first()

  output:

  set sampleId, file('*.bam'), experimentId, biosample, factor, treatment, replicate, controlId into mappedReads
  file '*.flagstat.qc' into mappedReadsStats
  file('version_*.txt') into alignReadsVersions

  script:

  if (pairedEnd) {
    """
    module load python/3.7.x-anaconda
    module load bowtie2/gcc/2.3.4.3
    module load samtools/1.6
    python3 $baseDir/scripts/map_reads.py -f ${reads[0]} ${reads[1]} -r ${index}/genome.fa -s $sampleId -p
    """
  }
  else {
    """
    module load python/3.7.x-anaconda
    module load bowtie2/gcc/2.3.4.3
    module load samtools/1.6
    python3 $baseDir/scripts/map_reads.py -f $reads -r ${index}/genome.fa -s $sampleId
    """
  }

}


// Dedup reads using sambamba
process filterReads {

  queue '128GB,256GB,256GBv1'
  tag "$sampleId-$replicate"
  publishDir "$outDir/${task.process}/${sampleId}", mode: 'copy'

  input:

  set sampleId, mapped, experimentId, biosample, factor, treatment, replicate, controlId from mappedReads

  output:

  set sampleId, file('*.bam'), file('*.bai'), experimentId, biosample, factor, treatment, replicate, controlId into dedupReads
  set sampleId, file('*.bam'), experimentId, biosample, factor, treatment, replicate, controlId into convertReads
  file '*.flagstat.qc' into dedupReadsStats
  file '*.pbc.qc' into dedupReadsComplexity
  file '*.dedup.qc' into dupReads
  file('version_*.txt') into filterReadsVersions

  script:

  if (pairedEnd) {
    """
    module load python/3.7.x-anaconda
    module load samtools/1.6
    module load picard/2.10.3
    python3 $baseDir/scripts/map_qc.py -b $mapped -p
    """
  }
  else {
    """
    module load python/3.7.x-anaconda
    module load samtools/1.6
    module load picard/2.10.3
    python3 $baseDir/scripts/map_qc.py -b $mapped
    """
  }

}


