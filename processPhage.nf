#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Script Parameters
params.refdir = "$projectDir/references"

// cutadapt based quality filtering and QC
// Filter to only reads with primer detected and trim poly-g
process cutadapt {
    
	publishDir "${params.outdir}/$pair_id",
		saveAs: {filename ->
			if (filename.indexOf(".fastq.gz") > 0) "filtered/$filename"
			else if (filename.indexOf(".cutadapt.log") > 0) "logs/$filename"
			else filename
		}

	input:
	tuple val(pair_id), path(reads) 

	output:
	tuple val(pair_id), path("filtered_${pair_id}_R{1,2}.fastq.gz"), path("${pair_id}.cutadapt.log")

 	conda 'bioconda::cutadapt'

	script:
	"""
	cutadapt \
	    --action=trim \
        --nextseq-trim=20 \
        --discard-untrimmed \
        -g ^NNNGTGGTTGGTGCTGTAGGAGCA -g ^NNGTGGTTGGTGCTGTAGGAGCA -g ^NGTGGTTGGTGCTGTAGGAGCA -g ^GTGGTTGGTGCTGTAGGAGCA \
		-G ^NNNGAGGCCATGGCATATGCTTATCA -G ^NNGAGGCCATGGCATATGCTTATCA -G ^NGAGGCCATGGCATATGCTTATCA -G ^GAGGCCATGGCATATGCTTATCA \
        -o filtered_${pair_id}_R1.fastq.gz \
        -p filtered_${pair_id}_R2.fastq.gz \
        ${reads[0]} \
        ${reads[1]} \
        > ${pair_id}.cutadapt.log
	"""
}


// bwa based alignment
process bwa_align {

	tag "bwa aligning ${pair_id}"
	label 'big_mem'

	publishDir "${params.outdir}/$pair_id"

	input:
	tuple val(pair_id), path(reads), path(cutadapt_log) 
	path refdir

	output:
	tuple val(pair_id), path("${pair_id}.bam"), path("${pair_id}.bam.bai")

	conda 'bioconda::bwa bioconda::samtools'

	script:
	"""
	# alignment
	bwa mem -t 8 $refdir/all_falciparome_targets_no_primer.fasta ${reads} | samtools view -b -o temp.bam

	# sorting reads
	samtools sort -@ 8 -o ${pair_id}.bam temp.bam

	# index bam
	samtools index ${pair_id}.bam

	# remove intermediary bams
	rm temp.bam
	"""
}


// target mapping stats
process mapping_statistics {

	publishDir "${params.outdir}/$pair_id"

	input:
	tuple val(pair_id), path(bam), path(bam_index)

	output:
	path("${pair_id}_mapping.tab")
	path("${pair_id}_mapping.summary.tab")
	path("${pair_id}_mapping.unfiltered.tab")
	path("${pair_id}_mapping.unfiltered.summary.tab")
	path("${pair_id}_mapping.semifiltered.tab")
	path("${pair_id}_mapping.semifiltered.summary.tab")

 	conda 'bioconda::samtools'

	errorStrategy { task.exitStatus=1 ? 'ignore' : 'terminate' }
	//validExitStatus 0,1

	// require mapped in proper pair, and first in the pair

	script:
	"""
	( samtools view -F 2304 ${bam} | \
		cut -f 3,6  | \
		grep -e "[IDSH*]" -v >> ${pair_id}_mapping.tab && \
	cut -f 1 ${pair_id}_mapping.tab | uniq -c > ${pair_id}_mapping.summary.tab ) || \
	( echo "#${pair_id} had no good reads!" > ${pair_id}_mapping.tab && \
	echo "#${pair_id} had no good reads!" > ${pair_id}_mapping.summary.tab )
	
	( samtools view -F 2304 ${bam} | \
		cut -f 3,6 >> ${pair_id}_mapping.unfiltered.tab && \
	cut -f 1 ${pair_id}_mapping.unfiltered.tab | uniq -c > ${pair_id}_mapping.unfiltered.summary.tab ) || \
	( echo "#${pair_id} had no good reads!" > ${pair_id}_mapping.unfiltered.tab && \
	echo "#${pair_id} had no good reads!" > ${pair_id}_mapping.unfiltered.summary.tab )
	
	( samtools view -F 2304 ${bam} | \
		cut -f 3,6  | \
		grep -e "[ID*]" -v >> ${pair_id}_mapping.semifiltered.tab && \
	cut -f 1 ${pair_id}_mapping.tab | uniq -c > ${pair_id}_mapping.semifiltered.summary.tab ) || \
	( echo "#${pair_id} had no good reads!" > ${pair_id}_mapping.semifiltered.tab && \
	echo "#${pair_id} had no good reads!" > ${pair_id}_mapping.semifiltered.summary.tab )
	"""
}

workflow {
	/*
	Create 'read_pairs' channel that emits for each read pair a
	tuple containing 3 elements: pair_id, R1, R2
	*/
	Channel
        .fromFilePairs(params.reads, checkIfExists: true)
		.set{read_pairs_ch}
		//.ifEmpty{error "Cannot find any reads matching: ${params.reads}"}

    // cutadapt reads
    filtered_reads_ch = cutadapt(read_pairs_ch)

	// bwa align 
	bam_list_ch = bwa_align(filtered_reads_ch, params.refdir)

	// target mapping stats
	mapping_stats = mapping_statistics(bam_list_ch)

}

