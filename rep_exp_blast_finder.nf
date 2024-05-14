#!/usr/bin/env nextflow
  
if (params.fq_folder_path == null) {
  println "Please enter the path to the fastq files using the '--fq_folder_path' parameter"
  System.exit(0)
}

if (params.pod5_dir_path == null) {
  println "Please enter the path to the pod5 directory using the '--pod5_dir_path' parameter"
  System.exit(0)
}

if (params.dorado_model_path == null) {
  println "Please enter the path to the dorado model using the '--dorado_model_path' parameter"
  System.exit(0)
}

if (params.out_path == null) {
  println "Please enter the path to the directory in which store the output using the '--out_path' parameter"
  System.exit(0)
}

process GUNZIP_CAT {
	cpus params.cpus
	
	output:
	path fastq
	
	script:
  """
	ls $params.fq_folder_path/*.f*q.gz | parallel -j ${task.cpus} gunzip {}
	cat $params.fq_folder_path/*.f*q >fastq
	"""
}

process MAP_REGION {
  cpus params.cpus
  publishDir params.out_path, mode: 'copy'

  input:
  path fastq

  output:
  path 'rep_region_reads.fa'

  script:
  """
  ${params.minimap2_exec_path} -ax map-ont --MD -t ${task.cpus} ${params.ref_mmi_path} $fastq -o sam
	samtools view -h -@ ${task.cpus} -L ${params.bed_path} sam | samtools fasta -@ ${task.cpus} >rep_region_reads.fa
  """
}

process MOD_BASECALL {
  cpus params.cpus
	publishDir params.out_path, mode: 'copy'	

	input:
	path 'rep_region_reads.fa'

  output:
  path mod_fq

  script:
  """
	export CUDA_DEVICE_ORDER=PCI_BUS_ID
	grep '>' rep_region_reads.fa | cut -f 2 -d '>' >rep_region_read_ids.txt
  ${params.dorado_exec_path} basecaller -x 'cuda:1,2,3,4' --modified-bases '5mCG_5hmCG' -l rep_region_read_ids.txt ${params.dorado_model_path} ${params.pod5_dir_path} >mod_bam
	samtools fastq -@ ${task.cpus} -T '*' mod_bam >mod_fq
  """
}

process MOD_ALIGN {
	cpus params.cpus
  publishDir params.out_path, mode: 'copy'

	input:
	path mod_fq

	output:
	path mod_aln_bam

	script:
  """
	${params.minimap2_exec_path} -ay -x map-ont --MD -t ${task.cpus} ${params.ref_mmi_path} $mod_fq -o sam
	samtools sort -@ ${task.cpus} sam >mod_aln_bam
  """

}

process EXTRACT_MODS {
	cpus params.cpus
  publishDir params.out_path, mode: 'copy'

	input:
  path mod_aln_bam

  output:
  path 'per_read_mods.txt'

  script:
  """
	samtools index -@ ${task.cpus} mod_aln_bam
	modkit call-mods -t ${task.cpus} mod_aln_bam mod_aln_bam_mods_called
	modkit extract -t ${task.cpus} mod_aln_bam_mods_called per_read_mods.txt
  """

}

process BLAST {

	cpus params.cpus
  publishDir params.out_path, mode: 'copy'

	input:
	path 'rep_region_reads.fa'

	output:
	path 'rep_region_reads_outfmt6.blastn'

	script:
  """
	blastn -db ${params.blast_db_path} -query rep_region_reads.fa -num_threads ${task.cpus} -outfmt 6 >rep_region_reads_outfmt6.blastn
  """
}

process PYTHON {
	cpus params.cpus
  publishDir params.out_path, mode: 'copy'

	input:
	path 'rep_region_reads.fa'
	path 'rep_region_reads_outfmt6.blastn'

	output:
	path 'blast_rep_region_finder_output.txt'

	script:
  """
	blast_rep_region_finder.py rep_region_reads.fa rep_region_reads_outfmt6.blastn 5000 >blast_rep_region_finder_output.txt
  """
}

process PYTHON2 {
	cpus params.cpus
  publishDir params.out_path, mode: 'copy'

  input:
	path 'rep_region_reads.fa'
  path 'blast_rep_region_finder_output.txt'
  path 'per_read_mods.txt'

  output:
  path 'per_read_expansion_and_methylation_info.txt'

  script:
  """
	parse_methylation_and_expansion_data.py rep_region_reads.fa blast_rep_region_finder_output.txt per_read_mods.txt  >per_read_expansion_and_methylation_info.txt
  """

}

workflow{
	cat_ch = GUNZIP_CAT()
  map_ch = MAP_REGION( cat_ch )
	mbam_ch = MOD_BASECALL ( map_ch )
	maln_ch = MOD_ALIGN( mbam_ch )
	prm_ch = EXTRACT_MODS( maln_ch )
  blast_ch = BLAST( map_ch )
	py_ch = PYTHON( map_ch, blast_ch )
	PYTHON2( map_ch, py_ch, prm_ch )
}
