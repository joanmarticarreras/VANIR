#!/usr/bin/env nextflow
/*
 *
 ======================================================================
 HERPES: Viral amplicon/enrichment analysis and intrahost variant calling
 ======================================================================
 # Homepage / Documentation

 # Authors
 Joan Marti-Carreras <joan.marti@kuleuven.be>
 ---------------------------------------------------------------------
 */

params.pipelineVersion = "0.1.0"

def helpMessage(){
 log.info"""
 =========================================================
       PIPELINE ~ version ${params.pipelineVersion}
 =========================================================
   Usage:

   Command for running the pipeline is as follows:
  
   nextflow run PIPELINE.nf OPTIONS

   OPTIONS:

     Mandatory options:
        Folder to raw fast5 reads [PATH]:                                          --fast5
        Folder to flipflop base-called reads [PATH]:                               --fastq
        Illumina 1st pair-end reads (can be .gz) [FILE]:                           --read1
        Illumina 2nd pair-end reads (can be .gz) [FILE]:                           --read2
        Illumina adapters [FILE]:                                                  --illumina_adapters
        Sample prefix [FILE]:                                                      --prefix
        Reference genome [FILE]:                                                   --reference
        Number of CPUs asked [NUMERIC]:                                            --cpu
        GPU-enabled ('true' enabled // 'false' disabled) [Boolean]:                --gpu
        Guppy-basecalling mode ('precise' for hac // 'fast' for fast) [Boolean]:   --guppy_algorithm
        Expected genome size (STRING, eg 35k, 8g or 250m):                         --genome_size

    NextFlow options [OPTIONAL]:
	Produce an html report with useful metrics of the pipeline [FILE]          -with-report
        Produce a tabular file with tracings of each processes [FILE]              -with-trace
        Produce an html graphic of all process executed [FILE]                     -with-timeline
        Produce a graph-image (.dot/.html/.pdf/.png/.svg) of the pipeline [FILE]   -with-dag


"""}

/*
 * Set-up configuration variables
 */


// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

/*
 * Define the default parameters
 */ 

params.fast5 = "$baseDir/fast5"
params.fastq = "$baseDir/fastq_flipflop"
params.read1 = "$baseDir/R1.fastq.gz"
params.read2 = "$baseDir/R2.fastq.gz"
params.illumina_adapters = "$baseDir/adapters.fasta"
params.reference = "$baseDir/reference.fasta"
params.cpu = 1
params.gpu = 'true'
params.guppy_algorithm = 'precise'
params.prefix = 'sample'


println """\

	NANOPORE SEQUENCING v$params.pipelineVersion - NF PIPELINE
	================================

	Folder to raw fast5 reads:		$params.fast5 
	Folder to flipflop base-called reads:	$params.fastq
	Illumina 1st pair-end reads:		$params.read1
	Illumina 2nd pair-end reads:		$params.read2 
	Illumina adapters:			$params.illumina_adapters
	Sample prefix:				$params.prefix
	Reference genome:			$params.reference
	Number of CPUs asked:			$params.cpu
	GPU-enabled:				$params.gpu
	Guppy basecalling algorithm:		$params.guppy_algorithm	
	Expected genome size:			$params.genome_size
	"""
	.stripIndent()

nanopore5Channel = Channel.fromPath(params.fast5)
nanoporeQChannel = Channel.fromPath(params.fastq)
illumina1Channel = Channel.fromPath(params.read1)
illumina2Channel = Channel.fromPath(params.read2)
illuminaAdaptersChannel = Channel.fromPath(params.illumina_adapters)
referenceGenomeChannel = Channel.fromPath(params.reference)
referenceGenomeChannel.into{referenceBLAST;referenceIlluminaMapping;referenceIlluminaVC;referenceNanoporeMapping}


process nanopore_basecalling{
	publishDir "${params.prefix}/fastq/nanopore", pattern:"nanopore/*" mode:'copy'
	input:
		val(fast5_dir) from nanopore5Channel
		val(fastq_dir) from nanoporeQChannel
		val cpu from params.cpu
		val gpu from params.gpu
		val type from params.guppy_algorithm
		val prefix from params.prefix
	output:
		file("nanopore/${prefix}_pass.fastq") into nanoporePassFastq
		file("nanopore/${prefix}_total.fastq") into nanoporeTotalFastq
		file("nanopore/*") into nanoporeFastq
	shell:
		model=''
		if (type == 'fast')
			model="/binaries/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg"
		else
			model="/binaries/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg"
		
		if (gpu == 'true')
			"""
			guppy_basecaller -i "!{fast5_dir}" -s nanopore --device auto -c "!{model}" --qscore_filtering --disable_pings -r
			cat nanopore/pass/*.fastq > "nanopore/!{prefix}_pass.fastq"
			cat nanopore/fail/*.fastq > "nanopore/!{prefix}_fail.fastq"
			cat "nanopore/!{prefix}_pass.fastq" "nanopore/!{prefix}_fail.fastq" > "nanopore/!{prefix}_total.fastq"
			"""
			
		else
			"""
			guppy_basecaller -i "!{fast5_dir}" -s nanopore --num_callers 1 --cpu_threads_per_caller "!{cpu}" -c "!{MODEL}" --qscore_filtering --disable_pings -r
			cat nanopore/pass/*.fastq > "nanopore/!{prefix}_pass.fastq"
			cat nanopore/fail/*.fastq > "nanopore/!{prefix}_fail.fastq"
			cat "nanopore/!{prefix}_pass.fastq" "nanopore/!{prefix}_fail.fastq" > "nanopore/!{prefix}_total.fastq"
			"""
}


process nanopore_trim_total_porechop{
	publishDir "${params.prefix}/fastq", mode:'copy'	
	input:
		file(fastq_total) from nanoporeTotalFastq	
		val prefix from params.prefix
		val cpu from params.cpu
	output:
		file "${prefix}_total_trim.fastq" into nanoporeTotalFastqTrim
	shell:
		"""
		porechop -i "!{fastq_total}" --middle_threshold 80  --discard_middle --format fastq -t "!{cpu}" -v 0 > "!{prefix}_total_trim.fastq"
		"""
}

nanoporeTotalFastqTrim.into{nanoporeTotalFastqTrimCanu;nanoporeTotalFastqVCReference;nanoporeTotalFastqVCDenovo}

process nanopore_trim_pass_porechop{
	publishDir "${params.prefix}/$params.fastq", mode:'copy'	
	input:
		file(fastq_pass) from nanoporePassFastq
		val prefix from params.prefix
		val cpu from params.cpu
	output:
		file "${prefix}_pass_trim_racon.fastq.gz" into nanoporePassFastqTrimRacon
		file "${prefix}_pass_trim_medaka.fastq.gz" into nanoporePassFastqTrimMedaka
	shell:
		"""
		porechop -i "!{fastq_pass}" --middle_threshold 80  --discard_middle --format fastq.gz -t "!{cpu}" -v 0 > "!{prefix}_pass_trim_racon.fastq.gz"
		cp "!{prefix}_pass_trim_racon.fastq.gz" "!{prefix}_pass_trim_medaka.fastq.gz"
		"""
}


process illumina_trim_bbduk{
	publishDir "${params.prefix}/fastq/illumina", mode:'copy'	
	input:
		val cpu from params.cpu
		val prefix from params.prefix
		file illumina_1 from illumina1Channel
		file illumina_2 from illumina2Channel
		file adapters from illuminaAdaptersChannel
	output:
		file "${prefix}_trim_1.fastq.gz" into illuminaFirstTrim
		file "${prefix}_trim_2.fastq.gz" into illuminaSecondTrim

	shell:
		"""
		bbduk in="!{illumina_1}" in2="!{illumina_2}" ref="!{adapters}" \
		out="!{prefix}_trim_1.fastq.gz" out2="!{prefix}_trim_2.fastq.gz" \
		ecco=t threads="!{cpu}" qtrim=rl trimq=7 minlength=30 \
		minavgquality=15 -eoom
		"""
}

illuminaFirstTrim.into{illuminaFirstTrimPilon;illuminaFirstTrimReference;illuminaFirstTrimDenovo}
illuminaSecondTrim.into{illuminaSecondTrimPilon;illuminaSecondTrimReference;illuminaSecondTrimDenovo}



process assembly_canu{
	publishDir "${params.prefix}/assembly/canu", mode:'copy'
	input:
		file(fastq_total_trim) from nanoporeTotalFastqTrimCanu
		val prefix from params.prefix
		val cpu from params.cpu
		val genome_size from params.genome_size
		

	output:
		file "canu/${prefix}.contigs.fasta" into nanoporeContigs
		file ("canu/*") into nanoporeCanu

	shell:
		"""
		canu -p "!{prefix}" -d canu \
			 genomeSize="!{genome_size}" minReadLength=300 minOverlapLength=100 \
			 readSamplingCoverage=100 maxThreads="!{cpu}" -nanopore-raw "!{fastq_total_trim}"
		"""

}

process filter_contigs_blast{
	publishDir "${params.prefix}/assembly/filter", pattern: "*.out6", mode:'copy'
	publishDir "${params.prefix}/assembly/filter", pattern: "*.list", mode:'copy'
	publishDir "${params.prefix}/assembly/genome", pattern: "*.fasta", mode:'copy'
	input:
		file reference from referenceBLAST
		file contigs from nanoporeContigs
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_blast_vs_reference.out6" into blastOut
		file "${prefix}_blast_vs_reference_best.list" into listHits
		file "${prefix}_genome.fasta" into targetGenome
	shell:
		"""
		blastn -query "!{contigs}" -subject "!{reference}" -outfmt 6 -evalue 1e-10 > "!{prefix}_blast_vs_reference.out6"
		awk '{print \$1}' "!{prefix}_blast_vs_reference.out6"  | sort -u > "!{prefix}_blast_vs_reference_best.list" 
		if [[ \$(awk 'NR==1{print \$9}' "!{prefix}_blast_vs_reference.out6") -lt  \$(awk 'NR==1{print \$10}' "!{prefix}_blast_vs_reference.out6") ]]
		then
		       seqtk subseq "!{contigs}" "!{prefix}_blast_vs_reference_best.list"  > "!{prefix}_genome.fasta" 
		else
		       seqtk subseq "!{contigs}" "!{prefix}_blast_vs_reference_best.list" | seqtk seq -r -A -  > "${prefix}_genome.fasta"
		fi
		"""
}

process nanopore_polishing_racon{
	publishDir "${params.prefix}/assembly/polishing/racon", pattern: "*.fasta", mode:'copy'
	publishDir "${params.prefix}/assembly/mapping/racon", pattern: "*.sam", mode:'copy'
	input:
		file fastq_pass_polishing_racon from nanoporePassFastqTrimRacon
		file target_genome from targetGenome
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_targetGenome_racon4.fasta" into targetGenomeNanoporeRaconPolished
		file "*.sam" into raconMappings
		file "*.fasta" into raconConsenus
	shell:
		"""
		minimap2 -t "!{cpu}" -L -a -x map-ont "!{target_genome}" "!{fastq_pass_polishing_racon}" > "!{prefix}_nanopore_vs_targetGenome.sam"
		racon -t "!{cpu}" -q 12 -e 0.05 "!{fastq_pass_polishing_racon}" "!{prefix}"_nanopore_vs_targetGenome.sam "!{target_genome}" > "!{prefix}_targetGenome_racon1.fasta"
		for i in {1..3}
		do
			ii=\$(( \$i + 1 ))
			minimap2 -t "!{cpu}" -L -a -x map-ont "!{prefix}_targetGenome_racon\$i.fasta" "!{fastq_pass_polishing_racon}" > "!{prefix}_nanopore_vs_targetGenome_racon\$i.sam"
			racon -t "!{cpu}" -q 12 -e 0.05  "!{fastq_pass_polishing_racon}" "!{prefix}_nanopore_vs_targetGenome_racon\$i.sam" "!{prefix}_targetGenome_racon\$i.fasta" > "!{prefix}_targetGenome_racon\$ii.fasta"
		done
		"""
}

process nanopore_polishing_medaka{
	publishDir "${params.prefix}/assembly/polishing/medaka", mode:'copy'
	input:
		file targetGenome_racon4 from targetGenomeNanoporeRaconPolished
		file fastq_pass_polishing_medaka from nanoporePassFastqTrimMedaka
		val cpu from params.cpu
		val TYPE from params.guppy_algorithm
	output:
		file "medaka/consensus.fasta" into targetGenomeNanoporeRaconMedakaPolished
		file "medaka/*" into medakaOutput
	shell:
		MODEL=''
		if(TYPE == 'fast')
			MODEL="r941_min_fast"
		else
			MODEL="r941_min_high"
		"""
		medaka_consensus -i  "!{fastq_pass_polishing_medaka}" -d "!{targetGenome_racon4}" -o medaka -m !{MODEL} -t "!{cpu}"
		"""
}

process illumina_polishing_pilon{
	publishDir "${params.prefix}/assembly/polishing/pilon", mode:'copy'
	input:
		file targetGenome_racon4_medaka from targetGenomeNanoporeRaconMedakaPolished
		file illumina_first_trim_pilon from illuminaFirstTrimPilon 
		file illumina_second_trim_pilon from illuminaSecondTrimPilon
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "pilon/${prefix}.fasta" into targetGenomeNanoporeRaconMedakaIlluminaPilonPolished
		file "pilon/*" into pilonOutput
		file "*am" into pilonMapping
	shell:
		"""
		bwa index "!{targetGenome_racon4_medaka}"
		bwa mem -t "!{cpu}" "!{targetGenome_racon4_medaka}" "!{illumina_first_trim_pilon}" "!{illumina_second_trim_pilon}" > "!{prefix}_illumina_vs_targetGenomeNanoporeRaconMedakaPolished.sam"
		samtools view -@ "!{cpu}" -b "!{prefix}_illumina_vs_targetGenomeNanoporeRaconMedakaPolished.sam" -T "!{targetGenome_racon4_medaka}" | samtools sort -@ "!{cpu}" > "!{prefix}_illumina_vs_targetGenomeNanoporeRaconMedakaPolished.sort.bam"
		samtools index "!{prefix}_illumina_vs_targetGenomeNanoporeRaconMedakaPolished.sort.bam"
		java -Xmx16G -jar /binaries/pilon-1.23.jar --genome "!{targetGenome_racon4_medaka}" --bam "!{prefix}_illumina_vs_targetGenomeNanoporeRaconMedakaPolished.sort.bam" --outdir pilon --output "!{prefix}" --changes --fix all --threads "!{cpu}" --verbose
		"""

}

targetGenomeNanoporeRaconMedakaIlluminaPilonPolished.into{denovoIlluminaMapping;denovoIlluminaVC;denovoNanoporeMapping;denovoNanoporeVC}

process illumina_mapping_reference{
	publishDir "${params.prefix}/variant_calling/reference/", mode:'copy'
	input:
		file reference_illumina_mapping from referenceIlluminaMapping
		file illumina_first_trim_bwa_reference from illuminaFirstTrimReference
		file illumina_second_trim_bwa_reference from illuminaSecondTrimReference
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_illumina_vs_reference.sort.bam" into illuminaMappingReference
		file "*am"
	shell:
		"""
		bwa index "!{reference_illumina_mapping}"
		bwa mem -t "!{cpu}" "!{reference_illumina_mapping}" "!{illumina_first_trim_bwa_reference}" "!{illumina_second_trim_bwa_reference}" > "!{prefix}_illumina_vs_reference.sam"
		samtools view -@ "!{cpu}" -b "!{prefix}_illumina_vs_reference.sam" -T "!{reference_illumina_mapping}" | samtools sort -@ "!{cpu}" - > "!{prefix}_illumina_vs_reference.sort.bam"
		"""

}


process illumina_mapping_denovo{
	publishDir "${params.prefix}/variant_calling/denovo/", mode:'copy'
	input:
		file denovo_illumina_mapping from denovoIlluminaMapping
		file illumina_first_trim_bwa_denovo from illuminaFirstTrimDenovo
		file illumina_second_trim_bwa_denovo from illuminaSecondTrimDenovo
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_illumina_vs_denovo.sort.bam" into illuminaMappingDenovo
		file "*am"
	shell:
		"""
		bwa index "!{denovo_illumina_mapping}"
		bwa mem -t "!{cpu}" "!{denovo_illumina_mapping}" "!{illumina_first_trim_bwa_denovo}" "!{illumina_second_trim_bwa_denovo}" > "!{prefix}_illumina_vs_denovo.sam"
		samtools view -@ "!{cpu}" -b "!{prefix}_illumina_vs_denovo.sam" -T "!{denovo_illumina_mapping}" | samtools sort -@ "!{cpu}" - > "!{prefix}_illumina_vs_denovo.sort.bam"
		"""
}

process illumina_VC_reference{
	publishDir "${params.prefix}/variant_calling/reference", mode:'copy'
	input:
		file reference_Illumina_VC from referenceIlluminaVC
		file illumina_mapping_reference from illuminaMappingReference
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_illumina_SNV_reference.vcf" into illuminaReferenceVCVCF
	shell:
		"""
		lofreq faidx "!{reference_Illumina_VC}"
		lofreq viterbi -f "!{reference_Illumina_VC}" "!{illumina_mapping_reference}" | samtools sort -@ "!{cpu}"  - > "!{prefix}_illumina_vs_reference.sort.vit.bam"
		lofreq alnqual -b -r "!{prefix}_illumina_vs_reference.sort.vit.bam" "!{reference_Illumina_VC}" | samtools sort -@ "!{cpu}" - > "!{prefix}_illumina_vs_reference.sort.vit.alnqual.bam"
		lofreq index "!{prefix}_illumina_vs_reference.sort.vit.alnqual.bam"
		lofreq2_call_pparallel --pp-threads "!{cpu}" -f "!{reference_Illumina_VC}" -o "!{prefix}_illumina_SNV_reference.vcf" --call-indels -b dynamic -s -a 0.001 --use-orphan "!{prefix}_illumina_vs_reference.sort.vit.alnqual.bam"
		"""

}

process illumina_VC_denovo{
	publishDir "${params.prefix}/variant_calling/denovo", mode:'copy'
	input:
		file denovo_Illumina_VC from denovoIlluminaVC
		file illumina_mapping_denovo from illuminaMappingDenovo
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_illumina_SNV_denovo.vcf" into illuminaDenovoVCVCF
	shell:
		"""
		lofreq faidx "!{denovo_Illumina_VC}"
		lofreq viterbi -f "!{denovo_Illumina_VC}" "!{illumina_mapping_denovo}" | samtools sort -@ "!{cpu}"  - > "!{prefix}_illumina_vs_denovo.sort.vit.bam"
		lofreq alnqual -b -r "!{prefix}_illumina_vs_denovo.sort.vit.bam" "!{denovo_Illumina_VC}" | samtools sort -@ "!{cpu}" - > "!{prefix}_illumina_vs_denovo.sort.vit.alnqual.bam"
		lofreq index "!{prefix}_illumina_vs_denovo.sort.vit.alnqual.bam"
		lofreq2_call_pparallel --pp-threads "!{cpu}" -f "!{denovo_Illumina_VC}" -o "!{prefix}_illumina_SNV_denovo.vcf" --call-indels -b dynamic -s -a 0.001 --use-orphan "!{prefix}_illumina_vs_denovo.sort.vit.alnqual.bam"
		"""
}

process nanopore_mapping_reference{
	publishDir "${params.prefix}/variant_calling/reference", mode:'copy'
	input:
		file reference_nanopore_mapping from referenceNanoporeMapping
		file nanopore_total_fastq_mapping_reference from nanoporeTotalFastqVCReference
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_nanopore_vs_reference.sort.bam" into nanoporeMappingReference
		file "*am"
	shell:
		"""
		ngmlr --bam-fix -x ont -t "!{cpu}" -r "!{reference_nanopore_mapping}" -q "!{nanopore_total_fastq_mapping_reference}" >  "!{prefix}_nanopore_vs_reference.sam"
		samtools view -@ "!{cpu}" -b "!{prefix}_nanopore_vs_reference.sam" -T "!{reference_nanopore_mapping}" | samtools sort -@ "!{cpu}" - > "!{prefix}_nanopore_vs_reference.sort.bam"
		"""
}

process nanopore_mapping_denovo{
	publishDir "${params.prefix}/variant_calling/denovo", mode:'copy'
	input:
		file denovo_nanopore_mapping from denovoNanoporeMapping
		file nanopore_total_fastq_mapping_denovo from nanoporeTotalFastqVCDenovo
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_nanopore_vs_denovo.sort.bam" into nanoporeMappingDenovo
		file "*am"
	shell:
		"""
		ngmlr --bam-fix -x ont -t "!{cpu}" -r "!{denovo_nanopore_mapping}" -q "!{nanopore_total_fastq_mapping_denovo}" >  "!{prefix}_nanopore_vs_denovo.sam"
		samtools view -@ "!{cpu}" -b "!{prefix}_nanopore_vs_denovo.sam" -T "!{denovo_nanopore_mapping}" | samtools sort -@ "!{cpu}" - > "!{prefix}_nanopore_vs_denovo.sort.bam"
		"""
}

process nanopore_VC_reference{
	publishDir "${params.prefix}/variant_calling/reference", mode:'copy'
	input:
		file nanopore_mapping_reference from nanoporeMappingReference
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_nanopore_SV_reference.vcf" into nanoporeReferenceVCVCF

	shell:
		"""
		sniffles -t "!{cpu}" -m "!{nanopore_mapping_reference}" -v "!{prefix}_nanopore_SV_reference.vcf" -s 2 --genotype --cluster --report_read_strands
		"""

}

process nanopore_VC_denovo{
	publishDir "${params.prefix}/variant_calling/denovo", mode:'copy'
	input:
		file nanopore_mapping_denovo from nanoporeMappingDenovo
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_nanopore_SV_denovo.vcf" into nanoporeDenovoVCVCF
	shell:
		"""
		sniffles -t "!{cpu}" -m "!{nanopore_mapping_denovo}" -v "!{prefix}_nanopore_SV_denovo.vcf" -s 2 --genotype --cluster --report_read_strands
		"""
}
