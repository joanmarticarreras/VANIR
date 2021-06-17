#!/usr/bin/env nextflow

/*
 *
 ================================================================================
 VANIR: Virus Variant calling and de novo Assembly of Nanopore and Illumine Reads
 ================================================================================
 # Homepage / Documentation
 https://github.com/joanmarticarreras/VANIR
 # Authors
 Joan Marti-Carreras <joan.marti.carreras@gmail.com> <joan.marti@kuleuven.be>
 --------------------------------------------------------------------------------
 *
 */

params.pipelineVersion = "0.1.0"

def helpMessage(){
 log.info"""
 =========================================================
      VANIR ~ version ${params.pipelineVersion}
 =========================================================
   Usage:

   Command for running the pipeline is as follows:
  
   nextflow run VANIR.nf OPTIONS

   OPTIONS:

     Mandatory options:
        Folder to raw fast5 reads [PATH]:                                          --fast5
        Illumina 1st pair-end reads (can be .gz) [FILE]:                           --read1
        Illumina 2nd pair-end reads (can be .gz) [FILE]:                           --read2
        Illumina adapters [FILE]:                                                  --illumina_adapters
        Sample prefix [FILE]:                                                      --prefix
        Reference genome (genbank) [FILE]:                                         --reference
        Number of CPUs asked [NUMERIC]:                                            --cpu
        GPU-enabled ('ON' enabled // 'OFF' disabled) [Boolean]:                    --gpu
        Guppy-basecalling mode ('precise' for hac // 'fast' for fast) [Boolean]:   --guppy_algorithm

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
params.read1 = "$baseDir/R1.fastq.gz"
params.read2 = "$baseDir/R2.fastq.gz"
params.illumina_adapters = "$baseDir/adapters.fasta"
params.reference = "$baseDir/reference.gb"
params.cpu = 1
params.gpu = 'ON'
params.guppy_algorithm = 'precise'
params.prefix = 'sample'

println """\

	VANIR v$params.pipelineVersion - NF PIPELINE
	================================

	Folder to raw fast5 reads:		$params.fast5 
	Illumina 1st pair-end reads:		$params.read1
	Illumina 2nd pair-end reads:		$params.read2 
	Illumina adapters:			$params.illumina_adapters
	Sample prefix:				$params.prefix
	Reference genome (genbank):		$params.reference
	Number of CPUs asked:			$params.cpu
	GPU-enabled:				$params.gpu
	Guppy basecalling algorithm:		$params.guppy_algorithm	

	"""
	.stripIndent()

nanopore5Channel = Channel.fromPath(params.fast5)
illumina1Channel = Channel.fromPath(params.read1)
illumina2Channel = Channel.fromPath(params.read2)
illuminaAdaptersChannel = Channel.fromPath(params.illumina_adapters)
referenceGenomeChannel = Channel.fromPath(params.reference)


process genbank_decompose{
	input:
		val gbk from referenceGenomeChannel
		val prefix from params.prefix
	output:
		file("${prefix}_reference.fasta") into reference_fasta
		file("${prefix}_reference.gff") into gff
	script:
		"""
		seqret -feature $gbk -fformat1 gb -offormat2 gff -ofname2 ${prefix}_reference.gff -ofdirectory2 ./ -osformat2 fasta -osdirectory2 ./ -osname2 ${prefix}_reference -auto
		"""
}

reference_fasta.into{referenceBLAST;referenceIlluminaMapping;referenceIlluminaVC;referenceNanoporeMapping;referenceFilter;referenceAnnotation;referenceMedaka}

gff.into{reference_gff;transport_gff}

process nanopore_basecalling{
	publishDir "${params.prefix}/fastq/nanopore/", pattern:"*.*", mode:'copy'
	input:
		val(fast5_dir) from nanopore5Channel
		val cpu from params.cpu
		val gpu from params.gpu
		val type from params.guppy_algorithm
		val prefix from params.prefix
	output:
		file("${prefix}/${prefix}_total.fastq") into nanoporeFastq
		file("${prefix}/fastq/sequencing_summary.txt") into sequencingSummary

	script:
		model=''
		if(type == 'fast')
			model="/binaries/ont-guppy/data/template_r9.4.1_450bps_fast.jsn"
		else
			model="/binaries/ont-guppy/data/template_r9.4.1_450bps_hac.jsn"
		
		if (gpu == "ON")
			"""
			guppy_basecaller -i $fast5_dir -s ${prefix}/fastq --device auto -m $model --qscore_filtering --disable_pings -r --chunk_size 500 --chunks_per_runner 768 --num_callers 14 --gpu_runners_per_device 20
			cat ${prefix}/fastq/pass/*.fastq > ${prefix}/fastq/${prefix}_pass.fastq
			cp ${prefix}/fastq/${prefix}_pass.fastq  ${prefix}/${prefix}_pass.fastq
			cat ${prefix}/fastq/fail/*.fastq > ${prefix}/fastq/${prefix}_fail.fastq
			cat ${prefix}/fastq/${prefix}_pass.fastq ${prefix}/fastq/${prefix}_fail.fastq > ${prefix}/fastq/${prefix}_total.fastq
			cp ${prefix}/fastq/${prefix}_total.fastq  ${prefix}/${prefix}_total.fastq
			"""
			
		else if (gpu == "OFF")
			"""
			guppy_basecaller -i ${fast5_dir} -s ${prefix}/fastq -m $model --qscore_filtering --disable_pings -r --chunk_size 1000 --num_callers 14 
			cat ${prefix}/fastq/pass/*.fastq > ${prefix}/fastq/${prefix}_pass.fastq
			cp ${prefix}/fastq/${prefix}_pass.fastq  ${prefix}/${prefix}_pass.fastq
                        cat ${prefix}/fastq/fail/*.fastq > ${prefix}/fastq/${prefix}_fail.fastq
                        cat ${prefix}/fastq/${prefix}_pass.fastq ${prefix}/fastq/${prefix}_fail.fastq > ${prefix}/fastq/${prefix}_total.fastq
			cp ${prefix}/fastq/${prefix}_total.fastq  ${prefix}/${prefix}_total.fastq
			"""
		else
			error "Invalid resource mode: --gpu ${gpu}"

}

nanoporeFastq.into{nanoporeQC;nanoporeTrim;nanoporeRacon}

process nanopore_QC_run_nanoplot{
	publishDir "${params.prefix}/QC/QC_nanopore_run", pattern:"*.*", mode:'copy'
	input:
		val prefix from params.prefix
		val raw_reads from nanoporeQC
	output:
		file "${prefix}_runQC.html"

	script:

		"""
		nanoQC $raw_reads
		mv nanoQC.html ${prefix}_runQC.html
		"""
}



process nanopore_trim_filter_reads{
	publishDir "${params.prefix}/fastq/nanopore", mode:'copy'	
	input:
		file(fastq_pass) from nanoporeTrim
		val prefix from params.prefix
		val cpu from params.cpu
		file reference from referenceFilter
	output:
		file "${prefix}_pass_trim_single_alvis.fastq" into nanoporePassMAPQ60Fastq
		file "${prefix}_pass_porechop.fastq" into nanoporePassFastq
	script:
		"""
		cp $fastq_pass ${prefix}_pass_porechop.fastq
		filtlong --min_mean_q 10 --min_length 250 -a $reference --trim --split 500 ${prefix}_pass_porechop.fastq > ${prefix}_pass_trim.fastq
		minimap2 -t $cpu -x map-ont -a -L $reference ${prefix}_pass_trim.fastq > reference_vs_reads.sam
		samtools view -q 20 -b reference_vs_reads.sam > reference_vs_reads.bam
		samtools fastq -s out reference_vs_reads.bam > ${prefix}_pass_trim_single.fastq
		minimap2 -t $cpu -x map-ont -L $reference ${prefix}_pass_trim_single.fastq > reference_vs_filter_reads.paf
		java -jar /binaries/alvis/dist/Alvis.jar -type contigAlignment -inputfmt paf -outputfmt svg -in reference_vs_filter_reads.paf -outdir alvis -out alvis -filter -chimeras -printChimeras -chimeraPositions -minChimeraCoveragePC 5 -minChimeraAlignmentPC 5
		cut -f1 alvis/chimeras.txt > alvis/chimeras.list
		cat ${prefix}_pass_trim_single.fastq | paste - - - - | grep -v -f alvis/chimeras.list | tr "\t" "\n" >  ${prefix}_pass_trim_single_alvis.fastq
		"""
}


nanoporePassFastq.into{nanoporePassFastqVCReference;nanoporePassFastqVCDenovo}


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

	script:
		"""
		bbduk in=$illumina_1 in2=$illumina_2 ref=$adapters \
		out=${prefix}_trim_1.fastq.gz out2=${prefix}_trim_2.fastq.gz \
		ecco=t threads=$cpu qtrim=rl trimq=7 minlength=30 \
		minavgquality=15 -eoom
		"""
}

illuminaFirstTrim.into{illuminaFirstTrimPilon;illuminaFirstTrimReference;illuminaFirstTrimDenovo;illuminaFirstTrimQCRun}

illuminaSecondTrim.into{illuminaSecondTrimPilon;illuminaSecondTrimReference;illuminaSecondTrimDenovo;illuminaSecondTrimQCRun}

process illumina_QC_run_fastqc{
	publishDir "${params.prefix}/QC/QC_illumina_run/", pattern:"*.*", mode:'copy'
	input:
		file read1 from illuminaFirstTrimQCRun
		file read2 from illuminaSecondTrimQCRun
		val prefix from params.prefix
		val cpu from params.cpu
	output:
		file "${prefix}_trim_1_fastqc.html"
		file "${prefix}_trim_2_fastqc.html"
		
	script:
		"""
		/binaries/FastQC/fastqc -t $cpu -o . -f fastq $read1 $read2
		"""
}


process assembly_flye{
	publishDir "${params.prefix}/assembly/flye", pattern:"*.*", mode:'copy'
	input:
		file fastq_pass_trim_subsampling from nanoporePassMAPQ60Fastq
		val prefix from params.prefix
		val cpu from params.cpu
	output:
		file "${prefix}_assembly.fasta" into nanoporeContigs


	script:
			"""
			flye -t $cpu -o flye -i 5 --nano-raw $fastq_pass_trim_subsampling
			cat flye/assembly.fasta flye/00-assembly/draft_assembly.fasta > ${prefix}_assembly.fasta
			"""	
}

process filter_structure{
	publishDir "${params.prefix}/assembly/filter", pattern:"*.out6", mode:'copy'
	publishDir "${params.prefix}/assembly/filter", pattern:"*.list", mode:'copy'
	publishDir "${params.prefix}/assembly/genome", pattern:"*.fasta", mode:'copy'
	input:
		file reference from referenceBLAST
		file contigs from nanoporeContigs
		val cpu from params.cpu
		val prefix from params.prefix
		file gff from reference_gff
	output:
		file "${prefix}_assembled_genome.fasta" into targetGenome
	script:
		"""
		blastn -query $contigs -subject $reference -evalue 1e-8 -outfmt 6 | sort -gk12 | tail -n 1 | cut -f 1  > best_contig.list
		seqtk subseq $contigs best_contig.list > best_contig.fasta
		minimap2 -x asm20 -L $reference best_contig.fasta  > corrected.paf
		seqtk subseq best_contig.fasta <( awk '{print \$1"\t"\$3"\t"\$4"\t"\$5}' corrected.paf  ) > corrected.fasta

		nucmer $reference corrected.fasta
		show-tiling -i 70.0 -V 0.0 -p ${prefix}_pseudo.fasta -v 0 out.delta

		grep -v "#" $gff  | awk '{print \$3}' | sort -u > ${prefix}_reference_features.list
		liftoff -p $cpu -f ${prefix}_reference_features.list -a 0.8 -s 0.9 -copies -flank 1 -g $gff -u ${prefix}_untransferred.gff -o ${prefix}_genome.gff ${prefix}_pseudo.fasta $reference
		length=\$(seqtk comp ${prefix}_genome.fasta | awk '{print \$2}')
		awk '{print \$1" "\$4" "\$5" "\$5-\$4}' ${prefix}_genome.gff | sort -k4,4nr | head -n 1 | awk -v len=\$length '{print \$1" "\$2" "\$3" "len}' > ${prefix}_genome.coordinates
		faidx ${prefix}_pseudo.fasta \$(awk '{print \$1":"1"-"\$2}' ${prefix}_genome.coordinates) > pre_${prefix}_genome.fasta
		faidx ${prefix}_pseudo.fasta \$(awk '{print \$1":"\$2"-"\$3}' ${prefix}_genome.coordinates) > in_${prefix}_genome.fasta
		faidx ${prefix}_pseudo.fasta \$(awk '{print \$1":"\$3"-"\$4}' ${prefix}_genome.coordinates) > post_${prefix}_genome.fasta

		if [[ \$(seqtk comp pre_${prefix}_genome.fasta | cut -f2) -lt 10 ]]
		then
			cp in_${prefix}_genome.fasta pre_${prefix}_genome.fasta
		fi

		if [[ \$(seqtk comp post_${prefix}_genome.fasta | cut -f2) -lt 10  ]]
		then
			cp in_${prefix}_genome.fasta post_${prefix}_genome.fasta
		fi

		megamerger -asequence in_${prefix}_genome.fasta -bsequence pre_${prefix}_genome.fasta -wordsize 10 -auto -stdout > in_pre.fasta
		megamerger -asequence in_pre.fasta -bsequence post_${prefix}_genome.fasta -wordsize 10 -auto -stdout > end.fasta
		cat end.fasta > ${prefix}_assembled_genome.fasta


		"""
}


process nanopore_polishing_racon{
	publishDir "${params.prefix}/assembly/polishing", pattern:"*.fasta", mode:'copy'
	input:
		file fastq_pass_polishing_racon from nanoporeRacon
		file target_genome from targetGenome
		val cpu from params.cpu
		val gpu from params.gpu
		val prefix from params.prefix
	output:
		file "${prefix}_targetGenome_racon.fasta" into targetGenomeNanoporeRaconPolished
		file "cleaner_reads.fastq" into nanoporePassFastqTrimMedaka
	script:
		if ( gpu == 'ON' )

			"""
			cat $fastq_pass_polishing_racon > cleaner_reads.fastq	
			minimap2 -t $cpu -L --score-N 0 -x map-ont $target_genome $fastq_pass_polishing_racon  > ${prefix}_nanopore_vs_targetGenome_1.paf
			racon --no-trimming -c 30 -t $cpu -q 12 -e 0.05 $fastq_pass_polishing_racon ${prefix}_nanopore_vs_targetGenome_1.paf $target_genome > ${prefix}_targetGenome_racon.fasta
			sed -i "s/>.*/>${prefix}_racon/g" ${prefix}_targetGenome_racon.fasta
		"""

		else
			"""
			cat $fastq_pass_polishing_racon > cleaner_reads.fastq
			minimap2 -t $cpu --score-N 0 -L -x map-ont $target_genome $fastq_pass_polishing_racon > ${prefix}_nanopore_vs_targetGenome.paf
			racon --no-trimming -c 30 -t $cpu -q 12 -e 0.05 $fastq_pass_polishing_racon ${prefix}_nanopore_vs_targetGenome.paf $target_genome > ${prefix}_targetGenome_racon.fasta
			sed -i "s/N//g" ${prefix}_targetGenome_racon.fasta
			sed -i "s/>.*/>${prefix}_racon/g" ${prefix}_targetGenome_racon.fasta
			"""

}

process nanopore_polishing_medaka{
	publishDir "${params.prefix}/assembly/polishing/", pattern:"*.*", mode:'copy'
	input:
		file targetGenome_racon1 from targetGenomeNanoporeRaconPolished
		file fastq_pass_polishing_medaka from nanoporePassFastqTrimMedaka
		file reference from referenceMedaka
		val cpu from params.cpu
		val type from params.guppy_algorithm
		val prefix from params.prefix
	output:
		file "${prefix}_target_genome_medaka.fasta" into targetGenomeNanoporeRaconMedakaPolished
	script:
		model=''
		if(type == 'fast')
			model="r941_min_fast_g351"
		else
			model="r941_min_high_g351"

		"""
		minimap2 --score-N 0 --MD -a -t $cpu -L -x map-ont $targetGenome_racon1 $fastq_pass_polishing_medaka cleaner_reads.fastq | samtools view -@ $cpu -b -F 2308 - | samtools sort -@ $cpu - > ${prefix}_target_genome_racon.bam
		samtools index ${prefix}_target_genome_racon.bam
		medaka consensus --threads $cpu --model $model ${prefix}_target_genome_racon.bam ${prefix}_target_genome_medaka.consensus
		medaka stitch --jobs $cpu ${prefix}_target_genome_medaka.consensus ${prefix}_pre_target_genome_medaka.fasta
		nucmer $reference ${prefix}_pre_target_genome_medaka.fasta
		show-tiling -v 0 -V 0.0 -i 70 -p ${prefix}_target_genome_medaka.fasta out.delta
		sed -i  "s/>.*/>${prefix}_medaka/g" ${prefix}_target_genome_medaka.fasta
		"""
}

process illumina_polishing_pilon{
	publishDir "${params.prefix}/assembly/polishing/pilon", pattern:"*.*", mode:'copy'
	input:
		file targetGenome_racon4_medaka from targetGenomeNanoporeRaconMedakaPolished
		file illumina_first_trim_pilon from illuminaFirstTrimPilon 
		file illumina_second_trim_pilon from illuminaSecondTrimPilon
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "pilon/${prefix}.fasta" into targetGenomeNanoporeRaconMedakaIlluminaPilonPolished
		file "${prefix}_targetGenome_racon_medaka_pilon.fasta"
	script:
		"""
		bwa index $targetGenome_racon4_medaka
		bwa mem -t $cpu $targetGenome_racon4_medaka $illumina_first_trim_pilon $illumina_second_trim_pilon > ${prefix}_illumina_vs_targetGenomeNanoporeRaconMedakaPolished.sam
		minimap2 --score-N 0 -a -L -x sr -t $cpu $targetGenome_racon4_medaka $illumina_first_trim_pilon $illumina_second_trim_pilon > ${prefix}_illumina_vs_targetGenomeNanoporeRaconMedakaPolished.sam
		samtools view -@ $cpu -b ${prefix}_illumina_vs_targetGenomeNanoporeRaconMedakaPolished.sam -T $targetGenome_racon4_medaka | samtools sort -@ $cpu > ${prefix}_illumina_vs_targetGenomeNanoporeRaconMedakaPolished.sort.bam
		samtools index ${prefix}_illumina_vs_targetGenomeNanoporeRaconMedakaPolished.sort.bam
		java -Xmx64G -jar /binaries/pilon-1.23.jar --genome $targetGenome_racon4_medaka --bam ${prefix}_illumina_vs_targetGenomeNanoporeRaconMedakaPolished.sort.bam --outdir pilon --output ${prefix} --changes --fix all --threads $cpu --verbose
		cp pilon/${prefix}.fasta ${prefix}_pre_targetGenome_racon_medaka_pilon.fasta
		sed  "s/N//g" ${prefix}_pre_targetGenome_racon_medaka_pilon.fasta > ${prefix}_targetGenome_racon_medaka_pilon.fasta
		sed -i "s/>.*/>${prefix}_pilon/g" ${prefix}_targetGenome_racon_medaka_pilon.fasta
		"""
}

targetGenomeNanoporeRaconMedakaIlluminaPilonPolished.into{denovoIlluminaMapping;denovoIlluminaVC;denovoNanoporeMapping;denovoNanoporeVC;QCIlluminaTargetGenome;AnnotationPolishedGenome}


process genome_annotation{
	publishDir "${params.prefix}/genome/", pattern:"*.*", mode:'copy'
	input:
		file gff from transport_gff
		file genome from AnnotationPolishedGenome
		file reference from referenceAnnotation
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_untransferred.gff"
		file "${prefix}.gff"
		file "${prefix}_genome.fasta"
		
	script:
		"""
		grep -v "#" $gff  | awk '{print \$3}' | sort -u > ${prefix}_reference_features.list
		liftoff -p $cpu -f ${prefix}_reference_features.list -a 0.8 -s 0.9 -copies -flank 1 -g $gff -u ${prefix}_untransferred.gff -o ${prefix}.gff $genome $reference
		sed -i 's/Liftoff/VANIR/g' ${prefix}_untransferred.gff
		sed -i 's/Liftoff/VANIR/g' ${prefix}.gff
		cp $genome ${prefix}_genome.fasta
		sed -i "s/>.*/>${prefix}/g" ${prefix}_genome.fasta
		"""

}


process illumina_mapping_reference{
	input:
		file reference_illumina_mapping from referenceIlluminaMapping
		file illumina_first_trim_bwa_reference from illuminaFirstTrimReference
		file illumina_second_trim_bwa_reference from illuminaSecondTrimReference
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_illumina_vs_reference.sort.bam" into illuminaMappingReference
	script:
		"""
		bwa index $reference_illumina_mapping
		bwa mem -t $cpu $reference_illumina_mapping $illumina_first_trim_bwa_reference $illumina_second_trim_bwa_reference > ${prefix}_illumina_vs_reference.sam
		samtools view -@ $cpu -b ${prefix}_illumina_vs_reference.sam -T $reference_illumina_mapping | samtools sort -@ $cpu - > ${prefix}_illumina_vs_reference.sort.bam
		"""
}


process illumina_mapping_denovo{
	input:
		file denovo_illumina_mapping from denovoIlluminaMapping
		file illumina_first_trim_bwa_denovo from illuminaFirstTrimDenovo
		file illumina_second_trim_bwa_denovo from illuminaSecondTrimDenovo
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_illumina_vs_denovo.sort.bam" into illuminaMappingDenovoVC
	script:
		"""
		bwa index $denovo_illumina_mapping
		bwa mem -t $cpu $denovo_illumina_mapping $illumina_first_trim_bwa_denovo $illumina_second_trim_bwa_denovo > ${prefix}_illumina_vs_denovo.sam
		samtools view -@ $cpu -b ${prefix}_illumina_vs_denovo.sam -T $denovo_illumina_mapping | samtools sort -@ $cpu - > ${prefix}_illumina_vs_denovo.sort.bam
		"""
}


process illumina_VC_reference{
	publishDir "${params.prefix}/variant_calling/reference", pattern:"*.*", mode:'copy'
	input:
		file reference_Illumina_VC from referenceIlluminaVC
		file illumina_mapping_reference from illuminaMappingReference
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_illumina_SNV_reference.vcf" into illuminaReferenceVCVCF
	script:
		"""
		lofreq faidx $reference_Illumina_VC
		lofreq viterbi -f $reference_Illumina_VC $illumina_mapping_reference | samtools sort -@ $cpu  - > ${prefix}_illumina_vs_reference.sort.vit.bam
		lofreq alnqual -b -r ${prefix}_illumina_vs_reference.sort.vit.bam $reference_Illumina_VC | samtools sort -@ $cpu - > ${prefix}_illumina_vs_reference.sort.vit.alnqual.bam
		lofreq index ${prefix}_illumina_vs_reference.sort.vit.alnqual.bam
		lofreq2_call_pparallel --pp-threads $cpu -f $reference_Illumina_VC -o ${prefix}_illumina_SNV_reference.vcf --call-indels -b dynamic -s -a 0.001 --use-orphan ${prefix}_illumina_vs_reference.sort.vit.alnqual.bam
		"""
}

process illumina_VC_denovo{
	publishDir "${params.prefix}/variant_calling/denovo", pattern:"*.*", mode:'copy'
	input:
		file denovo_Illumina_VC from denovoIlluminaVC
		file illumina_mapping_denovo from illuminaMappingDenovoVC
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_illumina_SNV_denovo.vcf" into illuminaDenovoVCVCF
	script:
		"""
		lofreq faidx $denovo_Illumina_VC
		lofreq viterbi -f $denovo_Illumina_VC $illumina_mapping_denovo | samtools sort -@ $cpu  - > ${prefix}_illumina_vs_denovo.sort.vit.bam
		lofreq alnqual -b -r ${prefix}_illumina_vs_denovo.sort.vit.bam $denovo_Illumina_VC | samtools sort -@ $cpu - > ${prefix}_illumina_vs_denovo.sort.vit.alnqual.bam
		lofreq index ${prefix}_illumina_vs_denovo.sort.vit.alnqual.bam
		lofreq2_call_pparallel --pp-threads $cpu -f $denovo_Illumina_VC -o ${prefix}_illumina_SNV_denovo.vcf --call-indels -b dynamic -s -a 0.001 --use-orphan ${prefix}_illumina_vs_denovo.sort.vit.alnqual.bam
		"""
}

process nanopore_mapping_reference{
	input:
		file reference_nanopore_mapping from referenceNanoporeMapping
		file nanopore_total_fastq_mapping_reference from nanoporePassFastqVCReference
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_nanopore_vs_reference.sort.bam" into nanoporeMappingReference
	script:
		"""
		ngmlr -i 0.9 --bam-fix -x ont -t $cpu -r $reference_nanopore_mapping -q $nanopore_total_fastq_mapping_reference >  ${prefix}_nanopore_vs_reference.sam
		samtools view -@ $cpu -b ${prefix}_nanopore_vs_reference.sam -T $reference_nanopore_mapping | samtools sort -@ $cpu - > ${prefix}_nanopore_vs_reference.sort.bam
		"""
}

process nanopore_mapping_denovo{
	input:
		file denovo_nanopore_mapping from denovoNanoporeMapping
		file nanopore_total_fastq_mapping_denovo from nanoporePassFastqVCDenovo
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_nanopore_vs_denovo.sort.bam" into nanoporeMappingDenovoVC
	script:
		"""
		ngmlr -i 0.9 --bam-fix -x ont -t $cpu -r $denovo_nanopore_mapping -q $nanopore_total_fastq_mapping_denovo >  ${prefix}_nanopore_vs_denovo.sam
		samtools view -@ $cpu -b ${prefix}_nanopore_vs_denovo.sam -T $denovo_nanopore_mapping | samtools sort -@ $cpu - > ${prefix}_nanopore_vs_denovo.sort.bam
		"""
}


process nanopore_VC_reference{
	publishDir "${params.prefix}/variant_calling/reference", pattern:"*.*", mode:'copy'
	input:
		file nanopore_mapping_reference from nanoporeMappingReference
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_nanopore_SV_reference.vcf" into nanoporeReferenceVCVCF

	script:
		"""
		sniffles -t $cpu -m $nanopore_mapping_reference -v ${prefix}_nanopore_SV_reference.vcf -s 2 -q 30 --report_seq -l 10 --cluster --report_read_strands
		"""
}

process nanopore_VC_denovo{
	publishDir "${params.prefix}/variant_calling/denovo", pattern:"*.*", mode:'copy'
	input:
		file(nanopore_mapping_denovo) from nanoporeMappingDenovoVC
		val cpu from params.cpu
		val prefix from params.prefix
	output:
		file "${prefix}_nanopore_SV_denovo.vcf" into nanoporeDenovoVCVCF
	script:
		"""
		sniffles -t $cpu -m $nanopore_mapping_denovo -v ${prefix}_nanopore_SV_denovo.vcf -s 2 -q 30 --report_seq -l 10 --cluster --report_read_strands
		"""
}

