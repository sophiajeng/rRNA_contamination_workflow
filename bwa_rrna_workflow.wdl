import "trimmomatic.wdl" as trim

workflow rrna_workflow {
    String fastq_file
    String samples_file

    #where bwa_fastqs.txt contains fastq files in the form sample1_r1.fastq.gz sample2_r2.fastq.gz (one row per sample)
    Array[Array[File]] input_fastqs = read_tsv(fastq_file)

    #this is simply a file that contains sample names for each row in bwa_fastqs.txt
    Array[String] input_samples = read_lines(samples_file)

    Array[Pair[String, Array[File]]] sample_fastqs = zip(input_samples, input_fastqs)
    
    String bwa_index_dir = "/home/exacloud/gscratch/mcweeney_lab/resources/bwa/ncbi_rrna_9606"    
    String trimmomatic_exe = "/home/exacloud/gscratch/mcweeney_lab/resources/programs/Trimmomatic-0.39/trimmomatic-0.39.jar"
	
	scatter (fq in sample_fastqs) {
		call trim.Trimmomatic{
		    input:
                        trim_exe=trimmomatic_exe,
		        input_fastqs=fq.right,
		        cur_sample=fq.left
        	}
		call bwa_mem{
			input:
			    cur_sample=fq.left,
				forward_samp=Trimmomatic.out_forward_paired,
				reverse_samp=Trimmomatic.out_reverse_paired,
				bwa_index_dir=bwa_index_dir
		}
		call samtools_view {
			input:
				cur_sample=fq.left,
				aligned_sam=bwa_mem.out_sam
		}
		call samtools_flagstat {
			input:
				cur_sample=fq.left,	
				aligned_bam=samtools_view.out_bam
		}
		call unique_reads {
			input:
				cur_sample=fq.left,
				aligned_bam=samtools_view.out_bam
		}
	}

	output{
        Array[File] flagstat_out = samtools_flagstat.out_txt
	Array[File] unique_reads_out = unique_reads.out_txt
    }
}

task bwa_mem{
	String cur_sample
    #need to figure out how to parameterize fastqs
	String forward_samp
	String reverse_samp
	String bwa_index_dir
	
	command {
		set -e
		bwa mem -t 8 "${bwa_index_dir}/rrna_9606.fasta" ${forward_samp} ${reverse_samp} > Sample_${cur_sample}_bwa_rrna.sam	
	}
	runtime {
        runtime_minutes: 600
        requested_memory_per_core: "30G"
        cpus: 9
        maxRetries: 1
    }
	output {
		File out_sam = "Sample_" + "${cur_sample}" + "_bwa_rrna.sam"
	}
}
task samtools_view {
	String cur_sample
	String aligned_sam
    command {
        set -e
        samtools view -bS ${aligned_sam} | samtools sort -o Sample_${cur_sample}_bwa_rrna_sorted.bam
	samtools index Sample_${cur_sample}_bwa_rrna_sorted.bam
    }
    runtime {
        runtime_minutes: 600
        requested_memory_per_core: "10G"
        cpus: 1
        maxRetries: 1
    }
    output {
        File out_bam = "Sample_" + "${cur_sample}" + "_bwa_rrna_sorted.bam"
	File out_bai = "Sample_" + "${cur_sample}" + "_bwa_rrna_sorted.bam.bai"
    }
}
task samtools_flagstat {
    String cur_sample
    String aligned_bam
    command {
        set -e
        samtools flagstat ${aligned_bam} > Sample_${cur_sample}_bwa_rrna.out
    }
    runtime {
        runtime_minutes: 600
        requested_memory_per_core: "10G"
        cpus: 1
        maxRetries: 1
    }
    output {
       File out_txt = "Sample_" + "${cur_sample}" + "_bwa_rrna.out"
    }

}

task unique_reads {
    String cur_sample
    String aligned_bam
    command {
        set -e
    	sh /home/exacloud/gscratch/mcweeney_lab/jengs/wdl_sequencing/unique_reads_chr.sh ${aligned_bam} Sample_${cur_sample}_bwa_rrna_unique_reads.txt
	}
    runtime {
        runtime_minutes: 600
        requested_memory_per_core: "10G"
        cpus: 1
        maxRetries: 1
    }
    output {
       File out_txt = "Sample_" + "${cur_sample}" + "_bwa_rrna_unique_reads.txt"
    }

}

