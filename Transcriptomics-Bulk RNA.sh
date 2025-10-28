#! /bin/bash
#PBS -l nodes=1:ppn=4:centos7,cput=24:00:00,walltime=24:00:00
#PBS -N A
#PBS -d /export/biostuds/2981631o/assignment/aim
#PBS -m abe
#PBS -M 2981631O@student.gla.ac.uk
#PBS -q bioinf-stud

# Resource Files
adapter="/export/projects/polyomics/biostuds/data/illumina_adapter.fa"
hs2index="/export/projects/polyomics/Genome/Mus_musculus/mm10/Hisat2Index/chr2"
gtf="/export/projects/polyomics/Genome/Mus_musculus/mm10/annotations/chr2.gtf"
data="/export/biostuds/2981631o/assignment/data"


# Make directories for storing outputs 
hisat_dir="/export/biostuds/2981631o/assignment/hisat2"
stringtie_dir="/export/biostuds/2981631o/assignment/stringtie"
mkdir -p ${hisat_dir}
mkdir -p ${stringtie_dir}

gtflist="list.gtf.txt" # Store paths of all generated GTF file
rm -f ${gtflist}       # Remove the file if it exists 



# Loop over each sample
for sample in s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12  
 
do
     	
    raw="${data}/${sample}.c2.fq"			        # Raw files
    trim1="${data}/${sample}_tr1.fq"			    # First trimmed output (adapter trimming)
    trim2="${data}/${sample}_tr2.fq"			    # Second trimmed output (quality filtering)
    bam="${hisat_dir}/${sample}.bam"			    # BAM file
    sam="${hisat_dir}/${sample}.sam"			    # SAM file
    sorted_bam="${hisat_dir}/${sample}.sort.bam"	# Final sorted BAM file

    scythe -q sanger -a ${adapter} -o ${trim1} ${raw}					            # Adapter trimming
    sickle se -f ${trim1} -t sanger -o ${trim2} -q 10 -l 50                         # Quality filtering
    hisat2 -p 4 --phred33 --rna-strandness R -x ${hs2index} -U ${trim2} -S ${sam} 	# Alignment (stranded protocol)
    samtools view -b -o ${bam} ${sam}							                    # Covert SAM to BAM 
    samtools sort -o ${sorted_bam} ${bam}						                    # Sort BAM file
    rm ${sam} ${bam}						
    rm  ${trim1} ${trim2}

    
    str_smp_dir="${stringtie_dir}/${sample}"	# Create a sample-specific directory for ouptput						
    mkdir -p ${str_smp_dir}
    
    sample_tr_gtf="${str_smp_dir}/${sample}_transcripts.gtf" 			        # Define output GTF file 
    stringtie -p 4 --rf -t -e -B -G ${gtf} -o ${sample_tr_gtf} ${sorted_bam} 	# Assembly(stranded protocol) 

    gtfline="${sample} ${sample_tr_gtf}" # Append GTF file path
    echo ${gtfline} >> ${gtflist}

done

prepDE.py -i ${gtflist} # Generate count matrice 


