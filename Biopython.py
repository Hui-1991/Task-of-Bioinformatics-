# Import libraries
import argparse
from Bio.Seq import Seq, MutableSeq
import gffutils
import logging
import math
import matplotlib.pyplot as plt
import os 
import pandas as pd 
import seaborn as sns
import vcf

# Create logger
logFile = 'Log.txt'
logger = logging.getLogger() 
logger.setLevel(logging.INFO)  # Set logging level to INFO
fh = logging.FileHandler(logFile) # Set logging level to INFO
fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s')) # Define format 
logger.addHandler(fh)

# Parse command-line arguments for VCF, GFF, and FASTA files
parser = argparse.ArgumentParser(description='This is BCPy Assignment', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf', required= True, help='includes information about genetic variants')
parser.add_argument('--gff', required= True, help='provides coordinate and general feature of genes')
parser.add_argument('--fasta', required= True, help='contains nucleotide sequence')

args = parser.parse_args()  # Parse the arguments

############################# Create Functions #############################################
#### Function 1: log and check if each file exist or not ####  
def logck_file(filePath, fileType):
    #filePath: Path to the file to check
    #fileType: Type of the file 
    try:
        if not os.path.isfile(filePath):
            raise FileNotFoundError
    except FileNotFoundError:
        logging.error(f"{fileType} file can not find:{filePath}")
        raise SystemExit(1)
    
#### Function 2: calcuate the location and sequence of variants within transcript ####
def cal_LocSeq(feature, transcript, pos, fastaFile): 
    # feature: CDS feature where the variant is located 
    # transcript: the transcript object from gfftils, which is the parents of feature
    # pos: The position of the variant. 
    # fastaFile: The path to the FASTA file.

    Loc = 0   # The location of the variant in the sequence(Seq)
    seq = ''  # The concatenated sequence 
    
    if transcript.strand == '+': # Forward strand
        # process for the '+' strand 
        for cds in db.children(transcript, featuretype='CDS', order_by='start'):
            seq += cds.sequence(fastaFile, use_strand=True )
            if cds == feature : # varient is in coding region 
                Loc += int(pos) - cds.start +1  
                break 
                                            
            else : 
                Loc += cds.end - cds.start + 1 
    
    elif transcript.strand == '-': # Reverse strand
        # process for the '-' strand 
        for cds in db.children(transcript, featuretype='CDS', order_by='start', reverse=True):
            seq += cds.sequence(fastaFile, use_strand=True) 
            if cds == feature : #  if varients in coding region 
                Loc += cds.end - int(pos) + 1
                break
            
            else : 
                Loc += cds.end - cds.start + 1 
            
    return Loc, seq
###############################################################################################
    
# Call the function 1 for validation of input files 
logck_file(args.vcf, "vcf")
logck_file(args.gff, "gff")
logck_file(args.fasta, "fasta")

#log : confirmation of the filenames given at command line 
logging.info(f"Input vcf file: {args.vcf}")
logging.info(f"Input gff file: {args.gff}")
logging.info(f"Input fasta file: {args.fasta}")

# Database setup
dbFile ='BCPy_assessment.db'

if os.path.isfile(dbFile): # Check if the database already exists
    print(f'{dbFile} already exists.Connecting\n')
    db = gffutils.FeatureDB(dbFile, keep_order=True )   

else: # Create the database if it doesn't exist
    print(f'create {dbFile}\n')
    db = gffutils.create_db(args.gff, dbfn=dbFile, keep_order=True)
    print(f'{dbFile} created successfully')

# Objects to count the number of each variant type for output
non_coding = 0
synonymous = 0
non_synonymous = 0

count_lowqual = 0 # Count lower quality record (record.QUAL<=20)   
highqual_var = [] # List to store high-quality variants (record.QUAL > 20)  

try:   
    vcfReader = vcf.Reader(filename=args.vcf)
    for record in vcfReader:   # Iterate over each variant in the VCF
        if record.QUAL <= 20 : # Filter out low-quality variants
            count_lowqual += 1   
            continue 
        # store the information of variants 
        chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alt = ",".join(str(base) for base in record.ALT) # if there are one more alts 
                
        # Determind the varients is in a coding region. 
        cds_features = list(db.region(seqid=record.CHROM, start=record.POS, end=record.POS, featuretype='CDS'))
        if cds_features:
            for feature in cds_features:                 
                for transcript in db.parents(feature, featuretype= 'mRNA'):  
                    transcript_id = transcript.id
                    
                    # Call the function 2 to calculate the location and sequence
                    varientLoc, seq = cal_LocSeq(feature, transcript, pos, args.fasta)   
                    # VarientLoc: he location of the variant in the sequence(Seq)
                    # seq: The concatenated sequence       
                    if len(seq) % 3 != 0: # Ensure the sequence length is a multiple of 3
                        seq = seq[:len(seq)-len(seq)%3] # Trimming the sequence length                                      
                                                 
                    # Calculate protein location and sequence of origin sequence. 
                    varientLoc_protein = math.ceil(varientLoc/3) # Protein location
                    protein_seq = Seq(seq).translate() # Translate to protein sequence                  
                    
                    # Create a mutable sequence for editing 
                    Mutableseq = MutableSeq(seq)                  
                    # Define a dictionary for complementary base mapping.
                    # This maps each nucleotide to its complementary base.
                    base = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                    try:
                        ## NOTE ##  location of varient -1: 0-based index adjustment for python 
                        # Process each alternate allele
                        for alt in record.ALT: 
                            alt = str(alt)
                            # Check if the transcript strand is reverse ('-').
                            # If so, transform the alternate allele into its complementary base.
                            if transcript.strand == '-': 
                                alt = base.get(alt) # Retrieve the complementary base using the dictionary.
                                Mutableseq[varientLoc-1] = alt # Apply the alternative allele
                            else: 
                                Mutableseq[varientLoc-1] = alt # Apply the alternative allele
                            # Translate the mutable sequence 
                            Mutable_Proteinseq= Mutableseq.translate() 

                            ## Compare the original and mutated protein sequences
                            # varient_type: types of varient (synonymous, non-synonymous, non-coding)
                            # ProteinLoc: variant's position in protein sequence
                            # RefAA: Reference amino acid at the variant's location in the protein sequence.
                            # AltAA: Amino acid produced when the alternative allele is applied.
                            if protein_seq[varientLoc_protein-1] == Mutable_Proteinseq[varientLoc_protein-1] : 
                                synonymous += 1 
                                varient_type = 'Synonymous'
                                ProteinLoc = varientLoc_protein
                                RefAA = protein_seq[varientLoc_protein-1]
                                AltAA = 'NA'
                                
                            else : 
                                non_synonymous += 1 
                                varient_type = 'Non-Synonymous'
                                ProteinLoc = varientLoc_protein
                                RefAA = protein_seq[varientLoc_protein-1]
                                AltAA = Mutable_Proteinseq[varientLoc_protein -1]
                            
                            #Add the variant to the list of high-quality 'Synonymous' and 'Non-synonymous variants 
                            highqual_var.append([chrom,pos,ref,alt,varient_type,transcript_id,ProteinLoc,RefAA,AltAA])
                        
                    except IndexError as e : # to catch index error
                        logger.error (f"Index varientLoc-1:{varientLoc-1} is out of range Mutableseq:{len(Mutableseq)} ")  
                        continue             
                   
        else: # Non-coding variants
            non_coding += 1 
            varient_type = 'Non-coding'
            transcript_id = 'NA'
            ProteinLoc = 'NA'
            RefAA = 'NA'
            AltAA = 'NA'

            # Add the variant to the list of high-quality 'Non-coding' variants 
            highqual_var.append([chrom,pos,ref,alt,varient_type,transcript_id,ProteinLoc,RefAA,AltAA])

    #### Write high-quality variants to a TSV file ####
    tsvFile = "highqual_variants.tsv"
    headers = ["Chrom", "Pos", "Ref", "Alt", "Type", "Transcript", "Protein Location", "Ref AA", "Alt AA"]
    with open(tsvFile, "w") as tsv:
        tsv.write("\t".join(headers) + "\n")
        for var in highqual_var:
            tsv.write("\t".join(map(str,var)) + "\n")
            

    #### Create a bar plot for the proportions of each variant type ####
    varCriteria= {'Synonymous':synonymous,'Non_Synonymous':non_synonymous,'Non_Coding':non_coding}
    total_var = sum(varCriteria.values()) # Calculate the total number of variants
    
    # Create a dataFrame for bar plot
    df = pd.DataFrame({
        'Variant Type': list(varCriteria.keys()),
        'Proportion': [round(val / total_var * 100, 1) for val in varCriteria.values()],
        'Raw Count': list(varCriteria.values())
    })
    
    # Generate the bar plot 
    plot = sns.barplot(data=df, x='Variant Type', y='Proportion') # create bar plot
    plot.bar_label(plot.containers[0], fmt='{:.1f}%') # display each proportion value 
    plot.bar_label(plot.containers[0], labels=[val for val in df['Raw Count']], label_type='center') # display raw count numbers
    plt.ylim(0,50)
    plt.ylabel("Proportion of Varients(%)")
    plt.savefig('Bar plot.png',format='png')
    plt.close()

    # log the outputs     
    logger.info(f'The number of low quality of variants is {count_lowqual}')  # The count of the variants where QUAL <= 20.
    logger.info(f'The location of dbfile is {os.path.abspath(dbFile)}')       # Location of dbFile    
    logger.info(f'The location of a log file is {os.path.abspath(logFile)}')  # Location of log file  
    logger.info(f'The location of tsv file is {os.path.abspath(tsvFile)}')    # Location of tsv file   
    logger.info(f'The barplot saved in {os.path.abspath("Bar plot.png")}')    # Location of barplot file  
    logger.info(f"This script runs successfully")
    print(f'\nThis script runs successfully')
# Handling unexpected errors 
except Exception as e:
    logger.error(f"An error occured: {e}")