# import necessary modules
import argparse                 # For parsing command-line arguments
import create_db                # Module to handle database creation
from load_csv import Load_csv   # Class to load and process CSV files
from load_tsv import Load_tsv   # Class to load and process TSV files
import logging                  # For logging messages and handling potential errors
import os                       # For file system operations
from query import Queries       # Class to handle database queries


# handling the proteintial errors
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(ch)

# set up argparse and define arguments
parser = argparse.ArgumentParser(description='This is PD CW2 ', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--createdb', required= False, action='store_true', help='create the database')
parser.add_argument('--loaddb', required= False, action= 'store_true', help='load files and insert data into the database')
parser.add_argument('--querydb', required= False, type= int, help= 'execute the specified query (1-9)')
parser.add_argument('db_file', help= 'This is SQLite datbase file ')

args = parser.parse_args()

# The block for --createdb
if args.createdb: 
   # Check if args.db_file exits in current directory
   if os.path.exists(args.db_file): 
      print(f"The database file: {args.db_file} alreay exists. please check and try again")
   # if not, create database using the create_db module
   else: 
      create_db.Create_DB(args.db_file)
      print(f'Database {args.db_file} created successfully')


# The block for --loaddb
# Assume that provided files are in current directory 
if args.loaddb:       
   try:
      # Load the csv files (Subject and Metabolome_Annotation)
      csv = Load_csv(args.db_file)  # call Load_csv class
      csv.load_subject('Subject.csv')
      csv.load_annotation('HMP_metabolome_annotation.csv')

      # Load the tsv files ( omics files)
      tsv = Load_tsv(args.db_file) # call Load_tsv class 
      tsv.load_transcripts('HMP_transcriptome_abundance.tsv')
      tsv.load_proteins('HMP_proteome_abundance.tsv' )
      tsv.load_peaks('HMP_metabolome_abundance.tsv')
   except FileNotFoundError as nofile: 
      logger.error(f"FilenotFoundError: {nofile}. Please check file's path or this file exist in current directory")
   else: 
      # Log a success message after loading data
      print(f'It completes to load and insert the data into the database')


# The block for --querydb(1-9)
if args.querydb: 
   # Initialize the Queries class with the database file
   q = Queries(args.db_file) # call Queries class 
   
   # Execute the query corresponding to the specified query number   
   if args.querydb == 1 : 
      q.query1()

   elif args.querydb == 2 : 
      q.query2()

   elif args.querydb == 3 : 
      q.query3()

   elif args.querydb == 4 :  
      q.query4()

   elif args.querydb == 5 :  
      q.query5()

   elif args.querydb == 6 : 
      q.query6()

   elif args.querydb == 7 : 
      q.query7()

   elif args.querydb == 8 : 
      q.query8()
      
   elif args.querydb == 9 : 
      q.query9() # scatter plot will be saved in current directory 
   
   # if a query number is not vaild  
   else: 
      print(f'Query number is not vaild. please check and try again between 1 and 9')
     
  
   
