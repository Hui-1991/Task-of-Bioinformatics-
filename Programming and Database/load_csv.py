import logging  # For logging messages and handling errors
import re       # For string pattern matching and manipulation
import sqlite3  # For interacting with the SQLite database

# handling the proteintial errors
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(ch)

# create class to load csv files (e.g. subject, annotation)
class Load_csv:
    # Call constructor
    def __init__(self, db_file):
       self._db_file = db_file 
    
    # Do not allow to change the attribute: db_file
    @property
    def db_file(self): 
        return self._db_file
        
    # create the mathod to load and insert the data from subject file
    def load_subject(self,subject_file):
        # Connection to database
        con = sqlite3.Connection(self.db_file)
        cur = con.cursor()
        #Load the subject.csv file 
        with open(subject_file) as csv_file: 
            # skip the header of the file
            csv_file.readline()
            for line in csv_file: 
                # 'NA' or 'Unknown' changes to 'NULL'
                line = line.replace('NA','NULL')
                line = line.replace('Unknown','NULL')
                # Assign the data to individual object for inserting into database schema afterwards 
                # Each object indicates each column in Subject table
                line = line.strip().split(',')
                SubjectID = line[0]
                Sex = line[2]
                Age = line[3]
                BMI = line[4]
                IRIS = line[6]
            
                # insert data to db_file 
                # handle the potential errors 
                try:
                    cur.execute('''
                        INSERT INTO Subject (SubjectID, Sex, Age, BMI, IR_IS_classification)
                        VALUES (?,?,?,?,?) 
                        ''', (SubjectID,Sex,Age,BMI,IRIS) )
                except sqlite3.OperationalError as o:
                    logger.error(f"OperationError: {o}. please check the Subject table or column exist in the database")
                except sqlite3.DataError as d: 
                    logger.error(f'DataError: {d}. please validate the data types and try again in Subject table')
        
        # close data connection        
        con.commit()
        con.close()

    # create the mathod to load and insert the data from metabolomics_annotation file
    def load_annotation(self, annot_file):
        # Connection to database 
        con = sqlite3.connect(self.db_file)
        cur = con.cursor()

        with open(annot_file) as csv_file: 
            # skip the header of the file
            csv_file.readline()
           
            for line in csv_file: 
                # parse the file 
                line = line.rstrip().split(',')
                # Assign the data to individual object for inserting into database schema afterwards 
                # Each object indicates each column in Metabolomics_Annotations table
                PeakID = line[0]
                Metabolite = line[1]
                # To merge the same metabolites by removing the suffix 
                # (e.g. ‘2,3-Dihydroxyvaleric acid(1)’ and ‘2,3-Dihydroxyvaleric acid(2)’)  
                match = re.search(r'\(\d\)$', Metabolite)
                if match: 
                    suffix = match.group()
                    Metabolite = Metabolite.replace(suffix,'')

                # handle many linkings between peaks and metabolites        
                # divide the data that has two KEGG or two metablolites 
                if Metabolite != '' and '|' in Metabolite : 
                    Metabolite = Metabolite.split('|')
                    Metabolite = f'{Metabolite[0]}\n{Metabolite[1]}'
                
                KEGG = line[2]
                if KEGG != '' and '|' in KEGG : 
                    KEGG = KEGG.split('|')
                    KEGG = f'{KEGG[0]}\n{KEGG[1]}'
                
                Pathway = line[5] 
                            
                # insert data to db_file 
                # handle the potential errors
                try:
                    cur.execute('''
                        INSERT INTO Metabolomics_Annotations (PeakID, Metabolite, KEGG, Pathway)
                        VALUES (?,?,?,?) 
                        ''', (PeakID,str(Metabolite),str(KEGG),Pathway) )
                except sqlite3.OperationalError as o:
                    logger.error(f"OperationError: {o}. please check the Annotation table or column exist in the database")
                except sqlite3.DataError as d: 
                    logger.error(f'DataError: {d}. please validate the data types and try again in Annotation table')   
        
        # close data connection  
        con.commit()
        con.close()

       
            