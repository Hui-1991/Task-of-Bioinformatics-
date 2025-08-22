import logging          # For logging messages and handling errors
import sqlite3          # For interacting with the SQLite database


# handling the proteintial errors
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(ch)


# create class to load tsv files (e.g. omics data)
class Load_tsv:
    
    # Call constructor 
    def __init__(self, db_file):
       self._db_file = db_file 
    
    # Do not allow to change the attribute: de_file
    @property
    def db_file(self): 
        return self._db_file
    
    # create the mathod to load and insert the data from transcriptomics file 
    def load_transcripts(self, transcripts_file) :
        # Connection to database
        con = sqlite3.connect(self.db_file)
        cur = con.cursor()

        # Open the transcriptomics.tsv File 
        with open(transcripts_file) as tsv: 
            # read the header and extract TranscriptIDs
            header = tsv.readline().rstrip().split('\t')
            transcriptIDs = header[1:]

            #Iterate through each row in the TSV file     
            for row in tsv: 
                row = row.rstrip().split('\t')
                SampleID = row[0]
                SampleID = SampleID.split('-')
                # Each object indicates each column in Transcripomics table
                # Split SampleID into SubjectID and VisitID
                SubjectID = SampleID[0]
                VisitID = SampleID[1]
                #Extract Abundances values 
                Abundances = row[1:]
                # Loop through PeakIDs and their corresponding Abundances
                for num in range(len(transcriptIDs)): 
                    TranscriptID = transcriptIDs[num]
                    Abundance = Abundances[num]


                    # Insert transcriptomics data into the database            
                    try:
                        cur.execute('''
                            INSERT INTO Transcriptomics (TranscriptID, Abundance, VisitID, SubjectID)
                            VALUES (?,?,?,?)
                        ''', (TranscriptID, Abundance, VisitID, SubjectID))
                    # handle the errors           
                    except sqlite3.OperationalError as o:
                        logger.error(f"OperationError: {o}. please check the transcripts table or column exist in the database")
                    except sqlite3.DataError as d: 
                        logger.error(f'DataError: {d}. please validate the data types and try again')

        # close data connection 
        con.commit()
        con.close()
        
    # create the mathod to load and insert the data from proteomics file 
    def load_proteins (self, proteins_file):
        # Connection to database
        con = sqlite3.connect(self.db_file)
        cur = con.cursor()
         # Open the proteomics.tsv File 
        with open(proteins_file) as tsv: 
            # read the header and extract ProteinIDs
            header = tsv.readline().rstrip().split('\t')
            ProteinIDs = header[1:]

            #Iterate through each row in the TSV file     
            for row in tsv: 
                row = row.rstrip().split('\t')
                SampleID = row[0]
                SampleID = SampleID.split('-')
                # Each object indicates each column in Proteomics table
                # Split SampleID into SubjectID and VisitID
                SubjectID = SampleID[0]
                VisitID = SampleID[1]
                #Extract Abundances values 
                Abundances = row[1:]
                # Loop through PeakIDs and their corresponding Abundances
                for num in range(len(ProteinIDs)): 
                    ProteinID = ProteinIDs[num]
                    Abundance = Abundances[num]

                    # Insert data into database
                    try:
                        cur.execute('''
                            INSERT INTO Proteomics (ProteinID, Abundance, VisitID, SubjectID)
                            VALUES (?,?,?,?)
                        ''', (ProteinID, Abundance, VisitID, SubjectID))
                    # handle the errors 
                    except sqlite3.OperationalError as o: 
                        logger.error(f"OperationError: {o}. please check the proteins table or column exist in the database")
                    except sqlite3.DataError as d: 
                        logger.error(f'DatabError: {d}. please validate the data types and try again')

        # close data connection 
        con.commit()
        con.close()
        
      
    # create the mathod to load and insert the data from metabolomics file 
    def load_peaks (self, metabolite_file):
        # Connection to database
        con = sqlite3.connect(self.db_file)
        cur = con.cursor()
        # Open the TSV File 
        with open(metabolite_file) as tsv: 
            # read the header and extract PeakIDs
            header = tsv.readline().rstrip().split('\t')
            PeakIDs = header[1:]

            #Iterate through each row in the TSV file     
            for row in tsv: 
                row = row.rstrip().split('\t')
                SampleID = row[0]
                SampleID = SampleID.split('-')
                # Each object indicates each column in Metabolomics table
                # Split SampleID into SubjectID and VisitID
                SubjectID = SampleID[0]
                VisitID = SampleID[1]
                #Extract Abundances values 
                Abundances = row[1:]
                # Loop through PeakIDs and their corresponding Abundances
                for num in range(len(PeakIDs)): 
                    PeakID = PeakIDs[num]
                    Abundance = Abundances[num]
                    
                    # Insert data into database
                    try:
                        cur.execute('''
                            INSERT INTO Metabolomics (PeakID, Abundance, VisitID, SubjectID)
                            VALUES (?,?,?,?)
                        ''', (PeakID, Abundance, VisitID, SubjectID))
                    # handle the errors 
                    except sqlite3.OperationalError as o: 
                        logger.error(f"OperationError: {o}. please check the peaks table or column exist in the database")
                        SystemExit
                    except sqlite3.DataError as d: 
                        logger.error(f'DataError: {d}. please validate the data types and try again')

        # close data connection 
        con.commit()
        con.close()
        
     