import sqlite3  # For interacting with the SQLite database


# Module to create the database 
def Create_DB (db_file): 
    # create connection to the dabase 'db_file' and connect a cursor
    con = sqlite3.connect(db_file)
    cursor = con.cursor()
    
    # Note: In the omics tables (Transcriptomics, Proteomics, Metabolomics), 
    # The primary key is each entityID (transcriptID, ProteinID and PeakID).
    # However, in this DDL script,the primary key isn't defined for these tables . 
    # Because the omics data need to be transformed from a wide to a long format
    # when inserting data, which allows EntityID to appear multiple times (duplicate)
    # across different samples (e.g. across VisitIDs or SubjectIDs). 
    
    # load and run SQL DDL script using the executescript method in sqlite3 .
    try:
        cursor.executescript('''
            CREATE TABLE Subject
        (
            SubjectID VARCHAR(50) NOT NULL,
            Sex CHAR(1) NOT NULL, 
            Age INT, 
            BMI INT,
            IR_IS_classification TEXT,
            PRIMARY KEY (subjectID)	
        );

        CREATE TABLE Transcriptomics
        (
            TranscriptID VARCHAR(50) NOT NULL, 
            Abundance DECIMAL,
            VisitID TEXT,
            SubjectID VARCHAR(50),  
            FOREIGN KEY (SubjectID) REFERENCES Subject(SubjectID) 
        );

        CREATE TABLE Proteomics
        (
            ProteinID VARCHAR(50) NOT NULL, 
            Abundance DECIMAL,
            VisitID TEXT,
            SubjectID VARCHAR(50),  
            FOREIGN KEY (SubjectID) REFERENCES Subject(SubjectID) 
        );

        CREATE TABLE Metabolomics
        (
            PeakID VARCHAR(50) NOT NULL, 
            Abundance DECIMAL,
            VisitID TEXT,
            SubjectID VARCHAR(50),  
            FOREIGN KEY (SubjectID) REFERENCES Subject(SubjectID) 
        );

        CREATE TABLE Metabolomics_Annotations
        (
            PeakID VARCHAR(50) NOT NULL,
            Metabolite VARCHAR(100),
            KEGG VARCHAR(50), 
            Pathway VARCHAR (200), 
            PRIMARY KEY (PeakID)
            FOREIGN KEY (PeakID) REFERENCES Metabolomics(PeakID)
        );

        ''')
    except sqlite3.OperationalError as o : 
        print(f'OperationalError: {o}, please check and try again')
        
    # close data connection 
    con.commit()
    con.close()

