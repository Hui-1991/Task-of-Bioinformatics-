import logging                          # For logging messages and handling errors
from matplotlib import pyplot as plt    # For creating plots and visualizations (query 9)
import pandas as pd                     # For data manipulation and analysis (query 9)
import seaborn as sns                   # For advanced data visualization (query 9)
import sqlite3                          # For interacting with the SQLite database

# handling the proteintial errors
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(ch)

# create class to retrieve data for each query
class Queries: 
    # Call constructor
    def __init__(self, db_file): 
        self._db_file = db_file 

    # Do not allow to change the attribute: de_file
    @property
    def db_file(self): 
        return self._db_file
    
    # query_db method is needed to execute the following methods 
    def query_db(self,sql):
        try:        
            # Connection to database
            connection = sqlite3.connect(self.db_file)
            cur = connection.cursor()
            cur.execute(sql)
            # save the returned all resulting rows 
            result = cur.fetchall()
            
            # cloes the connection to database
            cur.close()
            connection.close()
            return result
        except sqlite3.Error as error: 
            logger.error(f'[error]: {error} The error occured in the database')

    def query1 (self):
        try: 
            # DML scirpts for query1 
            sql = "SELECT SubjectID, Age FROM Subject WHERE Age > 70 AND Age != 'NULL';" 
            # call query_db method to get result for sql
            result = self.query_db(sql)
            
            # iterate each row from result, assign it to object and print them 
            # Each object indicates retrieved data 
            for row in result: 
                subjectid = row[0]
                age = row[1]
                print(f'{subjectid}\t{age}')
        except Exception as error: 
            logger.error(f'{error} query1 failed to execute')
            
    def query2 (self):
        try:    
            # DML scirpts for query2
            sql= '''SELECT SubjectID FROM Subject 
                    WHERE Sex = 'F' AND BMI >= 18.5 AND BMI <= 24.9 
                    ORDER BY SubjectID DESC ;'''   
            # call query_db method to get result for sql
            result = self.query_db(sql)

            # iterate each row from result, assign it to object and print them 
            # Each object indicates retrieved data 
            for row in result: 
                subjectid = row[0]
                print(f'{subjectid}')
        except Exception as error: 
            logger.error(f'{error} query2 failed to execute')

    def query3 (self):
        try:
            # DML scirpts for query3
            sql = '''SELECT VisitID FROM Proteomics WHERE SubjectID = 'ZNQOVZV' 
                    UNION SELECT VisitID FROM Metabolomics WHERE SubjectID = 'ZNQOVZV'
                    UNION SELECT VisitID FROM Transcriptomics WHERE SubjectID = 'ZNQOVZV';'''
            # call query_db method to get result for sql
            result = self.query_db(sql)

            # iterate each row from result, assign it to object and print them 
            # Each object indicates retrieved data 
            for row in result: 
                visitid= row[0]
                print(f'{visitid}')
        except Exception as error: 
            logger.error(f'{error} query3 failed to execute')
    
    def query4 (self):
        try:
            # DML scirpts for query4
            sql = '''SELECT DISTINCT Metabolomics.SubjectID
                    FROM Metabolomics
                    INNER JOIN Subject
                    ON Subject.SubjectID = Metabolomics.SubjectID
                    AND Subject.IR_IS_classification = 'IR';'''
            
            # call query_db method to get result for sql
            result = self.query_db(sql)

            # iterate each row from result, assign it to object and print them 
            # Each object indicates retrieved data 
            for row in result: 
                distinct_subjectid= row[0]
                print(f'{distinct_subjectid}')
        except Exception as error: 
            logger.error(f'{error} query4 failed to execute')
    
    def query5 (self):
        try:
            # DML scirpts for query5
            sql = '''SELECT DISTINCT KEGG
                    FROM Metabolomics_Annotations
                    WHERE Metabolomics_Annotations.PeakID 
                    IN ('nHILIC_121.0505_3.5','nHILIC_130.0872_6.3','nHILIC_133.0506_2.3','nHILIC_133.0506_4.4')'''
            # call query_db method to get result for sql
            result = self.query_db(sql)

            # iterate each row from result, assign it to object and print them 
            # Each object indicates retrieved data 
            for row in result: 
                distinct_kegg = row[0]
                print(distinct_kegg)
        except Exception as error: 
            logger.error(f'{error} query5 failed to execute')

    def query6 (self):
        try:
            # DML scirpts for query6
            sql =  '''SELECT MIN(AGE), MAX(AGE), AVG(AGE)
                        FROM Subject
                        WHERE AGE != "NULL"'''
            # call query_db method to get result for sql
            result = self.query_db(sql)

            # iterate each row from result, assign it to object and print them 
            # Each object indicates retrieved data 
            for row in result: 
                min = row[0]
                max = row[1]
                avg = round(row[2],2) # round() to show the second decimal place 
                print(f'Minimum Age: {min}\nMaximim Age: {max}\nAverage Age: {avg}')
        except Exception as error: 
            logger.error(f'{error} query6 failed to execute')
               
    def query7 (self):
        try:
            # DML scirpts for query7
            sql = '''SELECT Pathway, COUNT(*)
                    FROM Metabolomics_Annotations
                    GROUP BY Pathway
                    HAVING COUNT(*) >= 10 AND Pathway != ''
                    ORDER BY  COUNT(*) DESC;
                    '''
            # call query_db method to get result for sql
            result = self.query_db(sql)

            # iterate each row from result, assign it to object and print them 
            # Each object indicates retrieved data 
            for row in result: 
                pathwy = row[0]
                number = row[1]
                print(f'{pathwy:>43}\t{number:>3}') # :> to arrange the results to right aligned and improve readability
        except Exception as error: 
            logger.error(f'{error} query7 failed to execute')

    def query8 (self):
        try:
            # DML scirpts for query8
            sql =  '''SELECT MAX(Transcriptomics.Abundance)
                        FROM Transcriptomics
                        WHERE TranscriptID = 'A1BG' AND SubjectID = 'ZOZOW1T';
                        '''
            # call query_db method to get result for sql      
            result = self.query_db(sql)

            # iterate each row from result, assign it to object and print them 
            # Each object indicates retrieved data 
            for row in result: 
                max_abundance = row[0]
                print(max_abundance)
        except Exception as error: 
            logger.error(f'{error} query8 failed to execute')


    def query9 (self): 
        try:
            # DML scirpts for query9
            sql =  '''SELECT Age, BMI
                    FROM Subject
                    WHERE Age != 'NULL' AND BMI != 'NULL' ;'''
            
            # call query_db method to get result for sql
            result = self.query_db(sql)
            
            # to save each data and create dataframe, which is for making a plot
            AGE = []
            BMI = []

            # iterate each row from result, assign it to object and print them 
            # Each object indicates retrieved data 
            print(f'Age\tBMI')
            for row in result: 
                age = row[0]
                AGE.append(age)
                bmi = row[1]
                BMI.append(bmi)
                print(f'{age}\t{bmi}')

            # create a dataframe 
            df = pd.DataFrame ({
                                'Age': AGE,
                                'BMI': BMI
                                })  

            # make a scatter plot of age vs BMI 
            sns.scatterplot(data=df, x ='Age', y='BMI')
            plt.title('Age vs BMI')   
            # store the resulting plot as a png image file 
            plt.savefig('age_bmi_scatterplot.png')
        except Exception as error: 
            logger.error(f'{error} query9 failed to execute')




