from flask import render_template
from app import app
from flask import request
from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, create_database
import pandas as pd
import psycopg2
from retrieve_snps import retrieve_snps
from werkzeug import secure_filename
import pandas as pd
import numpy as np
import numpy.random as rand
from psycopg2.extensions import AsIs
from sklearn.ensemble import RandomForestClassifier
import warnings
from predict_phenotype import predict_phenotype
import pickle

warnings.filterwarnings('ignore')    

db = 'openSNP_genotypes'
username = 'postgres'

## 'engine' is a connection to a database
## Here, we're using postgres, but sqlalchemy can connect to other things too.
engine = create_engine('postgresql://%s:%s@localhost/%s'%(username,11235813,db))
print engine.url

## create a database (if it doesn't exist)
if not database_exists(engine.url):
    create_database(engine.url)
print(database_exists(engine.url))

##get genotype and phenotype data from database

# connect:
con = None
con = psycopg2.connect(database = db, user = username, password = "11235813", host = "/run/postgresql/")

# query:
sql_query = """
SELECT * FROM "openSNP_data"
"""
df = pd.read_sql_query(sql_query,con)

# query:
sql_query = """
SELECT * FROM "openSNP_data_hair_color"
"""
df2 = pd.read_sql_query(sql_query,con)

hairsnps = pd.read_csv('hair_color_SNPs.csv')
hairsnps = hairsnps[[g in df2.columns for g in hairsnps['ID']]] #only keep those snps that are in our training set

eyesnps = pd.read_csv('Eye_color_SNPs.csv')
eyesnps = eyesnps[[g in df.columns for g in eyesnps['ID']]] #only keep those snps that are in our training set

@app.route('/snps')
def show_snps():
    snps  = pd.read_csv('Eye_color_SNPs.csv')
    snps_list = []
    for i in range(0,len(snps.index)):
        snps_list.append(dict(ID=snps.iloc[i]['ID'], Chr=snps.iloc[i]['Chr'], 
		Position=snps.iloc[i]['Position'],Reference=snps.iloc[i]['Ref allele'],
		Alternate=snps.iloc[i]['Alt allele'],Genotyped=snps.iloc[i]['Genotyped by 23andMe?'] ))

    return render_template('snps.html',snps=snps_list)

@app.route('/')
def upload_file():
    return render_template("input.html")

@app.route('/output', methods = ['GET', 'POST'])
def upload_file_():
      if request.method == 'POST':
	  fXX = request.files['fileXX']
	  fXX.save(secure_filename(fXX.filename))
          fXY = request.files['fileXY']
          fXY.save(secure_filename(fXY.filename))

	  fXXeyesnps = retrieve_snps(fXX.filename, eyesnps)
	  fXYeyesnps = retrieve_snps(fXY.filename, eyesnps)

	  fXXeyesnps.columns = ['genotypeXX', 'refXX', 'altXX']
	  fXYeyesnps.columns = ['genotypeXY', 'refXY', 'altXY']

	  alleyesnps = pd.merge(fXXeyesnps, fXYeyesnps, how='inner', left_index=True, right_index=True)

	  #load random forest model

	  eye_color_model = pickle.load(open( 'eye_color_model_4colors.pkl', 'rb' ))	  

	  #fill in missing values
  	  eyemedians = pd.read_csv('eye_color_medians.csv')
	  eyemedians.index = eyemedians['Unnamed: 0']

	  fXXeyesnps[['altXX']][fXXeyesnps['genotypeXX'].isnull()] = eyemedians[fXXeyesnps['genotypeXX'].isnull()]
	  fXYeyesnps[['altXY']][fXYeyesnps['genotypeXY'].isnull()] = eyemedians[fXYeyesnps['genotypeXY'].isnull()]
	  eyeprobabilities = np.round(100*predict_phenotype(alleyesnps, eye_color_model)).astype(int) 

	  print 'Eye color probabilities:', eyeprobabilities
	  print 'Eye color classes:', eye_color_model.classes_


          fXXhairsnps = retrieve_snps(fXX.filename, hairsnps)
          fXYhairsnps = retrieve_snps(fXY.filename, hairsnps)

          fXXhairsnps.columns = ['genotypeXX', 'refXX', 'altXX']
          fXYhairsnps.columns = ['genotypeXY', 'refXY', 'altXY']

          allhairsnps = pd.merge(fXXhairsnps, fXYhairsnps, how='inner', left_index=True, right_index=True)

	  print allhairsnps
          #load random forest model

          hair_color_model = pickle.load(open( 'hair_color_model_5colors.pkl', 'rb' ))

          #fill in missing values
          hairmedians = pd.read_csv('hair_color_medians.csv')
          hairmedians.index = hairmedians['Unnamed: 0']

          fXXhairsnps[['altXX']][fXXhairsnps['genotypeXX'].isnull()] = hairmedians[fXXhairsnps['genotypeXX'].isnull()]
          fXYhairsnps[['altXY']][fXYhairsnps['genotypeXY'].isnull()] = hairmedians[fXYhairsnps['genotypeXY'].isnull()]
          hairprobabilities = np.round(100*predict_phenotype(allhairsnps, hair_color_model)).astype(int)

          print 'Hair color probabilities', hairprobabilities
          print 'Hair color classes', hair_color_model.classes_

	  return render_template("prediction.html", probabilities = eyeprobabilities, probabilities2 = hairprobabilities)

