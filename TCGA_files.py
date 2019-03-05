import requests as rq
import json
import pandas as pd
import numpy as np


cases_endpt = 'https://api.gdc.cancer.gov/cases'

fields = [
    "submitter_id",
    "case_id",
    "primary_site",
    "disease_type",
    "diagnoses.tumor_stage",
    "diagnoses.tumor_grade",
    "diagnoses.primary_diagnosis",
    "diagnoses.classification_of_tumor",
    "annotations.classification",
    "samples.tumor_code"
    ]

fields = ','.join(fields)

def queryFile(idFile):
    filters = {
    "op": "in",
    "content":{
        "field": "files.file_name",
        "value": [idFile]
        }
    }
    params = {
    "fields": fields,
    "filters": json.dumps(filters),
    "format": "TSV",
    "size": "1"
    }
    print("quering...%s"%idFile)
    response = rq.get(cases_endpt, params = params)
    #print(response.content.decode('utf-8'))
    r = response.content.decode("utf-8").split('\r')
    data = np.array(r[1].replace('\n','').split('\t'))
    data = data.reshape(1,len(data))
    return pd.DataFrame(data=data, columns=r[0].split('\t'), index=[0])
