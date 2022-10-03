# ============================================================================================
    # AGPL v3.0 NOTICE
    # PLEASE INCLUDE AT START OF ALL SOURCE CODE FILES PERTAINING TO THIS PROGRAM

    # PROGRAM TITLE: PathReporter
    # Copyright (C) 2022  Hansen Lam, M.D., Freddy Nguyen, M.D., Ph.D., Karl Eichorn, Dan Wild
    # Disclaimer: This software is intended for research purposes only. This software is not
    # for clinical use. The authors and contributors shall not be liable for any damages or
    # losses resulting from the use of this software. The authors and contributors shall not 
    # be liable for any actions or decisions resulting rom the use of this software.

    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU Affero General Public License as published
    # by the Free Software Foundation, either version 3 of the License, or
    # any later version.

    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU Affero General Public License for more details.

    # You should have received a copy of the GNU Affero General Public License
    # along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ============================================================================================

import os
import glob
import multiprocessing as mp
from extraction_style_methods import *
from multiprocessing import Queue, Process
import pandas as pd
import numpy as np
import re
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

# SCENARIO 1
result_all = []
xml_file_list = [<xml file path>]
for xml_file in xml_file_list:
    results = report_extract(xml_file)
    result_all = result_all + results
df = known_results_to_frame(result_all, extracted_type='entities',specimen_type='ain')
df.info()
df = df.drop(['ENTITIES','KEY_TERMS','HEDGE_ENTITIES','DEFINITE_ENTITIES','NEGATIVE_ENTITIES','SCORE_LIST','AVG_SCORE','CATEGORY','KEY_TERMS_STRING','DEFINITE_TERMS_STRING','INTERPRETATION'],axis=1)
df.to_excel('SCENARIO1_EXTRACTED_new.xlsx')
#t = XMLDictionary(<xml file path>)

# SCENARIO 2
result_all = []
xml_file_list = [<xml file path>]
for xml_file in xml_file_list:
    results = report_extract_feats(xml_file)
    result_all = result_all + results
df_feat = known_results_to_frame(result_all, extracted_type='feats',specimen_type='ibd')
df_feat
df_feat_deduplicated = (
    df_feat.drop_duplicates(subset=['MRN','CASE','PART'],keep='last')
)
df_feat_deduplicated[df_feat_deduplicated['MRN'].isin(df_feat[df_feat['ENTITY']=='ADENOCARCINOMA']['MRN'])]
df_feat_deduplicated['DYSPLASIA'] = (
    np.where(
        df_feat_deduplicated['ENTITY']=='ADENOCARCINOMA',
        'PRESENT',
        df_feat_deduplicated['DYSPLASIA']
    )
)
df_feat_deduplicated['GRADE'] = (
    np.where(
        df_feat_deduplicated['ENTITY']=='ADENOCARCINOMA',
        df_feat_deduplicated['ENTITY'],
        df_feat_deduplicated['GRADE']
    )
)

def specimen_bin(string):
    s = string
    s = re.sub(r'X\d','',s,re.IGNORECASE|re.MULTILINE)
    s = re.sub(r'#\d','',s,re.IGNORECASE|re.MULTILINE)
    digit_regex = ([
        r'\d+(?=\s*CM)',
        r'(?<=@)\s*\d+',
        r'(?<=COLON\s)\d+'
    ])
    dict_locations = ({
        'RECTUM/SIGMOID':[r'RECT|SIG',20],
        'DESCENDING/SF':[r'DESC|LEFT|SF|SPLEN',60],
        'TRANSVERSE':[r'TRAN|TC',100],
        'ASCENDING/HF':[r'ASC|RIGHT|AC|HEP|HF',140],
        'CECUM':[r'ILEO|CEC|VALV',150]
    })

    regex_found = []
    loc_interp = s
    for reg in digit_regex:
        if(len(regex_found)==0):
            if(any(re.findall(reg,s,re.IGNORECASE|re.MULTILINE))):
                if(re.findall(reg,s,re.IGNORECASE|re.MULTILINE)[0].strip().isnumeric()):
                    num = int(re.findall(reg,s,re.IGNORECASE|re.MULTILINE)[0].strip())
                    num_loc = None
                    for location, reg in dict_locations.items():
                        if((num_loc is None) & (num<=reg[1])):
                            num_loc = location
                            regex_found.append(num_loc)
    if(len(regex_found)==0):
        for location, reg in dict_locations.items():
            if(len(regex_found)==0):
                if(any(re.findall(reg[0],s,re.IGNORECASE|re.MULTILINE))):
                    regex_found.append(location)
    if(len(regex_found)>0):
        loc_interp = regex_found[0].strip()
    else:
        loc_interp = 'NON-COLON/OTHER'
    return loc_interp

df_feat_deduplicated['SPECIMEN'].head(60)
df_feat_deduplicated['SPECIMEN_BIN'] = (
    df_feat_deduplicated['SPECIMEN'].apply(lambda x: specimen_bin(x))
)
df_feat_deduplicated

df_feat_deduplicated.to_excel('SCENARIO2_EXTRACTED_binned_new.xlsx')
df_feat.info()


# SCENARIO 3
xml_file_list = [<xml file path>]
prostate_list = ([
    'TYPE', 'PRIMARY GLEASON GRADE', 'SECONDARY GLEASON GRADE', 'TOTAL GLEASON SCORE',
    'GRADE GROUP', 'TOTAL # OF CORES IDENTIFIED', 'TOTAL # OF CORES WITH CARCINOMA',
    'PERCENTAGE OF TISSUE WITH CARCINOMA', 'LINEAR AMOUNT OF TISSUE WITH CARCINOMA',
    '% OF GLEASON GRADE 4 AND/OR 5'
])
prostate_all = pd.DataFrame()
for xml_file in xml_file_list:
    prostate = report_extract_prostate(xml_file)
    prostate_all = pd.concat([prostate_all, prostate], ignore_index=True)

prostate_all['DIAG_FORMAT'] = prostate_all['Diagnosis'].map(prostate_format_class)
for col in prostate_list:
    prostate_all[col] = (
        np.where(
            prostate_all['DIAG_FORMAT'] == 'free',
            prostate_all['Diagnosis'].apply(lambda x: prostate_free_format_extract(x, col)),
            None
        )
    )
for col in prostate_list:
    prostate_all[col] = (
        np.where(
            prostate_all['DIAG_FORMAT'] == 'std',
            prostate_all['Diagnosis'].apply(lambda x: prostate_standard_format_extract(x, col)),
            prostate_all[col]
        )
    )
prostate_all.to_excel('SCENARIO3_EXTRACTED_11-24.xlsx')
t = XMLDictionary(<xml file path>)
