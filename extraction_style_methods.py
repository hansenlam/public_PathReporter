# ============================================================================================
    # AGPL v3.0 NOTICE
    # PLEASE INCLUDE AT START OF ALL SOURCE CODE FILES PERTAINING TO THIS PROGRAM

    # PROGRAM TITLE: PathReporter
    # Copyright (C) 2022  
    # Hansen Lam, M.D., Freddy Nguyen, M.D., Ph.D., Aryeh Stock, M.D., Xintong Wang, M.D., Volha Lenskaya, M.D.
    # Alexandros Polydorides, M.D., Qiusheng Si, M.D., John Kim, M.D., Karl Eichorn, Dan Wild
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

from lxml import etree
from lxml.etree import tostring, parse
import re
import pandas as pd
import ujson
import numpy as np
#import report_tabulate
#from report_tabulate import *
from class_XMLDictionary import *
from class_DiagnosisBlock import *
import multiprocessing as mp
from multiprocessing import Queue, Process

interesting_cases_17_18 = [195]
interesting_cases_15_16 = [20, 61, 100, 200, 12, 163, 227, 232, 410, 251]
interesting_cases_14_15 = [232, 230, 396]
interesting_cases_13_14 = [194, 383]

def report_extract_prostate(file):
    t = XMLDictionary(file)#, special_sigs = special_signatures)
    parser = etree.XMLParser()
    tree = etree.parse(file, parser)
    root = tree.getroot()
    t.partitionCases()
    cases = t.getDataFrames(return_item='case_table')
    specimens = t.getDataFrames(return_item='specimen_table')
    merged = cases.merge(specimens, right_on='Case', left_on='patientID_accession_no', how='right')
    return merged


def report_extract(file):
    with open('diagnostic_word_dictionary.json') as dict_file:
        test_dict = ujson.load(dict_file)
    with open('classification_list.json') as class_file:
        test_classes = ujson.load(class_file)

    report_map = []
    t = XMLDictionary(file)
    for case in t.getIdIndices().keys():
        print(case)
        mrn = t.getCaseDict(case)['patientID_med_rec_no']
        dob = t.getCaseDict(case)['patientID_patient_birth_date']
        sex = t.getCaseDict(case)['patientID_patient_sex']
        date = re.sub('  +','',t.getCaseDict(case)['patientID_specimen_received_date'])
        for part, text in t.getCaseDict(case)['fragments'].items():#, item='diags').items():
            f = DiagnosisBlock(case, part, text[0], ' '.join(text[1:]), dict=test_dict)
            print(part, f.diagnostic_summary())
            report_map.append([mrn, dob, sex, date, f.get_specimen()] + f.diagnostic_summary())
    return report_map


def row_create(term_dict, XMLDictOBJ, xml_case):
    mrn = XMLDictOBJ.getCaseDict(xml_case)['patientID_med_rec_no']
    dob = XMLDictOBJ.getCaseDict(xml_case)['patientID_patient_birth_date']
    sex = XMLDictOBJ.getCaseDict(xml_case)['patientID_patient_sex']
    date = re.sub('  +','',XMLDictOBJ.getCaseDict(xml_case)['patientID_specimen_received_date'])
    rows_all = []
    base = [mrn, dob, sex, date]
    for part, text in XMLDictOBJ.getCaseDict(xml_case, item='diags').items():
        f = DiagnosisBlock(xml_case, part, text[0], ' '.join(text[1:]), dict=term_dict)
        for item in f.diagnostic_summary(return_item='entity_feats'):
            row = base + [f.get_specimen()] + item
            rows_all.append(row)
    return rows_all


def report_extract_feats(file):
    with open('diagnostic_word_dictionary.json') as dict_file:
        test_dict = ujson.load(dict_file)
    with open('classification_list.json') as class_file:
        test_classes = ujson.load(class_file)

    report_map = []
    t = XMLDictionary(file)
    for case in t.getIdIndices().keys():
        print(case)
        dict_case = t.getCaseDict(case)
        mrn = dict_case['patientID_med_rec_no']
        dob = dict_case['patientID_patient_birth_date']
        sex = dict_case['patientID_patient_sex']
        date = re.sub('  +',' ',dict_case['patientID_specimen_received_date'])
        base = [mrn, dob, sex, date]
        for part, text in dict_case['fragments'].items():#, item='diags').items():
            f = DiagnosisBlock(case, part, text[0], text[1], dict=test_dict)#' '.join(text[1:]), dict=test_dict)
            print('report_extra_feats_call:',text[1])
            for item in f.diagnostic_summary(return_item='entity_feats'):
                row = base + [f.get_specimen()] + item
                report_map.append(row)
    '''
    PROCESSES = 3
    report_map = []
    t = XMLDictionary(file)
    with mp.Pool(PROCESSES) as pool:
        results = [pool.apply_async(row_create, args=(test_dict, t, case,)) for case in t.getIdIndices().keys()]
        for r in results:
            print(r.get())
            #report_map += r.get()
        pool.close()
        pool.join()
    '''
    return report_map


def known_results_to_frame(extracted_results, *args, **kwargs):
    specimen_type = kwargs.get('specimen_type','default')
    extracted_type = kwargs.get('extracted_type','default')

    def list_extract(list_input):
        term = ''
        if(len(list_input)>0):
            for item in list_input:
                if(list_input.index(item)<(len(list_input)-1)):
                    term = term + item + ', '
                else:
                    term = term + item
        return term

    def key_term_string_interpret(string):
        interpretation = ''
        lsil_list = ([
            'CONDYLOMA',
            'AIN1',
            'AIN I',
            'LOW GRADE'
        ])
        hsil_list = ([
            'HIGH GRADE',
            'AIN2',
            'AIN3',
            'AIN23',
            'AIN II',
            'AIN III',
            'AIN II-III',
            'IN SITU',
            'IN-SITU'
        ])
        scc_list = ([
            'INVASIVE'
        ])
        ain_dict = {}
        for lsil_item in lsil_list:
            ain_dict[lsil_item] = 'LSIL'
        for hsil_item in hsil_list:
            ain_dict[hsil_item] = 'HSIL'
        for scc_item in scc_list:
            ain_dict[scc_item] = 'SCC'
        for key in ain_dict.keys():
            if ((key in string) & ('NEGATIVE' not in string) & ('NO' not in string)):
                interpretation = ain_dict[key]
        return interpretation

    def definite_string_interpret(string):
        interpretation = None
        lsil_list = ([
            'CONDYLOMA',
            'AIN1',
            'AIN I',
            'LOW GRADE'
        ])
        hsil_list = ([
            'HIGH GRADE',
            'AIN2',
            'AIN3',
            'AIN23',
            'AIN II',
            'AIN III',
            'AIN II-III',
            'IN-SITU',
            'IN SITU'
        ])
        scc_list = ([
            'INVASIVE'
        ])
        carcinoma_dict = ({
            'SQUAMOUS':'SCC',
            'ADENO':'ADENOCA'
        })
        ain_dict = {}
        for lsil_item in lsil_list:
            ain_dict[lsil_item] = 'LOW GRADE'
        for hsil_item in hsil_list:
            ain_dict[hsil_item] = 'HIGH GRADE'
        for scc_item in scc_list:
            ain_dict[scc_item] = 'SCC'

        for key in ain_dict.keys():
            if key in string:
                interpretation = ain_dict[key]
        if interpretation is None:
            if any(re.findall('CARCINOMA',string,re.IGNORECASE|re.MULTILINE)):
                for key in carcinoma_dict.keys():
                    if any(re.findall(key,string,re.IGNORECASE|re.MULTILINE)):
                        interpretation = carcinoma_dict[key]
            else:
                interpretation = ''
        return interpretation

    def dysplasia_present(string, *args, **kwargs):
        interpretation = 'NO COMMENT'
        input_type = kwargs.get('input_type',None)
        if(type(string) is list):
            input_type = 'list'
        if(input_type is None):
            if(('DYSPLASIA' in string)&('NEGATIVE' not in string)):
                interpretation = 'PRESENT'
            elif('NEGATIVE DYSPLASIA' in string):
                interpretation = 'NOT PRESENT'
            elif('INDEFINITE DYSPLASIA'in string):
                interpretation = 'INDEFINITE'
            else:
                interpretation = 'NO COMMENT'
        elif(input_type == 'list'):
            if(len(string)>0):
                interpretation_list = []
                for elem in string:
                    if(('DYSPLASIA' in elem)&('NEGATIVE' not in elem)):
                        elem_interpretation = 'PRESENT'
                        interpretation = 'PRESENT'
                    elif('NEGATIVE DYSPLASIA' in elem):
                        elem_interpretation = 'NOT PRESENT'
                    else:
                        elem_interpretation = 'NO COMMENT'
                    interpretation_list.append(elem_interpretation)
            else:
                interpretation = 'NO COMMENT'
        return interpretation

    def dysplasia_grade(string, *args, **kwargs):
        interpretation = 'NO COMMENT'
        input_type = kwargs.get('input_type',None)
        if(type(string) is list):
            input_type = 'list'
        if(input_type is None):
            if(any(re.findall(r'LOW\-*',string,re.MULTILINE|re.IGNORECASE))):
                interpretation = 'LOW GRADE'
            elif(any(re.findall(r'HIGH\-*',string,re.MULTILINE|re.IGNORECASE))):
                interpretation = 'HIGH GRADE'
            elif(any(re.findall(r'ADENOCARCINOMA',string,re.MULTILINE|re.IGNORECASE))):
                interpretation = 'ADENOCARCINOMA'
            else:
                interpretation = 'NO COMMENT'
        elif(input_type == 'list'):
            if(len(string)>0):
                interpretation_list = []
                for elem in string:
                    if(any(re.findall(r'LOW\-*',elem,re.MULTILINE|re.IGNORECASE))):
                        elem_interpretation = 'LOW GRADE'
                        if(any(re.findall(r'NEG',elem,re.MULTILINE|re.IGNORECASE))):
                            elem_interpretation = 'NEGATIVE'
                    elif(any(re.findall(r'HIGH\-*',elem,re.MULTILINE|re.IGNORECASE))):
                        elem_interpretation = 'HIGH GRADE'
                        if(any(re.findall(r'NEG',elem,re.MULTILINE|re.IGNORECASE))):
                            elem_interpretation = 'NEGATIVE'
                    elif(any(re.findall(r'IND',elem,re.MULTILINE|re.IGNORECASE))):
                        elem_interpretation = 'INDEFINITE'
                    else:
                        elem_interpretation = 'NO COMMENT'
                    interpretation_list.append(elem_interpretation)
                for elem in interpretation_list:
                    if(elem=='INDEFINITE'):
                        interpretation = elem
                    if(elem=='LOW GRADE'):
                        interpretation = elem
                    elif(elem=='HIGH GRADE'):
                        interpretation = elem
            else:
                interpretation = 'NO COMMENT'
        return interpretation

    if(extracted_type=='entities'):
        cols = ([
            'MRN', 'DOB', 'SEX', 'DATE', 'SPECIMEN',
            'CASE', 'PART', 'TEXT', 'ENTITIES', 'KEY_TERMS',
            'HEDGE_ENTITIES', 'DEFINITE_ENTITIES', 'NEGATIVE_ENTITIES',
            'SCORE_LIST', 'AVG_SCORE', 'CATEGORY'
        ])
        df_test = pd.DataFrame(extracted_results, columns=cols)
        df_test_known = df_test
        df_test_known['KEY_TERMS_STRING'] = df_test_known['KEY_TERMS'].map(list_extract)
        df_test_known['DEFINITE_TERMS_STRING'] = df_test_known['DEFINITE_ENTITIES'].map(list_extract)
        if specimen_type == 'ain':
            #df_test_known['INTERPRETATION'] = df_test_known['KEY_TERMS_STRING'].map(key_term_string_interpret)
            df_test_known['INTERPRETATION'] = df_test_known['DEFINITE_TERMS_STRING'].map(definite_string_interpret)
            df_test_known['INTERPRETATION'] = (
                np.where(
                    (df_test_known['INTERPRETATION'].isna())&(df_test_known['TEXT'].str.contains('SQUAMOUS')),
                    'SCC',
                    df_test_known['INTERPRETATION']
                )
            )
            return df_test_known
        elif(specimen_type == 'default'):
            return df_test_known
    elif(extracted_type=='feats'):
        cols = ([
            'MRN', 'DOB', 'SEX', 'DATE', 'SPECIMEN', 'CASE', 'PART', 'TEXT', 'ENTITY', 'FEATURES'
        ])
        df_feat = pd.DataFrame(extracted_results, columns = cols)
        if specimen_type == 'ibd':
            df_feat['FEAT_STRING'] = df_feat['FEATURES'].map(list_extract)
            df_feat['DYSPLASIA'] = df_feat['FEATURES'].map(dysplasia_present)
            df_feat['GRADE'] = df_feat['FEATURES'].map(dysplasia_grade)
            return df_feat
        elif(specimen_type=='default'):
            return df_feat
    else:
        return 'extracted_type must be specified'


#determine format of prostate reports
def prostate_format_class(text):
    format = 'free'
    standard_format_regex = r'PROSTATE CANCER GRADING|TUMOR QUANTIFICATION'
    if any(re.findall(standard_format_regex, text, re.IGNORECASE|re.MULTILINE)):
        format = 'std'
    return format


# Standard format extractor
def prostate_standard_format_extract(text, section):
    dict_diagnostic_terms = {}
    #dict_diagnostic_terms['CASE'] = case
    #dict_diagnostic_terms['PART'] = part
    standard_format_regex = r'PROSTATE CANCER GRADING|TUMOR QUANTIFICATION'
    if any(re.findall(standard_format_regex, text, re.IGNORECASE|re.MULTILINE)):
        text_split = text.splitlines()
        for line in text_split:
            line = line.strip()
            if('ADENOCARCINOMA:' in line):
                dict_diagnostic_terms['TYPE'] = 'ADENOCARCINOMA'
            elif(len(line.split(':'))>1):
                if(len(line.split(':')[1])>0):
                    if('GRADE GROUP' in line):
                        if(any(re.findall(r'[0-9]', line))):
                            grade_group = re.findall(r'[0-9]', line)[0]
                        else:
                            grade_group = None
                        dict_diagnostic_terms[line.split(':')[0]] = grade_group
                    else:
                        line = re.sub(r' IS:',':',line)
                        line = re.sub(r'^THE ','',line)
                        if(any(re.findall(r'[0-9]\.*[0-9]*', line))):
                            number = re.findall(r'[0-9]\.*[0-9]*', line)[0]
                        else:
                            number = None
                        dict_diagnostic_terms[line.split(':')[0]] = number
    if('GLEASON SCORE' in dict_diagnostic_terms.keys()):
        dict_diagnostic_terms['TOTAL GLEASON SCORE'] = dict_diagnostic_terms.pop('GLEASON SCORE')
    if(section in dict_diagnostic_terms.keys()):
        return dict_diagnostic_terms[section]
    else:
        return None


#other formats
def prostate_free_format_extract(text, section):
    dict_diagnostic_terms = {}
    #dict_diagnostic_terms['CASE'] = case
    #dict_diagnostic_terms['PART'] = part
    dict_gleason = ({
        ('3', '3') : '1',
        ('3', '4') : '2',
        ('4', '3') : '3',
        ('3', '5') : '4',
        ('4', '4') : '4',
        ('5', '3') : '4',
        ('4', '5') : '5',
        ('5', '4') : '5',
        ('5', '5') : '5'
    })
    percentage = None
    percent_high_grade = None
    text = re.sub(r'\nMM',' MM',text)
    text_split = text.splitlines()
    for line in text_split:
        carc_split = line.split(r'[\.\,]')
        for carc_line in carc_split:
            if(any(re.findall('ADENOCARCINOMA',carc_line,re.IGNORECASE))&(not (any(re.findall(r'IMMUNOSTAIN|PIN4|INDEFINITE|NEGATIVE',carc_line,re.IGNORECASE))))):
                dict_diagnostic_terms['TYPE'] = 'ADENOCARCINOMA'
        if('+' in line):
            text_copy = line
            if(any(re.findall(r'([0-9]\s*\+\s*[0-9])', text_copy))):
                score_tuple = tuple(re.findall(r'([0-9]\s*\+\s*[0-9])', text_copy)[0].replace(' ','').split('+'))
            else:
                score_tuple = (0,0)
            dict_diagnostic_terms['PRIMARY GLEASON GRADE'] = score_tuple[0]
            dict_diagnostic_terms['SECONDARY GLEASON GRADE'] = score_tuple[1]
            dict_diagnostic_terms['TOTAL GLEASON SCORE'] = str(int(dict_diagnostic_terms['PRIMARY GLEASON GRADE']) + int(dict_diagnostic_terms['SECONDARY GLEASON GRADE']))
            if(score_tuple in dict_gleason.keys()):
                dict_diagnostic_terms['GRADE GROUP'] = dict_gleason[score_tuple]
            else:
                dict_diagnostic_terms['GRADE GROUP'] = 'SCORE ERROR ' + str(score_tuple[0]) + ',' + str(score_tuple[1])
        if(any(re.findall(r'\bCORE\b', line))|any(re.findall(r'\bCORES\b', line))|any(re.findall(r'CORES[/]CORE', line))|any(re.findall(r'TISSUE', line))):
            #print(line)
            dict_words_to_num = ({
                'ONE': '1',
                'TWO': '2',
                'THREE': '3',
                'FOUR': '4',
                'FIVE': '5',
                '1':'1',
                '2':'2',
                '3':'3',
                '4':'4',
                '5':'5'
            })
            cores = []
            quant_split = line.split(',')
            for quant_line in quant_split:
                quant_line = re.sub(r'TISSUE CORE', 'CORE', quant_line, re.MULTILINE)
                quant_line = re.sub(r'CORES[/]CORE', 'CORE', quant_line, re.MULTILINE)
                if(any(re.findall(r'GRADE|PATTERN',quant_line,re.MULTILINE))):
                    if(not (any(re.findall(r'TISSUE|INVOLV',quant_line,re.MULTILINE)))):
                        if(any(re.findall(r'(?<![(])([0-9]+)\s*\%(?![)])',quant_line,re.MULTILINE))):
                            percent_high_grade = re.findall(r'(?<![(])([0-9]+)\s*\%(?![)])',quant_line)[0]
                            dict_diagnostic_terms['% OF GLEASON GRADE 4 AND/OR 5'] = percent_high_grade
                if(any(re.findall(r'(\w+)[\s](?=CORE)',quant_line))|any(re.findall(r'(\w+)[\s](?=CORES\.*)',quant_line))):
                    if((percentage is None)*any(re.findall(r'\%',quant_line))*any(re.findall(r'TISSUE|INVOLV',quant_line))):
                        if(not (any(re.findall(r'GRADE|PATTERN',quant_line)))):
                            if(any(re.findall(r'([0-9]+)(?=%)',quant_line))):
                                percentage = re.findall(r'([0-9]+)(?=%)',quant_line)[0]
                            else:
                                precentage = '% ERROR'
                            dict_diagnostic_terms['PERCENTAGE OF TISSUE WITH CARCINOMA'] = percentage
                    if((re.findall(r'(\w+)[\s](?=CORE)',quant_line)[0] in ['THE','SINGLE'])|('% OF ONE CORE SUBMITTED' in quant_line)):
                        dict_diagnostic_terms['TOTAL # OF CORES IDENTIFIED'] = '1'
                        dict_diagnostic_terms['TOTAL # OF CORES WITH CARCINOMA'] = '1'
                    else:
                        if(any(re.findall(r'(\w+\sOF\s\w+)',quant_line))):
                            cores = re.findall(r'(\w+\sOF\s\w+)',quant_line)[0].split(' OF ')
                            if(len(cores)==2):
                                if((cores[0] in dict_words_to_num.keys())&(cores[1] in dict_words_to_num.keys())):
                                    cores_mapped = [core for core in map(dict_words_to_num.get, cores, cores)]
                                    dict_diagnostic_terms['TOTAL # OF CORES IDENTIFIED'] = cores_mapped[1]
                                    dict_diagnostic_terms['TOTAL # OF CORES WITH CARCINOMA'] = cores_mapped[0]
                        elif(any(re.findall(r'([0-9][/][0-9])',quant_line))):
                            cores = re.findall(r'([0-9][/][0-9])',quant_line)[0].split('/')
                            if(len(cores)==2):
                                if((cores[0] in dict_words_to_num.keys())&(cores[1] in dict_words_to_num.keys())):
                                    cores_mapped = [core for core in map(dict_words_to_num.get, cores, cores)]
                                    dict_diagnostic_terms['TOTAL # OF CORES IDENTIFIED'] = cores_mapped[1]
                                    dict_diagnostic_terms['TOTAL # OF CORES WITH CARCINOMA'] = cores_mapped[0]
                elif(any(re.findall(r'TISSUE',quant_line))):
                    if(any(re.findall(r'(\w+\sOF\s\w+)',quant_line))):
                        cores = re.findall(r'(\w+\sOF\s\w+)',quant_line)[0].split(' OF ')
                        if(len(cores)==2):
                            if((cores[0] in dict_words_to_num.keys())&(cores[1] in dict_words_to_num.keys())):
                                cores_mapped = [core for core in map(dict_words_to_num.get, cores, cores)]
                                dict_diagnostic_terms['TOTAL # OF CORES IDENTIFIED'] = cores_mapped[0]
                                dict_diagnostic_terms['TOTAL # OF CORES WITH CARCINOMA'] = cores_mapped[1]

        if((percentage is None)*any(re.findall(r'\%',line))*any(re.findall(r'TISSUE|INVOLV',line))):
            if(not (any(re.findall(r'GRADE|PATTERN',line)))):
                if(any(re.findall(r'([0-9]+)(?=%)',line))):
                    percentage = re.findall(r'([0-9]+)(?=%)',line)[0]
                else:
                    precentage = '% ERROR'
                dict_diagnostic_terms['PERCENTAGE OF TISSUE WITH CARCINOMA'] = percentage
        if(any(re.findall(r'MM',line,re.MULTILINE))*any(re.findall(r'(?<!BIOPSY )(LINEAR|LENGTH|MEASURING)',line,re.MULTILINE))):
            if(any(re.findall(r'\%',quant_line))*any(re.findall(r'TISSUE|INVOLV',quant_line))):
                if(not (any(re.findall(r'GRADE|PATTERN',quant_line)))):
                    if(any(re.findall(r'([0-9]+)(?=%)',quant_line))):
                        percentage = re.findall(r'([0-9]+)(?=%)',quant_line)[0]
                    else:
                        precentage = '% ERROR'
                    dict_diagnostic_terms['PERCENTAGE OF TISSUE WITH CARCINOMA'] = percentage
            linear_mark = None
            milli_mark = None
            millimeters = None
            linear_markers = [iter.start() for iter in re.finditer(r'(?<!BIOPSY )(LINEAR|LENGTH|MEASURING)',line,re.MULTILINE)]
            millimeter_list = []
            if(len(linear_markers)>0):
                linear_mark = linear_markers[0]
            milli_markers = [iter.start() for iter in re.finditer(r'(?<=[0-9])\s*(MM)\,*',line,re.MULTILINE)]
            if(len(milli_markers)>0):
                milli_mark = milli_markers[0]
            if((linear_mark is not None)&(milli_mark is not None)):
                if((linear_mark<milli_mark)):
                    millimeter_list = re.findall(r'[0-9]+\.*[0-9]*',line[linear_mark:milli_mark],re.MULTILINE)
                    #print(re.findall(r'[0-9]+\.*[0-9]*',line[linear_mark:milli_mark],re.MULTILINE))
                    #print(linear_mark, milli_mark, text[linear_mark:milli_mark])
                else:
                    millimeter_list = [re.findall(r'[0-9]+\.*[0-9]*(?=\s*MM\,*)',line,re.MULTILINE)[0]]
                    #print([re.findall(r'[0-9]+\.*[0-9]*(?=\s*MM)',line,re.MULTILINE)[0]])
            #else:
                #print(linear_mark, milli_mark)
            if(len(millimeter_list)>0):
                millimeters = 0
                for milli in millimeter_list:
                    millimeters += float(milli)
                millimeters = str(millimeters)
            dict_diagnostic_terms['LINEAR AMOUNT OF TISSUE WITH CARCINOMA'] = millimeters
        if((percent_high_grade is None)*any(re.findall(r'GRADE|PATTERN',line,re.MULTILINE))):
            #if(not (any(re.findall(r'TISSUE|INVOLV',line)))):
            quant_split = re.split(r'[,.\n]',line)
            for quant_line in quant_split:
                if(not (any(re.findall(r'INVOLV|OCCUPY',quant_line)))):
                #print(quant_line)
                #if(any(re.findall(r'%',quant_line))):
                    if(any(re.findall(r'(?<![(])([0-9]+)\s*\%(?![)])',quant_line))):
                        #percent_high_grade = re.findall(r'([0-9]+)\s*\%',quant_line)[0]
                        percent_high_grade = re.findall(r'(?<![(])([0-9]+)\s*\%(?![)])',quant_line)[0]
                        dict_diagnostic_terms['% OF GLEASON GRADE 4 AND/OR 5'] = percent_high_grade
                    #else:
                        #dict_diagnostic_terms['% OF GLEASON GRADE 4 AND/OR 5'] = '% ERROR' + ' ' + quant_line
    if(section in dict_diagnostic_terms.keys()):
        return dict_diagnostic_terms[section]
    else:
        return None
