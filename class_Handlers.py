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

import re
from collections import defaultdict
from case_dictionary_methods import *

class XML_to_Dict_Handler():
    def __init__(self, test_bx):
        self.test_bx = test_bx
        self.dict_sections = {}
        section = None
        self.consult = False
        text = []
        self.diagnosis_markers = []
        self.addendum_markers = []
        self.case_diagnosis = ['']
        self.case_addendums = ['']
        self.consult_cases = []
        self.dict_part_specimen = {}
        for key in self.test_bx.keys():
            if('name' in key):
                if(any(re.findall('DIAGNOSIS', self.test_bx[key], re.IGNORECASE))):
                    self.diagnosis_markers += [self.test_bx[key]]
                if(('RS' in self.test_bx['patientID_accession_no'])|('BS' in self.test_bx['patientID_accession_no'])):
                    if(any(re.findall('COMMENT', self.test_bx[key], re.IGNORECASE))):
                        self.diagnosis_markers += [self.test_bx[key]]
                        #print(diagnosis_markers)
                if('RP' in self.test_bx['patientID_accession_no']):
                    if(any(re.findall('NOTE', self.test_bx[key], re.IGNORECASE))):
                        #self.diagnosis_markers += [self.test_bx[key]]
                        self.addendum_markers += [self.test_bx[key]]
                if(any(re.findall('ADDENDUM', self.test_bx[key], re.IGNORECASE))):
                    self.addendum_markers += [self.test_bx[key]]
                if(any(re.findall('MATERIALS', self.test_bx[key], re.IGNORECASE))):
                    self.consult = True
                section = self.test_bx[key]
                text=[]
            elif('finding' in key):
                text += [self.test_bx[key]]
                #print(text)
            elif('description' in key):
                if(any(re.findall('SUBMITTED SLIDE', self.test_bx[key], re.IGNORECASE))|any(re.findall('CONSULT', self.test_bx[key], re.IGNORECASE))):
                    self.consult = True
                section = 'description'
                text += [self.test_bx[key].strip()]
            if(section is not None):
                self.dict_sections[section] = text
        #print(self.dict_sections)
        if self.consult:
            if ('MATERIALS RECEIVED:' in self.dict_sections):
                material_clean_regex = ({
                    r'(?<=[A-Z])(-)(?=[0-9])':'',
                    r'(?<=[0-9])(-)(?=[0-9])':'',
                    r'[;,]':' '
                })
                materials = self.dict_sections['MATERIALS RECEIVED:'][0]
                for regex, sub in material_clean_regex.items():
                    materials = re.sub(regex, sub, materials)
                self.consult_cases = re.findall(r'([A-Z]+[0-9]+)',materials,flags=re.IGNORECASE|re.MULTILINE)
        self.part_names = self.dict_sections['description']
        for diag in self.diagnosis_markers:
            self.case_diagnosis[0] += self.dict_sections[diag][0]
        for addendums in self.addendum_markers:
            self.case_addendums[0] += self.dict_sections[addendums][0]
    def getdiagmarkers(self):
        return self.diagnosis_markers
    def getaddmarkers(self):
        return self.addendum_markers
    def getconsult(self):
        return self.consult
    def getconsultcases(self):
        return self.consult_cases
    def getPartNames(self):
        return self.part_names
    def getdict(self):
        return self.dict_sections
    def getCaseDiagnosis(self):
        return self.case_diagnosis
    def getCaseAddendums(self):
        return self.case_addendums
    def setDictPartSpecimen(self, dict):
        self.dict_part_specimen = dict.copy()
    def getDictPartSpecimen(self):
        return self.dict_part_specimen
    def setCaseDiagnosis(self, text_list):
        self.case_diagnosis = text_list
    def setCaseAddendums(self, text_list):
        self.case_addendums = text_list


class TextBlock_Handler():
    def __init__(self, textblock, part_names, diagosis_text_reference, consult, consult_cases):
        self.key = textblock
        self.is_diagnosis = (textblock == diagosis_text_reference)
        self.part_names = part_names
        self.consult = consult
        self.consult_cases = consult_cases

        self.part_list = []
        self.grouped_letter_index_tuples = []
        self.dict_part_specimen = {}

        self.dict_parts_text = {}
        self.synoptic_found = [False]
        self.synoptic_text = []

        new_key = textblock_deconstruct(self.key, self.part_names, self.part_list, self.grouped_letter_index_tuples, self.consult, self.consult_cases, self.is_diagnosis, self.synoptic_found, self.synoptic_text)
        #self.dict_part_specimen = dict(zip(self.part_list, self.part_names))
        #for part in self.part_list:
        #    self.dict_parts_text[part] = []
        self.key = new_key
        for part in self.part_list:
            self.dict_parts_text[part] = []
        textblock_to_dict(self.key, self.part_list, self.grouped_letter_index_tuples, self.dict_parts_text)
        print(self.synoptic_found, self.synoptic_text)

    def getDictPartsText(self):
        return self.dict_parts_text
    def getDictPartSpecimen(self):
        return self.dict_part_specimen
    def hasSynoptic(self):
        return self.synoptic_found
    def getSynoptic(self):
        return self.synoptic_text
    def getIsDiagnosis(self):
        return self.is_diagnosis
    def getTextIndex(self,*args,**kwargs):
        item = kwargs.get('var', 'default')
        dict_return_items = ({
            'key': self.key,
            'part_names':self.part_names,
            'part_list':self.part_list,
            'grouped_letter_index_tuples':self.grouped_letter_index_tuples,
        })
        if item == 'default':
            return dict_return_items
        elif item in dict_return_items:
            return dict_return_items[item]
        else:
            return 0


def case_partition_oop(root, case_num):
    test_bx={}
    case = case_num
    test_bx = walker_bx(root[case], test_bx)
    xml_dict_handle = XML_to_Dict_Handler(test_bx)
    diagnosis_markers = xml_dict_handle.getdiagmarkers()
    addendum_markers = xml_dict_handle.getaddmarkers()
    consult = xml_dict_handle.getconsult()
    consult_cases = xml_dict_handle.getconsultcases()
    dict_sections = xml_dict_handle.getdict()
    part_names = xml_dict_handle.getPartNames()
    diag_add_blocks = ([
        xml_dict_handle.getCaseDiagnosis()[0],
        xml_dict_handle.getCaseAddendums()[0]
    ])
    part_list = []
    diag_block_deident = [textblock_deidentify(xml_dict_handle.getCaseDiagnosis()[0])]
    add_block_deident = [textblock_deidentify(xml_dict_handle.getCaseAddendums()[0])]
    xml_dict_handle.setCaseDiagnosis(diag_block_deident)
    xml_dict_handle.setCaseAddendums(add_block_deident)

    dict_diag_add_fragment = ({
        xml_dict_handle.getCaseDiagnosis()[0]:{},
        xml_dict_handle.getCaseAddendums()[0]:{}
    })
    dict_part_specimen = {}
    synoptic_summary = {}
    synoptic_text = ''
    synoptic = [False]
    for key in dict_diag_add_fragment.keys():
        #print(key)
        diagnosis = (key == xml_dict_handle.getCaseDiagnosis()[0])
        #print(diagnosis)
        if(len(key)>0):
            TextBlockHandle = TextBlock_Handler(key, part_names, xml_dict_handle.getCaseDiagnosis()[0], consult, consult_cases)
            synoptic = TextBlockHandle.hasSynoptic()
            if synoptic[0]:
                synoptic_text = TextBlockHandle.getSynoptic()
                test_bx['synoptic'] = synoptic_text
            dict_diag_add_fragment[key] = TextBlockHandle.getDictPartsText()
            #print(dict_diag_add_fragment[key])
            #if(TextBlockHandle.getIsDiagnosis()):
                #xml_dict_handle.setDictPartSpecimen(TextBlockHandle.getDictPartSpecimen())
                #print(xml_dict_handle.getDictPartSpecimen())
        #dict_part_specimen = xml_dict_handle.getDictPartSpecimen()
            for letter in dict_diag_add_fragment[key]:
                if(letter not in part_list):
                    part_list.append(letter)
                    #print(part_list)
            part_list.sort()
    list_name_diff = []
    if(len(part_list)>1):
        if(len(part_list)>len(part_names)):
            list_name_diff = part_list[len(part_names):]
            del part_list[len(part_names):]

    dict_part_specimen = dict(zip(part_list, part_names))
    dict_diag_with_add = defaultdict(list)

    for k in dict_part_specimen.keys():
        dict_diag_with_add[k] = [dict_part_specimen[k]]
    for key in dict_diag_add_fragment.keys():
        for k, v in dict_diag_add_fragment[key].items():
            if(k in list_name_diff):
                dict_diag_with_add['A'] += v
            else:
                dict_diag_with_add[k] += v
            #dict_diag_with_add[k].append(v)
    dict_diag_with_add = dict(dict_diag_with_add)
    test_bx['fragments'] = dict_diag_with_add
    print(test_bx['fragments'])
    #print(test_bx['fragments'])
    return test_bx
