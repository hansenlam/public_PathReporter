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

from lxml import etree
from class_Handlers import *
from xml_methods import *
from case_dictionary_methods import *

# CREATING XML to DICTIONARY CLASS:_____________________________________________
class XMLDictionary():
    #CHANGE CASE_PARTITION_NEW BACK TO CASE_PARTITION
    def __init__(self, file_name, *args, **kwargs):
        self.parser = etree.XMLParser()
        self.tree = etree.parse(file_name, self.parser)
        self.root = self.tree.getroot()
        range_list = list(kwargs.get('range', range(len(self.root)-1)))
        self.interfering_signatures = kwargs.get('special_sigs',[])
        self.case_ids_index_dict = rootID_Indices(self.root, range_list)
        self.case_nums = list(self.case_ids_index_dict.values())
        self.case_index_ids_dict = dict(zip(self.case_nums,self.case_ids_index_dict.items()))
    def partitionCases(self):
        self.case_bxs = {}
        for case in self.case_nums:
            print(self.case_index_ids_dict[case])
            self.case_bxs[case] = case_partition_oop(self.root, case) #, special_sigs=self.interfering_signatures)
            #self.case_bxs[case] = case_partition_test(self.root, case)
        self.dict_cases_diag = {}
        for index in self.case_bxs.keys():
            self.case_bxs_index = self.case_bxs[index]
            self.accession_no = self.case_bxs_index['patientID_accession_no']
            self.diag_phrases = self.case_bxs_index['fragments']
            self.dict_cases_diag[self.accession_no] = self.diag_phrases
    def getDataFrames(self, return_item):
        # return_item can either be 'case_table' or 'specimen_table'
        for case in self.case_nums:
            self.case_bxs[case] = case_partition_oop(self.root, case)#, special_sigs=self.interfering_signatures)
            #self.case_bxs[case] = case_partition_test(self.root, case)
        self.df = partition_dict_to_df(self.case_bxs, return_item=return_item)
        return self.df
    def getCaseDict(self, id, *args, **kwargs):
        return_item = kwargs.get('item', 'case')
        selected_case = self.case_ids_index_dict[id]
        partition = case_partition_oop(self.root, selected_case)#, special_sigs=self.interfering_signatures)
        #partition = case_partition_test(self.root, selected_case)
        if return_item == 'case':
            return partition
        if return_item == 'diags':
            return partition['fragments']
    def FindCaseIndex(self, id):
        return self.case_ids_index_dict[id]
    def getCaseIdList(self):
        return list(self.case_ids_index_dict.keys())
    def getIdIndices(self):
        return self.case_ids_index_dict
    def getCaseList(self):
        return self.case_nums
    def getXMLDict(self):
        self.partitionCases()
        return self.case_bxs
    def getDiagDict(self):
        self.partitionCases()
        return self.dict_cases_diag
