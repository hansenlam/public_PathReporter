# ============================================================================================
    # AGPL v3.0 NOTICE
    # PLEASE INCLUDE AT START OF ALL SOURCE CODE FILES PERTAINING TO THIS PROGRAM

    # PROGRAM TITLE: PathReporter
    # Copyright (C) 2022  
    # Hansen Lam, M.D., Freddy Nguyen, M.D., Ph.D., Aryeh Stock, M.D., John Kim, M.D., Karl Eichorn, Dan Wild
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
import ujson
import statistics
from case_dictionary_methods import *

# CREATING DIAGNOSISBLOCK CLASS:________________________________________________
class DiagnosisBlock():
    def __init__(self, case, part, specimen, text, *args, **kwargs):
        self.phrase = text
        self.case = case
        self.part = part
        self.specimen = specimen
        self.summary_options = ['DEFINITE', 'HEDGE', 'DESCRIPTIVE']
        self.summary_categories = ['MALIGNANT', 'PREMALIG', 'REACTIVE', 'BENIGN']
        self.dict_terms = kwargs.get('terms', None)
        self.test_dict = kwargs.get('dict', {})
        self.diagnosis_summary = None
        self.diagnosis_category = None
        self.assigned = True
        if(self.dict_terms is not None):
            for term in self.dict_terms:
                self.test_dict[term] = []
        self.term_map = diagnosis_mapper(self.phrase, self.test_dict)
        self.simple_term_map = self.simplify_map()
        self.stop_term_map = self.stop_map()
        self.init_pointer()
    def simplify_map(self):
        term_list = []
        map_new = []
        for term in self.term_map:
            term_list.append(term[0])
        num_terms = max(term_list)
        for counter in range(num_terms):
            term_new = None
            for term in self.term_map:
                if term[0] == (counter+1):
                    if term_new is None:
                        term_new = term[1]
                    else:
                        term_new = ' '.join([term_new, term[1]])
                    entity = term[3]
            map_new.append([counter+1, term_new, entity])
        return map_new
    def stop_map(self):
        stop_map = []
        stop_contents = []
        stop_position = 0
        for counter in range(len(self.simple_term_map)):
            term_meta = self.simple_term_map[counter]
            if term_meta[2] == 'STOP':
                stop_map.append(stop_contents)
                stop_position += 1
                stop_contents = []
            else:
                if((counter==(len(self.simple_term_map)-1)) & (len(stop_contents)>0)):
                    stop_map.append(stop_contents)
                    stop_contents = []
                else:
                    stop_contents.append(term_meta)
        return stop_map
    def diagnostic_summary(self, *args, **kwargs):
        entity_num = 0
        preceding_hedge = 0
        score = 0
        avg_score = 0
        feat_num = 0
        category = None
        score_list = []
        key_term = None
        key_term_list = []
        negative_list = []
        negative_term = None
        entity_dict = ({
            'hedge':[],
            'definite':[],
            'negative':[]
        })
        entity_feat_dict = ({
            'SEPARATE':[]
        })
        feat_list = []
        return_label = kwargs.get('return_item', 'default')
        negative_pointer = None
        for term_collection in self.stop_term_map:
            #print(entity_dict)
            print(term_collection)
            negative_term = None
            entity_feats = None
            entity_pointer = None
            negative_list = []
            key_term = None
            for term in term_collection:
                entity_list = entity_dict['definite'] + entity_dict['hedge']
                entity_reference = False
                if term[2] == 'HEDGE':
                    if key_term is None:
                        key_term = term[1]
                    else:
                        key_term = ' '.join([key_term, term[1]])
                    preceding_hedge += 1
                if term[2] == 'NEGATIVE':
                    negative_term = term[1]
                    if ((term_collection.index(term)-1) >= 0):
                        negative_term = (
                            term_collection[term_collection.index(term)-1][1]
                            + ' ' + negative_term
                        )
                    if ((term_collection.index(term)+1) < len(term_collection)):
                        negative_term = negative_term + ' ' + term_collection[term_collection.index(term)+1][1]
                    negative_list.append(negative_term)
                    print(negative_list)
                    negative_pointer = self.stop_term_map.index(term_collection)
                if term[2] == 'FEATURE':
                    if(entity_pointer is not None):
                        entity_feat_dict[entity_pointer].append(term[1])
                    else:
                        if((negative_term is not None)&(self.stop_term_map.index(term_collection)==negative_pointer)):
                            #print(negative_term)
                            #if any(re.findall(term[1], negative_term, flags=re.IGNORECASE)):
                            #    feat_list.append('NEGATIVE' + ' ' + term[1])
                            #else:
                            feat_list.append('NEGATIVE' + ' ' + term[1])
                        else:
                            feat_list.append(term[1])
                        print('feat_list: ', feat_list)
                    feat_num += 1
                if ((term[2]=='ENTITY')|(term[2]=='ENTITY REFERENCE')):
                    if term[2] == 'ENTITY':
                        print(term)
                        vague_entity = False
                        vague_identity = None
                        vague_term = None
                        if(entity_feat_dict.get(term[1]) is None):
                            for entity in entity_list:
                                if(entity in term[1]):
                                    vague_entity = True
                                    vague_identity = term[1]
                                    vague_term = entity
                                else:
                                    if(term[1] in entity):
                                        vague_entity = True
                                        vague_identity = entity
                            if vague_entity:
                                entity_pointer = vague_identity
                                if(vague_term is not None):
                                    entity_feat_dict[entity_pointer] = entity_feat_dict.pop(vague_term)
                                    for key in entity_dict.keys():
                                        if vague_term in entity_dict[key]:
                                            entity_dict[key] = [term.replace(vague_term, vague_identity) for term in entity_dict[key]]
                                            #print(entity_dict)
                            else:
                                entity_feat_dict[term[1]] = []
                                entity_pointer = term[1]
                        else:
                            entity_pointer = term[1]
                        if(len(feat_list)>0):
                            entity_feat_dict[entity_pointer] = entity_feat_dict[entity_pointer] + feat_list
                            feat_list = []
                        neg_match = False
                        neg_phrase = None
                        neg_ref = None
                        if key_term is None:
                            key_term = term[1]
                        else:
                            key_term = ' '.join([key_term, term[1]])
                        if(term_collection.index(term)+1<(len(term_collection))):
                            if(term_collection[term_collection.index(term)+1][2]=='HEDGE'):
                                key_term = ' '.join([key_term, term_collection[term_collection.index(term)+1][1]])
                                preceding_hedge += 1
                        entity_num += 1
                        if(len(negative_list)>0):
                            #print(negative_list)
                            for neg in negative_list:
                                for t in key_term.split(' '):
                                    if t in neg:
                                        neg_match = True
                                        neg_phrase = neg
                                        if(neg_ref is None):
                                            neg_ref = t
                                        else:
                                            neg_ref = neg_ref + ' ' + t
                                        score = 1
                            if not neg_match:
                                print('NOT NEG_MATCH')
                                print('neg_ref: ',neg_ref)
                                neg_match = True
                                neg_phrase = negative_list[len(negative_list)-1] + ' ' + key_term
                                if(neg_ref is None):
                                    neg_ref = key_term
                                else:
                                    neg_ref = neg_ref + ' ' + key_term
                                score = 1
                                print('neg_ref assigned: ', neg_ref)
                        if neg_match:
                            entity_dict['negative'].append(neg_ref)
                            print('entity_dict[negative] assigned: ', neg_ref)
                            key_term_list.append(neg_phrase)
                        if not neg_match:
                            key_term_list.append(key_term)
                            score = 1/(1+preceding_hedge)
                            if(0<score<=0.5):
                                entity_dict['hedge'].append(entity_pointer)
                            else:
                                entity_dict['definite'].append(entity_pointer)
                                print('entity_dict[definite] assigned: ', entity_pointer)
                        score_list.append(score)
                        #key_term = None
                        preceding_hedge = 0
                    if ((term[2]=='ENTITY REFERENCE')&(len(entity_list)>0)):
                        if(len(entity_list)==1):
                            entity_reference = True
                            entity_assigned = entity_list[0]
                        elif(len(entity_dict['definite'])==1):
                            entity_reference = True
                            entity_assigned = entity_list[0]
                        else:
                            for entity in entity_list:
                                if(term[1] in entity):
                                    entity_reference = True
                                    entity_assigned = entity
                            if not entity_reference:
                                if ((term_collection.index(term)-1) >= 0):
                                    preceding_term = term_collection[term_collection.index(term)-1][1]
                                    for entity in entity_list:
                                        if(preceding_term in entity):
                                            entity_reference = True
                                            entity_assigned = entity
                        if entity_reference:
                            entity_pointer = entity_assigned
                        else:
                            entity_pointer = None
            if(len(feat_list)>0):
                entity_feat_dict['SEPARATE'] = entity_feat_dict['SEPARATE'] + feat_list
                feat_list = []
        if(entity_num==0):
            if(feat_num==0):
                score_list.append(-1)
            else:
                score_list.append(0)
        avg_score = statistics.mean(score_list)
        if (avg_score>0.5):
            category = 'DEFINITE'
        elif (0<avg_score<=0.5):
            category = 'HEDGE'
        elif (avg_score==0):
            category = 'DESCRIPTIVE'
        else:
            category = 'UNKNOWN'
        return_base = [self.case, self.part]
        if(return_label=='default'):
            return_add = [self.phrase, entity_num, key_term_list, entity_dict['hedge'], entity_dict['definite'], entity_dict['negative'], score_list, avg_score, category]
            return_item = return_base + return_add
        elif(return_label=='entity_feats'):
            return_item = []
            return_row = return_base
            for key, value in entity_feat_dict.items():
                return_row = return_base + [self.phrase] + [key, value]
                return_item = return_item + [return_row]
        else:
            return_item = None
        return return_item

    def init_pointer(self):
        self.assigned = True
        self.pointer_index = 0
        while ((self.assigned) & (self.pointer_index < len(self.term_map))):
            if self.term_map[self.pointer_index][3] is None:
                self.assigned = False
            else:
                self.pointer_index += 1
        self.for_linker_index = self.pointer_index
        self.back_linker_index = self.pointer_index
        if(self.assigned):
            self.selected_term_meta = []
        else:
            self.selected_term_meta = self.term_map[self.back_linker_index:self.for_linker_index+1]
    def pointer_skip(self):
        if ((self.pointer_index+1) <= len(self.term_map)):
            skip_pointer = self.pointer_index+1
            assigned = True
            while assigned:
                if self.term_map[skip_pointer][3] is None:
                    assigned = False
                else:
                    skip_pointer += 1
            self.pointer_index = skip_pointer
            self.for_linker_index = self.pointer_index
            self.back_linker_index = self.pointer_index
            self.selected_term_meta = self.term_map[self.back_linker_index:self.for_linker_index+1]
    def link_term(self, linker, direction):
        if linker == 'forward':
            if direction == 'for':
                if (self.for_linker_index <= len(self.term_map)):
                    self.for_linker_index += 1
            if direction == 'back':
                if ((self.for_linker_index-1) >= self.pointer_index):
                    self.for_linker_index += -1
        if linker == 'backward':
            if direction == 'for':
                if (self.back_linker_index < self.pointer_index):
                    self.back_linker_index += 1
            if direction == 'back':
                if ((self.back_linker_index-1) > 0):
                    self.back_linker_index += -1
        if direction == 'reset':
            self.for_linker_index = self.pointer_index
            self.back_linker_index = self.pointer_index
        self.selected_term_meta = self.term_map[self.back_linker_index:self.for_linker_index+1]
    def get_selected_term_meta(self):
        return self.selected_term_meta
    def get_selected_pointer_range(self):
        return range(self.pointer_index, self.linker_index+1)
    def set_diag_summary(self, summary):
        if(summary in self.summary_options):
            self.diagnosis_summary = summary
    def set_diag_category(self, cat):
        if(cat in self.summary_categories):
            self.diagnosis_category = cat
    def set_dict(self, dict):
        self.test_dict = dict
    def set_diag_class(self, value):
        self.diagnosis_summary = value
    def set_diag_category(self, value):
        self.diagnosis_category = value
    def update_map(self):
        self.term_map = diagnosis_mapper(self.phrase, self.test_dict)
    def get_text(self):
        return self.phrase
    def get_dict(self):
        return self.test_dict
    def get_map(self):
        return self.term_map
    def get_case(self):
        return self.case
    def get_part(self):
        return self.part
    def get_specimen(self):
        return self.specimen
    def get_diag_class(self):
        return self.diagnosis_summary
    def get_diag_category(self):
        return self.diagnosis_category
    def get_assigned(self):
        return self.assigned
