from collections import defaultdict
import re
import pandas as pd
import string
import itertools
from collections import defaultdict
from xml_methods import *

def partition_dict_to_df(partition_dict, *args, **kwargs):
    return_item = kwargs.get('return_item', None)
    cases_bx = partition_dict
    df_specimen_diag_all = pd.DataFrame()
    df_case_table_all = pd.DataFrame()
    for index in cases_bx.keys():
        cases_bx_index = cases_bx[index]
        dict_case_table = {}
        accession_no = cases_bx_index['patientID_accession_no']
        specimens = []
        for key in cases_bx_index.keys():
            if('patientID' in key):
                dict_case_table[key] = cases_bx_index[key]
        df_specimen_diag = (
            pd.DataFrame.from_dict(cases_bx_index['fragments'], orient='index')
                        .reset_index()
                        .rename({'index':'Part',0:'Specimen',1:'Diagnosis',2:'Addendum'},axis=1)
        )
        df_specimen_diag['Case'] = accession_no
        df_case_table = (
            pd.DataFrame(dict_case_table, index=[0])
        )
        df_specimen_diag_all = pd.concat([df_specimen_diag_all, df_specimen_diag], ignore_index=True)
        df_case_table_all = pd.concat([df_case_table_all, df_case_table], ignore_index=True)
    #df_specimen_diag_all
    #df_case_table_all
    if return_item is None:
        return (df_case_table_all, df_specimen_diag_all)
    elif return_item=='case_table':
        return df_case_table_all
    elif return_item=='specimen_table':
        return df_specimen_diag_all
    else:
        return None


def string_clean(string, trash_list):
    phrase_clean = string.upper()
    clean_regex = ({
        r'\n':'',
        r' +':' ',
        r'(?<=[N])(.)(?=[0-9])':'',
        r'(?<=[0-9])(-)(?=[0-9])':'',
        r'(?<=Low)(.)(?=grade)':' ',
        r'(?<=High)(.)(?=grade)':' ',
        r'(?<=Intermediate)(.)(?=grade)':' ',
        r'(?<=Well)(.)(?=differentiated)':' ',
        r'(?<=Moderately)(.)(?=differentiated)':' ',
        r'(?<=Poorly)(.)(?=differentiated)':' '
    })
    for cleaner in clean_regex.keys():
        phrase_clean = re.sub(cleaner, clean_regex[cleaner], phrase_clean, flags=re.IGNORECASE|re.MULTILINE)
    for trash in trash_list:
        trash = r'\b'+trash.upper()+r'\b'
        if (len(re.findall(trash, phrase_clean))>0):
            phrase_clean = re.sub(trash, '', phrase_clean, flags=re.IGNORECASE|re.MULTILINE)
    phrase_clean = re.sub(r' +', ' ', phrase_clean, flags=re.IGNORECASE|re.MULTILINE)
    #phrase_clean = ' '.join([x.group() for x in re.finditer(r'\b\w+\b', phrase_clean, flags=re.IGNORECASE|re.MULTILINE)])
    return phrase_clean


def string_translate(string, thesaurus, *args, **kwargs):
    phrase_fragments = []
    start_key = kwargs.get('start', None)
    thesaurus_keys = list(thesaurus.keys())
    test_dict = thesaurus
    phrase = string.upper()
    key_order = [start_key] + [key for key in thesaurus_keys if key not in start_key]
    for key in key_order:
        for value in test_dict[key]:
            term = r'\b'+value.upper()+r'\b'
            if(len(re.findall(term, phrase, flags=re.IGNORECASE|re.MULTILINE))>0):
                iter = re.finditer(term, phrase, flags=re.IGNORECASE|re.MULTILINE)
                previous = [t for t in phrase_fragments]
                for m in iter:
                    assigned = False
                    if(len(previous)>0):
                        for i in previous:
                            if(((m.start()>=i[1][0]) & (m.end()<=i[1][1]))&(key!=i[2])):
                                assigned = True
                    if(not assigned):
                        word = value.upper()
                        span = m.span()
                        classification = key
                        phrase_fragments.append((word, span, classification))
    return phrase_fragments


# MAPPING DIAGNOSIS COMMENT SECTIONS INTO POS IN SENTENCE|WORD|SPAN|BOOL TRANSLATED
def string_map(phrase, assigned_fragments):
    phrase_clean = phrase
    phrase_fragments = assigned_fragments
    if(len(phrase_clean)>0):
        if(re.match(r'[A-Z]|\n',phrase_clean[len(phrase_clean)-1])):
            phrase_clean = phrase_clean + '.'
    term_iter = re.finditer(r"(?<=\s)\d*\.?\d+(?!=\s*CM)|(?<!\S)\d*\.?\d+(?!\S)|[\w']+|(?<!\d)\.(?!\d)|(?<![A-Z])\-(?![A-Z])|[,;]", phrase_clean, flags=re.IGNORECASE|re.MULTILINE)
    term_map = [[index, term.group(), term.span(), None] for index, term in enumerate(term_iter)]
    term_position = 0
    for word_pos in range(len(term_map)):
        word = term_map[word_pos][1]
        word_span = term_map[word_pos][2]
        assigned = False
        found = False
        frag_pos = 0
        if(re.match(r'(?<!\d)\.(?!\d)|(?<![A-Z])\-(?![A-Z])|[,;]', word)):
            term_map[word_pos][3] = 'STOP'
        else:
            while((not found) & (frag_pos < len(phrase_fragments))):
                frag = phrase_fragments[frag_pos]
                if((word_span[0]==frag[1][0])|(word_span[1]==frag[1][1])|((word_span[0]>frag[1][0]) & (word_span[1]<frag[1][1]))):
                    found = True
                    term_map[word_pos][3] = frag[2]
                    if (((word_span[0]==frag[1][0]) & (word_span[1] < frag[1][1]))
                        | ((word_span[0]>frag[1][0]) & (word_span[1]<frag[1][1]))
                        | ((word_span[0]>frag[1][0]) & (word_span[1]==frag[1][1]))):
                        assigned = True
                        if((word_span[0]==frag[1][0])):
                            term_position += 1
                else:
                    frag_pos += 1
        if((found & (not assigned))|((not found) & (not assigned))):
            term_position += 1
        term_map[word_pos][0] = term_position
    return term_map


def diagnosis_mapper(diagnosis, dict):
    phrase = diagnosis
    test_dict = dict
    remaining_classes = test_dict.copy()
    phrase_fragments = []
    phrase_clean = string_clean(phrase, remaining_classes['TRASH'])
    #del remaining_classes['TRASH']
    phrase_fragments = string_translate(phrase_clean, remaining_classes, start='ENTITY')
    term_map = string_map(phrase_clean, phrase_fragments)
    return term_map


def textblock_deidentify(textblock, *args, **kwargs):
    special_signatures = kwargs.get('special_sigs',[])
    physician_identifiers_regex = ([
        r'(?=([A-Z]+\s[A-Z]+,*\sM\.*D\.*[/]PHD))[A-Z]+',
        r'(?=([A-Z]+\s[A-Z]+,*\sM\.*D\.*))[A-Z]+',
        r'(?=([A-Z]+\.*\s[A-Z]+,*\sM\.*D\.*))[A-Z]+',
        r'(?=DR)([A-Z]+\.*\s*[A-Z]+\.*\s*[A-Z]+)',
        r'P\.S\.'
    ])
    block = textblock
    for regex in special_signatures:
        block = re.sub(regex, '', block, flags=re.IGNORECASE|re.MULTILINE)
    block = re.sub('  +', ' ', block)
    remove = string.punctuation
    remove = re.sub(r'[\*\-\([{})\]]', '', remove)
    identifiers = []
    text = block
    for regex in physician_identifiers_regex:
        for sub in re.findall(regex, text, flags=re.IGNORECASE|re.MULTILINE):
            if(sub not in identifiers):
                identifiers.append(sub)
    for person_verb in identifiers:
        for verb in ['review', 'agree', 'concur']:
            if verb in person_verb:
                for w in re.finditer(r'\b\w+\b', person_verb, flags=re.IGNORECASE|re.MULTILINE):
                    if verb in w.group():
                        word = w.group()
                        #print(word)
                person_only = re.findall(r'(?=DR\.*\s*[A-Z]+\s{0})(DR\.*\s*[A-Z]+)'.format(word), person_verb, flags=re.IGNORECASE|re.MULTILINE)
                #print(person)
                identifiers[identifiers.index(person_verb)] = person_only[0]
    for identity in identifiers:
        text = re.sub(identity, ' PHYSICIAN ', text, flags=re.IGNORECASE|re.MULTILINE)
    #print(identifiers)
    for part in re.findall(r'(?=[A-Z]+\-*\s*[A-Z]+\b,*\s*M\.*D\.*)([A-Z]+\-*\s*[A-Z]+)', text, flags=re.IGNORECASE|re.MULTILINE):
        text = re.sub(part, '', text, flags=re.IGNORECASE|re.MULTILINE)
    text = re.sub(' ,', ' ', text, flags=re.IGNORECASE|re.MULTILINE)
    counter = 0
    for part_iter in re.finditer(r'(?=[A-Z]+\.\s*[A-Z]*,*\sM\.*D\.*)([A-Z]+\.)', text, flags=re.IGNORECASE|re.MULTILINE):
        text = text[:part_iter.span()[0]-(counter*2)] + text[part_iter.span()[1]-(counter*2):]
        counter += 1
    for md in re.findall(r'\sM\.D\.\s|\bM\.D\.\b|\sMD\.*\b', text, flags=re.IGNORECASE|re.MULTILINE):
        text = re.sub(md, ' PHYSICIAN ', text, flags=re.IGNORECASE|re.MULTILINE)
    for extra_space in re.findall(r'(?<=\s)(\s[{}])'.format(remove), text, flags=re.IGNORECASE|re.MULTILINE):
        punct = re.findall(r'[{}]'.format(remove), extra_space)
        text = re.sub(extra_space, punct[0], text, flags=re.IGNORECASE|re.MULTILINE)
    text = re.sub('  +', ' ', text, flags=re.IGNORECASE|re.MULTILINE)
    for case in re.findall(r'([A-Z]{2}\-*\d{2}\-\d{5})', ' '+text+' ', re.IGNORECASE|re.MULTILINE):
        text = re.sub(case, 'CASE_ID', text, re.IGNORECASE|re.MULTILINE)
    if len(text) > 0:
        text = text.splitlines()
        series = pd.Series([line.upper() for line in text])
        series = series.drop(series[series.str.contains(r'PATIENT:|MRN:|DOB:|BIRTH:|NEW YORK')].index)
        text = '\n'.join(list(series))
    block = text
    return block


def text_decon_core(key, grouped_letter_index_tuples, synoptic_found, synoptic_text):
    #key = key.encode('utf-8')
    key = txt = re.sub(u'\u2013', '-', key)
    if(('MICROSCOPIC DESCRIPTIONS:' in key)|('GROSS DESCRIPTIONS:' in key)|('CONSULTATION DURING SURGERY' in key)|('FROZEN SECTION DIAGNOSIS:' in key)|(('DICT:' in key)&('PATHOLOGIST' in key))|('ATTENDING PATHOLOGIST' in key)):
        microscopic_start = [iter.start() for iter in re.finditer('MICROSCOPIC DESCRIPTIONS:',key,re.IGNORECASE|re.MULTILINE)]
        gross_start = [iter.start() for iter in re.finditer('MICROSCOPIC DESCRIPTIONS:',key,re.IGNORECASE|re.MULTILINE)]
        intraop_start = [iter.start() for iter in re.finditer('CONSULTATION DURING SURGERY',key,re.IGNORECASE|re.MULTILINE)]
        frozen_start = [iter.start() for iter in re.finditer('FROZEN SECTION DIAGNOSIS:',key,re.IGNORECASE|re.MULTILINE)]
        signature_start = [iter.start() for iter in re.finditer('PATHOLOGIST',key,re.IGNORECASE|re.MULTILINE)] + [iter.start() for iter in re.finditer('ATTENDING PATHOLOGIST',key,re.IGNORECASE|re.MULTILINE)]
        frozen_end = [iter.start() for iter in re.finditer('DICT:',key,re.IGNORECASE|re.MULTILINE)]
        list_index = []
        cut_list = microscopic_start + gross_start + intraop_start + frozen_start + signature_start + frozen_end
        for item in cut_list:
            #if(len(item)>0):
                list_index.append(item)
        cut_start = min(list_index)
        #print(cut_start)
        key = key[0:cut_start]
    synoptic_found[0] = False
    synoptic_start_indices = []
    synoptic_end_indices = []
    phd_index = []
    neuropath_start_indices = []
    neuropath_end_indices = []
    note_start_indices = []
    note_end_indices = []
    note_text = []
    #print('Pre_synoptic_found_key', key)
    part_text_line_split = key.splitlines()
    for line in part_text_line_split:
        if(any(re.findall(r'\,*\s*PH\.D\.',line,flags=re.IGNORECASE))):
            phd_index.append(part_text_line_split.index(line))
            #print('PH.D. index:',part_text_line_split.index(line))
    if(len(phd_index)>0):
        for x in range(len(phd_index)):
            del part_text_line_split[phd_index[x]-x]
            #print(part_text_line_split)
        part_text = '\n'.join(part_text_line_split)
        key = part_text
    for line in part_text_line_split:
        if(any(re.findall(r'NEUROPATHOLOGY DIAGNOSIS:',line,flags=re.IGNORECASE))):
            neuropath_start_indices.append(part_text_line_split.index(line))
            neuropath_start = min(neuropath_start_indices)
        elif(any(re.findall(r'PATHOLOGIST|ATTENDING|PHYSICIAN|Division of Neuropathology',line,flags=re.IGNORECASE))):
            neuropath_end_indices.append(part_text_line_split.index(line))
            neuropath_end = max(neuropath_end_indices)
    if((len(neuropath_start_indices)>0)&(len(neuropath_end_indices)>0)):
        del part_text_line_split[neuropath_start:neuropath_end]
        part_text = '\n'.join(part_text_line_split)
        key = part_text
    for line in part_text_line_split:
        if(re.search(r'SPECIMENS:|SPECIMEN:|PROCEDURE:|PROGNOSTIC INDICATORS|^SYNOPTIC REPORT',line, flags=re.IGNORECASE)):
            synoptic_start_indices.append(part_text_line_split.index(line))
        if(re.search(r'The part/slide number for a normal block is:|NORMAL BLOCK IS:|AJCC STAGE|ADDITIONAL PATHOLOGIC FINDINGS:|VENOUS/LYMPHATIC (LARGE/SMALL VESSEL) INVASION (V/L):', line, flags=re.IGNORECASE)):
            synoptic_end_indices.append(part_text_line_split.index(line))
    if((len(synoptic_start_indices)>0)&(len(synoptic_end_indices)>0)):
        if(min(synoptic_start_indices)<max(synoptic_end_indices)):
            synoptic_found[0] = True
            synoptic_start_line = min(synoptic_start_indices)
            synoptic_end_line = max(synoptic_end_indices)
            synoptic_line_split = part_text_line_split[synoptic_start_line:synoptic_end_line+1]
            for line in synoptic_line_split:
                synoptic_text.append(line)
            del part_text_line_split[synoptic_start_line:synoptic_end_line+1]
            part_text = '\n'.join(part_text_line_split)
            #print('synotpic_start:',synoptic_start_line, '\tsynoptic_end:',synoptic_end_line+1)
            #print('synoptic subtracted:',part_text)
            key = part_text
            #print('synoptic_found',synoptic_text)
    for line in part_text_line_split:
        if(re.findall(r'^NOTE\s*[\-\:]',line,flags=re.IGNORECASE)):
            note_start_indices.append(part_text_line_split.index(line))
        for reg in [r'(?<!PART)(?<!NEOPLASIA\s)^\s*?([A-Z])[.:)-](?!\s*\w*\sPHYSICIAN)(?!100)(?!\d)(?!CELL\t)',r'(?<!AIN)(?<!NEOPLASIA)(?<!GRADE)[,\s-](\n*\b^[A-Z])(,\s?)([A-Z][,.]\s?)*([A-Z][.])',r'PHYSICIAN']:
            if(re.findall(reg,line,flags=re.IGNORECASE)):
                note_end_indices.append(part_text_line_split.index(line))
    if((len(note_start_indices)>0)&(len(note_end_indices)>0)):
        note_start = max(note_start_indices)
        print('NOTE START',note_start)
        note_end_indices.append(note_start)
        note_end_indices.sort()
        note_end = None
        if(note_end_indices.index(note_start)<(len(note_end_indices)-1)):
            note_end = note_end_indices[note_end_indices.index(note_start)+1]
        else:
            note_end=note_start
        note_line_split = part_text_line_split[note_start:note_end+1]
        for line in note_line_split:
            note_text.append(line)
            print('NOTE TEXT: ','\s'.join(note_text))
        del part_text_line_split[note_start:note_end+1]
        part_text = '\n'.join(part_text_line_split)
        key = part_text

    list_letter_index_tuples = []
    single_part_list = []
    single_start_index = []
    rep_part_list = []
    rep_start_index = []
    single_group_counter = 0
    for exp in [r'(?<!PART)(?<!NEOPLASIA\s)^\s*?([A-Z])[.:)-](?!\s*\w*\sPHYSICIAN)(?!100)(?!\d)(?!CELL\t)',r'(?<=^PART\s)([A-Z])[.:)-](?!\s*\w*\sPHYSICIAN)']:
        for letter in re.finditer(exp, key, flags=re.MULTILINE|re.IGNORECASE):
            single_group_counter += 1
            for char in re.finditer(r'\b[A-Z]\b',letter.group()):
                single_part_list.append(char.group())
                single_start_index.append(letter.start())
                list_letter_index_tuples.append([char.group(),letter.start()])
    if(any(re.findall(r'(?<!AIN)(?<!NEOPLASIA)(?<!GRADE)[,\s-](\n*\b^[A-Z])(,\s?)([A-Z][,.]\s?)*([A-Z][.])',key,flags=re.IGNORECASE|re.MULTILINE))):
        rep_group_counter = 0
        for rep in re.finditer(r'(?<!AIN)(?<!NEOPLASIA)(?<!GRADE)[,\s-](\n*\b^[A-Z])(,\s?)([A-Z][,.]\s?)*([A-Z][.])',key,flags=re.IGNORECASE|re.MULTILINE):
            rep_group_counter += 1
            for char in re.finditer(r'\b[A-Z]\b',rep.group()):
                rep_part_list.append(char.group())
                rep_start_index.append(rep.start())
                list_letter_index_tuples.append([char.group(),rep.start()])

    list_letter_index_tuples.sort(key=lambda x:x[1])
    all_start_index = single_start_index + rep_start_index
    all_start_index.sort()
    #print('all_start_index:')
    #print(all_start_index)
    deduplicated_start_index = list(set(all_start_index))
    deduplicated_start_index.sort()
    deduplicated_start_index += [len(key)]
    del deduplicated_start_index[0]
    deduplicated_start_index
    grouped_letter_index_tuples += [list(g) for _,g in itertools.groupby(list_letter_index_tuples, key=lambda x:x[1])]
    for collect in grouped_letter_index_tuples:
        end_index = deduplicated_start_index[grouped_letter_index_tuples.index(collect)]
        for tuple in collect:
            tuple.append(end_index)
    return key

def textblock_deconstruct(key, part_names, part_list, grouped_letter_index_tuples, consult, consult_cases, is_diagnosis, synoptic_found, synoptic_text):
    consult_case_count = len(consult_cases)
    multi_consult = consult & (consult_case_count>1)
    if multi_consult:
        part_list.append('A')
        grouped_letter_index_tuples += [[['A',0,12]]]
        if(len(part_names)>1):
            part_names = [part_names[0]]
    else:
        key = text_decon_core(key, grouped_letter_index_tuples, synoptic_found, synoptic_text)
        part_count = 0
        for group in grouped_letter_index_tuples:
            for letter in group:
                if(letter[0] not in part_list):
                    part_list.append(letter[0])
                    part_count += 1
        if(part_count==0):
            grouped_letter_index_tuples += [[['A',0,len(key)]]]
            part_count = 1
            part_list.append('A')
        if(is_diagnosis):
            if((len(part_names) < part_count)&(len(part_names)==1)):
                for miss_count in range(len(part_list)-len(part_names)):
                    part_names += [part_names[0]]
            elif((len(part_names) < part_count)):
                for miss_count in range(len(part_list)-len(part_names)):
                    part_names += [part_names[len(part_names)-1]]

    return key


def textblock_to_dict(key, part_list, grouped_letter_index_tuples, dict_parts_text):
    #print(grouped_letter_index_tuples)
    for group in grouped_letter_index_tuples:
        for letter in group:
            if(len(dict_parts_text[letter[0]])>0):
                old_string = dict_parts_text[letter[0]][0]
                string_add = key[letter[1]:letter[2]]
                new_string = old_string + ' ' + string_add
                dict_parts_text[letter[0]] = [new_string]
            else:
                dict_parts_text[letter[0]] += [key[letter[1]:letter[2]]]
    #print(dict_parts_text)
    return 0
