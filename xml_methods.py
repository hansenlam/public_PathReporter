from lxml import etree
from lxml.etree import tostring, parse
import re
import string
import ujson

def num_child(elem):
# XML parsing function
# function returns number of child elements associated with input element
# arguments: parent element
# returns: number of child elements
    counter=0
    for child in elem:
        counter+=1
    return counter

def walker_bx(elem, dict):
    return_dict=dict
    for sub_elem in elem:
        ns = re.findall(r'\{.*?\}', elem.tag)[0]
        if ('FieldName' in sub_elem.attrib.keys()):
            for end_elem in sub_elem:
                if (end_elem.tag==(ns+'FormattedValue')):
                    if('patientID' in sub_elem.attrib['FieldName']):
                            if('name' not in sub_elem.attrib['FieldName']):
                                dict_key = re.sub(r'{|}','',sub_elem.attrib['FieldName']).split(r'_history_by_',1)[1]
                                dict_key = re.sub(r'\.','_',dict_key)
                                return_dict[dict_key] = end_elem.text
                    elif(len(re.sub(r'{|}','',sub_elem.attrib['FieldName']).split('_case_',1))>1):
                        existing_keys = return_dict.keys()
                        current_attrib = re.sub(r'\.','_',re.sub(r'{|}','',sub_elem.attrib['FieldName']).split('_case_',1)[1])
                        key_occur = (
                            sum(current_attrib in key for key in existing_keys)
                        )
                        if key_occur >= 0:
                            current_attrib += str(key_occur+1)
                        return_dict[current_attrib] = end_elem.text
        if(num_child(sub_elem)>0):
            walker_bx(sub_elem, return_dict)
    return return_dict

def rootIndex_CaseID(elem, dict):
    return_dict=dict
    for sub_elem in elem:
        ns = re.findall(r'\{.*?\}', elem.tag)[0]
        if ('FieldName' in sub_elem.attrib.keys()):
            for end_elem in sub_elem:
                if (end_elem.tag==(ns+'FormattedValue')):
                    if('accession_no' in sub_elem.attrib['FieldName']):
                        dict_key = re.sub(r'{|}','',sub_elem.attrib['FieldName']).split(r'_history_by_',1)[1]
                        dict_key = re.sub(r'\.','_',dict_key)
                        return_dict[dict_key] = end_elem.text
        if(num_child(sub_elem)>0):
            rootIndex_CaseID(sub_elem, return_dict)
    return return_dict


def rootID_Indices(root_obj, range_list):
    root_dict = {}
    for index, case in enumerate(root_obj):
        if index in range_list:
            autopsy = False
            id_dict = rootIndex_CaseID(root_obj[index], {})
            if(len(id_dict.keys())>0):
                id = id_dict['patientID_accession_no']
                if any(autopsy_marker in id for autopsy_marker in ['MA','RA','SA','LA','BA']):
                    autopsy = True
                if not autopsy:
                    root_dict[id] = index
    return root_dict