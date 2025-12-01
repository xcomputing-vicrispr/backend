import subprocess
import os

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAR_PARENT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORK_DIR = os.path.join(PAR_PARENT_DIR, "Lindel")

DATA_DIR = os.path.join(PARENT_DIR, "data")
json_dir = os.path.join(DATA_DIR, "re.json")

import json, re
IUPAC_MAP = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "W": "[AT]",
    "R": "[AG]",
    "M": "[AC]",
    "Y": "[CT]",
    "S": "[GC]",
    "K": "[GT]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]"
}

def load_enzymes():
    with open(json_dir, "r", encoding="utf-8") as f:
        enzymes = json.load(f)
    return enzymes


def motif_to_regex(motif):
    regex = ''
    for c in motif:
        regex += IUPAC_MAP.get(c.upper(), c)
    return regex

def parse_recognition_site(recognition):
    ##(a/b)NNN^N(c/d)
    # Nếu có dấu ^

    flag5 = ''
    flag3 = ''
    if '^' in recognition:
        cut_pos_in_motif = recognition.index('^')
        motif = recognition.replace('^','')
        motif_regex = motif_to_regex(motif)
        return motif_regex, cut_pos_in_motif, cut_pos_in_motif, flag5, flag3

    #parse 2 cai ngoac truoc sau ra    
    offsets = re.findall(r'\(([^/]+)/([^\)]+)\)', recognition)
    
    motif = re.sub(r'\([^\)]+\)', '', recognition)
    motif_regex = motif_to_regex(motif)
    
    def parse_offset(s):
        return int(s) if s.strip() != '?' else None
    
    cut5vt = []
    cut3vt = []
    
    if len(offsets) >= 1:
        cut5vt.append(parse_offset(offsets[0][0]))  
        cut3vt.append(parse_offset(offsets[0][1]))
    if len(offsets) == 2:
        cut5vt.append(parse_offset(offsets[1][0]))
        cut3vt.append(parse_offset(offsets[1][1]))

    if len(cut5vt) == 2 and cut5vt[1] is not None:
        cut5 = cut5vt[1]
        flag5 = 't'
    else: 
        cut5 = cut5vt[0]
        flag5 = 'p'

    if len(cut3vt) == 2 and cut3vt[1] is not None:
        cut3 = cut3vt[1]
        flag3 = 't'
    else: 
        cut3 = cut3vt[0]
        flag3 = 'p'

    return motif_regex, cut5, cut3, flag5, flag3


def find_cut_positions(sequence):

    results = {}
    enzim_len = 4
    enzymes = load_enzymes()
    
    for name, info in enzymes.items():
        recognition = info.get('recognition', '')
        motif_regex, cut5_offset, cut3_offset, flag5, flag3 = parse_recognition_site(recognition)

        if recognition.startswith("("):
            recognition = recognition.split(")", 1)[1]   
        if recognition.endswith(")"):
            recognition = recognition.rsplit("(", 1)[0]  
        if len(recognition) < enzim_len:
            continue
        
        for match in re.finditer(motif_regex, sequence):
            start, end = match.start(), match.end()


            if flag5 == 't':
                cut5 = start - cut5_offset
            elif flag5 == 'p':
                cut5 = end + cut5_offset
            else:
                cut5 = start + cut5_offset

            if flag3 == 't':
                cut3 = start - cut3_offset
            elif flag3 == 'p':
                cut3 = end + cut3_offset
            else:
                cut3 = start + cut3_offset
                 
            results.setdefault(name, []).append({
                "motif": match.group(),
                "motif_start": start,
                "motif_end": end,
                "cut5": cut5,
                "cut3": cut3
            })
    
    return results
    



