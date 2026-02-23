import pandas as pd
from .export import count_permu_IUPAC

df = pd.read_csv("./app/data/cfdtable.csv")
df[["rna", "dna"]] = df["pair"].str.extract(r"r([ATGC]):d([ATGC])")
df["pos"] = df["pos"].astype(int)

def get_percent_active(pos, rna, dna):
    row = df[(df["pos"] == pos) & (df["rna"] == rna) & (df["dna"] == dna)]
    if not row.empty:
        return float(row.iloc[0]["pa"])
    else:
        return 1

def get_rna(seq):
    ax = {
        'A': 'U',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    return ax.get(seq, '?')

def get_cfd_score(line, pam_name: str, ll: int):
    parts = line.strip().split('\t')

    permu_num = count_permu_IUPAC(pam_name)

    x = int(int(parts[0]) / permu_num)
    if len(parts) < 8:
        return x, 0
    
    mismatch_info = parts[7]
    mismatch_info_parse = mismatch_info.split(",")
    positions = []
    for item in mismatch_info_parse[1:]:
        pos_str = item.split(":")[0]  
        positions.append(int(pos_str))
    if any(pos >= ll for pos in positions):
        return -1, 0

    items = mismatch_info.split(',')

    result = []
    posl = []
    for item in items:
        pos_part, change = item.split(':')
        dna, rna = change.split('>')
        result.append((int(pos_part), dna, get_rna(rna)))
        posl.append(int(pos_part))

    scr = 1
    l = len(posl)
    for item in result:
        k = get_percent_active(item[0], item[2], item[1])
        if k == 0:
            k = 0.0025
        scr = scr * k

    scr = scr / l

    return x, scr
