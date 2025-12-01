from urllib import response
from app.api.calLindel import calLindelScore
from app.api.calRE import find_cut_positions
from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import JSONResponse
from celery.result import AsyncResult
from pydantic import BaseModel
from typing import Union
import httpx
import subprocess, gffutils
import redis.asyncio as aioredis

from .test3 import calMicroScore
from .cmdprocess import extract_exon_by_gene, get_fasta_from_twobit, fold_rna
from Bio.Seq import Seq
import re, os, json, random, string, time, shutil
from app.configs import get_settings
import redis

settings = get_settings() 

redis_client_fq = aioredis.from_url(
    url=settings.REDIS_FQ_URL,
    decode_responses=True
)
MAX_CONCURRENT_TASKS = 2

GAP = 500
router = APIRouter()
scaf_default = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"

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

def pam_to_regex(pam):
    return "".join(IUPAC_MAP.get(base, base) for base in pam.upper())

class TraitToID(BaseModel):
    species: str
    tt: str

class GeneralSetting(BaseModel):
    sgRNA_len: int
    target: str
    flank: bool
    iso_union: str
    lpromoter: int
    rpromoter: int
    min_gc: int
    max_gc: int
    scaffoldSeq: str


class cas9Setting(BaseModel):
    name: str = "cas9"
    pam: str
    requi5: str
    off_target: bool #0 la ca vung protospacer, 1 la chi vung distal region
    mismatch_num: int

class Cas12aSetting(BaseModel):
    name: str = "cas12a"
EnzymeSetting = Union[cas9Setting, Cas12aSetting]


class GeneInfo:
    def __init__(self, id: str, strand: str, start: int, end: int, seq: str):
        self.id = id
        self.strand = strand
        self.start = start
        self.end = end
class Data(BaseModel):
    gene_name: str
    species: str

class CoordinateEntry(BaseModel):
    coordinate: str
    species: str

class vpcName(BaseModel):
    idfile: str

class FastaEntry(BaseModel):
    dna_seq: str
    species: str

class Primer(BaseModel):
    namefile: str
    idRow: str

class PrimerSetting(BaseModel):
    min_product_size: int
    max_product_size: int
    min_primer_size: int
    max_primer_size: int
    optimal_primer_size: int
    min_tm: int
    max_tm: int
    optimal_tm: int

class LindelRequest(BaseModel):
    idfile: str
    idRow: str
class RERequest(BaseModel):
    idfile: str
    idRow: str
class SingleBowtieDetailsRequest(BaseModel):
    idfile: str
    idRow: str
class ScoreDetailsRequest(BaseModel):
    idfile: str
    idRow: str

def set_enzyme(setting: EnzymeSetting):
    if isinstance(setting, Cas9Setting):
        return setting
    elif isinstance(setting, Cas12aSetting):
        return setting
    else:
        raise ValueError("Unknown enzyme setting")

def find_pam_positions(seq, pam):
    regex_mo = pam_to_regex(pam)
    pattern = re.compile(f"(?={regex_mo})")
    return [m.start() for m in pattern.finditer(seq.upper())]

def fastaParse(fullstr: str):
    res = []
    entries = fullstr.strip().split('>')

    for entry in entries:
        if not entry.strip():
            continue
        lines = entry.strip().split('\n')
        seq_name = lines[0].strip()
        seq = ''.join(lines[1:]).replace(" ", "").upper()
        res.append((seq_name, seq))  #day la tuple
    return res

def gop(a: str, b:str) -> str:
    return f"{a}:{b}"

def write_sgrna_to_fasta(sgrnas, filename="app/data/sgrna_output.fa"):
    with open(filename, "w") as f:
        index = 1
        for seq in sgrnas:
            prefix = seq[:-3]
            f.write(f">gRNA{index}\n{prefix}\n")
            index += 1

def write_sgrna_to_fasta2(sgrnas, filename="app/data/sgrna_output.fa"):
    with open(filename, "w") as f:
        index = 1
        for seq in sgrnas:
            prefix = seq[:-3]
            f.write(f">{index}\n{prefix}TGG\n")
            index += 1
            f.write(f">{index}\n{prefix}AGG\n")
            index += 1
            f.write(f">{index}\n{prefix}GGG\n")
            index += 1
            f.write(f">{index}\n{prefix}CGG\n")
            index += 1

def gc_content(seq: str, len_without_pam: int) -> float:
    seq = seq[0:len_without_pam]
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100

def random_string():
    chars = string.ascii_letters + string.digits
    return ''.join(random.choices(chars, k=10))

def save_sgRNA_list(idd: str, data: dict, gene_name: str, spec: str, pam: str, sgRNA_len: int, type_task: str,
                    q1: int, q2: int, q3: int, q4: int, q5: int, q6: int, q7:int, q8: int, queue_task_id = "",
                     status: str = "processing", log = "", stage = 1, gene_strand="+"):
    
    if stage == 0: #gan index computing
        filename = "vcp" + idd + ".json"
        path = os.path.join(DATA_DIR, filename)

        meta_obj = {"name": gene_name, "spec": spec, "pam": pam, "sgRNA_len": sgRNA_len, "gene_strand": gene_strand, 
                    "type_task": type_task,
                    "min_product_size": q1, "max_product_size": q2,
                    "min_primer_size": q3, "max_primer_size": q4, "optimal_primer_size": q5,
                    "min_tm": q6, "max_tm": q7, "optimal_tm": q8, "status": status, "log": log, "queue_task_id": queue_task_id}
        
        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)
            data.pop(0)
            data.insert(0, meta_obj)
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=4, ensure_ascii=False)
        return

    if idd == "unk":
        x = random_string()
        filename = "vcp" + x + ".json"
        path = os.path.join(DATA_DIR, filename)

        meta_obj = {"name": gene_name, "spec": spec, "pam": pam, "sgRNA_len": sgRNA_len, "gene_strand": gene_strand,
                    "type_task": type_task,
                    "min_product_size": q1, "max_product_size": q2,
                    "min_primer_size": q3, "max_primer_size": q4, "optimal_primer_size": q5,
                    "min_tm": q6, "max_tm": q7, "optimal_tm": q8, "status": status, "log": log, "queue_task_id": queue_task_id}
        data_with_meta = [meta_obj] + data
    elif stage == 1:
        x = idd
        filename = "vcp" + x + ".json"
        path = os.path.join(DATA_DIR, filename)

        meta_obj = {"name": gene_name, "spec": spec, "pam": pam, "sgRNA_len": sgRNA_len, "gene_strand": gene_strand,
                    "type_task": type_task,
                    "min_product_size": q1, "max_product_size": q2,
                    "min_primer_size": q3, "max_primer_size": q4, "optimal_primer_size": q5,
                    "min_tm": q6, "max_tm": q7, "optimal_tm": q8, "status": status, "log": log, "queue_task_id": queue_task_id}
        data_with_meta = [meta_obj] + data
    
    print(os.getcwd())
    with open(path, "w") as f:
        json.dump(data_with_meta, f, indent=2)
    return x

def getPrimerSeq(seq: str, i: int):
        start = i - 200
        end = i + 180

        if start < 0 or end > len(seq):
            return ""
        return seq[start:end]

def check_5require(seq: str, requi5: str) -> bool:
    seq = seq.upper()
    if requi5 == "GNNG":
        x = seq[:2]
        return x in ("GA", "GC", "GT", "GG", "AG", "CG", "TG", "GG")
    elif requi5 == "GG":
        return seq.startswith("GG")
    else:
        return True

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "data")

def utr35_from_computational(request: Data, generalSetting: GeneralSetting, casData: cas9Setting, geneid: str, filename: str):

    final_st_end = []
    keyw_rna = f"Parent={geneid}"
    command = f'grep "{keyw_rna}" {filename}'
    result2 = subprocess.run(
        command,
        shell=True,
        capture_output=True,
        text=True,
        check=True,
        cwd=DATA_DIR,
    )
    x_rna = result2.stdout
    lists_rna = x_rna.splitlines()
    rna_types = ["mRNA", "transcript", "primary_transcript", "lnc_RNA"]
    lists_rna = [line for line in lists_rna if any(rna in line for rna in rna_types)]
    tier3_lines = []
    idx  = 0
    for line in lists_rna:
        idx = idx + 1
        fields = line.strip().split('\t')
        id_rna = line.split('\t')[8].split(';')[0].split('=')[1]
        keyw_child = f"Parent={id_rna}"
        command = f'grep "{keyw_child}" {filename}'

        result3 = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
            check=True,
            cwd=DATA_DIR,
        )
        tier3_lines = result3.stdout
        if generalSetting.target == "5utr": 
            l = fields[3]
            r = -1
            x = [line for line in tier3_lines.splitlines() if "CDS" in line]
            if not x:
                final_st_end.append((fields[0], int(l), int(l), idx)) #đoạn này dung la phải dùng ATG để xác định
            else:
                for line in x:
                    fields = line.strip().split('\t')
                    r = min(r, int(fields[3]))
            final_st_end.append((fields[0], int(l), r, idx))
        elif generalSetting.target == "3utr":
            l = 1e15
            r = fields[4]
            x = [line for line in tier3_lines.splitlines() if "CDS" in line]
            if not x:
                final_st_end.append((fields[0], int(l), int(l), idx)) #đoạn này dung la phải dùng stop-codon UAA, UAG hoac UGA để xác định
            else:
                for line in x:
                    fields = line.strip().split('\t')
                    l = max(l, int(fields[4]))
            final_st_end.append((fields[0], l, int(r), idx))
    return final_st_end


def getAnnotationFile(spec: str):
    filename_gff3 = f"{spec}.gff3"
    filename_gtf = f"{spec}.gtf"
    filename_gff = f"{spec}.gff"

    candidates = [filename_gff3, filename_gtf, filename_gff]
    filename = None
    for fname in candidates:
        path = os.path.join(DATA_DIR, fname)
        if os.path.isfile(path):
            filename = fname
            break
    if filename is None:
        return None
    return filename
    
def updatePrimerConfig(seq_in_primer_template: str, gene_strand: str, idfile: str, x: str, min_product_size: int, max_product_size: int, 
                                min_primer_size: int, max_primer_size: int, optimal_primer_size: int,
                                min_tm: int, max_tm: int, optimal_tm: int):
    

    x = x.upper()
    seq_in_primer_template = seq_in_primer_template.upper()

    idx = x.find(seq_in_primer_template)

    if idx == -1:
        seq_in_primer_template = str(Seq(seq_in_primer_template).reverse_complement())
        idx = x.find(seq_in_primer_template)

    taskPrimerConfig = f"{idfile}_primerInp.txt"

    file_path = os.path.join(DATA_DIR, "primerInp.txt")
    taskPrimerConfig = os.path.join(DATA_DIR, taskPrimerConfig)
    
    if os.path.exists(file_path):
        shutil.copy(file_path, taskPrimerConfig)
        print(f"Đã copy file thành công: {taskPrimerConfig}")
    else:
        print(f"Lỗi: Không tìm thấy file gốc {file_path}")


    with open(taskPrimerConfig, 'r', encoding='utf-8') as f:
        content = f.read()

    updated_content = re.sub(r'SEQUENCE_TEMPLATE=.*', f'SEQUENCE_TEMPLATE={x}', content)
    updated_content = re.sub(r'PRIMER_PRODUCT_SIZE_RANGE=.*', f'PRIMER_PRODUCT_SIZE_RANGE={min_product_size}-{max_product_size}', updated_content)
    updated_content = re.sub(r'SEQUENCE_TARGET=.*', f'SEQUENCE_TARGET={idx-35},{len(seq_in_primer_template)+35}', updated_content)
    updated_content = re.sub(r'PRIMER_MIN_SIZE=.*', f'PRIMER_MIN_SIZE={min_primer_size}', updated_content)
    updated_content = re.sub(r'PRIMER_MAX_SIZE=.*', f'PRIMER_MAX_SIZE={max_primer_size}', updated_content)
    updated_content = re.sub(r'PRIMER_OPT_SIZE=.*', f'PRIMER_OPT_SIZE={optimal_primer_size}', updated_content)
    updated_content = re.sub(r'PRIMER_MIN_TM=.*', f'PRIMER_MIN_TM={min_tm}', updated_content)
    updated_content = re.sub(r'PRIMER_MAX_TM=.*', f'PRIMER_MAX_TM={max_tm}', updated_content)
    updated_content = re.sub(r'PRIMER_OPT_TM=.*', f'PRIMER_OPT_TM={optimal_tm}', updated_content)
                             

    with open(taskPrimerConfig, 'w', encoding='utf-8') as f:
        f.write(updated_content)

def createPrimer(gene_strand, idfile: str, seq_in_primer_template: str, sgrna_pos: str):

    taskPrimerConfig = f"{idfile}_primerInp.txt"
    command = f"primer3_core < {taskPrimerConfig}"
    result = subprocess.run(
        command,
        shell=True,
        capture_output=True,
        text=True,
        check=True,
        cwd=DATA_DIR,
    )
    x = result.stdout
    items = x.strip().split()
    data = {}


    for item in items:
        if '=' in item:
            key, value = item.split('=', 1)
            data[key] = value

    num_pairs = int(data.get('PRIMER_PAIR_NUM_RETURNED', 0))
    template_sequence = data.get('SEQUENCE_TEMPLATE')
    x = template_sequence.find(seq_in_primer_template)



    print(template_sequence)
    print(seq_in_primer_template)

    rev_in_seq = 0


    #sgRNA o chieu 3' - 5'
    if x == -1:
        rev_in_seq = 1
        seq_in_primer_template = str(Seq(seq_in_primer_template).reverse_complement())
        x = template_sequence.find(seq_in_primer_template)


    print(sgrna_pos)
    print(seq_in_primer_template)
    print(len(template_sequence) - int(x) - len(seq_in_primer_template) + 1)

    print(111111111)
    print(rev_in_seq)
    print(x)
    print(111111111)



    dd = sgrna_pos.split(':')[0]
    sgrna_pos = sgrna_pos.split(':')[1]
    print(sgrna_pos)
    sgrna_pos = int(sgrna_pos)
    print(gene_strand)

    primers = []
    # if gene_strand == '-':
    #     for i in range(min(5, num_pairs)):
    #         print(i)
    #         val_left = data.get(f'PRIMER_LEFT_{i}', '')
    #         start, length = map(int, val_left.split(','))
    #         left_location = f"{dd}:{sgrna_pos + (x - start) - length}-{sgrna_pos + (x - start)}"
    #         #5' 3' on gene trand "-"

    #         val_right = data.get(f'PRIMER_RIGHT_{i}', '')
    #         start, length = map(int, val_right.split(','))
    #         right_location = f"{dd}:{sgrna_pos + (x - start) - length}-{sgrna_pos + (x - start)}"
            
    #         primer_data = {
    #             "pair_index": i,
    #             "left_sequence": data.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''),
    #             "right_sequence": data.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
    #             "left_position": left_location,
    #             "right_position": right_location,
    #             "left_tm": float(data.get(f'PRIMER_LEFT_{i}_TM', 0)) if data.get(f'PRIMER_LEFT_{i}_TM') else 0.0,
    #             "right_tm": float(data.get(f'PRIMER_RIGHT_{i}_TM', 0)) if data.get(f'PRIMER_RIGHT_{i}_TM') else 0.0,
    #             "product_size": int(data.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0)) if data.get(
    #                 f'PRIMER_PAIR_{i}_PRODUCT_SIZE') else 0
    #         }
    #         primers.append(primer_data)
    # else:
    for i in range(min(5, num_pairs)):
        val_left = data.get(f'PRIMER_LEFT_{i}', '')
        start, length = map(int, val_left.split(','))
        left_location = f"{dd}:{sgrna_pos - (x - start)}-{sgrna_pos - (x - start) + length}"

        val_right = data.get(f'PRIMER_RIGHT_{i}', '')
        start, length = map(int, val_right.split(','))
        right_location = f"{dd}:{sgrna_pos + (start - x) - length}-{sgrna_pos + (start - x)}"
        
        primer_data = {
            "pair_index": i,
            "left_sequence": data.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''),
            "right_sequence": data.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
            "left_position": left_location,
            "right_position": right_location,
            "left_tm": float(data.get(f'PRIMER_LEFT_{i}_TM', 0)) if data.get(f'PRIMER_LEFT_{i}_TM') else 0.0,
            "right_tm": float(data.get(f'PRIMER_RIGHT_{i}_TM', 0)) if data.get(f'PRIMER_RIGHT_{i}_TM') else 0.0,
            "product_size": int(data.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0)) if data.get(
                f'PRIMER_PAIR_{i}_PRODUCT_SIZE') else 0
        }
        primers.append(primer_data)
    start_tem = sgrna_pos - x
    end_tem = sgrna_pos - x + len(template_sequence)
    return primers, start_tem, end_tem

@router.get("/getSingleBowtie/result/{task_id}")
async def get_single_bowtie_result(task_id: str):
    from .tasks import celery
    task_result = AsyncResult(task_id, app=celery)
    if task_result.ready():
        return {"status": "done", "result": task_result.result}
    return {"status": "pending"}

class indexComputingRequest(BaseModel):
    idfile: str

class GetFasta(BaseModel):
    twobitfile: str
    loca: str
    s1: int
    s2: int
@router.post("/getFastaData")
async def getFastaData(data: GetFasta):
    try:
        x = get_fasta_from_twobit(data.twobitfile, data.loca, data.s1 - 1, data.s2)
        x = x.upper()
        return {'real_seq': x}
    except Exception as e:
        return {"error": str(e)}



@router.post("/getPrimer")
async def getPrimer(data: Primer, primerSetting: PrimerSetting):
    try:
        filename = "vcp" + data.namefile + ".json"
        file_path = os.path.join(DATA_DIR, filename)
        with open(file_path, 'r', encoding='utf-8') as f:
            datafile = json.load(f)
        header = datafile[0]
        datafile = datafile[int(data.idRow) + 1]
        seq_in_primer_template = datafile['sequence']
        sgrna_pos = datafile['location']
        
        updatePrimerConfig(seq_in_primer_template, header['gene_strand'], data.namefile, datafile['Primer'], primerSetting.min_product_size, primerSetting.max_product_size,
                            primerSetting.min_primer_size, primerSetting.max_primer_size, primerSetting.optimal_primer_size,
                            primerSetting.min_tm, primerSetting.max_tm, primerSetting.optimal_tm)
        primers, start_tem, end_tem = createPrimer(header['gene_strand'], data.namefile, seq_in_primer_template, sgrna_pos)
        print(datafile)
        return {'first': primers, 'second': datafile['sequence'], 'third': datafile['name'].split(',')[1],
                'fourth': datafile['Primer'],
                'start_tem': start_tem,
                'end_tem': end_tem,
                }
    except Exception as e:
        return {"error": str(e)}

@router.post("/getLindelPre")
async def getLindelPre(data: LindelRequest):
    try:
        filename = "vcp" + data.idfile + ".json"
        file_path = os.path.join(DATA_DIR, filename)
        with open(file_path, 'r', encoding='utf-8') as f:
            datafile = json.load(f)
        datafile = datafile[int(data.idRow) + 1]
        query = datafile['lindel']
        response = calLindelScore([query])
        return {'data': response}
    except Exception as e:
        return {"error": str(e)}
    
@router.post("/getREData")
async def getREData(data: RERequest):
    try:
        filename = "vcp" + data.idfile + ".json"
        file_path = os.path.join(DATA_DIR, filename)
        with open(file_path, 'r', encoding='utf-8') as f:
            datafile = json.load(f)
        datafile = datafile[int(data.idRow) + 1]
        query = datafile['Primer']
        response = find_cut_positions(query)
        return {'data': response}
    except Exception as e:
        return {"error": str(e)}
    
def rev_comp_base(b):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return comp[b]

def reverse_complement(seq):
    return "".join(rev_comp_base(b) for b in reversed(seq))

def getMMsequence(original, bowtie_details):
    original = original.upper()
    details_list = [x.strip() for x in bowtie_details.split(";") if x.strip()]

    results = []

    for entry in details_list:
        # entry dạng: "NC_...,,0:T>G,1:G>C,6:C>A,0"
        parts = entry.split(",,")
        if len(parts) < 2:
            continue
        
        mutations = parts[1].split(",")
        mutations = [m for m in mutations if ":" in m]

        # Kiểm tra orientation dựa vào mismatch đầu tiên
        first = mutations[0]
        pos, change = first.split(":")
        expected, actual = change.split(">")

        pos = int(pos)

        is_reverse = (original[pos] != expected)

        if is_reverse:
            seq = reverse_complement(original)
        else:
            seq = original

        seq = list(seq)

        for m in mutations:
            pos, change = m.split(":")
            pos = int(pos)
            expected, actual = change.split(">")

            if is_reverse:
                new_pos = len(original) - 1 - pos
                seq[new_pos] = expected
            else:
                seq[pos] = expected

        results.append("".join(seq))

    return results

@router.post("/getSingleBowtieDetails")
async def getSingleBowtieDetails(data: SingleBowtieDetailsRequest):
    try:
        filename = "vcp" + data.idfile + ".json"
        file_path = os.path.join(DATA_DIR, filename)
        with open(file_path, 'r', encoding='utf-8') as f:
            datafile = json.load(f)
        header = datafile[0]
        datafile = datafile[int(data.idRow) + 1]
        data = datafile['bowtie_details']
        details = datafile['mismatch_region']

        original = datafile['sequence']
        mm_seq = getMMsequence(original, data)

        return {'data': data, 'mismatch_region': details, 'mm_sequence': mm_seq}
    except Exception as e:
        return {"error": str(e)}
    
@router.post("/getScoreDetails")
async def getScoreDetails(data: ScoreDetailsRequest):
    try:
        filename = "vcp" + data.idfile + ".json"
        file_path = os.path.join(DATA_DIR, filename)
        with open(file_path, 'r', encoding='utf-8') as f:
            datafile = json.load(f)
        datafile = datafile[int(data.idRow) + 1]
        cfdScore = datafile['cfdScore']
        microScore = datafile['microScore']
        mlScore = datafile['mlScore'] #rule set 2, doench 2016
        structure = datafile['Secondary structure with scaffold']
        rs3 = datafile['rs3']
        return {'cfd_score': cfdScore, 'micro_score': microScore, 'ml_score': mlScore, 'structure': structure, 'rs3': rs3}
    except Exception as e:
        return {"error": str(e)}

@router.post("/getsgRNAListFromFile")
async def getsgRNAListFromFile(data: vpcName):
    filename = "vcp" + data.idfile + ".json"
    #toan bo data
    path = os.path.join(DATA_DIR, filename)

    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail="File not found")

    max_retries = 3
    for attempt in range(max_retries):
        try:
            with open(path, "r") as f:
                file_data = json.load(f)
            return JSONResponse(content=file_data)
        except Exception as e:
            if attempt < max_retries - 1:
                await asyncio.sleep(1)
            else:
                raise HTTPException(status_code=500, detail=f"Some unknown error happened, please try after 5 minutes")

import math

def micro(seq: str) -> float:

    res = 0
    l = len(seq)
    for i in range(l):
        if seq[i] in ['A', 'T']:
            res += 1
        else:
            res += 2
    return res

def calMicroScore(seq1: str, seq2: str) -> float:

    out_frame_score = 0.0
    frame_score = 0.0
    l = len(seq1)
    sad = set()

    for i in range(l):
        for j in range(2, l - i):
            tmp = seq1[i:i + j]
            for k in range(l - j):
                if tmp == seq2[k:k + j]:

                    ext = j
                    while (i + ext < l and k + ext < l and seq1[i + ext] == seq2[k + ext]):
                        ext += 1
                    tmp = seq2[k:k+ext]

                    key = (i, k, tmp)
                    if key in sad:
                        continue
                    sad.add(key)

                    delta = k - i + l
                    mh_score = math.exp(-delta / 20) * micro(tmp) * 100
                    frame_score += mh_score
                    if (delta) % 3 != 0:
                        out_frame_score += mh_score
    if (frame_score == 0):
        return 0
    return out_frame_score / frame_score

def consensus(regions, mode):
    if not regions:
        return []
    print(regions)
    if len(regions) == 1:
        regions = list(regions)
        return [(regions[0][0], regions[0][1], regions[0][2])]    

    regions = sorted(list(regions), key=lambda x: (x[0], x[1], x[2], x[3]))

    if mode == "no":
        merged = []
        curr_chr, curr_start, curr_end = regions[0][0], regions[0][1], regions[0][2]

        for c, s, e, _ in regions[1:]:
            if c == curr_chr and s <= curr_end:  
                curr_end = max(curr_end, e)
            else:
                merged.append((curr_chr, curr_start, curr_end))
                curr_chr, curr_start, curr_end = c, s, e
        merged.append((curr_chr, curr_start, curr_end))
        return merged
    elif mode == "yes":
        from collections import defaultdict
        nhom = defaultdict(list)
        giao = []
        for c, s, e, i in regions:
            nhom[i].append((c, s, e))
        
        for c, s, e in nhom[1]:
            giao.append((c, s, e))
        for i in range(1, len(nhom)):
            new_giao = []
            for c1, s1, e1 in nhom[i+1]:
                for c, s, e in giao:
                    if not c1 == c:
                        continue
                    ss = max(s1, s)
                    ee = min(e1, e)
                    if ss < ee:
                        new_giao.append((c, ss, ee))
            giao = new_giao
        return giao
    else:
        raise ValueError("ko bt")

@router.post("/getDNAfromGeneName")
async def getDNAfromGeneName(request_fe: Request,
    request: Data, generalSetting: GeneralSetting, casData: cas9Setting, primerConfigData: PrimerSetting):


    ip = request_fe.client.host
    queue_name = 'lookUpComputing'
    redis_key = f"running:{queue_name}:{ip}"

    current = await redis_client_fq.incr(redis_key)
    await redis_client_fq.expire(redis_key, 3600*12)
    print(current)
    if current > MAX_CONCURRENT_TASKS:
        await redis_client_fq.decr(redis_key) 
        raise HTTPException(
            status_code=429, 
            detail="Quá nhiều request đồng thời cho IP này, vui lòng thử lại sau"
        )
    #----------------------------------------------------------------------------

    from .tasks import lookUpComputing_celery, redis_client, QUEUE_POSITION_PREFIX

    results = []
    gene_name = request.gene_name
    spec = request.species
    PAM = casData.pam
    sgRNA_len = generalSetting.sgRNA_len

    q1 = primerConfigData.min_product_size
    q2 = primerConfigData.max_product_size
    q3 = primerConfigData.min_primer_size
    q4 = primerConfigData.max_primer_size
    q5 = primerConfigData.optimal_primer_size
    q6 = primerConfigData.min_tm
    q7 = primerConfigData.max_tm
    q8 = primerConfigData.optimal_tm
    idd = save_sgRNA_list("unk", results, gene_name, spec, PAM, sgRNA_len, "gene_name",
                            q1, q2, q3, q4, q5, q6, q7, q8, status="pending")

    task_id = None    
    try:
        task = lookUpComputing_celery.apply_async(
            kwargs={
                'redis_key': redis_key,
                'idd': idd,
                'request': request.dict(),
                'generalSetting': generalSetting.dict(),
                'casData': casData.dict(),
                'primerConfigData': primerConfigData.dict(),
                'type': "gene_name"
            },
            queue=queue_name
        )
        task_id = task.id

        idd = save_sgRNA_list(idd, results, gene_name, spec, PAM, sgRNA_len, "gene_name",
                        q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id=task_id, status="queued")
    
        return {'first': idd, "task_id": task.id, "queue_name": queue_name}

    except Exception as e:
        await redis_client_async.decr(redis_key)
        save_sgRNA_list(idd, results, gene_name, spec, PAM, sgRNA_len, "gene_name",
                       q1, q2, q3, q4, q5, q6, q7, q8,
                       status="failed", log=f"Queue error: {str(e)}")
        raise HTTPException(
            status_code=503, 
            detail=f"Không thể xếp hàng task: {e}"
        )



@router.post("/getDNAfromCoordinate")
async def getDNAfromCoordinate(request_fe: Request,
    request: CoordinateEntry, generalSetting: GeneralSetting, casData: cas9Setting, primerConfigData: PrimerSetting):
    
    from .tasks import lookUpComputing_celery
    
    ip = request_fe.client.host
    queue_name = 'lookUpComputing'
    redis_key = f"running:{queue_name}:{ip}"

    current = await redis_client_fq.incr(redis_key)
    await redis_client_fq.expire(redis_key, 3600*12)

    if current > MAX_CONCURRENT_TASKS:
        await redis_client_fq.decr(redis_key) 
        raise HTTPException(
            status_code=429, 
            detail="Quá nhiều request đồng thời cho IP này, vui lòng thử lại sau"
        )
    #----------------------------------------------------------------------------
    results = []
    gene_name = request.coordinate
    spec = request.species
    PAM = casData.pam
    sgRNA_len = generalSetting.sgRNA_len

    q1 = primerConfigData.min_product_size
    q2 = primerConfigData.max_product_size
    q3 = primerConfigData.min_primer_size
    q4 = primerConfigData.max_primer_size
    q5 = primerConfigData.optimal_primer_size
    q6 = primerConfigData.min_tm
    q7 = primerConfigData.max_tm
    q8 = primerConfigData.optimal_tm
    idd = save_sgRNA_list("unk", results, gene_name, spec, PAM, sgRNA_len, "coordinate",
                            q1, q2, q3, q4, q5, q6, q7, q8, status="pending")
    
    task_id = None    
    try:
        task = lookUpComputing_celery.apply_async(
            kwargs={
                'redis_key': redis_key,
                'idd': idd,
                'request': request.dict(),
                'generalSetting': generalSetting.dict(),
                'casData': casData.dict(),
                'primerConfigData': primerConfigData.dict(),
                'type': "coordinate"
            },
            queue=queue_name
        )
        task_id = task.id

        idd = save_sgRNA_list(idd, results, gene_name, spec, PAM, sgRNA_len, "coordinate",
                        q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id=task_id, status="queued")
    
        return {'first': idd, "task_id": task.id, "queue_name": queue_name}

    except Exception as e:
        await redis_client_async.decr(redis_key)
        save_sgRNA_list(idd, results, gene_name, spec, PAM, sgRNA_len, "coordinate",
                       q1, q2, q3, q4, q5, q6, q7, q8,
                       status="failed", log=f"Queue error: {str(e)}")
        raise HTTPException(
            status_code=503, 
            detail=f"Không thể xếp hàng task: {e}"
        )
    
    


@router.post("/getDNAfromFasta")
async def getDNAfromFasta(request_fe: Request,
    request: FastaEntry, generalSetting: GeneralSetting, casData: cas9Setting, primerConfigData: PrimerSetting):

    from .tasks import lookUpComputing_celery
    
    ip = request_fe.client.host
    queue_name = 'lookUpComputing'
    redis_key = f"running:{queue_name}:{ip}"

    current = await redis_client_fq.incr(redis_key)
    await redis_client_fq.expire(redis_key, 3600*12)

    if current > MAX_CONCURRENT_TASKS:
        await redis_client_fq.decr(redis_key)
        raise HTTPException(
            status_code=429, 
            detail="Quá nhiều request đồng thời cho IP này, vui lòng thử lại sau"
        )
    #----------------------------------------------------------------------------
    results = []
    gene_name = request.dna_seq
    spec = request.species
    PAM = casData.pam
    sgRNA_len = generalSetting.sgRNA_len

    q1 = primerConfigData.min_product_size
    q2 = primerConfigData.max_product_size
    q3 = primerConfigData.min_primer_size
    q4 = primerConfigData.max_primer_size
    q5 = primerConfigData.optimal_primer_size
    q6 = primerConfigData.min_tm
    q7 = primerConfigData.max_tm
    q8 = primerConfigData.optimal_tm
    idd = save_sgRNA_list("unk", results, gene_name, spec, PAM, sgRNA_len, "fasta",
                            q1, q2, q3, q4, q5, q6, q7, q8, status="pending")
    
    task_id = None    
    try:
        task = lookUpComputing_celery.apply_async(
            kwargs={
                'redis_key': redis_key,
                'idd': idd,
                'request': request.dict(),
                'generalSetting': generalSetting.dict(),
                'casData': casData.dict(),
                'primerConfigData': primerConfigData.dict(),
                'type': "fasta"
            },
            queue=queue_name
        )
        task_id = task.id

        idd = save_sgRNA_list(idd, results, gene_name, spec, PAM, sgRNA_len, "fasta",
                        q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id=task_id, status="queued")
    
        return {'first': idd, "task_id": task.id, "queue_name": queue_name}

    except Exception as e:
        await redis_client_async.decr(redis_key)
        save_sgRNA_list(idd, results, gene_name, spec, PAM, sgRNA_len, "fasta",
                       q1, q2, q3, q4, q5, q6, q7, q8,
                       status="failed", log=f"Queue error: {str(e)}")
        raise HTTPException(
            status_code=503, 
            detail=f"Không thể xếp hàng task: {e}"
        )
class CheckPosition(BaseModel):
    task_id: str
    queue_name: str
@router.post("/checkPosition")
def checkPosition(data: CheckPosition):
    """
    Kiểm tra ZSET (hàng đợi "ảo") trước, 
    sau đó kiểm tra AsyncResult (trạng thái "thật")
    để tránh lỗi "nhảy" (Race Condition).
    """

    task_id = data.task_id
    queue_name = data.queue_name
    
    from .tasks import celery, redis_client, QUEUE_POSITION_PREFIX
    # 1. Tạo tên "Bảng trắng" (ZSET)
    # Ví dụ: "vicrispr_queue:lookUpComputing"
    zset_key = f"{QUEUE_POSITION_PREFIX}{queue_name}"
    
    # 2. KIỂM TRA "PHÒNG CHỜ" (ZSET)
    #    Hỏi: "Task này có đang xếp hàng không?"
    position_index = redis_client.zrank(zset_key, task_id)
    
    if position_index is not None:
        # --> TÌM THẤY: Task vẫn đang PENDING (đang chờ)
        total_in_queue = redis_client.zcard(zset_key) # Lấy tổng số người đang chờ
        
        return {
            "status": "PENDING",
            "position": position_index + 1, # (index từ 0, nên + 1)
            "total": total_in_queue
        }

    # 3. KHÔNG TÌM THẤY TRONG ZSET
    #    Đây là lúc va chạm (Race Condition) có thể xảy ra.
    #    Lý do: Nó vừa bị Worker "vợt" đi (STARTED) hoặc đã chạy XONG.
    #    Hỏi: "Thế trạng thái THẬT của nó là gì?"
    
    result = AsyncResult(task_id, app=celery)
    state = result.state
    
    if state in ("STARTED", "PROGRESS"):
        # A-ha! Bắt được va chạm. 
        # Nó không có trong ZSET vì nó ĐANG CHẠY.
        return {"status": "STARTED"}
        
    elif state in ("SUCCESS", "FAILURE"):
        # Nó không có trong ZSET vì nó ĐÃ XONG.
        return {"status": state}
        
    else: 
        # Vẫn là PENDING (trường hợp hiếm, ví dụ lỗi lúc submit)
        # hoặc UNKNOWN (sai task_id)
        return {"status": state}