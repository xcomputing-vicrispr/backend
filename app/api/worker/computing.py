from email import header
from celery import Celery
import time
from pydantic import BaseModel
import httpx
import subprocess
import time, gffutils


import traceback 
import smtplib, numpy as np
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from collections import defaultdict
from app.models import EmailQueue

import re, os, json, random, string

from app.api.cmdprocess import get_fasta_from_twobit
from Bio.Seq import Seq
import re, os, json, random, string
from app.api.lookUpsgRNA import *
from app.api.cfdEffiencyScore import get_cfd_score
from app.api.mlEffiencyScore import get_ml_score, get_ml_score_azi3
from app.api.calLindel import calLindelScore
from app.api.export import sendMail, xuly, MailSession, sigmoid, gc_score, write_sgrna_to_fasta2, write_sgrna_to_fasta_with_NNAGAAW
from app.api.export import sender, password, OUTPUT_DIR, DATA_DIR, write_sgrna_to_fasta_with_IUPAC, count_permu_IUPAC
from app.database import get_db, SessionLocal
from sqlalchemy import func

class IndexComputingSession(BaseModel):
    idfile: str

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

primerDefault = "AAAAAAAAAAAAAAAAAA"
mlseqDefault = "CCACCAGGTGGTTGGTGATTTTGGCGGGGG"
lindelDefault = "CAATCATCGCCACCAGGTGGTTGGTGATTTTGGCGGGGGCAGAGAGGACGGTGGCCACCT"
CONST_GAP = 500

def iupac_combinations(seq: str):
    iupac_map = {
        'A': 1, 'C': 1, 'G': 1, 'T': 1,
        'R': 2, 'Y': 2, 'S': 2, 'W': 2, 'K': 2, 'M': 2,
        'B': 3, 'D': 3, 'H': 3, 'V': 3,
        'N': 4
    }

    seq = seq.upper()
    total = 1
    for base in seq:
        total *= iupac_map.get(base, 1)
    return total


def check_gene_match(target_name, attr_string):
    target = target_name.strip().upper()
    
    attrs = {}
    for item in attr_string.split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            attrs[k] = v 
    
    #1 name
    if 'Name' in attrs and attrs['Name'].upper() == target:
        return True

    #2 gene
    if 'gene' in attrs and attrs['gene'].upper() == target:
        return True

    #3: synonym
    if 'gene_synonym' in attrs:
        synonyms = [s.strip().upper() for s in attrs['gene_synonym'].split(',')]
        if target in synonyms:
            return True

    # 4 Dbxref
    if 'Dbxref' in attrs:
        dbxrefs = attrs['Dbxref'].split(',')
        for ref in dbxrefs:
            ref_upper = ref.upper()
            # lấy đoạn sau kiểu :671 hay kiểu đấy
            if ref_upper.endswith(f":{target}"):
                return True

    #5 ID
    if 'ID' in attrs:
        id_val = attrs['ID'].upper()
        main_id = id_val.split('-', 1)[-1]
        if target == id_val or target == main_id:
            return True

    return False

def GeneNameComputing(queue_task_id, idd, request, generalSetting, casData, primerConfigData):
    print("Da vao ham GeneNameComputing")
    GAP = CONST_GAP
    try:
        request = Data(**request)
        generalSetting = GeneralSetting(**generalSetting)
        casData = cas9Setting(**casData)
        primerConfigData = PrimerSetting(**primerConfigData)
        results = []
        auke = []
        gene_name = request.gene_name
        spec = request.species

        print(primerConfigData)
        q1 = primerConfigData.min_product_size
        q2 = primerConfigData.max_product_size
        q3 = primerConfigData.min_primer_size
        q4 = primerConfigData.max_primer_size
        q5 = primerConfigData.optimal_primer_size
        q6 = primerConfigData.min_tm
        q7 = primerConfigData.max_tm
        q8 = primerConfigData.optimal_tm

        PAM = casData.pam
        len_without_pam = generalSetting.sgRNA_len
        REV_PAM = Seq(PAM)
        REV_PAM = str(REV_PAM.reverse_complement())
        scaffold = generalSetting.scaffoldSeq
        
        idd = save_sgRNA_list_dbv(idd, results, gene_name, spec, PAM, len_without_pam, "gene_name",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="Finding", log="Finding sgRNA candidates")

        filename = getAnnotationFile(spec)

        if filename is None:
            raise FileNotFoundError(
                f"Not found {spec}.gff3 / {spec}.gtf / {spec}.gff trong CSDL"
            )

        command = f'grep "{gene_name}" {filename}'
        name = gene_name + "," + spec + "," + PAM
        resolved_paths = resolve_spec(spec)
        twobit_file = resolved_paths["twobit"]
        cutting_sites = []

        db_path = resolved_paths["gffdb"]
        db = gffutils.FeatureDB(db_path, keep_order=True)
        tier2_types = ["mRNA", "transcript", "primary_transcript", "lnc_RNA", "snoRNA", "miRNA", "rRNA", "tRNA"]
        final_st_end = []
        gene = GeneInfo("", "", 0, 0, "")
        idx = 0
        gene_strand = "+"
        if generalSetting.target != "promoter":  # exon hoac CDS
            query = f"""
            SELECT *
            FROM features
            WHERE featuretype IN ('gene', 'pseudogene')
            AND attributes LIKE '%{request.gene_name}%'
            """
            query2 = f"""
            SELECT *
            FROM features
            WHERE featuretype IN ('mRNA', 'transcript', 'primary_transcript', 'lnc_RNA', 'snoRNA', 'miRNA', 'rRNA', 'tRNA')
            AND attributes LIKE '%{request.gene_name}%'
            """
            flag = 0 
            geneid = None
            for row in db.execute(query):
                stra = row['attributes']

                if not check_gene_match(request.gene_name, row['attributes']):
                    continue

                invalid_pattern = r'[a-zA-Z0-9]' + re.escape(request.gene_name) + r'|' + re.escape(request.gene_name) + r'[a-zA-Z0-9]'

                if re.search(invalid_pattern, stra):
                    continue 
                print(row['seqid'], row['start'], row['end'], row['strand'])
                gene.id = row['id']
                geneid = row['id']
                gene.start = row['start']
                gene.end = row['end']
                gene.strand = row['strand']
                print("gene:", gene.id, gene.start, gene.end, gene.strand)
                gene_strand = gene.strand
                flag = 1

                for tier2 in db.children(gene.id, featuretype=tier2_types):
                    idx = idx + 1
                    flag = 2
                    print("tier 2:", tier2.id, tier2.start, tier2.end, tier2.strand)

                    if generalSetting.target == "exon":
                        for exon in db.children(tier2.id, featuretype="exon"):
                            print("exon:", exon.id, exon.start, exon.end)
                            final_st_end.append((row['seqid'], int(exon.start), int(exon.end), idx))

                    elif generalSetting.target == "CDS": #CDS
                            for cds in db.children(tier2.id, featuretype="CDS"):
                                print("CDS:", cds.id, cds.start, cds.end)
                                final_st_end.append((row['seqid'], int(cds.start), int(cds.end), idx))

                    elif generalSetting.target == "3utr": #3utr
                            for utr in db.children(tier2.id, featuretype="three_prime_UTR"):
                                print("3utr:", utr.id, utr.start, utr.end)
                                final_st_end.append((row['seqid'], int(utr.start), int(utr.end), idx))
                    else: #5utr
                            for utr in db.children(tier2.id, featuretype="five_prime_UTR"):
                                print("5utr:", utr.id, utr.start, utr.end)
                                final_st_end.append((row['seqid'], int(utr.start), int(utr.end), idx))

                if flag == 1:
                    if generalSetting.target == "exon":
                        print("was gare")
                        for exon in db.children(gene.id, featuretype="exon"):
                            print("exon:", exon.id, exon.start, exon.end)
                            final_st_end.append((row['seqid'], int(exon.start), int(exon.end), idx))
                    elif generalSetting.target == "CDS": #CDS
                        for cds in db.children(gene.id, featuretype="CDS"):
                            print("CDS:", cds.id, cds.start, cds.end)
                            final_st_end.append((row['seqid'], int(cds.start), int(cds.end), idx))
                
                if (generalSetting.target == "5utr" or generalSetting.target == "3utr") and len(final_st_end) == 0: #5utr, 3utr in case no UTR annotation
                    final_st_end = utr35_from_computational(request, generalSetting, casData, geneid, filename)

                print(final_st_end, 232323)
                final_st_end = consensus(set(final_st_end), generalSetting.iso_union)
                print(final_st_end, 56565656)
                break
            if flag == 0:
                #Tien hanh tim trong tier2
                print("da tim trong tier 1 roi nhung ko co")
                for row in db.execute(query2):
                    stra = row['attributes']
                    invalid_pattern = r'[a-zA-Z0-9]' + re.escape(request.gene_name) + r'|' + re.escape(request.gene_name) + r'[a-zA-Z0-9]'
                    if re.search(invalid_pattern, stra):
                        continue 
                    print("tier2:", row['seqid'], row['start'], row['end'], row['strand'])
                    gene.id = row['id']
                    gene.start = row['start']
                    gene.end = row['end']
                    gene.strand = row['strand']
                    gene_strand = gene.strand
                    idx = idx + 1

                    if generalSetting.target == "exon":
                        for exon in db.children(gene.id, featuretype="exon"):
                            print("exon:", exon.id, exon.start, exon.end)
                            final_st_end.append((row['seqid'], int(exon.start), int(exon.end), idx))

                    elif generalSetting.target == "CDS": #CDS
                            for cds in db.children(gene.id, featuretype="CDS"):
                                print("CDS:", cds.id, cds.start, cds.end)
                                final_st_end.append((row['seqid'], int(cds.start), int(cds.end), idx))

                    elif generalSetting.target == "3utr": #3utr
                            for utr in db.children(gene.id, featuretype="three_prime_UTR"):
                                print("3utr:", utr.id, utr.start, utr.end)
                                final_st_end.append((row['seqid'], int(utr.start), int(utr.end), idx))
                    else: #5utr
                            for utr in db.children(gene.id, featuretype="five_prime_UTR"):
                                print("5utr:", utr.id, utr.start, utr.end)
                                final_st_end.append((row['seqid'], int(utr.start), int(utr.end), idx))
                    
                    if (generalSetting.target == "5utr" or generalSetting.target == "3utr") and len(final_st_end) == 0: #5utr, 3utr in case no UTR annotation
                        final_st_end = utr35_from_computational(request, generalSetting, casData, geneid, filename)

                    final_st_end = consensus(set(final_st_end), generalSetting.iso_union)
                    break

        else: #promoter only
            result = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=DATA_DIR,
            )
            x = result.stdout
            lists = x.splitlines()

            raw_gene_lines = [line for line in lists if "gene" in line] # da bao gom pseudogene
            gene_lines = raw_gene_lines[0]
            for raw_gene in raw_gene_lines:
                fields = raw_gene.strip().split('\t')
                stra = fields[8]
                gene_strand = fields[6]
                invalid_pattern = r'[a-zA-Z0-9]' + re.escape(request.gene_name) + r'|' + re.escape(request.gene_name) + r'[a-zA-Z0-9]'
                if re.search(invalid_pattern, stra):
                    continue
                else:
                    gene_lines = raw_gene
                    break
            
            fields = gene_lines.strip().split('\t')
            final_st_end.append((fields[0], int(fields[3]) - generalSetting.lpromoter, int(fields[3]) + generalSetting.rpromoter, 0))
        print(flag)

        cds_lines = []
        pam_size = len(PAM)

        #Tim kiem tren vung flank
        tmp_final = final_st_end.copy()
        if(generalSetting.flank):
            final_st_end = []
            for s in tmp_final:
                s = list(s)
                s[1] = max(0, int(s[1]) - generalSetting.sgRNA_len - len(PAM) + 1)
                s[2] = int(s[2]) + generalSetting.sgRNA_len + len(PAM) - 1
                final_st_end.append((s[0], s[1], s[2], 0))

        print("--------------------")
        print("dong 312", final_st_end)
        print("--------------------")

        aval = 0

        for s in final_st_end:
            seq = get_fasta_from_twobit(twobit_file, s[0], str(int(s[1]) - 1), s[2])
            l = len(seq)        

            check_id = int(s[1]) - 1 - 500
            if check_id < 0:
                check_id = 0
                aval = 1
                GAP = 0

            template = get_fasta_from_twobit(twobit_file, s[0], str(check_id), str(int(s[2]) - 1 + 1000))
            ## xu ly pam nguoc
            x = find_pam_positions(seq, REV_PAM)
            for id in x:
                if id + int(s[1]) not in cutting_sites:         # chong PAM 2 dau giong nhau
                    cutting_sites.append(id + int(s[1]))    # id + pam_size + 3 la vi tri cua PAM
                else:
                    continue

                pam_seq = seq[id: min(l - 1, id + pam_size + len_without_pam)]        # i la vi tri dau tien cua chuoi PAM
                if len(pam_seq) != (pam_size + len_without_pam):                  # neu ko du ky tu thi bo qua
                    continue     

                pam_seq = Seq(pam_seq)
                rev_comp = str(pam_seq.reverse_complement())


                if not check_5require(rev_comp, casData.requi5):
                    continue
                gcc = gc_content(rev_comp, len_without_pam)
                if gcc < generalSetting.min_gc or gcc > generalSetting.max_gc:
                    continue

                (ss, mfe) = fold_rna(rev_comp)
                (ss2, mfe2) = fold_rna(str(rev_comp + scaffold))
                microScore = primer = mlseq = lindel = None
                clea1 = clea2 = ""
                try:
                    idl = max(0, id + GAP + pam_size + 3 - len_without_pam - 10)
                    clea1 = template[idl : id + GAP + pam_size + 3]
                    clea2 = template[id + GAP + pam_size + 3 : id + GAP + pam_size + 3 + len_without_pam + 10]
                    microScore = calMicroScore(clea1, clea2)

                    spe_seq = template[id + GAP - 3 : id + GAP + 27]
                    spe_seq = Seq(spe_seq)
                    spe_seq = str(spe_seq.reverse_complement())

                    lindel = Seq(template[id + GAP + pam_size + 3 - 30 : id + GAP + pam_size + 3 + 30])
                    lindel = str(lindel.reverse_complement())

                    primer = template[id + GAP - 280: id + GAP + 261]

                except Exception as e:
                    microScore = -111111
                    spe_seq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault
                
                if len(spe_seq) != 30:
                    microScore = -999999
                    spe_seq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault

                if id + GAP + pam_size + 3 - 30 < 0 or id + GAP + pam_size + 3 + 30 > len(template):
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault

                if primer == "":
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault



                results.append({
                    "sequence": rev_comp,
                    "location": gop(s[0], str(id + int(s[1]))),
                    "strand": "-",
                    "GC Content": gcc,
                    "Self-complementary": round(mfe, 2),
                    "Primer": primer,
                    "mlseq": spe_seq.upper(),
                    "mm0": -1,
                    "mm1": 0,
                    "mm2": 0,
                    "mm3": 0,
                    "cfdScore": 0,
                    "mlScore": 0,
                    "microScore": microScore,
                    "mmejpre": clea1 + "," + clea2,
                    "Secondary structure with scaffold": f"{ss2}, ({round(mfe2, 2)} kcal/mol)",
                    "name": name,
                    "bowtie_details": "",
                    "mismatch_region":"",
                    "lindel": lindel,
                    "rs3": ""
                })
                auke.append(rev_comp)


            ## XU LY PAM XUOI
            x = find_pam_positions(seq, PAM)
            for id in x:
                if id - len_without_pam + int(s[1]) not in cutting_sites:    
                    cutting_sites.append(id - len_without_pam + int(s[1]))     # chong PAM 2 dau giong nhau
                else:
                    continue

                pam_seq = seq[max(0, id - len_without_pam):id + pam_size]         # i la vi tri dau tien cua chuoi PAM
                if len(pam_seq) != (pam_size + len_without_pam):                  # neu ko du ky tu thi bo qua
                    continue     

                if not check_5require(pam_seq, casData.requi5):
                    continue
                gcc = gc_content(pam_seq, len_without_pam)
                if gcc < generalSetting.min_gc or gcc > generalSetting.max_gc:
                    continue

                (ss, mfe) = fold_rna(pam_seq)
                (ss2, mfe2) = fold_rna(pam_seq + scaffold)
                microScore = primer = mlseq = lindel = None
                clea1 = clea2 = ""
                try:
                    idl = max(id + GAP - len_without_pam - 3 - 10, 0)
                    clea1 = template[idl: id + GAP - 3]
                    clea2 = template[id + GAP - 3: id + GAP - 3 + len_without_pam + 10]
                    microScore = calMicroScore(clea1, clea2)

                    primer = template[id + GAP - 280: id + GAP + 261]
                    mlseq = template[id + GAP - 20 - 4: id + GAP - 19 + 25]

                    lindel = template[id + GAP - 3 - 30: id + GAP + pam_size + 30]

                except Exception as e:
                    microScore = -111111
                    mlseq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault
                if len(spe_seq) != 30:
                    microScore = -999999
                    mlseq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault

                if id + GAP - 3 - 30 < 0 or id + GAP + pam_size + 30 > len(template):
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault
                if primer == "":
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault

                results.append({
                        "sequence": pam_seq,
                        "location": gop(s[0], str(id - len_without_pam + int(s[1]))),
                        "strand": "+",
                        "GC Content": gcc,
                        "Self-complementary": round(mfe,2),
                        "Primer": primer,
                        "mlseq": mlseq.upper(),
                        "mm0": -1,
                        "mm1": 0,
                        "mm2": 0,
                        "mm3": 0,
                        "cfdScore": 0,
                        "mlScore": 0,
                        "microScore": microScore,
                        "mmejpre": clea1 + "," + clea2,
                        "Secondary structure with scaffold": f"{ss2}, ({round(mfe2, 2)} kcal/mol)",
                        "name": name,
                        "bowtie_details": "",
                        "mismatch_region":"",
                        "lindel": lindel,
                        "rs3": "",
                    })
                auke.append(pam_seq)
        print("truoc truoc", len(results))
        
        idd = save_sgRNA_list_dbv(idd, results, gene_name, spec, PAM, len_without_pam, "gene_name",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="calculating-index_processing", log="indexing", gene_strand=gene_strand)


        if len(results) == 0:
            idd = save_sgRNA_list_dbv(idd, results, gene_name, spec, PAM, len_without_pam, "gene_name",
                    q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="no_result", log="No result available, check your gene name or region, it must match with the genome", stage=0)
            
            return
        print("truoc", len(results))
        indexComputing_dbv(idd, casData.off_target, casData.mismatch_num)
                
        idd = save_sgRNA_list_dbv(idd, results, gene_name, spec, PAM, len_without_pam, "gene_name",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="success", log="done", stage=0, gene_strand=gene_strand)
        
        print("sau", len(results))
    except Exception as e:
        print(f"Error in GeneNameComputing: {str(e)}")

        save_sgRNA_list_dbv(
            idd, [], request.gene_name, request.species, 
            casData.pam, generalSetting.sgRNA_len, "gene_name",
            primerConfigData.min_product_size, primerConfigData.max_product_size,
            primerConfigData.min_primer_size, primerConfigData.max_primer_size,
            primerConfigData.optimal_primer_size, primerConfigData.min_tm,
            primerConfigData.max_tm, primerConfigData.optimal_tm,
            queue_task_id, status="failed", log=f"Processing error: {str(e)}"
        )
        raise       
    return
def CoordinateComputing(queue_task_id, idd, request, generalSetting, casData, primerConfigData):
    print("Da vao ham CoordinateComputing")
    GAP = CONST_GAP
    try:
        request = CoordinateEntry(**request)
        generalSetting = GeneralSetting(**generalSetting)
        casData = cas9Setting(**casData)
        primerConfigData = PrimerSetting(**primerConfigData)

        print(primerConfigData)
        q1 = primerConfigData.min_product_size
        q2 = primerConfigData.max_product_size
        q3 = primerConfigData.min_primer_size
        q4 = primerConfigData.max_primer_size
        q5 = primerConfigData.optimal_primer_size
        q6 = primerConfigData.min_tm
        q7 = primerConfigData.max_tm
        q8 = primerConfigData.optimal_tm

        results = []
        auke = []
        cutting_sites = []
        query = request.coordinate
        spec = request.species

        PAM = casData.pam
        len_without_pam = generalSetting.sgRNA_len
        REV_PAM = Seq(PAM)
        REV_PAM = str(REV_PAM.reverse_complement())
        scaffold = generalSetting.scaffoldSeq
        request_strand = "+"
        
        idd = save_sgRNA_list_dbv(idd, results, query, spec, PAM, len_without_pam, "coordinate",
                        q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="Finding", log="Finding sgRNA candidates")
                                

        resolved_paths = resolve_spec(spec)
        twobit_file = resolved_paths["twobit"]

        final_st_end = query.split(',')
        final_st_end = [x.split(':') for x in final_st_end]
        final_st_end = [(x[0], int(x[1].split('-')[0]), int(x[1].split('-')[1]) + 1, idx + 1) for idx, x in enumerate(final_st_end)]
        pam_size = len(PAM)
        print(final_st_end)
        
        tmp_final = final_st_end.copy()
        if(generalSetting.flank):
            final_st_end = []
            for s in tmp_final:
                s = list(s)
                s[1] = max(0, int(s[1]) - generalSetting.sgRNA_len - len(PAM) + 1)
                s[2] = int(s[2]) + generalSetting.sgRNA_len + len(PAM) - 1
                final_st_end.append((s[0], s[1], s[2], 0))

        for s in final_st_end:
            print(s[0])
            print(s[1])
            print(s[2])

            if (int(s[2]) - int(s[1]) > 30000):
                raise Exception(f"The maximum bp allowed per query is 30000 bp")

            seq = get_fasta_from_twobit(twobit_file, s[0], (int(s[1]) - 1), int(s[2]))                
            l = len(seq)        
            
            template = None
            try:
                check_id = int(s[1]) - 1 - 500
                if check_id < 0:
                    check_id = 0
                    GAP = 0
                template = get_fasta_from_twobit(twobit_file, s[0], str(check_id), str(int(s[2]) - 1 + 1000))
                print("do dai template la", len(template))
            except Exception as e:
                print("checkpoint", "out of template")
                template = seq
                GAP = 0

            ## xu ly pam nguoc
            x = find_pam_positions(seq, REV_PAM)
            for id in x:
                if id + int(s[1]) not in cutting_sites:         # chong PAM 2 dau giong nhau
                    cutting_sites.append(id + int(s[1]))    # id + pam_size + 3 la vi tri cua PAM
                else:
                    continue

                pam_seq = seq[id: min(l - 1, id + pam_size + len_without_pam)]        # i la vi tri dau tien cua chuoi PAM
                if len(pam_seq) != (pam_size + len_without_pam):                  # neu ko du ky tu thi bo qua
                    continue     

                pam_seq = Seq(pam_seq)
                rev_comp = str(pam_seq.reverse_complement())


                if not check_5require(rev_comp, casData.requi5):
                    continue
                gcc = gc_content(rev_comp, len_without_pam)
                if gcc < generalSetting.min_gc or gcc > generalSetting.max_gc:
                    continue

                (ss, mfe) = fold_rna(rev_comp)
                (ss2, mfe2) = fold_rna(str(rev_comp + scaffold))

                microScore = primer = mlseq = lindel = None
                clea1 = clea2 = ""
                try:
                    idl = id + GAP + pam_size + 3 - len_without_pam - 10
                    clea1 = template[idl : id + GAP + pam_size + 3]
                    clea2 = template[id + GAP + pam_size + 3 : id + GAP + pam_size + 3 + len_without_pam + 10]
                    microScore = calMicroScore(clea1, clea2)

                    spe_seq = template[id + GAP - 3 : id + GAP + 27]
                    spe_seq = Seq(spe_seq)
                    spe_seq = str(spe_seq.reverse_complement())

                    lindel = Seq(template[id + GAP + pam_size + 3 - 30 : id + GAP + pam_size + 3 + 30])
                    lindel = str(lindel.reverse_complement())

                    primer = template[id + GAP - 260 : id + GAP + 281]

                except Exception as e:
                    traceback.print_exc()
                    print("checkpoint", e, id + GAP + 281)
                    microScore = -999999
                    spe_seq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault
                
                if len(spe_seq) != 30:
                    print("checkpoint", pam_seq)
                    microScore = -999999
                    spe_seq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault
                if id + GAP + pam_size + 3 - 30 < 0 or id + GAP + pam_size + 3 + 30 > len(template):
                    print("checkpoint", pam_seq)
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault

                if primer == "":
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault
                

                results.append({
                    "sequence": rev_comp,
                    "location": gop(s[0], str(id + int(s[1]))),
                    "strand": "-",
                    "GC Content": gcc,
                    "Self-complementary": round(mfe, 2),
                    "Primer": primer,
                    "mlseq": spe_seq.upper(),
                    "mm0": -1,
                    "mm1": 0,
                    "mm2": 0,
                    "mm3": 0,
                    "cfdScore": 0,
                    "mlScore": 0,
                    "microScore": microScore,
                    "mmejpre": clea1 + "," + clea2,
                    "Secondary structure with scaffold": f"{ss2}, ({round(mfe2, 2)} kcal/mol)",
                    "name": request.coordinate + "," + spec + "," + PAM,
                    "bowtie_details": "",
                    "mismatch_region":"",
                    "lindel": lindel,
                    "rs3": ""
                })
                auke.append(rev_comp)


            ## XU LY PAM XUOI
            x = find_pam_positions(seq, PAM)
            for id in x:
                if id - len_without_pam + int(s[1]) not in cutting_sites:    
                    cutting_sites.append(id - len_without_pam + int(s[1]))     # chong PAM 2 dau giong nhau
                else:
                    continue

                pam_seq = seq[max(0, id - len_without_pam):id + pam_size]         # i la vi tri dau tien cua chuoi PAM
                if len(pam_seq) != (pam_size + len_without_pam):                  # neu ko du ky tu thi bo qua
                    continue     

                if not check_5require(pam_seq, casData.requi5):
                    continue
                gcc = gc_content(pam_seq, len_without_pam)
                if gcc < generalSetting.min_gc or gcc > generalSetting.max_gc:
                    continue

                (ss, mfe) = fold_rna(pam_seq)
                (ss2, mfe2) = fold_rna(pam_seq + scaffold)
                microScore = primer = mlseq = lindel = None
                clea1 = clea2 = ""
                try:
                    idl = id + GAP - len_without_pam - 3 - 10
                    clea1 = template[idl: id + GAP - 3]
                    clea2 = template[id + GAP - 3: id + GAP - 3 + len_without_pam + 10]
                    microScore = calMicroScore(clea1, clea2)

                    primer = template[id + GAP - 280: id + GAP + 261]
                    mlseq = template[id + GAP - 20 - 4: id + GAP - 19 + 25]

                    lindel = template[id + GAP - 3 - 30: id + GAP + pam_size + 30]

                except Exception as e:
                    traceback.print_exc()
                    print("checkpoint", e, id + GAP + 261)
                    microScore = -999999
                    mlseq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault
                if len(mlseq) != 30:
                    print("checkpoint", mlseq)
                    microScore = -999999
                    mlseq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault
                if id + GAP - 3 - 30 < 0 or id + GAP + pam_size + 30 > len(template):
                    print("checkpoint", id + GAP)
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault

                if primer == "":
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault

                results.append({
                        "sequence": pam_seq,
                        "location": gop(s[0], str(id - len_without_pam + int(s[1]))),
                        "strand": "+",
                        "GC Content": gcc,
                        "Self-complementary": round(mfe,2),
                        "Primer": primer,
                        "mlseq": mlseq.upper(),
                        "mm0": -1,
                        "mm1": 0,
                        "mm2": 0,
                        "mm3": 0,
                        "cfdScore": 0,
                        "mlScore": 0,
                        "microScore": microScore,
                        "mmejpre": clea1 + "," + clea2,
                        "Secondary structure with scaffold": f"{ss2}, ({round(mfe2, 2)} kcal/mol)",
                        "name": request.coordinate + "," + spec + "," + PAM,
                        "bowtie_details": "",
                        "mismatch_region":"",
                        "lindel": lindel,
                        "rs3": "",
                    })
                auke.append(pam_seq)

        idd = save_sgRNA_list_dbv(idd, results, query, spec, PAM, len_without_pam, "coordinate",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="calculating-index_processing", log="indexing")

        if len(results) == 0:
            idd = save_sgRNA_list_dbv(idd, results, query, spec, PAM, len_without_pam, "coordinate",
                    q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="no_result", log="No result available, check your gene name or region", stage=0)

            return

        indexComputing_dbv(idd, casData.off_target, casData.mismatch_num)
        
        idd = save_sgRNA_list_dbv(idd, results, query, spec, PAM, len_without_pam, "coordinate",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="success", log="done", stage=0)

    except Exception as e:
        print(f"Error in CoordinateComputing: {str(e)}")
        save_sgRNA_list_dbv(
            idd, [], request.coordinate, request.species, 
            casData.pam, generalSetting.sgRNA_len, "coordinate",
            primerConfigData.min_product_size, primerConfigData.max_product_size,
            primerConfigData.min_primer_size, primerConfigData.max_primer_size,
            primerConfigData.optimal_primer_size, primerConfigData.min_tm,
            primerConfigData.max_tm, primerConfigData.optimal_tm,
            queue_task_id, status="failed", log=f"Processing error: {str(e)}"
        )
        raise 
    return
    

def FastaComputing(queue_task_id, idd, request, generalSetting, casData, primerConfigData):
    print("Da vao ham FastaComputing")
    GAP = CONST_GAP
    seq = ""
    try:
        request = FastaEntry(**request)
        generalSetting = GeneralSetting(**generalSetting)
        casData = cas9Setting(**casData)
        primerConfigData = PrimerSetting(**primerConfigData)

        print(primerConfigData)
        q1 = primerConfigData.min_product_size
        q2 = primerConfigData.max_product_size
        q3 = primerConfigData.min_primer_size
        q4 = primerConfigData.max_primer_size
        q5 = primerConfigData.optimal_primer_size
        q6 = primerConfigData.min_tm
        q7 = primerConfigData.max_tm
        q8 = primerConfigData.optimal_tm

        results = []
        auke = []
        didang = []
        cutting_sites = []
        spec = request.species

        records = request.dna_seq.split('>')

        header = records[1].split('\n')[0]
        header = header.split(',')[0]
        header = header.split(' ')[0]
        seq = re.sub(r"\s+", "", records[1].split('\n')[1]).upper()
        

        print(seq)

        PAM = casData.pam
        len_without_pam = generalSetting.sgRNA_len
        REV_PAM = Seq(PAM)
        REV_PAM = str(REV_PAM.reverse_complement())
        scaffold = generalSetting.scaffoldSeq
        

        pam_size = len(PAM)
        l = len(seq)        
        template = seq
        
        idd = save_sgRNA_list_dbv(idd, results, header, spec, PAM, len_without_pam, "fasta",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="Finding", log="Finding sgRNA candidates")
        
        x = find_pam_positions(seq, REV_PAM)
        for id in x:
            if id not in cutting_sites:         # chong PAM 2 dau giong nhau
                cutting_sites.append(id) 
            else:
                continue

            pam_seq = seq[id: min(l - 1, id + pam_size + len_without_pam)]        # i la vi tri dau tien cua chuoi PAM
            if len(pam_seq) != (pam_size + len_without_pam):                  # neu ko du ky tu thi bo qua
                continue     

            pam_seq = Seq(pam_seq)
            rev_comp = str(pam_seq.reverse_complement())


            if not check_5require(rev_comp, casData.requi5):
                continue
            gcc = gc_content(rev_comp, len_without_pam)
            if gcc < generalSetting.min_gc or gcc > generalSetting.max_gc:
                continue

            (ss, mfe) = fold_rna(rev_comp)
            (ss2, mfe2) = fold_rna(str(rev_comp + scaffold))

            clea1 = clea2 = ""
            try:
                idl = id + pam_size + 3 - len_without_pam - 10
                clea1 = template[idl: id + pam_size + 3]
                clea2 = template[id + pam_size + 3 : id + pam_size + 3 + len_without_pam + 10]
                microScore = calMicroScore(clea1, clea2)

                spe_seq = template[id - 3 : id + 27]
                spe_seq = Seq(spe_seq)
                spe_seq = str(spe_seq.reverse_complement())

                lindel = Seq(template[id + pam_size + 3 - 30 : id + pam_size + 3 + 30])
                lindel = str(lindel.reverse_complement())

                primer = template[id - 260 : id + 281]

            except Exception as e:
                microScore = -999999
                spe_seq = mlseqDefault
                lindel = lindelDefault
                primer = primerDefault
            
            if len(spe_seq) != 30:
                microScore = -999999
                spe_seq = mlseqDefault
                lindel = lindelDefault
                primer = primerDefault
            

            results.append({
                "sequence": rev_comp,
                "location": gop(header, str(id)),
                "strand": "-",
                "GC Content": gcc,
                "Self-complementary": round(mfe, 2),
                "Primer": primer,
                "mlseq": spe_seq.upper(),
                "mm0": -1,
                "mm1": 0,
                "mm2": 0,
                "mm3": 0,
                "cfdScore": 0,
                "mlScore": 0,
                "microScore": microScore,
                "mmejpre": clea1 + "," + clea2,
                "Secondary structure with scaffold": f"{ss2}, ({round(mfe2, 2)} kcal/mol)",
                "name": header + "," + spec + "," + PAM,
                "bowtie_details": "",
                "mismatch_region":"",
                "lindel": lindel,
                "rs3": ""
            })
            auke.append(rev_comp)


        ## XU LY PAM XUOI
        x = find_pam_positions(seq, PAM)
        for id in x:
            if id - len_without_pam not in cutting_sites:    
                cutting_sites.append(id - len_without_pam)     # chong PAM 2 dau giong nhau
            else:
                continue

            pam_seq = seq[max(0, id - len_without_pam):id + pam_size]         # i la vi tri dau tien cua chuoi PAM
            if len(pam_seq) != (pam_size + len_without_pam):                  # neu ko du ky tu thi bo qua
                continue     

            if not check_5require(pam_seq, casData.requi5):
                continue
            gcc = gc_content(pam_seq, len_without_pam)
            if gcc < generalSetting.min_gc or gcc > generalSetting.max_gc:
                continue

            (ss, mfe) = fold_rna(pam_seq)
            (ss2, mfe2) = fold_rna(pam_seq + scaffold)


            clea1 = clea2 = microScore = primer = mlseq = lindel = ""
            try:
                idl = id - len_without_pam - 3 - 10
                clea1 = template[idl: id - 3]
                clea2 = template[id - 3: id - 3 + len_without_pam + 10]
                microScore = calMicroScore(clea1, clea2)

                primer = template[id - 280: id + 261]
                mlseq = template[id - 20 - 4: id - 19 + 25]

                lindel = template[id - 3 - 30: id + pam_size + 30]

            except Exception as e:
                microScore = -999999
                mlseq = mlseqDefault
                lindel = lindelDefault
                primer = primerDefault
            if len(mlseq) != 30:
                microScore = -999999
                mlseq = mlseqDefault
                lindel = lindelDefault
                primer = primerDefault
            results.append({
                    "sequence": pam_seq,
                    "location": gop(header, str(id - len_without_pam)),
                    "strand": "+",
                    "GC Content": gcc,
                    "Self-complementary": round(mfe,2),
                    "Primer": primer,
                    "mlseq": mlseq,
                    "mm0": -1,
                    "mm1": 0,
                    "mm2": 0,
                    "mm3": 0,
                    "cfdScore": 0,
                    "mlScore": 0,
                    "microScore": microScore,
                    "mmejpre": clea1 + "," + clea2,
                    "Secondary structure with scaffold": f"{ss2}, ({round(mfe2, 2)} kcal/mol)",
                    "name": header + "," + spec + "," + PAM,
                    "bowtie_details": "",
                    "mismatch_region":"",
                    "lindel": lindel,
                    "rs3": "",
                })
            auke.append(pam_seq)

        print(primerConfigData)
        
        idd = save_sgRNA_list_dbv(idd, results, header, spec, PAM, len_without_pam, "fasta",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="calculating-index_processing", log="indexing", gene_strand="+")

        if len(results) == 0:            
            idd = save_sgRNA_list_dbv(idd, results, header, spec, PAM, len_without_pam, "fasta",
                    q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="no_result", log="No result available, check your gene name or region", stage=0)

            return
        indexComputing_dbv(idd, casData.off_target, casData.mismatch_num)
        idd = save_sgRNA_list_dbv(idd, results, header, spec, PAM, len_without_pam, "fasta",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="success", log="done", stage=0, gene_strand="+")

    
    except Exception as e:
        print(f"Error in GeneNameComputing: {str(e)}")

        save_sgRNA_list_dbv(
            idd, [], header, request.species, 
            casData.pam, generalSetting.sgRNA_len, "fasta",
            primerConfigData.min_product_size, primerConfigData.max_product_size,
            primerConfigData.min_primer_size, primerConfigData.max_primer_size,
            primerConfigData.optimal_primer_size, primerConfigData.min_tm,
            primerConfigData.max_tm, primerConfigData.optimal_tm,
            queue_task_id, status="failed", log=f"Processing error: {str(e)}"
        )
        raise 
    return


def getMMDT_dbv(name: str, idfile: str, pos_list: list, bowtiedata: list):
    """
    Gộp getMMRegion + getMMDetails:
    - Gán nhãn exon/intron/intergenic cho mismatch regions
    - Trả về list[(stt, mismatch_region)] để bulk insert vào DB
    
    Args:
        name: Tên genome (để tìm file .gff3/.gff/.gtf)
        idfile: ID của task
        pos_list: List[(chrom, start, end, idseq)] từ bowtie output
    
    Returns:
        List[(stt, mismatch_region)] - dùng để bulk insert
    """    

    if name.startswith("hash:"):
        _, fasta_h, gff3_h = name.split(":")
        # Thư mục chứa annotation: D:\UET\ViCRISPR\backend\app\data\anno_7de6b...
        current_anno_dir = BASE_DATA_DIR / f"anno_{gff3_h}"
        # Tên file gốc bên trong thư mục đó: 7de6b...
        file_base_name = gff3_h
    else:
        current_anno_dir = BASE_DATA_DIR
        file_base_name = name

    NEW_DATA_DIR = str(current_anno_dir)

    raw_bed_file = os.path.join(NEW_DATA_DIR, f"{idfile}_raw.bed")
    bed_file = os.path.join(NEW_DATA_DIR, f"{idfile}_sorted.bed")
    
    possible_ext = [".gff3", ".gff", ".gtf"]
    found_ext = None
    
    for ext in possible_ext:
        file_path = os.path.join(NEW_DATA_DIR, f"{file_base_name}{ext}")
        if os.path.exists(file_path):
            found_ext = ext
            break
    
    if not found_ext:
        print(f"Không tìm thấy file annotation: {name}.gff3 / .gff / .gtf")
        return None
    
    name = file_base_name
    
    exon_file = os.path.join(NEW_DATA_DIR, f"{name}_exons.sorted{found_ext}")
    gene_file = os.path.join(NEW_DATA_DIR, f"{name}_genes.sorted{found_ext}")
    result_file = os.path.join(NEW_DATA_DIR, f"{idfile}_mm_annotation.bed")
    
    grouped = defaultdict(list)
    for chrom, start, end, idseq in pos_list:
        grouped[idseq].append((chrom, start, end))
    
    with open(raw_bed_file, "w") as f:
        for idseq, regions in grouped.items():
            for chrom, start, end in regions:
                f.write(f"{chrom}\t{start}\t{end}\t{idseq}\n")
    
    cmds = [
        f"sort -k1,1V -k2,2n {raw_bed_file} > {bed_file}",
        
        # A. Gán nhãn exon
        f"bedtools intersect -a {bed_file} -b {exon_file} -wa | "
        f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $1, $2, $3, $4, \"exon\"}}' > annotated.part1.exonic",
        
        # B. Lấy vùng còn lại
        f"bedtools intersect -a {bed_file} -b {exon_file} -v > remaining_regions.1.bed",
        
        # C. Gán nhãn intron
        f"bedtools intersect -a remaining_regions.1.bed -b {gene_file} -wa | "
        f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $1, $2, $3, $4, \"intron\"}}' > annotated.part2.intronic",
        
        # D. Lấy vùng còn lại (không exon, không intron)
        f"bedtools intersect -a remaining_regions.1.bed -b {gene_file} -v > remaining_regions.2.bed",
        
        # E. Gán nhãn intergenic
        f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $1, $2, $3, $4, \"intergenic\"}}' remaining_regions.2.bed > annotated.part3.intergenic",
        
        # F. Gộp & sắp xếp & dọn file tạm
        f"cat annotated.part1.exonic annotated.part2.intronic annotated.part3.intergenic | "
        f"sort -k4,4V -k2,2n -u > {result_file} && "
        f"rm remaining_regions.*.bed annotated.part* && "
        f"rm {bed_file}"
    ]
    
    for cmd in cmds:
        try:
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", cwd=NEW_DATA_DIR)
        except subprocess.CalledProcessError as e:
            print(f"Lỗi khi chạy: {cmd}\nChi tiết: {e}")
            return None
    
    mm_details = {}
    with open(result_file, 'r', encoding='utf-8') as bed_file_handle:
        for line in bed_file_handle:
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            chrom = fields[0]
            start = int(fields[1])
            region_type = fields[4]
            
            key = f"{chrom}:{start}"
            mm_details[key] = region_type
    
    
    result_list = []


    for entry in bowtiedata:
        
        if not entry or entry.strip() == "":
            result_list.append("")
            continue
        
        regions = []
        offtargets = entry.strip().rstrip(';').split(';')
        
        print("ggom", offtargets)
        for offtarget in offtargets:
            offtarget = offtarget.strip()
            if not offtarget:
                continue
            
            # Format: "NC_000017.11:79479867,,2:C>A,4:C>T,..."
            chrom_pos = offtarget.split(',,')[0].strip()
            region = mm_details.get(chrom_pos, "unknown")
            regions.append(region)
        
        mismatch_region = ";".join(regions)
        result_list.append(mismatch_region)
    
    print("checkpoint")

    try:
        os.remove(result_file)
        os.remove(raw_bed_file)
    except:
        pass
    
    print(f"Hoàn tất! Trả về {len(result_list)} records với mismatch regions")
    return result_list

def indexComputing_dbv(idfile: str, off_target: bool = 0, num_of_mismatches: int = 3):
    db = SessionLocal()
    try:
        task = db.query(TaskMetadata).filter(TaskMetadata.query_id == idfile).first()
        if not task:
            raise ValueError(f"Task with id {idfile} not found in database")
        
        sgrna_records = db.query(Sgrna).filter(Sgrna.query_id == idfile).all()
        if not sgrna_records:
            raise ValueError(f"No sgRNA records found for task {idfile}")
        
        datafile = []
        for record in sgrna_records:
            row = {
                "stt": record.stt,
                "sequence": record.sequence,
                "location": record.location,
                "strand": record.strand,
                "GC Content": record.gc_content,
                "Self-complementary": record.self_complementary,
                "Primer": record.primer,
                "mlseq": record.mlseq,
                "lindel": record.lindel,
                "mm0": record.mm0 or 0,
                "mm1": record.mm1 or 0,
                "mm2": record.mm2 or 0,
                "mm3": record.mm3 or 0,
                "cfdScore": record.cfd_score or 0,
                "mlScore": record.ml_score or None,
                "microScore": record.micro_score or None,
                "rs3": record.rs3_score or None,
                "mmejpre": record.mmej_pre,
                "Secondary structure with scaffold": record.sec_structure,
                "bowtie_details": ""
            }
            datafile.append(row)
        
        pam_name = task.pam
        resolved_paths = resolve_spec(task.spec)
        bowtie_index_file = resolved_paths["bowtie_index"]

        print(f"Resolved bowtie index file: {bowtie_index_file}")
        spec_name = task.spec
        sgRNA_len = task.sgrna_len
        
        seq_list = [row["sequence"] for row in datafile]
        seq_list_ml = [row["mlseq"] for row in datafile]
        
        write_sgrna_to_fasta_with_IUPAC(seq_list, pam_name, idfile)

        ml_score = [999] * len(datafile)
        rs3_score = [999] * len(datafile)
        
        if pam_name == "NGG":
            ml_score = get_ml_score(seq_list_ml)
            rs3_score = get_ml_score_azi3(seq_list_ml)
        
        pol = 0
        limit_num = 1000
        ss_map = {19: 5000, 20: 4000, 21: 3000, 22: 2000, 23: 1000}
        ss = ss_map.get(len(pam_name) + sgRNA_len, 700)
        pp = count_permu_IUPAC(pam_name)
        limit_num = max(20, ss / pp)
        sg_file = f"{idfile}_sgrna_output.fa"
        
        para_mm = 3
        if off_target == 0:
            para_mm = num_of_mismatches

        bowtie_index_file = os.path.relpath(bowtie_index_file, DATA_DIR)

        command = [
            "bowtie",
            "-v", str(para_mm),
            "-k", str(limit_num),
            "-f",
            "-x", bowtie_index_file,
            sg_file,
            "/dev/stdout"
        ]
        
        process = subprocess.Popen(
            command,
            cwd=DATA_DIR,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1
        )
        
        pos_list = []
        count_dict = {}
        max_hits_per_id = 5
        
        danh_dau = set()
        vcl = []
        m = defaultdict(int)
        ll = len(datafile[0]['sequence']) - len(pam_name)
        
        for line in process.stdout:
            pol = pol + 1
            line = line.strip()
            
            idseq, ttmm = xuly(line, datafile, pam_name, ll, off_target, num_of_mismatches)
            id_in_fa = int(int(line.strip().split('\t')[0]))
            m[id_in_fa] += 1
            
            if m[id_in_fa] >= limit_num:
                danh_dau.add(id_in_fa / pp)
                vcl.append(id_in_fa)
            
            if idseq == -1:
                continue
            
            datafile[idseq]["bowtie_details"] += ",".join(map(str, ttmm)) + "; "
            
            chrom = line.strip().split('\t')[2]
            start = int(line.strip().split('\t')[3])
            end = start + sgRNA_len - 1
            pos_list.append((chrom, start, end, idseq))
            
            x, scr = get_cfd_score(line, pam_name, ll)
            if x == -1:
                continue
            datafile[x]["cfdScore"] += scr
            
            print(str(x) + " " + str(datafile[x]["cfdScore"]))
        
        process.stdout.close()
        process.wait()
        
        # this mean it break the threshold, we will save in db that -a mean >= a
        for i, row in enumerate(datafile):
            if i in danh_dau:
                if str(row.get("mm3", 0)) != 0:
                    #row["mm3"] = f">={row['mm3']}"
                    row["mm3"] = -row["mm3"]
                elif str(row.get("mm2", 0)) != 0:
                    #row["mm2"] = f">={row['mm2']}"
                    row["mm2"] = -row["mm2"]
                elif str(row.get("mm1", 0)) != 0:
                    #row["mm1"] = f">={row['mm1']}"
                    row["mm1"] = -row["mm1"]
            print(
                f"{i}: {row.get('sequence', '')}, "
                f"location={row.get('location', '')}, "
                f"mm0={row.get('mm0', 0)}, "
                f"mm1={row.get('mm1', 0)}, "
                f"mm2={row.get('mm2', 0)}, "
                f"mm3={row.get('mm3', 0)}"
            )
        
        print(datafile)
        print(ml_score)

        for i in range(len(datafile)):
            datafile[i]["mlScore"] = ml_score[i]
            datafile[i]["rs3"] = rs3_score[i]
            datafile[i]["cfdScore"] = round(100 / (100 + datafile[i]["cfdScore"]), 2)
        
        grouped = defaultdict(list)
        for chrom, start, end, idseq in pos_list:
            grouped[idseq].append((chrom, start, end))
            
        
        parent_dir = os.path.dirname(resolved_paths["gff3"])
        raw_bed_dir = os.path.join(parent_dir, f"{idfile}_raw.bed")
        with open(raw_bed_dir, "w") as f:
            for idseq, regions in grouped.items():
                for chrom, start, end in regions:
                    f.write(f"{chrom}\t{start}\t{end}\t{idseq}\n")


        bowtie_data = []
        for idx, row in enumerate(datafile):
            bowtie_data.append(row.get("bowtie_details"))
        mm_results = getMMDT_dbv(spec_name, idfile, pos_list, bowtie_data)


        if mm_results:
            for idx, row in enumerate(datafile):
                row['mismatch_region'] = mm_results[idx]

        sgrna_updates = []



        for idx, row in enumerate(datafile):
            sgrna_updates.append({
                "stt": idx + 1,
                "query_id": idfile,
                "sequence": row.get("sequence"),
                "location": row.get("location"),
                "strand": row.get("strand"),
                "gc_content": row.get("GC Content"),
                "self_complementary": row.get("Self-complementary"),
                "primer": str(row.get("Primer")),
                "mlseq": row.get("mlseq"),
                "mm0": row.get("mm0"),
                "mm1": row.get("mm1"),
                "mm2": row.get("mm2"),
                "mm3": row.get("mm3"),
                "cfd_score": row.get("cfdScore"),
                "ml_score": row.get("mlScore"),
                "micro_score": row.get("microScore"),
                "rs3_score": safe_float(row.get("rs3")),
                "mmej_pre": str(row.get("mmejpre")),
                "sec_structure": str(row.get("Secondary structure with scaffold")),
                "lindel": str(row.get("lindel")),
                "mismatch_region": row.get("mismatch_region", ""),
                "bowtie_details": row.get("bowtie_details", "")
            })
        
        db.query(Sgrna).filter(Sgrna.query_id == idfile).delete()
        db.bulk_insert_mappings(Sgrna, sgrna_updates, return_defaults=False)
        
        db.commit()
        
        print("Da tinh toan xong")
        print(limit_num)
        print(pol)
        print(ll)
        
        checkAndSendMail(idfile)
        
        return
        
    except Exception as e:
        db.rollback()
        print(f"Error in indexComputing_dbv: {str(e)}")
        try:
            task = db.query(TaskMetadata).filter(TaskMetadata.query_id == idfile).first()
            if task:
                task.status = "error"
                task.log = str(e)
                task.completed_at = func.now()
                db.commit()
        except:
            pass
        raise
    finally:
        db.close()

def checkAndSendMail(idfile: str):
    db = SessionLocal()
    #trong ham tinh toan index
    #check xem co phai gui mail ko, co thi tim trong db idfile va gui toi mailist tuong ung 
    try:
        print(f"{idfile} starting check mail in queue")

        email_records = db.query(EmailQueue).filter(
            EmailQueue.idfile == idfile
        ).all()

        if not email_records:
            print(f"{idfile} not found any mails")
            return

        mail_list = [record.email for record in email_records]
        
        print(f"found {len(mail_list)} email in queue")

        success_count = 0

        sendMail(idfile, mail_list)

        delete_count = db.query(EmailQueue).filter(
            EmailQueue.idfile == idfile
        ).delete(synchronize_session='fetch')

        db.commit()
        
        print(f"successfully sent {idfile}:")
        print(f"mail num: {len(mail_list)}")
        print(f"mail successfully sent: {success_count}")

    except Exception as e:
        print(f"error {idfile}: {e}")

    finally:
        db.close()