from email import header
from celery import Celery
import time
from pydantic import BaseModel
import httpx
import subprocess
import time, gffutils

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

        idd = save_sgRNA_list(idd, results, gene_name, spec, PAM, len_without_pam, "gene_name",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="Finding", log="Finding sgRNA candidates")

        filename = getAnnotationFile(spec)

        if filename is None:
            raise FileNotFoundError(
                f"Không tìm thấy {spec}.gff3 / {spec}.gtf / {spec}.gff trong CSDL"
            )

        command = f'grep "{gene_name}" {filename}'
        name = gene_name + "," + spec + "," + PAM
        twobit_file = spec + ".2bit"
        cutting_sites = []

        db_path = os.path.join(DATA_DIR, f"{spec}.db")
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
                print("da tim trong day roi nhung ko co")
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

            template = get_fasta_from_twobit(twobit_file, s[0], str(check_id), str(int(s[2]) - 1 + 500))
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
        idd = save_sgRNA_list(idd, results, gene_name, spec, PAM, len_without_pam, "gene_name",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="calculating-index_processing", log="indexing", gene_strand=gene_strand)
        
        if len(results) == 0:
            idd = save_sgRNA_list(idd, results, gene_name, spec, PAM, len_without_pam, "gene_name",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="no_result", log="No result available, check your gene name or region", stage=0)
            return
        print("truoc", len(results))
        indexComputing(idd, casData.off_target, casData.mismatch_num)
        idd = save_sgRNA_list(idd, results, gene_name, spec, PAM, len_without_pam, "gene_name",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="success", log="done", stage=0, gene_strand=gene_strand)
        print("sau", len(results))
    except Exception as e:
        print(f"Error in GeneNameComputing: {str(e)}")
        save_sgRNA_list(
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
        
        idd = save_sgRNA_list(idd, results, query, spec, PAM, len_without_pam, "coordinate",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="Finding", log="Finding sgRNA candidates")


        twobit_file = spec + ".2bit"

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
        aval = 0
        for s in final_st_end:
            print(s[0])
            print(s[1])
            print(s[2])
            seq = get_fasta_from_twobit(twobit_file, s[0], (int(s[1]) - 1), int(s[2]))                
            l = len(seq)        
            
            template = None
            try:
                check_id = int(s[1]) - 1 - 500
                if check_id < 0:
                    check_id = 0
                    aval = 1
                template = get_fasta_from_twobit(twobit_file, s[0], str(check_id), str(int(s[2]) - 1 + 500))
            except Exception as e:
                template = seq

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
                    microScore = -999999
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
                    clea2 = template[id + GAP - 3: id + GAP - 3 + len_without_pam]
                    microScore = calMicroScore(clea1, clea2)

                    primer = template[id + GAP - 280: id + GAP + 261]
                    mlseq = template[id + GAP - 20 - 4: id + GAP - 19 + 25]

                    lindel = template[id + GAP - 3 - 30: id + GAP + pam_size + 30]

                except Exception as e:
                    microScore = -999999
                    mlseq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault
                if len(mlseq) != 30 or aval == 1:
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
                        "name": request.coordinate + "," + spec + "," + PAM,
                        "bowtie_details": "",
                        "mismatch_region":"",
                        "lindel": lindel,
                        "rs3": "",
                    })
                auke.append(pam_seq)

        idd = save_sgRNA_list(idd, results, query, spec, PAM, len_without_pam, "coordinate",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="calculating-index_processing", log="indexing")

        if len(results) == 0:
            idd = save_sgRNA_list(idd, results, query, spec, PAM, len_without_pam,"coordinate",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="no_result", log="No result available, check your gene name or region", stage=0)
            return

        indexComputing(idd, casData.off_target, casData.mismatch_num)
        idd = save_sgRNA_list(idd, results, query, spec, PAM, len_without_pam,"coordinate",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="success", log="done", stage=0)
    except Exception as e:
        print(f"Error in GeneNameComputing: {str(e)}")
        save_sgRNA_list(
            idd, [], query, request.species, 
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
        idd = save_sgRNA_list(idd, results, seq, spec, PAM, len_without_pam,"fasta",
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

        idd = save_sgRNA_list(idd, results, seq, spec, PAM, len_without_pam,"fasta",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="calculating-index_processing", log="indexing",gene_strand="+")
        
        if len(results) == 0:
            idd = save_sgRNA_list(idd, results, seq, spec, PAM, len_without_pam,"fasta",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="no_result", log="No result available, check your gene name or region", stage=0)
            return

        indexComputing(idd, casData.off_target, casData.mismatch_num)
        idd = save_sgRNA_list(idd, results, seq, spec, PAM, len_without_pam,"fasta",
                                q1, q2, q3, q4, q5, q6, q7, q8, queue_task_id, status="success", log="done", stage=0, gene_strand="+")
    except Exception as e:
        print(f"Error in GeneNameComputing: {str(e)}")
        save_sgRNA_list(
            idd, [], seq, request.species,
            casData.pam, generalSetting.sgRNA_len,"fasta",
            primerConfigData.min_product_size, primerConfigData.max_product_size,
            primerConfigData.min_primer_size, primerConfigData.max_primer_size,
            primerConfigData.optimal_primer_size, primerConfigData.min_tm,
            primerConfigData.max_tm, primerConfigData.optimal_tm,
            queue_task_id, status="failed", log=f"Processing error: {str(e)}"
        )
        raise 
    return

def getMMDetails(bed_file_path: str, json_final_file: str):
    """
    Chuyển file bed đã gán nhãn exon/intron/intergenic sang file json chi tiết.
    """
    mm_details = {}

    with open(bed_file_path, 'r', encoding='utf-8') as bed_file:
        for line in bed_file:
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            region_type = fields[4]

            key = f"{chrom}:{start}"
            mm_details[key] = region_type

    #Mo file json ra
    #bowtie detail co dang "bowtie_details": "NC_000017.11:79479867,,2:C>A,4:C>T,17:G>C,0; NC_000009.12:134635608,,3:T>C,8:T>A,17:C>G,0; NC_000008.11:109469324,,0:G>T,6:T>C,17:A>C,0; NC_000017.11:82854914,,7:A>C,17:C>G,18:T>C,0; NC_000009.12:78397677,,1:T>A,5:G>A,17:A>C,0; NC_000007.14:158190273,,0:A>T,17:G>C,18:A>G,0; NC_000012.12:119804485,,1:G>A,8:A>T,17:T>C,0; NC_000020.11:8839449,,2:G>T,5:A>T,17:A>G,0; NC_000003.12:100090496,,1:C>A,9:A>C,19:G>A,0; ",
    #duyet tat ca cac {} trong file json va lay ra bowtie detail
    #gan nhan lan luot, vi du kia la NC_000017.11:79479867, NC_000009.12:134635608
    #thi gan tuong ung la mm_details["NC_000017.11:79479867"] va mm_details["NC_000009.12:134635608"]
    #cong vao thanh xau roi luu lai trong mis match region "Secondary structure with scaffold": "........(((((((..((....))..)))))))........((((((...((((.(......).)))).)))))).......(((((((...)))))))..., (-24.1 kcal/mol)",
    with open(json_final_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    header = data[0]  # Phan chua thong tin chung
    data = data[1:]  # Chi lay phan chua thong tin sgRNA
    for entry in data:
        bowtie_details = entry.get('bowtie_details', '')
        
        if not bowtie_details or bowtie_details.strip() == "":
            entry['mismatch_region'] = ""
            continue
        
        regions = []
        print(bowtie_details)
        offtargets = bowtie_details.strip().rstrip(';').split(';')
        for offtarget in offtargets:
            offtarget = offtarget.strip()
            if not offtarget:
                continue
            
            chrom_pos = offtarget.split(',,')[0].strip()
            region = mm_details.get(chrom_pos, "unknown")
            regions.append(region)
        
        # Format: "exon, intron, intergenic"
        print(regions)
        entry['mismatch_region'] = ";".join(regions)


    final_data = [header] + data
    with open(json_final_file, 'w', encoding='utf-8') as f:
        json.dump(final_data, f, indent=4)
    

    return


def getMMRegion(name: str, idfile: str):
    """
    Gán nhãn exon/intron/intergenic cho my_regions.sorted.bed dựa vào annotation file của genome.
    Ưu tiên tìm {name}.gff3 > {name}.gff > {name}.gtf
    """

    raw_bed_file = os.path.join(DATA_DIR, f"{idfile}_raw.bed")
    bed_file = f"{idfile}_sorted.bed"    


    # ---  Xác định phần mở rộng thực tế ---
    possible_ext = [".gff3", ".gff", ".gtf"]
    found_ext = None

    for ext in possible_ext:
        file_path = os.path.join(DATA_DIR, f"{name}{ext}")
        if os.path.exists(file_path):
            found_ext = ext
            break

    if not found_ext:
        print(f"Không tìm thấy file annotation: {name}.gff3 / .gff / .gtf")
        return None

    # --- Tạo tên file exon & gene theo đúng phần mở rộng ---
    exon_file = f"{name}_exons.sorted{found_ext}"
    gene_file = f"{name}_genes.sorted{found_ext}"
    result_file = f"{idfile}_mm_annotation.bed"

    exon_file = os.path.join(DATA_DIR, exon_file)
    gene_file = os.path.join(DATA_DIR, gene_file)
    bed_file = os.path.join(DATA_DIR, bed_file)
    result_file = os.path.join(DATA_DIR, result_file)



    # --- lệnh pipeline ---
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
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        except subprocess.CalledProcessError as e:
            print(f" Lỗi khi chạy: {cmd}\nChi tiết: {e}")
            return None

    json_file_path = os.path.join(DATA_DIR, f"vcp{idfile}.json")
    getMMDetails(result_file, json_file_path)



    print(f"Hoàn tất! File kết quả: {result_file} (dựa trên {found_ext})")
    return

def indexComputing(idfile: str, off_target: bool = 0, num_of_mismatches: int = 3):
    try:
        seq_list = []
        seq_list_ml = []
        seq_list_lindel = []

        filename = "vcp" + idfile + ".json"
        file_path = os.path.join(OUTPUT_DIR, filename)
        with open(file_path, 'r', encoding='utf-8') as f:
            datafile = json.load(f)

        pam_name = datafile[0]["pam"]
        bowtie_index_file = datafile[0]["spec"] + "_index"
        spec_name = datafile[0]["spec"]
        sgRNA_len = datafile[0]["sgRNA_len"]

        tmp = datafile[0]
        datafile = datafile[1:]

        for i, row in enumerate(datafile):
            seq_list.append(row["sequence"])
            seq_list_ml.append(row["mlseq"])
            seq_list_lindel.append(row["lindel"])

        write_sgrna_to_fasta_with_IUPAC(seq_list, pam_name, idfile)
        lindel_scores = []

        if (pam_name == "NGG"):
            ml_score = get_ml_score(seq_list_ml)
            rs3_score = get_ml_score_azi3(seq_list_ml)
        else:
            ml_score = ["N/S"] * len(datafile)
            rs3_score = ["N/S"] * len(datafile)


        pol = 0
        limit_num = 1000
        ss_map = {19: 5000, 20: 4000, 21: 3000, 22: 2000, 23: 1000}
        ss = ss_map.get(len(pam_name) + sgRNA_len, 700)
        pp = count_permu_IUPAC(pam_name)
        limit_num = max(20, ss/pp)
        sg_file = f"{idfile}_sgrna_output.fa"

        para_mm = 3
        if off_target == 0:
            para_mm = num_of_mismatches
        command = [
        "bowtie",
        "-v", str(para_mm),   # số mismatch tối đa trên toàn bộ chuỗi
        # "-a",  
        "-k", str(limit_num),
        "-f",
        "-x", bowtie_index_file,
        sg_file,
        "/dev/stdout"
    ]
        
        # if off_target == 1:
        #     command = [
        #     "bowtie",
        #     "-n", str(0),   # số mismatch tối đa trong distal, seed region phai dc giu nguyen
        #     "--seedlen", "9",     
        #     # "-a",  
        #     "-k", str(limit_num),
        #     "-f",
        #     "-x", bowtie_index_file,
        #     sg_file,
        #     "/dev/stdout"
        #     ]

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
        m = defaultdict(int)  # mặc định giá trị là 0
        ll = len(datafile[1]['sequence']) - len(pam_name)
        for line in process.stdout:
            pol = pol + 1
            line = line.strip()
            #off_target là xem user chọn kiểu mismatches trên cả toàn bộ chuỗi hay chỉ distal regions
            #num_of_mismatches là số mismatches tối đa user chọn, chỉ phục vụ cho phần distal regions
            idseq, ttmm = xuly(line, datafile, pam_name, ll, off_target, num_of_mismatches)
            id_in_fa = int(int(line.strip().split('\t')[0]))
            m[id_in_fa] += 1
            if (m[id_in_fa] >= limit_num):
                danh_dau.add(id_in_fa/pp)
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

        for i, row in enumerate(datafile):
            if i in danh_dau:
                if str(row.get("mm3", "0")) != "0":
                    row["mm3"] = f">={row['mm3']}"
                elif str(row.get("mm2", "0")) != "0":
                    row["mm2"] = f">={row['mm2']}"
                elif str(row.get("mm2", "0")) != "0":
                    row["mm1"] = f">={row['mm1']}"
            print(
                f"{i}: {row.get('sequence', '')}, "
                f"location={row.get('location', '')}, "
                f"mm0={row.get('mm0', 0)}, "
                f"mm1={row.get('mm1', 0)}, "
                f"mm2={row.get('mm2', 0)}, "
                f"mm3={row.get('mm3', 0)}"
            )

        for i in range(len(datafile)):
            datafile[i]["mlScore"] = ml_score[i]
            datafile[i]["rs3"] = rs3_score[i]

            datafile[i]["cfdScore"] = round(100 / (100 + datafile[i]["cfdScore"]), 2)

        
        
        
        real_datafile = [tmp] + datafile
        with open(file_path, 'w', encoding='utf-8') as f:
            json.dump(real_datafile, f, ensure_ascii=False, indent=2)

        grouped = defaultdict(list)
        for chrom, start, end, idseq in pos_list:
            grouped[idseq].append((chrom, start, end))

        raw_bed_dir = os.path.join(DATA_DIR, f"{idfile}_raw.bed")
        with open(raw_bed_dir, "w") as f:
            for idseq, regions in grouped.items():
                for chrom, start, end in regions:
                    f.write(f"{chrom}\t{start}\t{end}\t{idseq}\n")

        getMMRegion(spec_name, idfile)

        print("Da tinh toan xong")
        print(limit_num)
        print(pol)
        print(ll)

        checkAndSendMail(idfile)



        return
    except Exception as e:
        print(f"Error in indexComputing: {str(e)}")
        raise
def checkAndSendMail(idfile: str):
    db = SessionLocal()
    #trong ham tinh toan index
    #check xem co phai gui mail ko, co thi tim trong db idfile va gui toi mailist tuong ung 
    try:
        print(f"🔔 ID File {idfile}: Đủ điều kiện gửi email. Bắt đầu truy vấn...")

        email_records = db.query(EmailQueue).filter(
            EmailQueue.idfile == idfile
        ).all()

        if not email_records:
            print(f"⚠️ ID File {idfile}: Không tìm thấy email nào trong hàng đợi.")
            return

        mail_list = [record.email for record in email_records]
        
        print(f"Đã tìm thấy {len(mail_list)} email trong hàng đợi.")

        # 3. Gửi mail tới từng địa chỉ
        success_count = 0

        sendMail(idfile, mail_list)

        delete_count = db.query(EmailQueue).filter(
            EmailQueue.idfile == idfile
        ).delete(synchronize_session='fetch')

        db.commit()
        
        print("---")
        print(f"Hoàn tất gửi mail cho ID File {idfile}:")
        print(f"- Tổng số email trong DB: {len(mail_list)}")
        print(f"- Số email gửi thành công: {success_count}")

    except Exception as e:
        print(f"Lỗi xảy ra khi xử lý ID File {idfile}: {e}")

    finally:
        db.close()