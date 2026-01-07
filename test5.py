import subprocess 
import os      

anno_file = "r570.gff3"  
DATA_DIR = "./app/data"
# Lấy (.gff, .gtf, .gff3)
ext = os.path.splitext(anno_file)[1]

pt = anno_file.split(".")[0]

# Tạo tên file đầu ra cho exon và gene
exon_out = f"{pt}_exons.sorted{ext}"
gene_out = f"{pt}_genes.sorted{ext}"


exon_out = os.path.join(DATA_DIR, exon_out)
gene_out = os.path.join(DATA_DIR, gene_out)

anno_file = os.path.join(DATA_DIR, anno_file)

# Lệnh trích xuất exon
cmd_exon = f"awk 'BEGIN{{OFS=\"\\t\"}} $3 == \"exon\"' {anno_file} | sort -k1,1V -k4,4n > {exon_out}"
# Lệnh trích xuất gene
cmd_gene = f"awk 'BEGIN{{OFS=\"\\t\"}} $3 == \"gene\"' {anno_file} | sort -k1,1V -k4,4n > {gene_out}"

try:
    subprocess.run(cmd_exon, shell=True, check=True, executable='/bin/bash')
    subprocess.run(cmd_gene, shell=True, check=True, executable='/bin/bash')
    print(f"Tạo thành công:\n  {exon_out}\n  {gene_out}")
except subprocess.CalledProcessError as e:
    print(" Lỗi khi xử lý:", e)