"""
Genome Deduplication - Two-Phase Hashing Module.

Phase 1 (fast):
    - File size ±5MB filter

Phase 2 (canonical, slow):
    - Fasta: seqkit sort + normalize line-wrap + lowercase to uppercase (-u) → MD5
    - GFF3: strip comments/blank lines/trailing whitespace + LC_ALL=C sort → MD5
    → Definitive comparison: handles different formatting/ordering.
"""

import hashlib
import json
import os
import subprocess
from sqlalchemy.orm import Session
from app.models import Genome

SIZE_TOLERANCE_BYTES = 5 * 1024 * 1024

def canonical_fasta_hash(filepath: str) -> str | None:
    try:
        cmd = f"seqkit sort -n {filepath} 2>/dev/null | seqkit seq -w 0 -u | md5sum"
        output = subprocess.check_output(cmd, shell=True).decode()
        return output.split()[0]
    except Exception as e:
        print(f"❌ canonical_fasta_hash error for {filepath}: {e}")
        return None

def canonical_gff3_hash(filepath: str) -> str | None:
    try:
        cmd = (
            f"export LC_ALL=C; "
            f"grep -v '^#' {filepath} | "
            f"grep -v '^[[:space:]]*$' | "
            f"sed 's/[[:space:]]*$//' | "
            f"sort -S 25% --parallel=4 | "
            f"md5sum"
        )
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode()
        return output.split()[0]
    except Exception as e:
        print(f"❌ canonical_gff3_hash error for {filepath}: {e}")
        return None

def find_existing_fasta(db: Session, filepath: str) -> tuple[str | None, str | None]:
    file_size = os.path.getsize(filepath)

    size_candidates = (
        db.query(Genome)
        .filter(
            Genome.fasta_size.between(file_size - SIZE_TOLERANCE_BYTES, file_size + SIZE_TOLERANCE_BYTES),
            Genome.id_use_for_us_fasta.isnot(None),
        )
        .all()
    )

    canonical = canonical_fasta_hash(filepath)
    if not canonical:
        return None, None

    if size_candidates:
        for c in size_candidates:
            if c.id_use_for_us_fasta == canonical:
                return canonical, canonical
                
    return None, canonical

def find_existing_gff3(db: Session, filepath: str) -> tuple[str | None, str | None]:
    file_size = os.path.getsize(filepath)

    size_candidates = (
        db.query(Genome)
        .filter(
            Genome.anno_size.between(file_size - SIZE_TOLERANCE_BYTES, file_size + SIZE_TOLERANCE_BYTES),
            Genome.id_use_for_us_gff3.isnot(None),
        )
        .all()
    )

    canonical = canonical_gff3_hash(filepath)
    if not canonical:
        return None, None

    if size_candidates:
        for c in size_candidates:
            if c.id_use_for_us_gff3 == canonical:
                return canonical, canonical

    return None, canonical

def generate_query_hash(params: dict) -> str:
    """
    Generate a SHA-256 hash from a dictionary of query parameters.
    Normalization: 
    - Sort keys alphabetically
    - Convert values to strings and handle casing (strip + uppercase)
    - Serialize to JSON
    """
    normalized = {}
    for k in sorted(params.keys()):
        val = params[k]
        if isinstance(val, str):
            normalized[k] = val.strip().upper()
        elif isinstance(val, dict):
            normalized[k] = generate_query_hash(val)
        else:
            normalized[k] = val
            
    serialized = json.dumps(normalized, sort_keys=True)
    return hashlib.sha256(serialized.encode()).hexdigest()
