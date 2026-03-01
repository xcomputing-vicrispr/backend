import numpy as np
from pydantic import BaseModel
from fastapi import APIRouter
from .nonModel import GenomeUpdate, update_genome_status
from app.database import SessionLocal
from app.models import Genome

router = APIRouter()

NUC2ONEHOT = {
    'A': 0b1000,
    'C': 0b0100,
    'G': 0b0010,
    'T': 0b0001
}

mapping = {
    'A': np.array([1,0,0,0], dtype=np.float32),
    'C': np.array([0,1,0,0], dtype=np.float32),
    'G': np.array([0,0,1,0], dtype=np.float32),
    'T': np.array([0,0,0,1], dtype=np.float32)
}

class BuildFaissIndexRequest(BaseModel):
    owner_id: int
    genome_name: str
    PAM: str = "NGG"
    sgRNA_length: int = 20

class QueryFaissIndexRequest(BaseModel):
    owner_id: int
    genome_name: str
    PAM: str = "NGG"
    sgRNA_length: int = 20
    seed_length: int = 9
    hamming_distance: int = 3
    flank_before: int = 100
    flank_after: int = 100
    maillist: list[str] = []

@router.post("/runFaissPipeline")
async def runFaissPipeline(data: QueryFaissIndexRequest):

    #neu nhu availible thi ok ko thi cut
    print("bao sam jkqA2")
    update_data = GenomeUpdate(gname=data.genome_name,owner_id=data.owner_id, gw_state="navailable")
    update_genome_status(update_data)

    db = SessionLocal()
    td = db.query(Genome).filter(Genome.gname == data.genome_name, Genome.owner_id == data.owner_id).first()

    x = td.kbstorage
    
    if x > 60956679:
        print("ko dc dau")
        return {"status": "Not Availible"}

    from .tasks import run_pipeline
    run_pipeline.delay(
        data.owner_id,
        data.genome_name,
        data.PAM,
        data.sgRNA_length,
        data.seed_length,
        data.hamming_distance,
        data.flank_before,
        data.flank_after,
        data.maillist
    )
    return {"status": "pending in queue"}