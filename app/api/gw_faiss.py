import numpy as np
from pydantic import BaseModel
from fastapi import APIRouter, HTTPException
from .nonModel import GenomeUpdate, update_genome_status
from app.database import SessionLocal
from app.models import Genome
from typing import Optional

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

class RunFaissPipelineRequest(BaseModel):
    display_id: str
    PAM: str = "NGG"
    sgRNA_length: int = 20
    seed_length: int = 9
    hamming_distance: int = 3
    flank_before: int = 100
    flank_after: int = 100
    maillist: list[str] = []

@router.post("/runFaissPipeline")
async def runFaissPipeline(data: RunFaissPipelineRequest):

    # Lookup genome by display_id
    db = SessionLocal()
    try:
        td = db.query(Genome).filter(Genome.id_for_user_display == data.display_id).first()
        if not td:
            raise HTTPException(status_code=404, detail="Genome not found for the given display_id")
        
        if td.status != 'success':
            raise HTTPException(status_code=400, detail=f"Genome is not ready (status: {td.status})")

        x = td.kbstorage
        if x > 60956679:
            print("ko dc dau")
            return {"status": "Not Available"}

        # Update genome-wide state
        update_data = GenomeUpdate(display_id=data.display_id, gw_state="navailable")
        update_genome_status(update_data)

        from .tasks import run_pipeline
        run_pipeline.delay(
            data.display_id,
            data.PAM,
            data.sgRNA_length,
            data.seed_length,
            data.hamming_distance,
            data.flank_before,
            data.flank_after,
            data.maillist
        )
    finally:
        db.close()

    return {"status": "pending in queue"}