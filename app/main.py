from fastapi import FastAPI, Depends
from pydantic import BaseModel
from app.api.cmdprocess import extract_exon_by_gene
import httpx
from starlette.requests import Request
from fastapi.middleware.cors import CORSMiddleware
import os, re
from .database import Base, engine, get_db
from . import models, database
from sqlalchemy.orm import Session
import asyncio
from app.models import Genome

app = FastAPI()

from app.api.lookUpsgRNA import router as LookUpSgRNArouter
from app.api.export import router as ExportRouter
from app.api.nonModel import router as NonModelRouter
from app.api.authen.auth import router as AuthRouter
from app.api.gw_faiss import router as FaissRouter
from app.configs import get_settings
from app import cron_jobs

settings = get_settings()

raw_ips = os.getenv("CORS_ALLOWED_IPS")
ip_list = [re.escape(ip.strip()) for ip in raw_ips.split(",")]
combined_regex = rf"^https?://({'|'.join(ip_list)})(:\d+)?$"

print(combined_regex)
app.add_middleware(
    CORSMiddleware,
    allow_origin_regex=combined_regex,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
class fastaEntry(BaseModel):
    dna_seq: str
    species: str


PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "app/data")

@app.get("/")
async def read_root(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})

class getdt(BaseModel):
    user_id: str
@app.post("/getData")
async def root(new: getdt, db: Session = Depends(get_db)):
    genomes = db.query(Genome).filter(Genome.owner_id == new.user_id).all()
    return {"first": genomes, "second": "zerooooo", "users": "nah"}

class Data(BaseModel):
    gen_name: str
class TraitToID(BaseModel):
    species: str
    tt: str

@app.post("/getGeneID")
async def root(dl: TraitToID):

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene"
    query = dl.species + "[Organism]" + "+" + dl.tt
    params = {'db': "gene", 'term': query, 'retmode': "json"}

    async with httpx.AsyncClient() as client:
        response = await client.get(url, params=params)

    if response.status_code == 200:
        try:
            res = response.json()
            return {'first': res["esearchresult"]["idlist"], 'second': 'LOI ROI HUHUHUHU'}
        except ValueError:
            return {'error': 'ko phai file json '}

    return {'first': "call api that bai", 'second': 'chua ho tro'}
@app.post("/getSeq")
async def root(dl: Data):

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {'db': "gene", 'id': dl.gen_name, 'retmode': "json"}

    async with httpx.AsyncClient() as client:
        response = await client.get(url, params=params)

    if response.status_code == 200:
        try:
           # res = response.json()
           # start = res["result"][dl.gen_name]["genomicinfo"][0]["chrstart"]
           # end = res["result"][dl.gen_name]["genomicinfo"][0]["chrstop"]
           # chromosome = res["result"][dl.gen_name]["genomicinfo"][0]["chraccver"]

            twobit_file = "hg38.2bit"

            #realSeq = get_fasta_from_twobit(twobit_file, chromosome, start, end)

            exonSeq = extract_exon_by_gene("EXT1")
            return {'first': exonSeq, 'second': 'LOI ROI HUHUHUHU'}
        except ValueError:
            return {'error': 'ko phai file json '}
    else:
        return {'first': 'ukelele', 'second': 'LOI ROI HUHUHUHU'}


class Test(BaseModel):
    genome: str
    chrom: str
    start: int
    end: int


app.include_router(LookUpSgRNArouter, prefix="/api", tags=["genes"])
app.include_router(ExportRouter, prefix="/export", tags=["export"])
app.include_router(NonModelRouter, prefix="/non_model", tags=["nonModel"])
app.include_router(AuthRouter, prefix="/auth", tags=["auth"])
app.include_router(FaissRouter, prefix="/faiss", tags=["faiss"])


@app.on_event("startup")
async def start_jobs():
    asyncio.create_task(cron_jobs.cleanup_files())

