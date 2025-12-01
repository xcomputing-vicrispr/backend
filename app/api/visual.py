from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from .cfdEffiencyScore import get_cfd_score
from .mlEffiencyScore import get_ml_score



import subprocess
import smtplib, numpy as np
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication

import re, os, json, random, string


router = APIRouter()

sender = "vicrispr@gmail.com"
password = "yghnybppmeaehpxs"


class Exo_Intro(BaseModel):
    gene_name: str
    spec_name: list


PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PARENT_DIR, "data")
OUTPUT_DIR = "app/data"


@router.post("/submitSendMail")
async def getExoIntroData(data: Exo_Intro):
    gene = data.gene_name
    spec = data.spec_name
