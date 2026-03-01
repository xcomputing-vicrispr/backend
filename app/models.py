from sqlalchemy import Column, Integer, String, BigInteger, ForeignKey, Float, Text, DateTime, func
from sqlalchemy.orm import relationship
from .database import Base

class User(Base):
    __tablename__ = "users"

    id = Column(Integer, primary_key=True, index=True)
    username = Column(String(50), unique=True, index=True, nullable=False)
    password = Column(String(255), nullable=False)
    email = Column(String(100), unique=True, index=True)

    genomes = relationship("Genome", back_populates="owner")


class Genome(Base):
    __tablename__ = "genome"  

    gname = Column(String(100), primary_key=True, index=True)
    kbstorage = Column(BigInteger)
    status = Column(String(100), server_default="unknown")
    log = Column(String(100), server_default="unknown")
    owner_id = Column(Integer, ForeignKey("users.id"), nullable=False, primary_key=True)
    task_queue_id = Column(String(100))
    gw_state = Column(String(100), server_default="available")
    upload_id = Column(String(100))
    

    owner = relationship("User", back_populates="genomes")

class EmailQueue(Base):
    __tablename__ = "emailqueue"  

    idfile = Column(String(100), primary_key=True, index=True)
    email = Column(String(100), primary_key=True, index=True)

class TaskMetadata(Base):
    __tablename__ = 'task_metadata'

    query_id = Column(String(100), primary_key=True, unique=True, nullable=False, autoincrement=False)
    query_name = Column(Text)
    spec = Column(String(100))
    pam = Column(String(50))
    sgrna_len = Column(Integer)
    gene_strand = Column(String(5))
    type_task = Column(String(100))
    min_product_size = Column(Integer)
    max_product_size = Column(Integer)
    min_primer_size = Column(Integer)
    max_primer_size = Column(Integer)
    optimal_primer_size = Column(Integer)
    min_tm = Column(Float)
    max_tm = Column(Float)
    optimal_tm = Column(Float)
    status = Column(String(50))

    queue_task_id = Column(String(100))
    log = Column(Text, nullable=True)

    created_at = Column(DateTime(timezone=True), server_default=func.now())

    sgrnas = relationship("Sgrna", back_populates="task", cascade="all, delete-orphan")

class Sgrna(Base):
    __tablename__ = 'sgrnas'

    query_id = Column(String(100), ForeignKey('task_metadata.query_id', ondelete='CASCADE'), primary_key=True, autoincrement=False)
    stt = Column(Integer, primary_key=True, autoincrement=False)
    sequence = Column(Text)
    location = Column(String(255))
    strand = Column(String(5))
    gc_content = Column(Float)
    self_complementary = Column(Float)
    primer = Column(Text)
    mlseq = Column(String(255))
    mm0 = Column(Integer)
    mm1 = Column(Integer)
    mm2 = Column(Integer)
    mm3 = Column(Integer)
    cfd_score = Column(Float)
    ml_score = Column(Float)
    micro_score = Column(Float)
    rs3_score = Column(Float)
    mmej_pre = Column(Text)
    sec_structure = Column(Text)
    lindel = Column(Text)
    bowtie_details = Column(Text)
    mismatch_region = Column(Text)

    task = relationship("TaskMetadata", back_populates="sgrnas")