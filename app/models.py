# app/models.py
from sqlalchemy import Column, Integer, String, BigInteger, ForeignKey
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
    status = Column(String(20), server_default="unknown")
    log = Column(String(20), server_default="unknown")
    owner_id = Column(Integer, ForeignKey("users.id"), nullable=False)
    task_queue_id = Column(String(100))
    gw_state = Column(String(20), server_default="available")
    upload_id = Column(String(100))
    

    owner = relationship("User", back_populates="genomes")

class EmailQueue(Base):
    __tablename__ = "emailqueue"  

    idfile = Column(String(100), primary_key=True, index=True)
    email = Column(String(100), primary_key=True, index=True)
