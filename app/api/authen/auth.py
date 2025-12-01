from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from app.database import get_db
from app.models import User
from app.api.authen.jwthandle import create_access_token, hash_password, verify_password, decode_access_token
from app.api.authen.deps import get_current_user
from pydantic import BaseModel

router = APIRouter()

class UserCreate(BaseModel):
    username: str
    password: str

class UserSignUp(BaseModel):
    username: str
    password: str
    email: str

@router.post("/signup")
def signup(user: UserSignUp, db: Session = Depends(get_db)):
    if db.query(User).filter(User.username == user.username).first():
        raise HTTPException(status_code=400, detail="Username already exists")
    new_user = User(username=user.username,
                    password=hash_password(user.password),
                    email=user.email)
    db.add(new_user)
    db.commit()
    return {"msg": "User created"}

@router.post("/login")
def login(user: UserCreate, db: Session = Depends(get_db)):
    db_user = db.query(User).filter(User.username == user.username).first()
    if not db_user or not verify_password(user.password, db_user.password):
        raise HTTPException(status_code=401, detail="Invalid credentials")
    token = create_access_token({"name": db_user.username, "id": db_user.id})
    return {"access_token": token, "token_type": "bearer"}


class getUserID(BaseModel):
    token: str

@router.post("/getUserID")
def getUserID(data: getUserID):
    print(data)
    payload = decode_access_token(data.token)
    if payload is None:
        raise HTTPException(status_code=401)
    print(payload)
    return payload
