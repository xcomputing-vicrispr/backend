import os
from functools import lru_cache
from dotenv import load_dotenv

# Load file .env
load_dotenv(os.path.join(os.path.dirname(__file__), "../.env"))

class Settings:
    # CORS
    CORS_ALLOW_ORIGINS: str = os.getenv("CORS_ALLOW_ORIGINS")

    # SMTP Email
    SMTP_SENDER: str = os.getenv("SMTP_SENDER")
    SMTP_PASSWORD: str = os.getenv("SMTP_PASSWORD")
    SMTP_HOST: str = os.getenv("SMTP_HOST", "smtp.gmail.com")
    SMTP_PORT: int = int(os.getenv("SMTP_PORT", 587))
    SMTP_TLS: bool = os.getenv("SMTP_TLS", "true").lower() == "true"

    # Redis
    REDIS_BACKEND_URL: str = os.getenv("REDIS_BACKEND_URL")
    REDIS_APP_URL: str = os.getenv("REDIS_APP_URL")
    REDIS_FQ_URL: str = os.getenv("REDIS_FQ_URL")

    # RabbitMQ / Celery
    RABBITMQ_BROKER_URL: str = os.getenv("RABBITMQ_BROKER_URL")

    # Database
    DB_HOST: str = os.getenv("DB_HOST")
    DB_PORT: int = int(os.getenv("DB_PORT", 5432))
    DB_USER: str = os.getenv("DB_USER")
    DB_PASSWORD: str = os.getenv("DB_PASSWORD")
    DB_NAME: str = os.getenv("DB_NAME")
    DATABASE_URL: str = os.getenv("DATABASE_URL")

@lru_cache()
def get_settings() -> Settings:
    return Settings()

DATABASE_URL = os.getenv("DATABASE_URL")
print(DATABASE_URL)
