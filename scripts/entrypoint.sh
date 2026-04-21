#!/bin/bash
set -e
eval "$(micromamba shell hook --shell bash)"

# ENV to run the backend (fix: choose correct one)
ENV_NAME=${ENV_NAME:-base}

echo ">>> Using environment: $ENV_NAME"

# Initialize conda
# source /opt/conda/etc/profile.d/conda.sh
# conda activate $ENV_NAME
micromamba activate $ENV_NAME

echo ">>> Waiting for PostgreSQL..."
until pg_isready -h $POSTGRES_HOST -p $POSTGRES_PORT; do
    echo "PostgreSQL is unavailable - sleeping"
    sleep 1
done

echo ">>> PostgreSQL is ready!"

echo ">>> Running DB initialization script..."
/app/scripts/init_db.sh || true

echo ">>> Starting FastAPI..."
exec uvicorn app.main:app --host 0.0.0.0 --port 8099