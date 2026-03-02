#!/bin/bash

echo ">>> Waiting for PostgreSQL..."

# Wait until Postgres is ready
until pg_isready -h db -p 5432 -U postgres >/dev/null 2>&1; do
  echo "PostgreSQL is unavailable - sleeping"
  sleep 2
done

echo ">>> PostgreSQL is ready!"

echo ">>> Running init_db.py..."
python /app/app/init_db.py
