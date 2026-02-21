#!/bin/bash
set -e
PORT="${PORT:-8080}"
echo "Starting on port $PORT"
exec gunicorn app:app --bind "0.0.0.0:$PORT" --workers 2 --timeout 120 --log-level info --access-logfile - --error-logfile -
