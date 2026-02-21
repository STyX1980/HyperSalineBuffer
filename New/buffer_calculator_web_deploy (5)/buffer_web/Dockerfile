FROM python:3.11-slim

# cache-bust: 2026-02-21
RUN apt-get update && apt-get install -y gcc g++ make cmake libxml2-dev && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN pip install --upgrade pip setuptools wheel

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Confirm key packages installed
RUN python -m flask --version
RUN python -c "import scipy; print('scipy ok')"
RUN python -c "import phreeqpython; print('phreeqpython ok')"

COPY . .

CMD gunicorn app:app --bind 0.0.0.0:$PORT --workers 2 --timeout 120 --log-level debug
