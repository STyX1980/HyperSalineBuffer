FROM python:3.11-slim

RUN apt-get update && apt-get install -y gcc g++ make cmake libxml2-dev && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN pip install --upgrade pip setuptools wheel

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Railway sets $PORT dynamically — must use shell form (not exec form) to expand it
CMD gunicorn app:app --bind 0.0.0.0:$PORT --workers 2 --timeout 120
