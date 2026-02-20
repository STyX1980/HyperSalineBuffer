FROM python:3.11-slim

RUN apt-get update && apt-get install -y gcc g++ make cmake libxml2-dev && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN pip install --upgrade pip setuptools wheel

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["gunicorn", "app:app", "--bind", "0.0.0.0:8080", "--workers", "4", "--timeout", "120"]
