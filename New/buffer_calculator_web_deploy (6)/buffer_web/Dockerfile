FROM python:3.11-slim

RUN apt-get update && apt-get install -y gcc g++ make cmake libxml2-dev && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN pip install --upgrade pip setuptools wheel

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Verify installs
RUN python -c "import flask, scipy, phreeqpython, gunicorn; print('all packages OK')"

COPY . .

RUN chmod +x start.sh

EXPOSE 8080
CMD ["bash", "start.sh"]
