FROM python:3.11-slim

# System deps needed to compile phreeqpython's C extension
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    cmake \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Upgrade pip first, then install wheel so binary packages are preferred
RUN pip install --upgrade pip setuptools wheel

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Verify phreeqpython installed correctly and can find pitzer.dat
RUN python -c "
from phreeqpython import PhreeqPython
import os, glob
pkg_dir = os.path.dirname(__import__('phreeqpython').__file__)
matches = glob.glob(os.path.join(pkg_dir, '**', 'pitzer.dat'), recursive=True)
if not matches:
    raise RuntimeError('pitzer.dat not found in phreeqpython package')
pp = PhreeqPython(database=matches[0])
print('phreeqpython OK, pitzer.dat found at:', matches[0])
"

# Copy application
COPY . .

# Use $PORT env var (Railway sets this automatically)
CMD ["gunicorn", "app:app", "--bind", "0.0.0.0:8080", "--workers", "4", "--timeout", "120"]
