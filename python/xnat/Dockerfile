FROM python:3.12-slim

# System dependencies needed by opencv and scikit-image
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgl1 \
    libglib2.0-0 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies first (cached layer)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code
COPY xnat_runner.py .
COPY functions/ functions/

ENTRYPOINT ["python", "xnat_runner.py"]
