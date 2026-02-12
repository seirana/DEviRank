# syntax=docker/dockerfile:1
FROM python:3.12-slim

# System deps (keep minimal; add others only if needed)
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    build-essential \
    ca-certificates \
    curl \
    wget \
    unzip \
 && rm -rf /var/lib/apt/lists/*

# Workdir inside container
WORKDIR /app

# Install Python requirements first (better caching)
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

# Copy repo code
COPY . /app

# Default command: show help / no-op (we will override in run scripts)
CMD ["python3", "-c", "print('Container ready. Use the provided run scripts.')"]

