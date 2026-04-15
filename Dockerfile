FROM python:3.10-slim

# System dependencies: hmmer for ANARCI, git for pip installs from GitHub
RUN apt-get update && apt-get install -y --no-install-recommends \
    hmmer \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy app code
COPY . .

# HuggingFace Spaces runs as non-root user 1000
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user PATH=/home/user/.local/bin:$PATH

EXPOSE 7860

CMD ["python", "app.py"]
