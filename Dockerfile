FROM python:3.10-slim

# System dependencies: hmmer for ANARCI, git for pip installs from GitHub
RUN apt-get update && apt-get install -y --no-install-recommends \
    hmmer \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy app code
COPY . .

# Download vanilla ProteinMPNN weights (not included in repo)
RUN mkdir -p /app/ProteinMPNN/vanilla_model_weights && \
    wget -q -O /app/ProteinMPNN/vanilla_model_weights/v_48_020.pt \
    https://media.githubusercontent.com/media/dauparas/ProteinMPNN/main/vanilla_model_weights/v_48_020.pt

# HuggingFace Spaces runs as non-root user 1000
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user PATH=/home/user/.local/bin:$PATH

EXPOSE 7860

CMD ["python", "app.py"]
