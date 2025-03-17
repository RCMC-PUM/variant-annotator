FROM ubuntu:latest

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        wget \
        curl \
        bcftools \
        python3.12 \
        python3.12-venv \
        python3-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
RUN ln -sf /usr/bin/python3.12 /usr/bin/python

# Create a working directory
WORKDIR /app

# Copy
COPY . .

# Install dependencies
RUN python -m pip install --break-system-package --no-cache-dir poetry poetry-plugin-export
RUN python -m poetry export --without-hashes --format=requirements.txt -o requirements.txt
RUN python -m pip install --break-system-package -r requirements.txt

# Change mode of bin files
RUN chmod +x bin/*.py

# Set the default command
CMD ["/bin/bash"]
