# Getting python from Docker Hub
FROM python:3.9.5

# Setting working directory
WORKDIR /app

# Install the dependencies
RUN apt-get update && \
        apt-get install -y --no-install-recommends \
        gcc \
        build-essential \
        libssl-dev \
        git

# Install ATaRVa
RUN git clone https://github.com/SowpatiLab/ATaRVa.git /app

RUN pip install --no-cache-dir . && \
        rm -rf /var/lib/apt/lists/*

# Define the default command to run your application
ENTRYPOINT ["python", "-m", "ATARVA.core"]
