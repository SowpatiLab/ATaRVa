# Getting python from Docker Hub
FROM python:3.9.18

# Setting working directory
WORKDIR /ATARVA

# Copy the requirements file into the container
COPY requirements.txt .

# Install the dependencies
RUN apt-get update && \
        apt-get install -y --no-install-recommends \
        gcc \
        build-essential \
        libssl-dev \
        git \
        && \
        pip install --no-cache-dir -r requirements.txt && \
        rm -rf /var/lib/apt/lists/*
        
# Install ATaRVa
RUN git clone https://github.com/SowpatiLab/ATaRVa.git /ATARVA


# Install SSW
RUN git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.git /ATARVA/Complete-Striped-Smith-Waterman-Library

# Set up build directory
WORKDIR /ATARVA/Complete-Striped-Smith-Waterman-Library/src

# build ssw
RUN gcc -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h

# Add shared library to LD_LIBRART_PATH
ENV LD_LIBRARY_PATH="/ATARVA/Complete-Striped-Smith-Waterman-Library/src:$LD_LIBRARY_PATH"
RUN ldconfig

WORKDIR /ATARVA

# Add src path of SSW to pythonpath
ENV PYTHONPATH="/ATARVA/Complete-Striped-Smith-Waterman-Library/src"

# Define the default command to run your application
ENTRYPOINT ["python", "core.py"]
