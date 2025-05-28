FROM python:3.13-slim

LABEL org.opencontainers.image.authors="mfoster11@mgh.harvard.edu" \
    org.opencontainers.image.source="https://github.com/mjfos2r/rst_caller" \
    org.opencontainers.image.description="Python3 application for the classification of ospC allele type in Borrelia burgdorferi genome assemblies" \
    org.opencontainers.image.version="1.0.0" \
    maintainer="mfoster11@mgh.harvard.edu"

RUN apt-get update && apt-get install -y \
    bash \
    && rm -rf /var/lib/apt/lists/*
# add UV
ADD https://astral.sh/uv/0.7.8/install.sh /uv-installer.sh
RUN sh /uv-installer.sh && rm /uv-installer.sh
ENV PATH="/root/.local/bin/:$PATH"
RUN mkdir -p /app/rst_caller /data
ADD --checksum=sha256:d8b160af8548639691a6f2ad153bff068c8047a3f99f61ea16ee50cf2e793dab https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz /tmp/
RUN tar -xzf /tmp/seqkit_linux_amd64.tar.gz -C /usr/local/bin --strip-components=1 --no-same-owner \
    && rm /tmp/seqkit_linux_amd64.tar.gz
WORKDIR /app/rst_caller
COPY . .
WORKDIR /app
RUN uv venv
ENV PATH="/app/.venv/bin:$PATH"
RUN uv pip install -e rst_caller/
RUN rst_caller --version
WORKDIR /data
ENTRYPOINT [ "/bin/bash" ]