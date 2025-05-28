FROM python:3.13-slim
COPY --from=ghcr.io/astral-sh/uv:0.7.8 /uv /uvx /bin/

LABEL org.opencontainers.image.authors="mfoster11@mgh.harvard.edu" \
    org.opencontainers.image.source="https://github.com/mjfos2r/rst_caller" \
    org.opencontainers.image.description="Python3 application for the classification of ospC allele type in Borrelia burgdorferi genome assemblies" \
    maintainer="mfoster11@mgh.harvard.edu"

RUN apt-get update && apt-get install -y \
    bash \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*


ENV PATH="/root/.local/bin/:$PATH"
RUN mkdir -p /app/rst_caller /data
RUN curl -LsS https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz >/tmp/seqkit_linux_amd64.tar.gz \
    && curl -LsS https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz.md5.txt >/tmp/seqkit_linux_amd64.tar.gz.md5.txt

RUN [ "$(md5sum '/tmp/seqkit_linux_amd64.tar.gz' | awk '{print $1}')" = "$(awk '{print $1}' /tmp/seqkit_linux_amd64.tar.gz.md5.txt)" ]\
    && { echo -e "\033[0;32mChecksum valid!\033[0m"; tar -xzf /tmp/seqkit_linux_amd64.tar.gz -C /usr/local/bin --no-same-owner; rm /tmp/seqkit_linux_amd64.tar.gz; } \
    || { echo -e "\033[0;31mChecksum invalid!\033[0m"; exit 1; }
WORKDIR /app
COPY . .
RUN uv sync
RUN uv run rst_caller --version
ENV PATH="/app/.venv/bin:$PATH"
WORKDIR /data
ENTRYPOINT [ "/bin/bash" ]
