FROM mambaorg/micromamba:1.5.8

USER root
# Java 17 + curl + certs + unzip
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre-headless curl ca-certificates unzip \
 && rm -rf /var/lib/apt/lists/*

# Install Temurin 11 JDK for Toxtree
RUN mkdir -p /opt/java && \
    curl -fsSL https://api.adoptium.net/v3/binary/latest/11/ga/linux/x64/jdk/hotspot/normal/eclipse \
      -o /tmp/temurin11.tar.gz && \
    tar -xzf /tmp/temurin11.tar.gz -C /opt/java && \
    rm /tmp/temurin11.tar.gz && \
    ln -s $(ls -d /opt/java/jdk-11* | head -n 1) /opt/java/temurin-11

# Create conda env
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

# App
WORKDIR /app
COPY . /app

# --- BioTransformer: download full distribution ZIP (contains config.json + KBs)
RUN mkdir -p /app/tools/biotransformer \
 && curl -fL -o /tmp/biotransformer.zip \
      "https://github.com/mrstephenb0/readacross/releases/download/assets/wishartlab-biotransformer3.0jar-6432cf887ed7.zip" \
 && unzip -q /tmp/biotransformer.zip -d /app/tools/biotransformer \
 && rm /tmp/biotransformer.zip

# Optional sanity checks in build logs:
RUN find /app/tools/biotransformer -maxdepth 3 -type f -name 'config.json' -print || true \
 && ls -lah /app/tools/biotransformer || true


# tools for downloading and unzipping
RUN apt-get update && apt-get install -y --no-install-recommends unzip curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# --- Toxtree: download full distribution ZIP (contains ext/ + index.properties)
RUN mkdir -p /app/tools/toxtree \
 && curl -fL -o /tmp/toxtree.zip \
      "https://github.com/mrstephenb0/readacross/releases/download/assets/Toxtree-v3.1.0.1851.zip" \
 && unzip -q /tmp/toxtree.zip -d /app/tools/toxtree \
 && rm /tmp/toxtree.zip

# Optional: quick sanity check in build logs
RUN ls -lah /app/tools/toxtree/Toxtree-v3.1.0.1851/Toxtree \
 && ls -lah /app/tools/toxtree/Toxtree-v3.1.0.1851/Toxtree/ext \
 && grep -H . /app/tools/toxtree/Toxtree-v3.1.0.1851/Toxtree/ext/index.properties || true

# Streamlit
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
# Toxtree uses Java 11
ENV JAVA_BIN_TOXTREE=/opt/java/temurin-11/bin/java

EXPOSE 8501
CMD ["/usr/local/bin/_entrypoint.sh","bash","-lc","streamlit run app.py --server.address=0.0.0.0 --server.port=$PORT"]
