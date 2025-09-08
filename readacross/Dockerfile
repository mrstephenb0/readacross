FROM mambaorg/micromamba:1.5.8

USER root
# Java 17 for BioTransformer + tools to fetch Temurin 11
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre-headless curl ca-certificates \
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

# Download real jars into the exact nested paths your code expects
RUN mkdir -p tools/biotransformer/wishartlab-biotransformer3.0jar-6432cf887ed7 \
 && curl -fL -o tools/biotransformer/wishartlab-biotransformer3.0jar-6432cf887ed7/BioTransformer3.0_20230525.jar \
      "https://github.com/mrstephenb0/readacross/releases/download/assets/BioTransformer3.0_20230525.jar" \
 && mkdir -p tools/toxtree/Toxtree-v3.1.0.1851/Toxtree \
 && curl -fL -o tools/toxtree/Toxtree-v3.1.0.1851/Toxtree/Toxtree-3.1.0.1851.jar \
      "https://github.com/mrstephenb0/readacross/releases/download/assets/Toxtree-3.1.0.1851.jar" \
 && test $(stat -c%s tools/biotransformer/wishartlab-biotransformer3.0jar-6432cf887ed7/BioTransformer3.0_20230525.jar) -ge 10000000 \
 && test $(stat -c%s tools/toxtree/Toxtree-v3.1.0.1851/Toxtree/Toxtree-3.1.0.1851.jar) -ge 5000000

# Streamlit
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
# Toxtree uses Java 11
ENV JAVA_BIN_TOXTREE=/opt/java/temurin-11/bin/java

EXPOSE 8501
CMD ["/usr/local/bin/_entrypoint.sh","bash","-lc","streamlit run app.py --server.address=0.0.0.0 --server.port=$PORT"]
