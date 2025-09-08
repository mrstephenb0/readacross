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

# >>> Your conda env creation stays as-is <<<
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

# App
WORKDIR /app
COPY . /app

# Streamlit
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

# Tell the app which Java to use for Toxtree
ENV JAVA_BIN_TOXTREE=/opt/java/temurin-11/bin/java

EXPOSE 8501
# Expand $PORT at runtime
CMD ["/usr/local/bin/_entrypoint.sh","bash","-lc","streamlit run app.py --server.address=0.0.0.0 --server.port=$PORT"]
