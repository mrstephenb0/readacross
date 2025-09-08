FROM mambaorg/micromamba:1.5.8

# Create conda env with RDKit + deps
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

# Install Java
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre-headless openjdk-11-jre-headless \
 && rm -rf /var/lib/apt/lists/*

# App files
WORKDIR /app
COPY . /app

# Streamlit config (headless)
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0
ENV JAVA_BIN_TOXTREE=/usr/lib/jvm/java-11-openjdk-amd64/bin/java

# Use micromamba's python to run Streamlit
CMD ["/usr/local/bin/_entrypoint.sh", "bash", "-lc", "streamlit run app.py --server.address=0.0.0.0 --server.port=$PORT"]
