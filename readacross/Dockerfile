FROM mambaorg/micromamba:1.5.8

# Create conda env with RDKit + deps
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

# Install Java
USER root
RUN apt-get update && apt-get install -y --no-install-recommends openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# App files
WORKDIR /app
COPY . /app

# Streamlit config (headless)
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0
# Some hosts (e.g., Render/Cloud Run) provide PORT env var:
ENV PORT=8501

EXPOSE 8501
# Use micromamba's python to run Streamlit
CMD ["/usr/local/bin/_entrypoint.sh", "streamlit", "run", "app.py", "--server.port=${PORT}"]
