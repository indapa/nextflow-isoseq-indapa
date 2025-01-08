# Start from the Python 3.11 slim image as a base
FROM python:3.11-slim

# Install prerequisites
RUN apt-get update && apt-get install -y --no-install-recommends build-essential curl gcc g++ git unzip zip sudo htop less gzip \
    curl ca-certificates gnupg \
    && rm -rf /var/lib/apt/lists/*



# Install Python packages
RUN pip install --no-cache-dir pandas numpy

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
        unzip awscliv2.zip && \
        sudo ./aws/install

# Create a non-root user (named 'vscode') to match devcontainer defaults
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID
RUN groupadd --gid $USER_GID $USERNAME && \
    useradd --uid $USER_UID --gid $USER_GID -m $USERNAME

# (Optional) Give the non-root user sudo access without password
RUN echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers


# Install Miniconda
RUN curl -fsSLo /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    chmod +x /tmp/miniconda.sh && \
    /tmp/miniconda.sh -b -p /opt/miniconda && \
    rm -f /tmp/miniconda.sh && \
    /opt/miniconda/bin/conda clean -afy

# Make conda available on the PATH
ENV PATH="/opt/miniconda/bin:${PATH}"

# Install nf-core from bioconda/conda-forge
RUN conda install -y -c conda-forge -c bioconda nf-core samtools && \
    conda clean -afy

USER $USERNAME


# Set up SDKMAN environment variables
ENV SDKMAN_DIR="/home/vscode/.sdkman"
ENV SDKMAN_NON_INTERACTIVE="true"

# Install SDKMAN
RUN curl -s "https://get.sdkman.io" | bash

# Install a specific version of Java with SDKMAN
RUN bash -c "source $SDKMAN_DIR/bin/sdkman-init.sh && sdk install java 17.0.8-tem"

# Make the installed Java available on the PATH in future layers
ENV PATH="$SDKMAN_DIR/candidates/java/current/bin:${PATH}"

WORKDIR /workspaces

RUN curl -s https://get.nextflow.io | bash 
RUN sudo mv nextflow /usr/local/bin/


WORKDIR /workspaces