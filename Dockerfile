FROM nfcore/base
LABEL authors="Urmo Võsa" \
      description="Docker image containing tools for MR and coloc comparisons"

COPY environment.yml /
RUN apt-get update && apt install -y libgmp-dev && apt install -y build-essential
RUN conda install -n base -c conda-forge mamba && \
    mamba env create -f environment.yml && \
    mamba clean -a
ENV PATH /opt/conda/envs/eqtlgenmrlink2/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
COPY temp/IGUtilityPackage_0.3.02.tar.gz /
RUN R -e "install.packages('IGUtilityPackage_0.3.02.tar.gz', repos = NULL, type = 'source', dependencies = TRUE)"
