FROM debian:buster
SHELL ["/bin/bash", "-c"]

# Install system libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends dstat moreutils zip unzip gzip less socat wget curl && \
    apt-get install -y --no-install-recommends sudo gdebi-core build-essential alien && \
    apt-get install -y --no-install-recommends apt-transport-https ca-certificates gnupg && \
    apt-get install -y --no-install-recommends libfontconfig1-dev libharfbuzz-dev libfribidi-dev && \
    apt-get install -y --no-install-recommends libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev     && \
    apt-get install -y --no-install-recommends libcairo2-dev zlib1g libssl-dev libhdf5-dev libxml2-dev && \
    apt-get install -y --no-install-recommends python3-pip

# Install R
RUN export R_VERSION=4.3.0 && \
    curl -O https://cdn.rstudio.com/r/debian-10/pkgs/r-${R_VERSION}_1_amd64.deb && \
    echo y | gdebi r-${R_VERSION}_1_amd64.deb                                   && \
    sudo ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R                       && \
    sudo ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# Install Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && apt-get update -y && apt-get install google-cloud-sdk -y

# Add bcl2fastq2
ADD bcl2fastq2-v2-20-0-linux-x86-64.zip /software/
RUN unzip -d /software/ /software/bcl2fastq2-v2-20-0-linux-x86-64.zip && \
    alien -i /software/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && \
    rm /software/bcl2fastq2-v2*

# Add Cell Ranger
ADD cellranger-8.0.0 /software/cellranger-8.0.0
ENV PATH=/software/cellranger-8.0.0:$PATH

# Install Julia
RUN wget -P /software https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz && \
    tar -xzvf /software/julia-1.8.5-linux-x86_64.tar.gz -C /software                                       && \
    rm /software/julia-1.8.5-linux-x86_64.tar.gz
RUN /software/julia-1.8.5/bin/julia -e 'using Pkg;               \
                                        Pkg.add("CSV");          \
                                        Pkg.add("GZip");         \
                                        Pkg.add("HDF5");         \
                                        Pkg.add("FASTX");        \
                                        Pkg.add("IterTools");    \
                                        Pkg.add("CodecZlib");    \
                                        Pkg.add("StatsBase");    \
                                        Pkg.add("DataFrames");   \
                                        Pkg.add("StringViews");  \
                                        Pkg.add("Combinatorics");\
                                        Pkg.add("LinearAlgebra")'

# Install Python packages
RUN pip3 install --no-cache-dir gspread dnspython

# Install R libraries
RUN R -e "if (!require('BiocManager',quietly=T)) {install.packages('BiocManager',repos='http://cran.rstudio.com/')}; \
          BiocManager::install('rhdf5')"
RUN R -e "install.packages(c('qs',          \
                             'glue',        \
                             'qpdf',        \
                             'ineq',        \
                             'purrr',       \
                             'furrr',       \
                             'dplyr',       \
                             'rlist',       \
                             'rdist',       \
                             'tidyr',       \
                             'hdf5r',       \
                             'e1071',       \
                             'Seurat',      \
                             'future',      \
                             'dbscan',      \
                             'ggrastr',     \
                             'stringr',     \
                             'viridis',     \
                             'cowplot',     \
                             'ggplot2',     \
                             'magrittr',    \
                             'parallel',    \
                             'jsonlite',    \
                             'gridExtra',   \
                             'stringdist',  \
                             'ggnewscale',  \
                             'matrixStats', \
                             'SeuratObject'), repos='http://cran.rstudio.com/')"
