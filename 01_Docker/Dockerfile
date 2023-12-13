FROM cumulusprod/cellranger:7.1.0
SHELL ["/bin/bash", "-c"]

# Install system libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends dstat moreutils zip less socat wget sudo curl gdebi-core && \
    apt-get install -y --no-install-recommends libfontconfig1-dev libharfbuzz-dev libfribidi-dev        && \
    apt-get install -y --no-install-recommends libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev     && \
    apt-get install -y --no-install-recommends libcairo2-dev zlib1g libssl-dev libhdf5-dev libxml2-dev

# Install R
RUN export R_VERSION=4.3.0 && \
    curl -O https://cdn.rstudio.com/r/debian-11/pkgs/r-${R_VERSION}_1_amd64.deb && \
    echo y | gdebi r-${R_VERSION}_1_amd64.deb                                   && \
    sudo ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R                       && \
    sudo ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# Install bcl2fastq (https://cumulus.readthedocs.io/en/latest/bcl2fastq.html)
ADD bcl2fastq.zip /software/
RUN cd /software                                            && \
    unzip -d /software/ /software/bcl2fastq.zip             && \
    tar -zxf /software/bcl2fastq2-v2.20.0.422-Source.tar.gz && \
    rm /software/bcl2fastq.zip /software/bcl2fastq2-v2.20.0.422-Source.tar.gz
ENV C_INCLUDE_PATH=/usr/include/x86_64-linux-gnu
ENV INSTALL_DIR=/usr/local/bcl2fastq
ENV SOURCE=/software/bcl2fastq
ENV BUILD=/software/bcl2fastq-build
RUN mkdir ${BUILD}                                            && \
    cd ${BUILD}                                               && \
    chmod ugo+x ${SOURCE}/src/configure                       && \
    chmod ugo+x ${SOURCE}/src/cmake/bootstrap/installCmake.sh && \
    ${SOURCE}/src/configure --prefix=${INSTALL_DIR}           && \
    make -j 30                                                && \
    make install                                              && \
    rm -rf /software/bcl2fastq-build
ENV PATH=$INSTALL_DIR/bin:$PATH

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
RUN pip install --no-cache-dir gspread dnspython

# Install R libraries
RUN R -e "if (!require('BiocManager',quietly=T)) {install.packages('BiocManager',repos='http://cran.rstudio.com/')}; \
          BiocManager::install('rhdf5')"
RUN R -e "install.packages(c('qs',          \
                             'glue',        \
                             'qpdf',        \
                             'purrr',       \
                             'furrr',       \
                             'dplyr',       \
                             'rlist',       \
                             'rdist',       \
                             'tidyr',       \
                             'hdf5r',       \
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


