FROM cumulusprod/cellranger:7.1.0
SHELL ["/bin/bash", "-c"]

RUN apt-get update                                                                                   && \
    apt-get install -y --no-install-recommends zlib1g-dev libcurl4-gnutls-dev libssl-dev libxml2-dev && \
    apt-get install -y --no-install-recommends libpng-dev libudunits2-dev libhdf5-dev                && \
    apt-get install -y --no-install-recommends liblapack-dev libblas-dev gfortran                    && \
    apt-get install -y --no-install-recommends dstat moreutils zip less socat wget                   && \
    apt-get install -y --no-install-recommends r-base

# zlib1g needed for bcl2fastq

# ==> FROM THIS

# https://github.com/lilab-bcb/cumulus/blob/707f52bb4a9834f3a8ef5bba03345ffd75f26dd1/docker/cellranger/7.1.0/Dockerfile
# Gives monitoring not using
# Can remove some apt/pip if wanted to, esp. stratocumulus

# https://cumulus.readthedocs.io/en/latest/bcl2fastq.html
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

# <== END from bcl2fastq

# already has gcloud installed, yay. Don't need to spend storage on AWS
# check storage
# apt update; apt install ncdu; ncdu /
RUN rm -rf /aws /usr/local/aws-cli

# For usage monitoring, added above
# RUN apt-get update && apt-get install -y --no-install-recommends dstat moreutils zip less socat

RUN apt-get -qq -y autoremove && \
    apt-get clean             && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

# need dnspython for hacks, see make_samplesheet.py
RUN pip install --no-cache-dir gspread dnspython

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

RUN R -e "install.packages(c('qs',        \
                             'glue',      \
                             'qpdf',      \
                             'purrr',     \
                             'dplyr',     \
                             'rlist',     \
                             'rdist',     \
                             'hdf5r',     \
                             'Seurat',    \
                             'dbscan',    \
                             'stringr',   \
                             'viridis',   \
                             'cowplot',   \
                             'ggplot2',   \
                             'magrittr',  \
                             'parallel',  \
                             'jsonlite',  \
                             'gridExtra', \
                             'matrixStats'), repos='http://cran.rstudio.com/')"

# in this directory:
# docker build -t slidetags_presquash .
# test out running docker
# docker run -it slidetags_presquash

# BUT TO PUSH:
# go into 01_for_squash and run there
