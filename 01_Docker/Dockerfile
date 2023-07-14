FROM cumulusprod/cellranger:7.1.0
SHELL ["/bin/bash", "-c"]

RUN apt-get update                                        && \
    apt-get install -y --no-install-recommends zlib1g-dev && \
    apt-get install -y --no-install-recommends dstat moreutils zip less socat

# zlib1g needed for bcl2fastq
# dstat..., needed for capturing the usage, etc. Not for bcl2fastq
# socat is optional but helpful for debugging with reverse shell

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

ADD ./make_samplesheet.py /make_samplesheet.py
ADD ./upload_for_google_key.json /upload_for_google_key.json

# need dnspython for hacks, see make_samplesheet.py
RUN pip install --no-cache-dir gspread dnspython


# in this directory:
# docker build -t slidetags_presquash .
# test out running docker
# docker run -it slidetags_presquash

# BUT TO PUSH:
# go into 01_for_squash and run there

