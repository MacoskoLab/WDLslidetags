version 1.0

# TODO timeout, or email if runs for more than a day (didn't die)
# kill 15?
# ps -e|grep -A500 tee|grep bash|head -n 1|awk '{print $1}'|xargs kill -9
# root@054a5a697e0f:/cromwell_root# ps -e
# PID TTY          TIME CMD
# 1 ?        00:00:00 bash
# 13 ?        00:00:00 tee
# 14 ?        00:00:00 tee
# 15 ?        00:00:00 bash
# 16 ?        00:00:00 socat
# 17 pts/0    00:00:00 bash
# 920 pts/0    00:00:00 ps

# TODO: upload resource logging and outs logging
# TODO make preemptible

# TODO: better matching, positioning, plots (.jl/.R), joinpath
# TODO remove unneeded files from the BCL (.tifs?)
# TODO record reference
# umi collapsing / chimerism
# check if the r worked

task read_sheet {
  input {
      String bcl
      Boolean run_mkfastq
      Boolean run_RNAcounts
      Boolean run_SBcounts
      Boolean run_spatial

      String bucket
      String docker
    }
    command <<<

      # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201    

      echo "true" > DONE

      # Indexes.csv, Transcriptomes.csv, Spatial.csv
      gsutil cp gs://~{bucket}/scripts/upload_for_google_key.json .
      gsutil cp gs://~{bucket}/scripts/read_sheet.py .
      python3 read_sheet.py ~{bcl} ~{bucket}
      # Checks on the file contents
      ~{run_mkfastq} && [[ "$(wc -l < Indexes.csv)" -lt 1 ]] && echo "ERROR: Indexes.csv is blank but run_mkfastq is true" && rm DONE
      ~{run_RNAcounts} && [[ "$(wc -l < RNAcounts.tsv)" -lt 1 ]] && echo "ERROR: RNAcounts.tsv is blank but run_RNAcounts is true" && rm DONE
      ~{run_SBcounts} && [[ "$(wc -l < SBcounts.tsv)" -lt 1 ]] && echo "ERROR: SBcounts.tsv is blank but run_SBcounts is true" && rm DONE
      ~{run_spatial} && [[ "$(wc -l < Spatial.tsv)" -lt 1 ]] && echo "ERROR: Spatial.tsv is blank but run_spatial is true" && rm DONE
      # Empty files if not needed
      ! ~{run_mkfastq} && echo -n > Indexes.csv
      ! ~{run_RNAcounts} && echo -n > RNAcounts.tsv
      ! ~{run_SBcounts} && echo -n > SBcounts.tsv
      ! ~{run_spatial} && echo -n > Spatial.tsv

      # Get the existence of important files
      gsutil ls gs://~{bucket}/01_BCLS/~{bcl} > /dev/null 2>&1
      exists_bcl=$?
      gsutil ls gs://~{bucket}/02_FASTQS/~{bcl} > /dev/null 2>&1
      exists_mkfastq=$?
      gsutil ls gs://~{bucket}/03_COUNTS/~{bcl} > /dev/null 2>&1
      exists_RNAcounts=$?
      gsutil ls gs://~{bucket}/04_SPATIAL/~{bcl} > /dev/null 2>&1
      exists_SBcounts=$?
      # Abort if the files are needed but don't exist
      ~{run_mkfastq} && [[ $exists_bcl -ne 0 ]] && echo "ERROR: RUNNING MKFASTQ BUT NO BCL" && rm DONE
      ~{run_RNAcounts} && ! ~{run_mkfastq} && [[ $exists_mkfastq -ne 0 ]] && echo "ERROR: RUNNING RNACOUNTS BUT NO MKFASTQ" && rm DONE
      ~{run_SBcounts} && ! ~{run_mkfastq} && [[ $exists_mkfastq -ne 0 ]] && echo "ERROR: RUNNING SBCOUNTS BUT NO MKFASTQ" && rm DONE
      ~{run_spatial} && ! ~{run_RNAcounts} && [[ $exists_RNAcounts -ne 0 ]] && echo "ERROR: RUNNING SPATIAL BUT NO RNACOUNTS" && rm DONE
      ~{run_spatial} && ! ~{run_SBcounts} && [[ $exists_SBcounts -ne 0 ]] && echo "ERROR: RUNNING SPATIAL BUT NO SBCOUNTS" && rm DONE
      # Abort if the files will be created but already exist (except spatial, easy to generate so will overwrite)
      ~{run_mkfastq} && [[ $exists_mkfastq -eq 0 ]] && echo "ERROR: MKFASTQ ALREADY DONE" && rm DONE
      ~{run_RNAcounts} && [[ $exists_counts -eq 0 ]] && echo "WARNING: COUNTS ALREADY EXISTS, WILL ONLY WRITE NEW LIBRARIES"

      # Get the size of the BCL, then make sure 3x size is below 6TB
      gsutil du -sc gs://~{bucket}/01_BCLS/~{bcl} | grep total | awk '{size=$1/1024/1024/1024 ; size=size*3 ; if (size<127) size=127 ; printf "%d\n", size+1}' > MKFASTQSIZE
      [[ $(cat MKFASTQSIZE) -gt 6000 ]] && echo "BCL is too large, please increase the size limit in the read_sheet task" && rm DONE

      echo 'END' # set the return code to 0

    >>>
    output {
      File Indexes                   = "Indexes.csv"
      Array[Array[String]] RNAcounts = read_tsv("RNAcounts.tsv")
      Array[Array[String]] SBcounts  = read_tsv("SBcounts.tsv")
      Array[Array[String]] Spatial   = read_tsv("Spatial.tsv")

      Int MKFASTQSIZE                = read_int("MKFASTQSIZE")
      Boolean DONE                   = read_boolean("DONE")
    }
    runtime {
      docker: docker
      memory: "5 GB"
      disks: "local-disk 10 HDD"
      cpu: 1
      preemptible: 0
    }
}

task mkfastq {
  input {
    String bcl
    File Indexes

    String bucket
    String docker
    Int MKFASTQSIZE
  }
  command <<<
    
    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

    export PATH="/software/cellranger-7.1.0/bin:$PATH"
    export PATH="/usr/local/bcl2fastq/bin:$PATH"
    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    echo "downloading BCL"
    mkdir BCL
    gcloud storage cp -r gs://~{bucket}/01_BCLS/~{bcl}/* ./BCL |& ts

    echo "running mkfastq"
    time stdbuf -oL -eL cellranger mkfastq                     \
      --run=BCL                                                \
      --id=mkfastq                                             \
      --csv=~{Indexes}                                         \
      --disable-ui |& ts | tee ./mkfastq.log

    echo "removing MAKE_FASTQS_CS"
    rm -rf ./mkfastq/MAKE_FASTQS_CS

    echo "checking for success"
    if [ -f mkfastq/outs/fastq_path/Reports/html/index.html ]
    then
      echo "success, uploading fastqs"
      gcloud storage cp -r mkfastq gs://~{bucket}/02_FASTQS/~{bcl}
      echo "true" > DONE
    else
      echo "FAILURE, CANNOT FIND: index.html"
    fi

    echo 'END' # set the return code to 0
  >>>
  output {
    Boolean DONEmkfastq = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "64 GB"
    disks: "local-disk ~{MKFASTQSIZE} LOCAL"
    cpu: 8
    preemptible: 0
  }
}

task compute_sizes {
  input {
      Array[Array[String]] array2D
      String bcl
      String bucket
      String docker
      Boolean DONEmkfastq
    }
    command <<<

      socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

      echo -n > SIZES

      while IFS=$'\t' read -r index rest_of_line
      do
        size=$(gsutil du -sce "*_I[0-9]_[0-9][0-9][0-9].fastq.gz" "gs://~{bucket}/02_FASTQS/~{bcl}/outs/fastq_path/*/$index_*_R*.fastq.gz" | grep total | awk '{print $1}')
        if [ "$size" -ne 0 ]; then
          echo $size | awk '{$1/1024/1024/1024 ; size=size*6+20 ; if (size<127) size=127 ; printf "%d\n", size+1}' >> SIZES
        else
          echo "CANNOT FIND THE FASTQS FOR $index"
          echo "0" >> SIZES
        fi
      done < ~{array2D}

    >>>
  output {
    Array[Int] SIZES  = read_lines("SIZES")
  }
  runtime {
    docker: docker
    memory: "5 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 0
  }
}





task spatial {
  input {
    String bcl
    File Spatial
    Int rowindex
    Boolean run_spatial

    String bucket
    String docker
    Int COUNTSSIZE
    Array[File] DONEcounts
  }
  command <<<

    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

    # OVERWRITE!

    if ~{run_spatial}
    then
      export PATH="/software/julia-1.8.5/bin:$PATH"

      gsutil cp gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/scripts/3M-february-2018.txt .
      gsutil cp gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/scripts/read_fastq.jl .
      gsutil cp gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/scripts/spatial.R .

      RNAIND=$(awk -v i="~{rowindex}" -F ',' 'NR == i {print $1}' ~{Spatial} | tr -d '\r\n')

      date > DONE

      julia read_fastq.jl ~{Spatial} ~{rowindex}

      echo "checking for .jl success"
      if [ -d $RNAIND ]
      then
        echo "uploading to 04_SPATIAL"
        gcloud storage cp -r $RNAIND gs://~{bucket}/04_SPATIAL/~{bcl}/$RNAIND
      else
        echo "FAILURE, CANNOT FIND: RNAINDEX/"
        rm DONE
      fi

      Rscript spatial.R

      echo "checking for .R success"
      if [ -f $RNAIND.qs ]
      then
        echo "uploading to 05_SEURATS"
        gcloud storage cp $RNAIND.qs gs://~{bucket}/05_SEURATS/~{bcl}/$RNAIND.qs
        gcloud storage cp $RNAIND/summary.pdf gs://~{bucket}/04_SPATIAL/~{bcl}/$RNAIND/summary.pdf
      else
        echo "FAILURE, CANNOT FIND: RNAINDEX.qs"
        rm DONE
      fi
    else
      echo "skipping spatial"
      date > DONE
    fi

    echo 'END' # set the return code to 0

  >>>
  output {
    File DONEspatial = "DONE"
  }
  runtime {
    docker: docker
    memory: "128 GB"
    disks: "local-disk ~{COUNTSSIZE} LOCAL"
    cpu: "16"
    preemptible: 0
  }
}




task counts {
  input {
    String index
    String transcriptome

    String bcl
    String bucket
    String docker
    Int COUNTSSIZE
  }
  command <<<

    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

    export PATH="/software/cellranger-7.1.0/bin:$PATH"
    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    echo "downloading FASTQs"
    mkdir -p ./mkfastq/outs/fastq_path/flowcell
    gcloud storage cp "gs://~{bucket}/02_FASTQS/~{bcl}/**/~{index}_*.fastq.gz" ./mkfastq/outs/fastq_path/flowcell |& ts

    echo "downloading reference"
    mkdir ./~{transcriptome}
    gcloud storage cp -r gs://~{bucket}/references/~{transcriptome}/* ./~{transcriptome}/ |& ts

    echo "running counts"
    time stdbuf -oL -eL cellranger count \
      --id=~{index}                      \
      --fastqs=mkfastq/outs/fastq_path   \
      --sample=~{index}                  \
      --transcriptome=~{transcriptome}   \
      --jobmode=local --disable-ui       \
      --nosecondary                      \
      --include-introns=true |& ts | tee -a ./counts.log
    echo "removing unnecessary files"
    rm -rf ~{index}/SC_RNA_COUNTER_CS
    echo "checking for success"
    if [ -f ~{index}/outs/metrics_summary.csv ]
    then
      echo "SUCCESS: uploading counts"
      gcloud storage cp -r $INDEX gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_COUNTS/~{bcl}/~{index}
      echo "true" > DONE
    else
      echo "FAILURE, CANNOT FIND: outs/metrics_summary.csv"
    fi

    echo 'END' # set the return code to 0

  >>>
  output {
    Boolean DONEcounts = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "64 GB"
    disks: "local-disk ~{COUNTSSIZE} LOCAL"
    cpu: "8"
    preemptible: 0
  }
}





