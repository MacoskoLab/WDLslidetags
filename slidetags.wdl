version 1.0

# TODO timeout, or email if runs for more than a day
# TODO make preemptible
# TODO remove unneeded files from the BCL (.tifs?)

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
      gsutil cp gs://~{bucket}/scripts/google_key.json .
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
      # Abort if the files will be created but already exist (except SBcounts and spatial, easy to generate so will overwrite)
      ~{run_mkfastq} && [[ $exists_mkfastq -eq 0 ]] && echo "ERROR: MKFASTQ ALREADY DONE" && rm DONE
      ~{run_RNAcounts} && [[ $exists_RNAcounts -eq 0 ]] && echo "WARNING: COUNTS ALREADY EXISTS, WILL ONLY WRITE NEW LIBRARIES"

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
    Int disksize
  }
  command <<<
    
    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

    dstat --time --cpu --disk --mem --io > mkfastq.usage 2>&1 &

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

    echo "uploading logs"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    ( echo; echo "CPU INFO:"; lscpu ) >> mkfastq.log
    gcloud storage cp mkfastq.log gs://~{bucket}/logs/~{bcl}/mkfastq.log
    gcloud storage cp mkfastq.usage gs://~{bucket}/logs/~{bcl}/mkfastq.usage

    echo 'END' # set the return code to 0
  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "64 GB"
    disks: "local-disk ~{disksize} LOCAL"
    cpu: 8
    preemptible: 0
  }
}

task compute_sizes {
  input {
      Array[Array[String]] preCOUNTS # contains a .tsv where the first column is the index

      String bcl
      String bucket
      String docker
      Boolean DONE
    }
    command <<<

      # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

      # ensure files are empty before appending
      echo -n > SIZES
      echo -n > COUNTS.tsv

      # get the size of the fastqs and only run if the fastqs exist
      while IFS=$'\t' read -r index rest_of_line
      do
        size=$(gsutil du -sce "*_I[0-9]_[0-9][0-9][0-9].fastq.gz" "gs://~{bucket}/02_FASTQS/~{bcl}/**/$index_*_R*.fastq.gz" | grep total | awk '{print $1}')
        if [ "$size" -gt 0 ]; then
          echo $size | awk '{size=$1/1024/1024/1024 ; printf "%d\n", size+1}' >> SIZES
          echo -e "$index\t$rest_of_line" >> COUNTS.tsv
        else
          echo "CANNOT FIND THE FASTQS FOR $index, will not run"
        fi
      done < ~{write_tsv(preCOUNTS)}

    >>>
  output {
    Array[Int] SIZES  = read_lines("SIZES")
    Array[Array[String]] COUNTS = read_tsv("COUNTS.tsv")
  }
  runtime {
    docker: docker
    memory: "5 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 0
  }
}

task RNAcounts {
  input {
    String index
    String transcriptome

    String bcl
    String bucket
    String docker
    Int disksize
  }
  command <<<

    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

    dstat --time --cpu --disk --mem --io > RNAcounts.usage 2>&1 &

    export PATH="/software/cellranger-7.1.0/bin:$PATH"
    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    echo "downloading FASTQs"
    mkdir -p ./mkfastq/outs/fastq_path/flowcell
    gcloud storage cp "gs://~{bucket}/02_FASTQS/~{bcl}/**/~{index}_*_R*.fastq.gz" ./mkfastq/outs/fastq_path/flowcell |& ts

    echo "downloading reference"
    mkdir ./~{transcriptome}
    gcloud storage cp -r gs://~{bucket}/references/~{transcriptome}/* ./~{transcriptome}/ |& ts

    echo "running counts"
    time stdbuf -oL -eL cellranger count \
      --id=~{index}                      \
      --fastqs=mkfastq/outs/fastq_path   \
      --sample=~{index}                  \
      --transcriptome=~{transcriptome}   \
      --create-bam=true                  \
      --jobmode=local --disable-ui       \
      --nosecondary                      \
      --include-introns=true |& ts | tee -a ./RNAcounts.log
    
    echo "removing unnecessary files"
    rm -rf ~{index}/SC_RNA_COUNTER_CS

    echo "checking for success"
    if [ -f ~{index}/outs/metrics_summary.csv ]
    then
      echo "SUCCESS: uploading counts"
      gcloud storage cp -r ~{index} gs://~{bucket}/03_COUNTS/~{bcl}/~{index}
      echo "true" > DONE
    else
      echo "FAILURE, CANNOT FIND: outs/metrics_summary.csv"
    fi

    echo "uploading logs"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    ( echo; echo "CPU INFO:"; lscpu ) >> RNAcounts.log
    gcloud storage cp RNAcounts.log gs://~{bucket}/logs/~{bcl}/RNAcounts-~{index}.log
    gcloud storage cp RNAcounts.usage gs://~{bucket}/logs/~{bcl}/RNAcounts-~{index}.usage

    echo 'END' # set the return code to 0

  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "64 GB"
    disks: "local-disk ~{disksize} LOCAL"
    cpu: "8"
    preemptible: 0
  }
}

task RNAcounts_FFPE {
  input {
    String index
    String transcriptome
    String probeset

    String bcl
    String bucket
    String docker
    Int disksize
  }
  command <<<

    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

    dstat --time --cpu --disk --mem --io > RNAcounts.usage 2>&1 &

    export PATH="/software/cellranger-7.1.0/bin:$PATH"
    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    echo "downloading FASTQS"
    mkdir ./FASTQS
    gcloud storage cp "gs://~{bucket}/02_FASTQS/~{bcl}/**/~{index}_*_R*.fastq.gz" ./FASTQS/ |& ts

    echo "downloading reference"
    mkdir ~{transcriptome}
    gcloud storage cp -r gs://~{bucket}/references/~{transcriptome}/* ~{transcriptome}/ |& ts

    echo "downloading probeset"
    gcloud storage cp gs://~{bucket}/probesets/~{probeset} . |& ts

    echo "creating the config file"
    echo -n > multiconfig.csv
    echo "[gene-expression]" >> multiconfig.csv
    echo "reference,$(readlink -f ~{transcriptome})" >> multiconfig.csv
    echo "probe-set,$(readlink -f ~{probeset})" >> multiconfig.csv
    echo "no-secondary,true" >> multiconfig.csv
    echo "no-bam,false" >> multiconfig.csv
    echo "include-introns,true" >> multiconfig.csv
    echo "" >> multiconfig.csv
    echo "[libraries]" >> multiconfig.csv
    echo "fastq_id,fastqs,feature_types" >> multiconfig.csv
    echo "~{index},$(readlink -f FASTQS),Gene Expression" >> multiconfig.csv

    echo "running counts"
    time stdbuf -oL -eL cellranger multi \
      --id=~{index} \
      --csv=multiconfig.csv \
      --jobmode=local --disable-ui |& ts | tee -a ./RNAcounts.log

    echo "removing unnecessary files"
    rm -rf ~{index}/SC_MULTI_CS

    echo "checking for success"
    if [ -f ~{index}/outs/per_sample_outs/~{index}/metrics_summary.csv ]
    then
      echo "SUCCESS: uploading counts"
      gcloud storage cp -r ~{index} gs://~{bucket}/03_COUNTS/~{bcl}/~{index}
      echo "true" > DONE
    else
      echo "FAILURE, CANNOT FIND: ~{index}/outs/per_sample_outs/~{index}/metrics_summary.csv"
    fi

    echo "uploading logs"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    ( echo; echo "CPU INFO:"; lscpu ) >> RNAcounts.log
    gcloud storage cp RNAcounts.log gs://~{bucket}/logs/~{bcl}/RNAcounts-~{index}.log
    gcloud storage cp RNAcounts.usage gs://~{bucket}/logs/~{bcl}/RNAcounts-~{index}.usage

    echo 'END' # set the return code to 0

  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "64 GB"
    disks: "local-disk ~{disksize} LOCAL"
    cpu: "8"
    preemptible: 0
  }
}

task SBcounts {
  input {
    String index
    String pucklist

    String bcl
    String bucket
    String docker
    Int disksize
  }
  command <<<

    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

    dstat --time --cpu --disk --mem --io > SBcounts.usage 2>&1 &
    
    export PATH="/software/julia-1.8.5/bin:$PATH"
    gcloud storage cp "gs://~{bucket}/scripts/3M-february-2018.txt" .
    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    echo "downloading FASTQs"
    mkdir -p ./fastqs
    gcloud storage cp "gs://~{bucket}/02_FASTQS/~{bcl}/**/~{index}_*_R*.fastq.gz" ./fastqs |& ts

    echo "downloading pucks"
    mkdir ./pucks
    IFS=',' read -ra pucklist <<< ~{pucklist}
    for puck in "${pucklist[@]}"; do
        gcloud storage cp gs://~{bucket}/pucks/$puck ./pucks |& ts
    done
    echo "$(ls -1 pucks | wc -l) puck(s) downloaded"

    echo "running read_fastq.jl"
    gsutil cp gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/scripts/read_fastq.jl .
    time stdbuf -oL -eL julia read_fastq.jl fastqs pucks |& ts | tee -a ./SBcounts.log
    
    echo "checking for success"
    if [ -f SBcounts.h5 ]
    then
      echo "SUCCESS: uploading SBcounts.h5"
      gcloud storage cp SBcounts.h5 gs://~{bucket}/04_SPATIAL/~{bcl}/~{index}/SBcounts.h5
      echo "true" > DONE
    else
      echo "FAILURE, CANNOT FIND: SBcounts.h5"
    fi

    echo "uploading logs"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    ( echo; echo "CPU INFO:"; lscpu ) >> SBcounts.log
    gcloud storage cp SBcounts.log gs://~{bucket}/logs/~{bcl}/SBcounts-~{index}.log
    gcloud storage cp SBcounts.usage gs://~{bucket}/logs/~{bcl}/SBcounts-~{index}.usage

    echo 'END' # set the return code to 0

  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "32 GB"
    disks: "local-disk ~{disksize} LOCAL"
    cpu: "1"
    preemptible: 0
  }
}

task spatial {
  input {
    String RNApath
    String SBpath

    String bcl
    String bucket
    String docker
    Int DONE_RNA_SB
  }
  command <<<

    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

    dstat --time --cpu --disk --mem --io > spatial.usage 2>&1 &

    gcloud storage cp "gs://~{bucket}/scripts/spatial.R" .

    time stdbuf -oL -eL Rscript spatial.R ~{RNApath} ~{SBpath}

    echo "checking for success"
    if [ -f seurat.qs ]
    then
      echo "SUCCESS: uploading seurat.qs, coords.csv, and summary.pdf"
      gcloud storage cp seurat.qs coords.csv summary.pdf gs://~{bucket}/05_SEURATS/~{bcl}/~{basename(RNApath)}/
      echo "true" > DONE
    else
      echo "FAILURE, CANNOT FIND: seurat.qs"
    fi

    echo "uploading log"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    gcloud storage cp spatial.usage gs://~{bucket}/logs/~{bcl}/spatial-~{basename(RNApath)}.usage

    echo 'END' # set the return code to 0

  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "64 GB"
    disks: "local-disk 100 HDD"
    cpu: "20"
    preemptible: 0
  }
}

# import "https://raw.githubusercontent.com/MacoskoLab/WDLslidetags/main/tasks.wdl" as tasks

workflow tagspipeline {
  String pipeline_version = "1.0.0"
  input {
    String bcl
    Boolean run1_mkfastq   = true
    Boolean run2_RNAcounts = true
    Boolean run3_SBcounts  = true
    Boolean run4_spatial   = true
    String zz_bucket = "fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d"
    String zz_docker = "us-central1-docker.pkg.dev/velina-208320/jonah-slidetags/img:latest"
  }

  parameter_meta {
  }

  call read_sheet as read_sheet {
    input:
      bcl = bcl,
      run_mkfastq   = run1_mkfastq,
      run_RNAcounts = run2_RNAcounts,
      run_SBcounts  = run3_SBcounts,
      run_spatial   = run4_spatial,
      
      bucket = zz_bucket,
      docker = zz_docker
  }

  ### RUN MKFASTQ ###
  if (run1_mkfastq) {
    call mkfastq as mkfastq {
      input:
        Indexes = read_sheet.Indexes,
        
        bcl = bcl,
        bucket = zz_bucket,
        docker = zz_docker,
        disksize = read_sheet.MKFASTQSIZE
    }
  }
  Boolean DONE_mkfastq = select_first([mkfastq.DONE, true])

  ### COMPUTE DISK SIZES ###
  if (run2_RNAcounts) {
    call compute_sizes as compute_RNAsizes {
      input:
        preCOUNTS = read_sheet.RNAcounts,
        DONE = DONE_mkfastq,

        bcl = bcl,
        bucket = zz_bucket,
        docker = zz_docker
    }
  }
  Array[Int] RNAsizes_array = select_first([compute_RNAsizes.SIZES, []])
  Array[Array[String]] RNAcounts_array = select_first([compute_RNAsizes.COUNTS, []])

  if (run3_SBcounts) {
    call compute_sizes as compute_SBsizes {
      input:
        preCOUNTS = read_sheet.SBcounts,
        DONE = DONE_mkfastq,

        bcl = bcl,
        bucket = zz_bucket,
        docker = zz_docker
    }
  }
  Array[Int] SBsizes_array = select_first([compute_SBsizes.SIZES, []])
  Array[Array[String]] SBcounts_array = select_first([compute_SBsizes.COUNTS, []])

  ### RUN RNACOUNTS ###
  scatter(pair in zip(RNAcounts_array, RNAsizes_array)) {
    Array[String] RNAparams = pair.left
    Int basernasize = pair.right * 6 + 20
    Int RNASIZE = if (basernasize < 64) then 64 else if (basernasize > 2048) then 2048 else basernasize
    if (length(RNAparams) == 2) {
      call RNAcounts as RNAcounts {
        input:
          index = RNAparams[0],
          transcriptome = RNAparams[1],

          bcl = bcl,
          bucket = zz_bucket,
          docker = zz_docker,
          disksize = RNASIZE
      }
    }
    if (length(RNAparams) == 3) {
      call RNAcounts_FFPE as RNAcounts_FFPE {
        input:
          index = RNAparams[0],
          transcriptome = RNAparams[1],
          probeset = RNAparams[2],

          bcl = bcl,
          bucket = zz_bucket,
          docker = zz_docker,
          disksize = RNASIZE
      }
    }
  }
  Array[Boolean] DONE_RNAcounts = flatten([select_all(RNAcounts.DONE), select_all(RNAcounts_FFPE.DONE)])

  ### RUN SBCOUNTS ###
  scatter(pair in zip(SBcounts_array, SBsizes_array)) {
    Array[String] SBparams = pair.left
    Int basesbsize = pair.right * 3
    Int SBSIZE = if (basesbsize < 64) then 64 else if (basesbsize > 2048) then 2048 else basesbsize
    call SBcounts as SBcounts {
      input:
        index = SBparams[0],
        pucklist = SBparams[1],

        bcl = bcl,
        bucket = zz_bucket,
        docker = zz_docker,
        disksize = SBSIZE
    }
  }
  Array[Boolean] DONE_SBcounts = select_all(SBcounts.DONE)

  Int DONE_RNA_SB = length(DONE_RNAcounts) + length(DONE_SBcounts)

  ### RUN SEURAT GENERATION ###
  if (run4_spatial) {
    scatter(row in read_sheet.Spatial) {
      call spatial as spatial {
        input:
          RNApath = row[0],
          SBpath = row[1],

          bcl = bcl,
          bucket = zz_bucket,
          docker = zz_docker,
          DONE_RNA_SB = DONE_RNA_SB
      }
    }
  }

  output {
  }
}
