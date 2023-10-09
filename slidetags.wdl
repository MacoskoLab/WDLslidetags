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

# TODO: make the sheet support worksheets other than Tags (read_sheet.py)
# TODO: other counts types - currently assumes that the input is Tags
# TODO: upload resource logging and outs logging
# TODO: better matching, positioning, plots (.jl/.R), joinpath
# TODO make preemptible
# TODO remove unneeded files from the BCL (.tifs?)
# TODO record reference
# TODO: check docker size: apt install ncdu; ncdu
# TODO: add support for commas in the spreadsheet
# TODO: "_" vs "_S" and/in file size calculation

# real TODOs:
# umi collapsing / chimerism
# check if the r worked

task read_sheet {
  input {
      String bcl
      Boolean run_mkfastq
      Boolean run_counts
      Boolean run_spatial

      String bucket
      String docker
    }
    command <<<

      # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201    

      # Indexes.csv, Transcriptomes.csv, Spatial.csv
      gsutil cp gs://~{bucket}/scripts/upload_for_google_key.json .
      gsutil cp gs://~{bucket}/scripts/read_sheet.py .
      python3 read_sheet.py ~{bcl} ~{bucket}

      date > DONE

      # Get the size of the BCL, then make sure 3x size is below 6TB
      gsutil du -sc gs://~{bucket}/01_BCLS/~{bcl} | grep total | awk '{size=$1/1024/1024/1024 ; size=size*3 ; if (size<127) size=127 ; printf "%d\n", size+1}' > MKFASTQSIZE
      awk 'NR==1 && $1>6000 { system("echo TOOBIG ; rm DONE") }' MKFASTQSIZE

      # Get the existence of important files
      gsutil ls gs://~{bucket}/01_BCLS/~{bcl} > /dev/null 2>&1
      exists_bcl=$?
      gsutil ls gs://~{bucket}/02_FASTQS/~{bcl} > /dev/null 2>&1
      exists_mkfastq=$?
      gsutil ls gs://~{bucket}/03_COUNTS/~{bcl} > /dev/null 2>&1
      exists_counts=$?
      # Remove files if not needed
      ! ~{run_counts} && > Transcriptomes.csv
      ! ~{run_spatial} && > Spatial.csv
      # Abort if the files are needed but don't exist
      ~{run_mkfastq} && [[ $exists_bcl -ne 0 ]] && echo "ERROR: NO BCL" && rm DONE
      ~{run_counts} && ! ~{run_mkfastq} && [[ $exists_mkfastq -ne 0 ]] && echo "ERROR: NO MKFASTQ" && rm DONE
      ~{run_spatial} && ! ~{run_mkfastq} && [[ $exists_mkfastq -ne 0 ]] && echo "ERROR: NO MKFASTQ" && rm DONE
      ~{run_spatial} && ! ~{run_counts} && [[ $exists_counts -ne 0 ]] && echo "ERROR: NO COUNTS" && rm DONE
      # Abort if the files will be created but already exist (except spatial, easy to generate so will overwrite)
      ~{run_mkfastq} && [[ $exists_mkfastq -eq 0 ]] && echo "ERROR: MKFASTQ ALREADY DONE" && rm DONE
      ~{run_counts} && [[ $exists_counts -eq 0 ]] && echo "ERROR: COUNTS ALREADY DONE" && rm DONE
      # Checks on the file contents
      ~{run_mkfastq} && [[ "$(wc -l < Indexes.csv)" -lt 2 ]] && echo "ERROR: Indexes.csv is blank but mkfastq is true" && rm DONE
      ~{run_counts} && [[ "$(wc -l < Transcriptomes.csv)" -lt 1 ]] && echo "ERROR: Transcriptomes.csv is blank but counts is true" && rm DONE
      ~{run_spatial} && [[ "$(wc -l < Spatial.csv)" -lt 1 ]] && echo "ERROR: Spatial.csv is blank but spatial is true" && rm DONE

      # Get files ready for counts and spatial
      awk 'BEGIN{a=0} {a++; print a}' Transcriptomes.csv > COUNTSROWS
      awk 'BEGIN{a=0} {a++; print a}' Spatial.csv > SPATIALROWS

      echo 'END' # set the return code to 0

    >>>
    output {
      File Indexes           = "Indexes.csv"
      File Counts            = "Transcriptomes.csv"
      File Spatial           = "Spatial.csv"

      Array[Int] COUNTSROWS  = read_lines("COUNTSROWS")
      Array[Int] SPATIALROWS = read_lines("SPATIALROWS")

      Int MKFASTQSIZE        = read_int("MKFASTQSIZE")
      File DONE              = "DONE"
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
    Boolean run_mkfastq

    String bucket
    String docker
    Int MKFASTQSIZE
  }
  command <<<
    
    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

    if ~{run_mkfastq}
    then
      export PATH="/software/cellranger-7.1.0/bin:$PATH"
      export PATH="/usr/local/bcl2fastq/bin:$PATH"
      gcloud config set storage/process_count 16
      gcloud config set storage/thread_count  2

      echo "downloading BCL"
      mkdir BCL
      gcloud storage cp -r gs://~{bucket}/01_BCLS/~{bcl}/* BCL |& ts

      echo "running mkfastq" 
      time stdbuf -oL -eL cellranger mkfastq                     \
        --run=BCL                                                \
        --id=mkfastq                                             \
        --csv=~{Indexes}                                         \
        --jobmode=local --disable-ui  |& ts | tee ./mkfastq.log

      echo "removing unnecessary files"
      rm -rf ./mkfastq/MAKE_FASTQS_CS

      echo "checking for success"
      if [ -f mkfastq/outs/fastq_path/Reports/html/index.html ]
      then
        echo "uploading fastqs"
        gcloud storage cp -r mkfastq gs://~{bucket}/02_FASTQS/~{bcl}
        date > DONE
      else
        echo "FAILURE, CANNOT FIND: index.html"
      fi
    else
      echo "skipping mkfastq"
      date > DONE
    fi

    # At this point, assert there are either pipeline-generated fastqs or user-input fastqs
    [ $(gsutil ls -r gs://~{bucket}/02_FASTQS/~{bcl} | grep -F ".fastq.gz" | grep -v "Undetermined" | wc -l) -eq 0 ] && echo "NOFASTQS" && rm DONE

    # Calculate the amount of disk space to use for counts/spatial
    gsutil du gs://~{bucket}/02_FASTQS/~{bcl} | grep -F ".fastq.gz" | grep -v "Undetermined" > SIZES
    Rscript -e "
      library(dplyr) ; library(purrr)
      df = read.table('SIZES', header=F, sep='', stringsAsFactors=F)
      df[[2]] = df[[2]] %>% basename %>% stringr::str_split('_S') %>% map_chr(pluck(1))
      df %>% group_by(V2) %>% summarise(total=sum(V1)) %>% pull(total) %>% max %>% cat(sep='\n')
    " | awk '{size=$1/1024/1024/1024 ; size=size*6+20 ; if (size<127) size=127 ; printf "%d\n", size+1}' > COUNTSSIZE
    
    # Assert that the disk size is below 1TB
    awk 'NR==1 && $1>1000 { system("echo TOOBIG ; rm DONE") }' COUNTSSIZE

    echo 'END' # set the return code to 0
  >>>
  output {
    Int COUNTSSIZE  = read_int("COUNTSSIZE")
    File DONEmkfastq = "DONE"
  }
  runtime {
    docker: docker
    memory: "64 GB"
    disks: "local-disk ~{MKFASTQSIZE} LOCAL"
    cpu: 8
    preemptible: 0
  }
}

task counts {
  input {
    String bcl
    File Counts
    Int rowindex
    Boolean run_counts

    String bucket
    String docker
    Int COUNTSSIZE
  }
  command <<<

    # socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9201

    if ~{run_counts}
    then
      export PATH="/software/cellranger-7.1.0/bin:$PATH"
      gcloud config set storage/process_count 16
      gcloud config set storage/thread_count  2

      echo "get the technique, index, and transcriptome from the sheet"
      line=$(awk -v i="~{rowindex}" -F ',' 'NR == i {print $1 ,$2, $3}' ~{Counts})
      read TECH INDEX REF <<< "$line"
      REF=$(echo $REF|tr -d '\r\n')

      echo "downloading FASTQs"
      mkdir -p ./mkfastq/outs/fastq_path/flowcell
      gcloud storage cp "gs://~{bucket}/02_FASTQS/~{bcl}/**/$INDEX*.fastq.gz" ./mkfastq/outs/fastq_path/flowcell |& ts

      echo "downloading reference"
      mkdir ./ref_folder
      gcloud storage cp -r gs://~{bucket}/references/$REF/* ./ref_folder/ |& ts

      date > DONE

      echo "running counts (with introns)"
      time stdbuf -oL -eL cellranger count \
        --id=$INDEX                        \
        --fastqs=mkfastq/outs/fastq_path   \
        --sample=$INDEX                    \
        --transcriptome=ref_folder         \
        --jobmode=local --disable-ui       \
        --nosecondary                      \
        --include-introns=true |& ts | tee -a ./counts.log
      echo "removing unnecessary files"
      rm -rf $INDEX/SC_RNA_COUNTER_CS
      echo "checking for success"
      if [ -f $INDEX/outs/metrics_summary.csv ]
      then
        echo "SUCCESS: uploading counts"
        gcloud storage cp -r $INDEX gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_COUNTS/~{bcl}/$INDEX
      else
        echo "FAILURE"
        rm DONE
      fi

      rm -rf $INDEX

      echo "running counts (without introns)"
      time stdbuf -oL -eL cellranger count \
        --id=$INDEX                        \
        --fastqs=mkfastq/outs/fastq_path   \
        --sample=$INDEX                    \
        --transcriptome=ref_folder         \
        --jobmode=local --disable-ui       \
        --nosecondary                      \
        --include-introns=false |& ts | tee -a ./counts.log
      echo "removing unnecessary files"
      rm -rf $INDEX/SC_RNA_COUNTER_CS
      echo "checking for success"
      if [ -f $INDEX/outs/metrics_summary.csv ]
      then
        echo "SUCCESS: uploading counts"
        gcloud storage cp -r $INDEX gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_COUNTS_NOINTRONS/~{bcl}/$INDEX
      else
        echo "FAILURE"
        rm DONE
      fi

    else
      echo "skipping counts"
      date > DONE
    fi

    echo 'END' # set the return code to 0

  >>>
  output {
    File DONEcounts = "DONE"
  }
  runtime {
    docker: docker
    memory: "64 GB"
    disks: "local-disk ~{COUNTSSIZE} LOCAL"
    cpu: "8"
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

workflow tagspipeline {
  String pipeline_version = "1.0.0"
  input {
    String bcl
    Boolean run_mkfastq = true
    Boolean run_counts  = true
    Boolean run_spatial = true
    String zz_bucket = "fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d"
    String zz_docker = "us-central1-docker.pkg.dev/velina-208320/jonah-slidetags/img:latest"
  }

  parameter_meta {
  }

  call read_sheet {
    input:
      bcl = bcl,
      run_mkfastq = run_mkfastq,
      run_counts  = run_counts,
      run_spatial = run_spatial,
      
      bucket = zz_bucket,
      docker = zz_docker
  }

  call mkfastq {
    input:
      bcl = bcl,
      Indexes = read_sheet.Indexes,
      run_mkfastq = run_mkfastq,
      
      bucket = zz_bucket,
      docker = zz_docker,
      MKFASTQSIZE = read_sheet.MKFASTQSIZE
  }

  scatter(rowindex in read_sheet.COUNTSROWS) {
    call counts {
      input:
        bcl = bcl,
        Counts = read_sheet.Counts,
        rowindex = rowindex,
        run_counts = run_counts,

        bucket = zz_bucket,
        docker = zz_docker,
        COUNTSSIZE = mkfastq.COUNTSSIZE
    }
  }

  scatter(rowindex in read_sheet.SPATIALROWS) {
    call spatial {
      input:
        bcl = bcl,
        Spatial = read_sheet.Spatial,
        rowindex = rowindex,
        run_spatial = run_spatial,

        bucket = zz_bucket,
        docker = zz_docker,
        COUNTSSIZE = mkfastq.COUNTSSIZE,
        DONEcounts = counts.DONEcounts
    }
  }
  output {
  }
}
