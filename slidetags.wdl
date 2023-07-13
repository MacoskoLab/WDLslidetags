version 1.0
# TODO instead of bcl_dir__noendslash being a full path, harcode prefix and just take folder. 
# TODO back up index file
# str replace ending slash when needed to be sure
# TODO timeout
# see what files in bcl need, inclue .tif?

# export FLOWCELL=230310_SL-NVN_0914_AHCWK5DSX5
# gsutil ls gs://macosko_data/slidetags/flowcells/$FLOWCELL > /tmp/t
# date > d
# gcloud storage cp d  gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/01_INBCLS/${FLOWCELL}/d
# cat /tmp/t|grep -v cellranger|xargs -I@ echo 'gcloud storage cp -r -n @ gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/01_INBCLS/${FLOWCELL}'|parallel -j 10

# TODO remove MAKE_FASTQS_CS/
# TODO gzip genes.gtf

# TODO to kill from inside
# TODO can use size and then set localizationOptional to false
# TODO, no use size as output and pass along, so just takes directory. Or read_float and write du -shc, so get recursive?
# TODO timeout manually

# TODO do mkfastq each lane separate, splay out, maybe only needed for big ones
# so meh

# TODO check if already have gcloud recently on docker
# and check gcs_transfer.sh for good flags

# TODO if less than 3 hours or 2   TB, then preemtible 2x?

# TODO cell ranger disable analysis


# TODO AHHH DIDN'T DIE
# https://app.terra.bio/#workspaces/testmybroad/Slide-tags/job_history/afad731d-e2dc-424a-9d6f-8e54233afe59

# TODO if running more than 5 hours, email me with curl
# https://api.firecloud.org/#/Submissions/listSubmissions -> see running -> get time, and email. Check every day, if more than 24 hours old!
# https://support.terra.bio/hc/en-us/articles/360042259232-Manage-data-automate-workflows-with-the-FISS-API



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


task getBCLSize {
    input {
        String bcl_dir__noendslash
        String docker
        String testbalah
    }

    command <<<
      echo ~{bcl_dir__noendslash}
      # https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md
      # Apparently terra says "GB" but means GiB, so need 1024
      gsutil du -sc ~{bcl_dir__noendslash} | grep total | awk '{print $1/1024/1024/1024}' > BCL_GiB_file

      # as a check so not expensive in next step, make sure less than 1.5 TB
      # which seems huuge
      awk 'NR==1 && $1<1500 { system("date > DONE_getBCLSize") }' BCL_GiB_file
    >>>

    output {
      Float GiB_size_bcl              = read_float("BCL_GiB_file")
      File  DONE_mkfastq             = "DONE_getBCLSize"
    }
    runtime {
      docker: docker
      memory: "5 GB"
      disks: "local-disk 10 HDD"
      cpu: "1"
      preemptible: 0
    }
}

task mkfastq {
    input {
        String bcl_dir__noendslash

        String local_workspace_bucket__noendslash
        String output_fastq_path

        Int zz_bcl_disksize
        Int zz_bcl_numCores
        Int zz_bcl_ram_gb
        String docker

    }

    command <<<
      socat exec:'bash -li',pty,stderr,setsid,sigint,sane tcp:167.172.130.57:9006

      touch usage.csv; dstat -t --cpu --mem --disk --io --freespace > usage.csv 2>&1 &
      cat /proc/cpuinfo > cpuinfo
      df -h . >> cpuinfo
      df -H . >> cpuinfo
      df -H

      gsutil -q stat '~{output_fastq_path}/*'
      output_gs_status=$?
      # 1 = does not exist
      if [[ $output_gs_status == 0 ]]; then
        echo "output fastq directory: $OUTPUT_GS already exists, exiting!"
        exit 1
      else
        echo "File does not exist, good"
      fi

      gsutil -q stat '~{bcl_dir__noendslash}/*'
      output_gs_status=$?
      # 1 = does not exist
      if [[ $output_gs_status == 0 ]]; then
        echo "Has bcl directory, good"
      else
        echo "bcl directory does not exist"
        exit 1
      fi

      # Seems like it could have weird consequences, don't mess with resolv.conf
      # maybe because includes google.com and is cached, but isn't supposed to
      # cp /etc/resolv.conf .
      # echo "nameserver 8.8.8.8"  > /etc/resolv.conf
      # echo "nameserver 8.8.4.4" >> /etc/resolv.conf
      # ...
      # cp resolv.conf /etc/resolv.conf

      BCL_BASENAME=`basename ~{bcl_dir__noendslash}`
      python /make_samplesheet.py $BCL_BASENAME
      cat Indexes.csv |awk -F, '{print $2}' | grep RNA | grep -v ^Sample$ | sort -u > justRNASampleNames

      cat Indexes.csv
      if grep -q Lane Indexes.csv; then
         echo "Created samplesheet happily"
      else
         echo "Failed creating samplesheet"
         exit 1
      fi


      gcloud config set storage/process_count 16
      gcloud config set storage/thread_count  2
      # could mess with this, depeneding on big vs small files
      # copy_chunk_size

      # for debugging
      # mkdir a;
      # # rm ~/.config/gcloud/surface_data/storage/tracker_files/*
      # gcloud storage cp --verbosity=debug -r 'gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/01_INBCLS/221217_VL00297_113_AACHKNMM5/Data/Intensities/BaseCalls/L001/C5*' a
      # gcloud storage cp -r 'gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/01_INBCLS/221217_VL00297_113_AACHKNMM5/Data/Intensities/BaseCalls/L001/C5*' a

      # gcloud config set storage/process_count 1
      # gcloud config set storage/thread_count 10
      # gcloud storage cp -r 'gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/01_INBCLS/221217_VL00297_113_AACHKNMM5/FocusModelGeneration/Cycle1/*' a


      # # speeds on smaller example
      # 588
      #  5/2 =443
      # 10/2 =639
      # 10/10=443
      # 16/1 =676,647,632
      # 16/2 =669,662,655,652
      # 16/3 =647,688,630
      # 16/5 =451,554
      # 25/2 =557
      # 16/5 =504

      mkdir ./bcl_folder
      # TODO don't need logs, probably not thumbnails either
      # Don't want nested folders, copy directly into bcl_folder
      gcloud storage cp -r '~{bcl_dir__noendslash}/*' ./bcl_folder/ |& ts

      df -H

      mkdir ./fastq_out

      # --localcores=~{zz_bcl_numCores} --localmem=~{zz_bcl_ram_gb} \

      # actual command :-)
      time stdbuf -oL -eL cellranger mkfastq                          \
          --jobmode=local --disable-ui                                \
          --id=fastq_out                                              \
          --run=bcl_folder                                            \
          --csv=Indexes.csv |& ts | tee ./mkfastq.log

      df -H

      # Simple check to make sure fastq_out worked
      if [ -f fastq_out/outs/fastq_path/Reports/html/index.html ] ; then date > DONE_mkfastq; fi

      # Don't need, often small. When errors can be big but still 99% don't want
      rm -rf fastq_out/MAKE_FASTQS_CS

      date > d
      gcloud storage cp d ~{output_fastq_path}/d
      # Big upload
      gcloud storage cp -r fastq_out/* ~{output_fastq_path} |& ts

      # Make extra sure dstat is dead
      touch usage.csv; pkill dstat; sleep 5s; pkill dstat; pkill dstat;
      gzip -9 usage.csv

      zip -j mkfastq_logs.zip ./mkfastq.log usage.csv.gz cpuinfo
    >>>

    output {
      # fails in order AFAIK, so put DONE_ last so still get usage
      File mkfastq_logs_zip            = "mkfastq_logs.zip"
      String output_fastq_path_copy    = output_fastq_path
      File Indexes                     = "Indexes.csv"
      Array[String] justRNASampleNames = read_lines("justRNASampleNames")
      File DONE_mkfastq                = "DONE_mkfastq"

    }
    runtime {
      # Short and fun, tested and almost perfectly fastq size
      # Takes no memory, might as well make preemtible and HDD
      docker: docker
      # TODO is it GB or GiB
      memory: "~{zz_bcl_ram_gb} GB"
      disks: "local-disk ~{zz_bcl_disksize} LOCAL"
      cpu: "~{zz_bcl_numCores}"
      preemptible: 0
    }
}

task runcounts {
    input {
      String fastq_path

      String ref_path__noendslash

      String output_count_path

      String index_name

      Int zz_count_disksize
      Int zz_count_numCores
      Int zz_count_ram_gb

      String docker
    }

    command <<<

      touch usage.csv; dstat -t --cpu --mem --disk --io --freespace > usage.csv 2>&1 &
      cat /proc/cpuinfo > cpuinfo
      df -h . >> cpuinfo
      df -H . >> cpuinfo
      df -H

      gsutil -q stat '~{fastq_path}/*'
      output_gs_status=$?
      if [[ $output_gs_status == 0 ]]; then
        echo "Has fastq input, good"
      else
        echo "Fastq does not exist"
        exit 1
      fi

      gsutil -q stat '~{output_count_path}/*'
      output_gs_status=$?
      if [[ $output_gs_status == 0 ]]; then
        echo "Already has counts!"
        exit 1
      else
        echo "Counts does not exist, good"
      fi

      justSamplePath='~{fastq_path}/outs/fastq_path/*/~{index_name}*'
      gsutil -q stat $justSamplePath
      output_gs_status=$?
      if [[ $output_gs_status == 0 ]]; then
        echo "Has this sample's fastq input, good"
      else
        echo "Fastq does not exist"
        exit 1
      fi


      gcloud config set storage/process_count 12
      gcloud config set storage/thread_count  2

      mkdir -p ./fastq_out/tmpsubfolder
      #gcloud storage cp -r '~{fastq_path}/*' fastq_out
      # Don't need -r becase just a glob of raw fastq.gz's
      gcloud storage cp $justSamplePath fastq_out/tmpsubfolder |& ts

      df -H

      mkdir ./ref_folder
      # Don't want tested folders, copy directly into ref_folder
      gcloud storage cp -r '~{ref_path__noendslash}/*' ./ref_folder/ |& ts

      mkdir count_out


      # --localcores=~{zz_count_numCores} --localmem=~{zz_count_ram_gb} \

      # actual command :-)
      time stdbuf -oL -eL cellranger count \
          --jobmode=local                  \
          --disable-ui                     \
          --nosecondary                    \
          --id=count_out                   \
          --fastqs=fastq_out/              \
          --sample=~{index_name}           \
          --transcriptome=ref_folder       \
          --include-introns=true |& ts | tee ./count.log

      # no-secondary means no clustering/etc which don't use

      df -H

      rm -rf count_out/SC_RNA_COUNTER_CS

      # Basic check to see that count finished all the way to the end
      if [ -f count_out/outs/metrics_summary.csv ] ; then date > DONE_counts; fi

      # Big upload
      date > d
      gcloud storage cp d ~{output_count_path}/d
      gcloud storage cp -r count_out/* ~{output_count_path} |& ts

      # Make extra sure dstat is dead
      touch usage.csv; pkill dstat; sleep 5s; pkill dstat; pkill dstat;
      gzip -9 usage.csv

      zip -j count_logs.zip ./count.log usage.csv.gz cpuinfo
    >>>

    output {
      File count_logs_zip = "count_logs.zip"
      # Make sure is last so if fails, back up other things
      File DONE_counts = "DONE_counts"
    }
    runtime {
      # Short and fun, tested and almost perfectly fastq size
      # Takes no memory, might as well make preemtible and HDD
      # docker: "us-central1-docker.pkg.dev/velina-208320/jonah-scsnv/img:latest"
      # docker: "scsnv_presquash"
      docker: docker
      # TODO how much memory do need for fastq runs?
      memory: "~{zz_count_ram_gb} GB"
      disks: "local-disk ~{zz_count_disksize} LOCAL"
      cpu: "~{zz_count_numCores}"
      preemptible: 0
    }
}


workflow tags_through_cellranger {
  String pipeline_version = "1.0.0"
  input {
    String bcl_dir__noendslash

    String local_workspace_bucket__noendslash = "gs://fc-7dd29fd0-8983-4611-b3cd-46123534add7"
    String zz_docker = "us-central1-docker.pkg.dev/velina-208320/jonah-scsnv/img:latest"

    String ref_path__noendslash

    Array[String] justRNASampleNames
  }

  parameter_meta {
    # fastqs: "List of fastqs for this library"
    # libraryname: "Library name to call this"
    # tar_reference_without_fasta: "tar of reference (with scsnv index + bwa index)"
  }

  call getBCLSize {
      input:
         bcl_dir__noendslash = bcl_dir__noendslash,
         docker = zz_docker,
  }

  # call mkfastq {
  #   input:
  #     bcl_dir__noendslash = bcl_dir__noendslash,

  #     local_workspace_bucket__noendslash = local_workspace_bucket__noendslash,
  #     output_fastq_path = local_workspace_bucket__noendslash+"/02_OUTFASTQS/"+basename(bcl_dir__noendslash),
  #     # TODO disksize, resize up for actual usages
  #     # zz_bcl_disksize = 180,
  #     # Made sure is decimal not binary GB (not GiB)
  #     # zz_bcl_disksize = 2 * getBCLSize.GiB_size_bcl,
  #     zz_bcl_disksize = 2500,
  #     # 2k ->
  #     # -h=2       2.2 TiB
  #     # -H=SI=1000 2.4 TB
  #     zz_bcl_numCores = 8,
  #     # zz_bcl_ram_gb   = 30,
  #     zz_bcl_ram_gb   = 60,
  #     docker = zz_docker
  # }

  # TODO in mkfastq, return the max size of any of the count chunks, if grouped together for the next step disk size

  # TODO put back
  scatter(thisSamplename in justRNASampleNames){
  # scatter(thisSamplename in mkfastq.justRNASampleNames){
    call runcounts {
      input:
        # fastq_path = mkfastq.output_fastq_path_copy,

        ref_path__noendslash = ref_path__noendslash,

        index_name = thisSamplename,

        output_count_path = local_workspace_bucket__noendslash+"/04_COUNTS/"+basename(bcl_dir__noendslash)+"/"+thisSamplename,

        zz_count_disksize = 180,
        zz_count_numCores = 8,
        zz_count_ram_gb   = 30,

        docker = zz_docker,

        # TODO put back
        fastq_path = local_workspace_bucket__noendslash+"/02_OUTFASTQS/"+basename(bcl_dir__noendslash),
    }
  }

  output {
    # File barcode_counts_totals_gz = scsnv_count.barcode_counts_totals_gz
    # File count_logs_zip           = scsnv_count.count_logs_zip

    # File map_logs_zip   = scsnv_map.map_logs_zip
    # File map_allout_tar = scsnv_map.map_allout_tar

    # File collapsed_outbam   = scsnv_collapse.collapsed_outbam
    # File collapse_logs_zip = scsnv_collapse.collapse_logs_zip

    # File pileup_allout_tar   = scsnv_pileup.pileup_allout_tar
    # File pileup_logs_zip = scsnv_pileup.pileup_logs_zip
  }
}

# Debug with tmux
# apt update && apt install -y tmux htop
# export TERM=linux
# stty rows 32 cols 130
# tmux
# tmux set prefix C-a

# curl -L 'https://github.com/aristocratos/btop/releases/download/v1.2.13/btop-x86_64-linux-musl.tbz' > tmp.tbz
# tar xavf tmp.tbz

# gsutil perfdiag -o out.json 'gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/02_OUTFASTQS/221217_VL00297_113_AACHKNMM5/'

# PATH=$PATH:/software/cellranger-7.1.0/bin/:/usr/local/bcl2fastq/bin/
# /software/cellranger-7.1.0/external/cellranger_tiny_fastq

# RUNNING THEM
# get the bcl flowcell
# gcloud storage cp -r gs://macosko_backup/flowcells/221217_VL00297_113_AACHKNMM5/ gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/01_INBCLS


# date > d
# gsutil cp d gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/01_INBCLS/230125_SL-NSA_0530_BHN2GJDRX2/d
# gsutil ls gs://macosko_data/slidetags/flowcells/230125_SL-NSA_0530_BHN2GJDRX2/ > t
# cat t|grep -v cellranger|xargs -I@ echo 'gcloud storage cp -r @ gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/01_INBCLS/230125_SL-NSA_0530_BHN2GJDRX2/'|parallel -j 5
