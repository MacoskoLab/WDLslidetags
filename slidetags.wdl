version 1.0

import "https://github.com/MacoskoLab/WDLslidetags/blob/main/tasks.wdl" as tasks

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

  call tasks.read_sheet {
    input:
      bcl = bcl,
      run_mkfastq = run_mkfastq,
      run_counts  = run_counts,
      run_spatial = run_spatial,
      
      bucket = zz_bucket,
      docker = zz_docker
  }

  call tasks.mkfastq {
    input:
      bcl = bcl,
      Indexes = read_sheet.Indexes,
      run_mkfastq = run_mkfastq,
      
      bucket = zz_bucket,
      docker = zz_docker,
      MKFASTQSIZE = read_sheet.MKFASTQSIZE
  }

  scatter(rowindex in read_sheet.COUNTSROWS) {
    call tasks.counts {
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
    call tasks.spatial {
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
