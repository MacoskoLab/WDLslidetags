version 1.0

import "https://github.com/MacoskoLab/WDLslidetags/blob/main/tasks.wdl" as tasks

workflow tagspipeline {
  String pipeline_version = "1.0.0"
  input {
    String bcl
    Boolean run1_mkfastq   = true
    Boolean run2_RNAcounts = true
    Boolean run2_SBcounts  = true
    Boolean run3_spatial   = true
    String zz_bucket = "fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d"
    String zz_docker = "us-central1-docker.pkg.dev/velina-208320/jonah-slidetags/img:latest"
  }

  parameter_meta {
  }

  call tasks.read_sheet as read_sheet {
    input:
      bcl = bcl,
      run_mkfastq   = run1_mkfastq,
      run_RNAcounts = run2_RNAcounts,
      run_SBcounts  = run2_SBcounts,
      run_spatial   = run3_spatial,
      
      bucket = zz_bucket,
      docker = zz_docker
  }

  if (run1_mkfastq) {
    call tasks.mkfastq as mkfastq {
      input:
        bcl = bcl,
        Indexes = read_sheet.Indexes,
        
        bucket = zz_bucket,
        docker = zz_docker,
        MKFASTQSIZE = read_sheet.MKFASTQSIZE
    }
  }
  Boolean DONEmkfastq = select_first([mkfastq.DONEmkfastq, true])

  call tasks.compute_sizes as RNAsizes {
    input:
      array2D = read_sheet.RNAcounts,
      bcl = bcl,
      bucket = zz_bucket,
      docker = zz_docker,
      DONEmkfastq = DONEmkfastq
  }

  output {
  }
}
