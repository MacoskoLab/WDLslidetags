version 1.0

import "https://raw.githubusercontent.com/MacoskoLab/WDLslidetags/main/tasks.wdl" as tasks

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
        Indexes = read_sheet.Indexes,
        
        bcl = bcl,
        bucket = zz_bucket,
        docker = zz_docker,
        SIZE = read_sheet.MKFASTQSIZE
    }
  }
  Boolean DONEmkfastq = select_first([mkfastq.DONE, true])

  if (run2_RNAcounts) {
    call tasks.compute_sizes as compute_RNAsizes {
      input:
        array2D = read_sheet.RNAcounts,
        DONEmkfastq = DONEmkfastq,

        bcl = bcl,
        bucket = zz_bucket,
        docker = zz_docker
    }
  }
  Array[Int] RNAsizes_arr = select_first([compute_RNAsizes.SIZES, []])
  Array[Array[String]] RNAcounts_arr = select_first([compute_RNAsizes.COUNTS, []])

  scatter(pair in zip(RNAcounts_arr, RNAsizes_arr)) {
    Array[String] RNAparams = pair.left
    Int RNASIZE = pair.right
    if (length(RNAparams) == 2) {
      call tasks.RNAcounts as RNAcounts {
        input:
          index = RNAparams[0],
          transcriptome = RNAparams[1],

          bcl = bcl,
          bucket = zz_bucket,
          docker = zz_docker,
          SIZE = RNASIZE
      }
    }
    if (length(RNAparams) == 3) {
      call tasks.RNAcounts_FFPE as RNAcounts_FFPE {
        input:
          index = RNAparams[0],
          transcriptome = RNAparams[1],
          probeset = RNAparams[2],

          bcl = bcl,
          bucket = zz_bucket,
          docker = zz_docker,
          SIZE = RNASIZE
      }
    }
  }
  Array[Boolean] DONERNAcounts = flatten([select_all(RNAcounts.DONE), select_all(RNAcounts_FFPE.DONE)])


  output {
  }
}
