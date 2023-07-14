# You need 2 files to build this docker

- gs://macosko_data/WDLslidetags_Input/upload_for_google_key.json

Gives the service key to modify the google sheet which hods sample metadata

- gs://macosko_data/WDLslidetags_Input/bcl2fastq.zip

bcl2fastq from Illumina. Can't distribute publically bc of license

# How to run

I did a hack to compress the docker file. First I build `slidetags_presquash`
which has all the docker layers. Then I build in `01_for_squash` which imports
it all at once. And push to `gcr` (`docker push     us-central1-docker.pkg.dev/velina-208320/jonah-slidetags/img`)

