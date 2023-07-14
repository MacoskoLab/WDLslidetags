FROM slidetags_presquash as initial
FROM scratch

SHELL ["/bin/bash", "-c"]


COPY --from=initial / /

# RUN apt update && apt install -y valgrind


# In this directory:
# docker build -t us-central1-docker.pkg.dev/velina-208320/jonah-slidetags/img .
# docker push     us-central1-docker.pkg.dev/velina-208320/jonah-slidetags/img
