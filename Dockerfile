FROM debian:stable-slim
LABEL maintainer <torrisimirko@yahoo.com>

# satisfy the requirements
RUN apt-get update && apt-get upgrade -y
RUN apt-get install git python3 python3-numpy hhsuite ncbi-blast+ -y
RUN apt-get autoremove -y && rm -rf /var/lib/apt/lists/*

# get Porter5
RUN git clone https://github.com/mircare/Porter5/
RUN git clone http://github.com/soedinglab/hh-suite

ENV HHLIB=/hh-suite
ENV PATH="$HHLIB/bin:$HHLIB/scripts:${PATH}"

# initialize Porter5
RUN echo "[DEFAULT]" >> Porter5/scripts/config.ini
RUN echo "psiblast = psiblast" >> Porter5/scripts/config.ini
RUN echo "uniref90 = /uniref90/uniref90_2018-03.fasta" >> Porter5/scripts/config.ini
RUN echo "hhblits = hhblits" >> Porter5/scripts/config.ini
RUN echo "uniprot20 = /uniprot20/uniprot20_2016_02" >> Porter5/scripts/config.ini

WORKDIR /Porter5/scripts/Predict_BRNN
RUN bash set_models.sh && chmod a+x Predict

WORKDIR /Porter5