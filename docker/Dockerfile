FROM continuumio/anaconda3:2019.07

LABEL maintainer='Filippo Valle filippo.valle@edu.unito.it'

USER root

RUN apt-get update --yes && apt-get upgrade --yes
RUN apt-get install sudo wget apt-utils --yes

ENV PATH="/opt/conda/bin:$PATH"

RUN useradd -m -g users filippo

RUN conda update conda --yes

RUN mkdir /home/filippo/tacos
RUN mkdir /home/filippo/tacos/files

COPY environment.yml /home/filippo/tacos

WORKDIR /home/filippo/tacos
RUN conda env create --name phd --file environment.yml

WORKDIR /home/filippo/tacos

RUN apt-get update && apt-get install git cmake libboost-all-dev g++ --yes
RUN git clone https://github.com/fvalle1/tacos.git
WORKDIR /home/filippo/tacos/tacos/
RUN mkdir build
WORKDIR /home/filippo/tacos/tacos/build
RUN cmake ..
RUN make data_mining
RUN ./data_mining
RUN mv /home/filippo/tacos/tacos/build/data_mining /home/filippo/tacos/.

WORKDIR /home/filippo/tacos

RUN chmod a+rwx -R /home/filippo/
COPY pipeline.py /home/filippo/tacos/
COPY addo.py /home/filippo/tacos/
COPY Pipeline.ipynb /home/filippo/tacos/
COPY hsbmpy.py /home/filippo/tacos/
COPY tableanalyser.py /home/filippo/tacos/
COPY tacos_plot.py /home/filippo/tacos/
COPY geneontology.py /home/filippo/tacos/
COPY gene_symbol.txt /home/filippo/tacos/
COPY genes.txt /home/filippo/tacos/
RUN chmod a+rwx -R /home/filippo/tacos/*


USER filippo
RUN /bin/bash -c "conda init bash"
RUN echo "source activate phd" >> $HOME/.bashrc
RUN /bin/bash -c "source $HOME/.bashrc"

ENTRYPOINT /bin/bash -c "/opt/conda/envs/phd/bin/python pipeline.py"
CMD [""]
