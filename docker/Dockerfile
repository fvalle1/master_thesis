FROM continuumio/anaconda3:2019.07

MAiNTAINER Filippo Valle filippo.valle@edu.unito.it

USER root

RUN apt-get update --yes && apt-get upgrade --yes
RUN apt-get install sudo wget apt-utils --yes

ENV PATH="/opt/conda/bin:$PATH"

RUN useradd -m -g users filippo

RUN conda update conda

COPY environment.yml /home/filippo/

WORKDIR /home/filippo
RUN conda env create --name phd --file environment.yml

RUN mkdir /home/filippo/files
RUN chmod a+x /home/filippo/files
COPY pipeline.py /home/filippo/files/
COPY Pipeline.ipynb /home/filippo/files/
COPY hsbmpy.py /home/filippo/files/
COPY geneontology.py /home/filippo/files/

USER filippo
RUN /bin/bash -c "conda init bash"
RUN echo "source activate phd" >> $HOME/.bashrc
RUN /bin/bash -c "source $HOME/.bashrc"

WORKDIR /home/filippo/files

ENTRYPOINT /bin/bash -c "python pipeline.py"
CMD [""]
