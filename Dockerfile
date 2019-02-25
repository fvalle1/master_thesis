FROM jupyter/scipy-notebook

USER root


USER $NB_UID
RUN mkdir /home/jovyan/thesis
COPY *.ipynb /home/jovyan/thesis/

WORKDIR /home/jovyan/
COPY spec-file.txt /home/jovyan/
RUN conda install
