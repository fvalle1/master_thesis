FROM jupyter/scipy-notebook

USER root

RUN sudo useradd -g wheel filippo

USER filippo
RUN mkdir /home/filippo/thesis
COPY *.ipynb /home/filippo/thesis/

WORKDIR /home/filippo/
COPY spec-file.txt /home/filippo/

ENTRYPOINT "jupyter notebook"
CMD [""]
