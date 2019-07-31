FROM tiagopeixoto/graph-tool

MAINTAINER Filippo Valle <filippo.valle@edu.unito.it>

RUN cd /home/user/
COPY hSBM.ipynb /home/user/
COPY sbmtm.py /home/user/
ADD run.py /home/user/

WORKDIR /home/user/

ENTRYPOINT ["/home/user/run.py"]
CMD ["python", "./run.py"]
