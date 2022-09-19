FROM ubuntu:latest
WORKDIR /usr/src/drugstone/

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update --no-install-recommends && apt-get upgrade -y && apt-get install -y supervisor nginx libgtk-3-dev wget
RUN apt-get autoclean -y && apt-get autoremove -y


ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh -O ~/miniconda.sh && /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN conda init bash

RUN pip install gunicorn

COPY ./requirements.txt /usr/src/drugstone/requirements.txt
RUN pip install -r /usr/src/drugstone/requirements.txt

COPY ./supervisord.conf /etc/supervisor/conf.d/supervisord.conf
RUN pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple nedrex==0.1.4

COPY . /usr/src/drugstone/
