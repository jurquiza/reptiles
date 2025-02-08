From python:3.10-slim
RUN apt update -y
RUN apt install build-essential -y
RUN pip install --upgrade pip
COPY requirements_reptiles.txt /
RUN pip install -r /requirements_reptiles.txt

#ADD . /simple_worker
WORKDIR /simple_worker

ENTRYPOINT celery -A tasks worker --loglevel=info