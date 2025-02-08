FROM python:3.10-slim
RUN apt update -y
RUN apt install build-essential -y
RUN pip install --upgrade pip
COPY requirements_reptiles.txt /
RUN pip install -r /requirements_reptiles.txt

#ADD . /flask_app
WORKDIR /flask_app

ENV FLASK_ENV=development
ENV FLASK_DEBUG=1

CMD flask run --host=0.0.0.0