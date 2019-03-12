FROM nvidia/cuda:10.1-base

RUN apt-get update && apt-get install -y python3 python3-pip

RUN pip3 install numpy pandas sklearn keras tensorflow matplotlib pillow argparse

COPY src /app

WORKDIR /app

ENTRYPOINT ["python3", "exp_gan.py"]
