FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y python3-pip git ffmpeg awscli && rm -rf /var/lib/apt/lists/*
WORKDIR /workspace
COPY requirements.txt ./
RUN pip3 install -r requirements.txt
COPY . /workspace/
ENV PYTHONPATH=/workspace
ENTRYPOINT ["bash"]
