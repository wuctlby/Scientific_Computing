FROM almalinux:9

ARG USERNAME=wuct_docker
ARG USER_UID=1000
ARG USER_GID=1000

RUN dnf update -y && \
    dnf install -y \
    sudo \
    git \
    wget \
    tar \
    gzip \
    python3 \
    python3-pip \
    gcc \
    gcc-c++ \
    make \
    cmake \
    which \
    gsl gsl-devel fftw fftw-devel \
    && dnf clean all \
    && rm -rf /var/cache/dnf

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

WORKDIR /workspace
USER $USERNAME
ENV LANG=C.UTF-8

ENV PATH="/home/$USERNAME/.local/bin:$PATH"
RUN pip3 install --user --upgrade pip
RUN pip3 install --user numpy scipy pytest