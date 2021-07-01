FROM ubuntu:20.04 as base

LABEL org.label-schema.name="LDServer"
LABEL org.label-schema.description="LDServer for calculating linkage disequilibrium of genetic variants"
LABEL org.label-schema.vendor="University of Michigan, Center for Statistical Genetics"
LABEL org.label-schema.url="https://github.com/statgen/LDServer"
LABEL org.label-schema.usage="https://github.com/statgen/LDServer#docker"
LABEL org.label-schema.vcs-url="https://github.com/statgen/LDServer"
LABEL org.label-schema.schema-version="1.0"

# Install required packages for LDServer to install.
ENV DEBIAN_FRONTEND="noninteractive"
RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential \
  curl \
  cmake \
  python3 \
  python3-dev \
  python3-distutils \
  python3-setuptools \
  python3-pip \
  python3-wheel \
  zlib1g-dev \
  liblzma-dev \
  libopenblas-base \
  libopenblas-dev \
  liblapack-dev \
  libarpack2 \
  libarpack2-dev \
  redis \
  locales \
  git \
  pkg-config \
  && rm -rf /var/lib/apt/lists/* \
  && locale-gen en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# Upgrade pip
RUN pip3 install --upgrade pip

# Install required python packages for building later packages
COPY rest/build.txt /
RUN pip3 install -r build.txt

# Install required python packages
COPY rest/requirements.txt /
RUN pip3 install -r requirements.txt

# Create a group and user to execute as, then drop root
ARG UID
ARG GID
RUN \
  if [ -n "$GID" ]; then \
    addgroup --gid $GID ldserver; \
  else \
    addgroup ldserver; \
  fi && \
  if [ -n "$UID" ]; then \
    adduser --gecos "User for running LDServer as non-root" --shell /bin/bash --disabled-password --uid $UID --ingroup ldserver ldserver; \
  else \
    adduser --gecos "User for running LDServer as non-root" --shell /bin/bash --disabled-password --ingroup ldserver ldserver; \
  fi

WORKDIR /home/ldserver
USER ldserver

# Install cpp dependencies
COPY --chown=ldserver:ldserver core/requirements.txt /home/ldserver/core/requirements.txt
COPY --chown=ldserver:ldserver core/*.cmake /home/ldserver/core/
ARG CMAKE_BUILD_PARALLEL_LEVEL
ARG MAKEFLAGS
RUN cget install -f core/requirements.txt

# Next stage: compiled server/binaries
FROM base as compile

# Copy source
COPY --chown=ldserver:ldserver . /home/ldserver/

# Compile ldserver cpp
ENV CGET_PREFIX="/home/ldserver/cget"
ENV INSTALL_PREFIX="/home/ldserver/cget"
RUN \
  mkdir build \
  && cd build \
  && cmake .. \
    -DCMAKE_TOOLCHAIN_FILE=${CGET_PREFIX}/cget/cget.cmake \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
    -DCMAKE_BUILD_TYPE=Release \
  && cmake --build . --target install

# Run test cases
RUN tox

# Frequently changing metadata here to avoid cache misses
ARG BUILD_DATE
ARG GIT_SHA
ARG LDSERVER_VERSION

LABEL org.label-schema.version=$LDSERVER_VERSION
LABEL org.label-schema.vcs-ref=$GIT_SHA
LABEL org.label-schema.build-date=$BUILD_DATE

# Set the default stage to be the base files + compiled binaries + test cases.
FROM compile