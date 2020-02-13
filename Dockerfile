FROM ubuntu:18.04

ARG BUILD_DATE
ARG GIT_SHA
ARG LDSERVER_VERSION

LABEL org.label-schema.name="LDServer"
LABEL org.label-schema.version=$LDSERVER_VERSION
LABEL org.label-schema.description="LDServer for calculating linkage disequilibrium of genetic variants"
LABEL org.label-schema.vendor="University of Michigan, Center for Statistical Genetics"
LABEL org.label-schema.url="https://github.com/statgen/LDServer"
LABEL org.label-schema.usage="https://github.com/statgen/LDServer#docker"
LABEL org.label-schema.vcs-url="https://github.com/statgen/LDServer"
LABEL org.label-schema.vcs-ref=$GIT_SHA
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.schema-version="1.0"

# Install required packages for swiss to install. Many of swiss' dependencies
# require compiling C/C++ code.
RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  python \
  python-pip \
  zlib1g-dev \
  liblzma-dev \
  libopenblas-base \
  libopenblas-dev \
  liblapack-dev \
  libarpack2 \
  libarpack2-dev \
  redis \
  locales \
  && rm -rf /var/lib/apt/lists/* \
  && locale-gen en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# Install necessary python packages (backports.lzma needed for cget to extract .xz archives)
RUN pip install backports.lzma cget pytest invoke tox

# Copy the python requirements for install
COPY rest/requirements.txt /
RUN pip install -r requirements.txt

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

# Copy source
COPY --chown=ldserver:ldserver . /home/ldserver/

# This executes a compile, then runs C++ and python tests
# We want the docker build to fail if the tests fail
RUN invoke test
