#!/bin/bash
LDSERVER_VERSION=`git describe --tags --abbrev=11 | sed 's/^v//' | sed 's/-g/-/'`
GIT_SHA=`git rev-parse HEAD`
BUILD_DATE=`date -u +'%Y-%m-%dT%H:%M:%SZ'`

# Build the base ldserver image.
docker build --pull -t ldserver:base \
  --build-arg MAKEFLAGS="-j 4" \
  --build-arg CMAKE_BUILD_PARALLEL_LEVEL=4 \
  --build-arg BUILD_DATE=${BUILD_DATE} \
  --build-arg GIT_SHA=${GIT_SHA} \
  --build-arg LDSERVER_VERSION=${LDSERVER_VERSION} \
  --target base \
  "$@" .

# Create the final compiled ldserver image.
docker build --pull -t ldserver:${LDSERVER_VERSION} \
  --build-arg MAKEFLAGS="-j 4" \
  --build-arg CMAKE_BUILD_PARALLEL_LEVEL=4 \
  --build-arg BUILD_DATE=${BUILD_DATE} \
  --build-arg GIT_SHA=${GIT_SHA} \
  --build-arg LDSERVER_VERSION=${LDSERVER_VERSION} \
  --target compile \
  "$@" .

# Tag final ldserver image as latest.
docker tag ldserver:${LDSERVER_VERSION} ldserver:latest

# Create development image for CLion.
# The - < is important, it pipes the Dockerfile-CLion contents into
# the docker build engine, which means it will not attempt to copy the local
# context. CLion will rsync source files into the container so we don't need the
# entire context in the image.
# A password should be set for CLion to SSH into the container (see below for example.)
docker build \
  --build-arg LDSERVER_SSH_PASSWORD=ldserver \
  -t ldserver:dev-clion \
  - < Dockerfile-CLion