#!/bin/bash
LDSERVER_VERSION=`git describe --tags --abbrev=11`
GIT_SHA=`git rev-parse HEAD`
BUILD_DATE=`date -u +'%Y-%m-%dT%H:%M:%SZ'`

if [ `git rev-parse --abbrev-ref HEAD` != 'master' ]; then
  echo "Must be on master branch to build docker image"
  exit 1
fi

docker build -t ldserver:${LDSERVER_VERSION} --build-arg BUILD_DATE=${BUILD_DATE} --build-arg GIT_SHA=${GIT_SHA} --build-arg LDSERVER_VERSION=${LDSERVER_VERSION} . && \
  docker tag ldserver:${LDSERVER_VERSION} ldserver:latest
