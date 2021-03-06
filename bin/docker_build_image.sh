#!/bin/bash
LDSERVER_VERSION=`git describe --tags --abbrev=11 | sed 's/^v//' | sed 's/-g/-/'`
GIT_SHA=`git rev-parse HEAD`
BUILD_DATE=`date -u +'%Y-%m-%dT%H:%M:%SZ'`

docker build --pull -t ldserver:${LDSERVER_VERSION} --build-arg BUILD_DATE=${BUILD_DATE} --build-arg GIT_SHA=${GIT_SHA} --build-arg LDSERVER_VERSION=${LDSERVER_VERSION} "$@" .
