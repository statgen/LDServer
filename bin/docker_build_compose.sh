#!/bin/bash
export LDSERVER_VERSION=`git describe --tags --abbrev=11 | sed 's/^v//' | sed 's/-g/-/'`
export GIT_SHA=`git rev-parse HEAD`
export BUILD_DATE=`date -u +'%Y-%m-%dT%H:%M:%SZ'`

docker-compose build --pull "$@"
