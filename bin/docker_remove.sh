#!/bin/bash
LDSERVER_VERSION=`git describe --tags --abbrev=11`
docker system prune -f
docker image rm ldserver:${LDSERVER_VERSION}
