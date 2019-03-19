#!/bin/bash

docker system prune -f
docker image rm -f `docker image ls --filter "label=org.label-schema.name=LDServer" -q | sort -u`
