#!/bin/bash

# This script should be sourced, otherwise docker-compose will continuously complain about
# missing environment variables.
export FLASK_APP="rest/raremetal"
export CONFIG_DATA="var/test.yaml"
export WORKERS=4

docker-compose pull --ignore-pull-failures && docker-compose up -d
