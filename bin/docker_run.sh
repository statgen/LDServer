#!/bin/bash
LDSERVER_VERSION=`git describe --tags --abbrev=11`
CONFIG_DATA="/home/ldserver/var/test.yaml"
WORKERS=6

# Run the docker image, passing various parameters into the container.
docker run -it \
  -v /mnt/data:/home/ldserver/var \
  -p 5000:5000 \
  -e FLASK_APP="rest/raremetal" \
  ldserver:${LDSERVER_VERSION} \
  /bin/bash -c "flask add-yaml ${CONFIG_DATA} && gunicorn -b 127.0.0.1:5000 -w ${WORKERS} -k gevent --pythonpath rest 'raremetal:create_app()'"
