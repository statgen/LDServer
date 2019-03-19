#!/bin/bash
LDSERVER_VERSION=`git describe --tags --abbrev=11`
RAREMETAL_CONFIG_DATA="/home/ldserver/var/test.yaml"
RAREMETAL_WORKERS=6

# Run the raremetal flask app.
# This is an example of how to run the container for the raremetal flask app directly,
# though you most likely want to use docker-compose in production.
docker run -it \
  -v /mnt/data:/home/ldserver/var \
  -v /mnt/data/config.py:/home/ldserver/rest/instance/config.py \
  -p 5000:5000 \
  -e FLASK_APP="rest/raremetal" \
  ldserver:${LDSERVER_VERSION} \
  /bin/bash -c "flask add-yaml ${RAREMETAL_CONFIG_DATA} && gunicorn -b 0.0.0.0:5000 -w ${RAREMETAL_WORKERS} -k gevent --pythonpath rest 'raremetal:create_app()'"
