# LDServer

[![Build Status](https://travis-ci.com/statgen/LDServer.svg?branch=master)](https://travis-ci.com/statgen/LDServer) [![Docker Pulls](https://img.shields.io/docker/pulls/statgen/ldserver.svg)](https://cloud.docker.com/u/statgen/repository/docker/statgen/ldserver)

LDServer is a fast implementation of various metrics of linkage disequilibrium (LD) between genetic variants.

Features:

* **Fast**: calculating LD in a typical window (500kb) and on typical VCFs (10,000's of samples) will return in under 3-4 seconds.
* **API**: serves statistics over a REST API, suitable for web applications (used by [LocusZoom](http://locuszoom.org)).
* **Caching**: statistics are cached by Redis; common queries will return in under 1 second.
* **Paging**: iterate over pages of results, rather than processing immediately all at once.

This project contains multiple components that work together to provide these features:

* LDServer C++ shared library, wrapped by boost python.
* Flask apps for calculating and serving the statistics in various ways:
  * A standard `rest` app for serving r2, D', and covariance between genetic markers. This app is typically used by web applications that display LD between variants.
  * A `playground` app to help developers get started with the various APIs.
  * A `raremetal` app for serving covariance in addition to score statistics for calculating aggregation tests of rare genetic variants.

## Documentation

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [LDServer](#ldserver)
  - [Documentation](#documentation)
  - [Installation](#installation)
    - [Docker](#docker)
      - [Getting started](#getting-started)
      - [.env](#env)
      - [docker-compose config](#docker-compose-config)
      - [Redis config](#redis-config)
      - [Flask config](#flask-config)
      - [Configuration for ldserver app](#configuration-for-ldserver-app)
      - [Configuration for raremetal app](#configuration-for-raremetal-app)
      - [Running the services](#running-the-services)
    - [Manual installation](#manual-installation)
  - [Updating](#updating)
  - [Configuring & running the flask apps](#configuring-running-the-flask-apps)
    - [ldserver app](#ldserver-app)
      - [Configuring the ldserver app](#configuring-the-ldserver-app)
        - [Add new reference](#add-new-reference)
        - [Add new (sub-)population to the reference](#add-new-sub-population-to-the-reference)
        - [List loaded references](#list-loaded-references)
        - [List loaded genotype files](#list-loaded-genotype-files)
        - [List loaded population samples](#list-loaded-population-samples)
      - [Running the rest app](#running-the-rest-app)
    - [playground app](#playground-app)
      - [Running the playground app](#running-the-playground-app)
    - [raremetal app](#raremetal-app)
      - [Required data and formats](#required-data-and-formats)
        - [Genotype data](#genotype-data)
        - [Phenotype data](#phenotype-data)
        - [Masks](#masks)
        - [Summary statistics](#summary-statistics)
      - [Configuring the raremetal app](#configuring-the-raremetal-app)
        - [YAML: Adding multiple datasets with one YAML config file](#yaml-adding-multiple-datasets-with-one-yaml-config-file)
          - [Genotype records](#genotype-records)
          - [Phenotype records](#phenotype-records)
          - [Mask records](#mask-records)
          - [Summary statistic records](#summary-statistic-records)
          - [Add all records in YAML](#add-all-records-in-yaml)
        - [CLI: Quickly add data for testing purposes](#cli-quickly-add-data-for-testing-purposes)
          - [Adding genotype datasets](#adding-genotype-datasets)
          - [Adding phenotype datasets](#adding-phenotype-datasets)
          - [Adding masks of genetic variants](#adding-masks-of-genetic-variants)
          - [Listing available datasets](#listing-available-datasets)
      - [Running the raremetal app](#running-the-raremetal-app)
      - [Exposed API endpoints](#exposed-api-endpoints)

<!-- /code_chunk_output -->

## Installation

### Docker

#### Getting started

To run LDServer in production, we recommend using [docker-compose](https://docs.docker.com/compose/install/). We provide a base configuration `docker-compose.yml` that specifies the core services.

To customize the config to your environment, either create a [docker-compose.override.yml](https://docs.docker.com/compose/extends/#example-use-case) file, or copy the base `docker-compose.yml` to a new file to edit, for example `docker-compose.prod.yml`.

Using an override file is convenient because you can continue using `docker-compose <command>` as you normally would. Using a new compose file, such as `docker-compose.prod.yml` requires an additional flag upon each invocation: `docker-compose -f docker-compose.prod.yml <command>`.

To begin, first checkout a copy of the repository with git clone:

```bash
git clone https://github.com/statgen/LDServer.git
cd LDServer
```

#### .env

Create a `.env` file to specify docker configuration settings:

```bash
LDSERVER_PORT=4546
LDSERVER_CONFIG_SCRIPT=/home/ldserver/startup.sh
LDSERVER_WORKERS=4
RAREMETAL_CONFIG_DATA=var/config.yaml
RAREMETAL_WORKERS=4
RAREMETAL_PORT=4545
```

* `LDSERVER_PORT` is the exposed port that the ldserver app will run on. You will likely want to put apache or nginx
  in front of the server as an HTTP proxy.
* `LDSERVER_CONFIG_SCRIPT` is the path *inside the container* to the startup script. Usually you will create this
  script on your local filesystem and mount it into the container with a volume directive (see the example
  `docker-compose.override.yml` file below.)
* `LDSERVER_WORKERS` is the number of asynchronous workers to be started to serve ldserver requests.
* `RAREMETAL_CONFIG_DATA` is the path *inside the container* to the raremetal app's data config yaml, which
  specifies all of the datasets to load. Similarly to the `LDSERVER_CONFIG_SCRIPT`, this should be mounted into the
  container in the `docker-compose.override.yml` file.
* `RAREMETAL_WORKERS` see `LDSERVER_WORKERS`
* `RAREMETAL_PORT` see `LDSERVER_PORT`


#### docker-compose config

Create a `docker-compose.override.yml` file where you can provide your own additional docker configuration. For
example:

```YAML
version: '3'
services:
  ldserver:
    build:
      args:
        UID: 1000
        GID: 1001
    volumes:
      - /opt/ld/ld_panels:/data/ld_panels
      - /opt/ld/logs:/data/logs
      - /opt/ld/config.py:/home/ldserver/rest/instance/config.py
      - /opt/ld/startup.sh:/home/ldserver/startup.sh
    restart: "unless-stopped"
    command: >
      /bin/bash -c "
      source $$LDSERVER_CONFIG_SCRIPT &&
      gunicorn -b 0.0.0.0 -w $$LDSERVER_WORKERS -k gevent \
        --access-logfile /data/logs/gunicorn.access.log \
        --error-logfile /data/logs/gunicorn.error.log \
        --pythonpath rest 'ldserver:create_app()'"

  raremetal:
    build:
      args:
        UID: 1000
        GID: 1001
    volumes:
      - /opt/raremetal:/home/ldserver/var
      - /opt/raremetal/config.py:/home/ldserver/rest/instance/config.py
    restart: "unless-stopped"

  redis:
    volumes:
      - /opt/ld/redis_cache:/data/redis/redis_cache
      - /opt/ld/logs:/data/redis/logs
      - /opt/ld/redis.conf:/usr/local/etc/redis/redis.conf
    command: "/usr/local/etc/redis/redis.conf"
    restart: "unless-stopped"
```

In the above example, we overrode a few pieces of the various services.

For the `ldserver` service:
  * We overrode the command that is executed when the container first runs to pass additional arguments to gunicorn, which starts the `ldserver` flask app. In this case we've added new logfile arguments `--access-logfile` and `--error-logfile`.
  * We specified a number of directories and files to mount within the container using `volumes`. The syntax is
  `HOST_PATH:CONTAINER_PATH`. For example, the file `/opt/ld/startup.sh` on the host machine will appear as
  `/home/ldserver/startup.sh` within the ldserver container when it runs.

For the `raremetal` service:
  * We're also mounting a directory `/opt/raremetal` and config file `/opt/raremetal/config.py` from the host when 
    the container runs. Volumes have the format `/path/to/data`:`/path/to/mount/in/container`.

For both services:
  * We've turned on `restart: "unless-stopped"`. This means the services will be restarted if they fail, or if the docker daemon is restarted. They will not restart if they are manually stopped.
  * We decided to override the UID and GID of the user inside the docker image. It is important to specify UID/GID
    for both the ldserver and raremetal images to avoid cache misses during the docker build process. **WARNING**:
    do not specify a very large UID or GID, or the build process will take a long time due to a [known bug in docker](https://github.com/moby/moby/issues/5419).

#### Redis config

You must supply a slightly modified `redis.conf` for running under docker. You can grab a copy of the default `redis.conf` to modify by going to https://redis.io/topics/config (we are currently using [redis 5.0](https://raw.githubusercontent.com/redis/redis/5.0/redis.conf)).

It is a good idea to customize redis to control memory usage, and path to where the periodic database saves are done. The `docker-compose.override.yml` file above shows to how to mount your own `redis.conf` into the redis service container when it runs.

The following settings should be set in your `redis.conf`:

```
bind 0.0.0.0
port 6379
logfile "/data/redis/logs/redis.log"
save 14400 1
dbfilename ld.rdb
dir /data/redis/redis_cache/
maxmemory 4g
maxmemory-policy allkeys-lru
```

`bind 0.0.0.0` is needed to allow redis to connect to other docker containers (such as the ldserver or raremetal services).

Change `maxmemory` to something sensible considering the amount of available RAM on your server.

#### Flask config

There is a default config file that can be used as a starting place in `rest/config/default.py`. It looks like:

```
SQLALCHEMY_DATABASE_URI = 'sqlite:///sql.db'
SQLALCHEMY_TRACK_MODIFICATIONS = False
PROXY_PASS = None # set if Apache mod_proxy is used e.g. http://my.host.com/prefix/
API_MAX_PAGE_SIZE = 100000
API_MAX_REGION_SIZE = 4000000
API_MAX_COV_REGION_SIZE = 1000000
SEGMENT_SIZE_BP = 1000
CACHE_ENABLED = True
CACHE_REDIS_HOSTNAME = '127.0.0.1'
CACHE_REDIS_PORT = 6379
GZIP_COMPRESSION = True # enable build-in response compression if for any reason it was not possible to enable it through Apache
API_BASE_URL = 'http://127.0.0.1:5000' # specify the base URL address of the LD server. It is used by the Playground.
SENTRY_DSN = None # include your Sentry DSN for error reporting to your own Sentry instance
SENTRY_ENV = None # name of your deployment - usually 'production' or 'staging' or 'travis'
```

You can copy this file and customize it. It should be placed in a location you have specified in your `docker-compose.override.yml` (see above), for example: `/opt/ld/config.py` or `/opt/raremetal/config.py`.

:warning: For docker, **you must change** `CACHE_REDIS_HOSTNAME` to `redis`.

If you have an Apache proxy forwarding requests, you should specify `PROXY_PASS`. For example, let's say your server has Apache setup to do:  

* https://myserver.com/ld/ -> http://localhost:4546
* https://myserver.com/raremetal/ -> http://localhost:4545

Then, in your ldserver's config.py, specify `PROXYPASS = 'https://myserver.com/ld/'`. And so forth for the raremetal server as well.

#### Configuration for ldserver app

Mount a configuration script into the container by changing:

```YAML
volumes:
  - XXX:/home/ldserver/startup.sh
```

The `XXX` is a path on your local filesystem to a script that contains commands for adding datasets to the ldserver. See the [Configuring the ldserver app](#configuring-the-ldserver-app) section for commands to use for adding reference genotype files. If desired, you can modify the name of this script, but you must also modify the `LDSERVER_CONFIG_SCRIPT` variable to refer to the new path (as it would appear inside the container.)

You will also need to mount your data into the container. As an example, the line `/opt/ld/ld_panels:/data/ld_panels
` mounts a directory of data `/opt/ld/ld_panels` on your local filesystem into the container at `/data/ld_panels`.

An example `startup.sh` looks like:

```bash
#!/bin/bash

# If the database already exists, then we don't want to run the commands below or they will fail.
# This can happen when the container is being restarted, instead of a clean down/up. The volume
# persists, and therefore so does the database. Remember to always do `docker-compose down && docker-compose up -d`.
if [[ ! -f "rest/ldserver/sql.db" ]]; then
  echo "Running startup flask add commands..."

  flask add-reference 1000G "1000 Genomes Project Phase 3" GRCh37 /data/ld_panels/sav/1000G_GRCh37/ALL.samples.txt /data/ld_panels/sav/1000G_GRCh37/1000G_phase3_GRCh37_ALL_chr*sav
  flask add-reference 1000G "1000 Genomes Project Phase 3" GRCh38 /data/ld_panels/sav/1000G_GRCh38/ALL.samples.txt /data/ld_panels/sav/1000G_GRCh38/1000G_phase3_GRCh38_ALL_chr*sav
  flask create-population GRCh37 1000G AFR /data/ld_panels/sav/1000G_GRCh37/AFR.samples.txt
  flask create-population GRCh37 1000G EUR /data/ld_panels/sav/1000G_GRCh37/EUR.samples.txt
else
  echo "Database file already existed (docker container restart?), skipping flask add commands..."
fi
```

#### Configuration for raremetal app

The raremetal app uses a YAML file for configuration. See the raremetal app section "[YAML: Adding multiple datasets with one YAML config file](#yaml-adding-multiple-datasets-with-one-yaml-config-file)" for more information. Map one into the container by modifying:

```YAML
volumes:
  - XXX:/home/ldserver/var
```

Remember that when specifying paths in the YAML file, you must specify paths **as they would appear in the container.** In the example above, we have mapped in our data files into the directory `/home/ldserver/var`, and therefore the paths in `var/test.yaml` should be `var/data_file.txt`, etc. For example:

```YAML
genotypes:
- id: 1
  name: "1000G"
  description: "1000G chr22 Testing VCF"
  filepath: "var/test.vcf.gz"
  genome_build: "GRCh37"
```

#### Running the services

Assuming you have created your own `docker-compose.override.yml` file for production, first build the docker images using:

```bash
bin/docker_build_compose.sh
```

You don't have to use the above script, and could instead directly run `docker-compose build`, but it does set a few extra labels on the resulting docker image that can be useful.

You can also customize your build with your own script, which can be useful for setting various settings, such as parallel build CPUs:

```bash
#!/bin/bash
LDSERVER_VERSION=`git describe --tags --abbrev=11 | sed 's/^v//' | sed 's/-g/-/'`
GIT_SHA=`git rev-parse HEAD`
BUILD_DATE=`date -u +'%Y-%m-%dT%H:%M:%SZ'`

docker-compose build --pull \
  --build-arg MAKEFLAGS="-j 9" \
  --build-arg CMAKE_BUILD_PARALLEL_LEVEL=9 \
  --build-arg BUILD_DATE=${BUILD_DATE} \
  --build-arg GIT_SHA=${GIT_SHA} \
  --build-arg LDSERVER_VERSION=${LDSERVER_VERSION} "$@"
```

In the example above, we added some extra args (`MAKEFLAGS` and `CMAKE_BUILD_PARALLEL_LEVEL`) to specify 9 CPUs should be used in parallel when compiling the C++ components.

:warning: **Stick with your build script consistently.** Removing build args like the ones above can cause a docker build cache miss, and require you to recompile far more than you would like.

Now, you can run all services once the build is finished using:

```bash
# -d runs services detached in background
docker-compose up -d
```

If you instead chose to use a copied production yaml file, you would do the following. Remember to include this `-f` flag in all docker-compose commands. We will not include it from here onwards.

```bash
docker-compose -f docker-compose.prod.yml up -d
```

To debug the container, the following command is very useful. It will start the raremetal service, using all of your configuration values from the docker-compose file, and drop you into a bash shell as root. From there, you can run the services manually to debug, and install any additional packages that could be useful. Obviously, you should not use this for production, as the container is not meant to be run as root.

```bash
docker-compose run -u root ldserver bash
```

If you've already started services, replace `run` with `exec` to jump into the already running container.

### Manual installation

Before installing, the following software/packages are required:

* CMake 3.14+
* gcc 5.4+
* Python 3.6+
* A BLAS implementation, such as OpenBLAS or Intel MKL

The following software is optional:

* Redis

Most Linux distributions provide these packages already. For example, in Ubuntu:

```bash
sudo apt install \
  build-essential \
  cmake \
  python \
  virtualenv \
  libopenblas-base \
  libopenblas-dev \
  liblapack-dev \
  libarpack2 \
  libarpack2-dev \
  redis
```

Now follow these steps to complete installation:

- Clone repository.

  ```bash
  git clone https://github.com/statgen/LDServer
  cd LDServer
  ```

- Compile and install.

  ```bash
  virtualenv -p python3 venv &&
    source venv/bin/activate &&
    pip install cget &&
    cget install core
  ```

- *Optional*: Test the LDServer C++ shared library. The '127.0.0.1' and '8888' arguments are the hostname and port for Redis server which will be started during unit tests. If you wish to run Redis server on another port, then change these values correspondingly.

  ```bash
  cd cget/test/ && ./testAll 127.0.0.1 8888 && cd ../../
  ```

- Install required python packages for the flask apps.

  ```bash
  # Activate the environment
  source venv/bin/activate

  # Install required packages
  pip install -r rest/requirements.txt
  ```

- *Optional*: Run the flask app test suite.

  ```bash
  cd rest
  python -m pytest
  ```

- Configure REST API by copying the config file `rest/config/default.py` to `rest/instance/config.py` and then modify values as needed. Modifying `default.py` is not recommended as it will be overwritten when pulling updates from github.

## Updating

Update your files using `git pull`, then recompile the core server by doing `cget install --update core`.

## Configuring & running the flask apps

### ldserver app

This app serves r2, D', and covariance of genetic variants over a REST API.

#### Configuring the ldserver app

##### Add new reference

- Execute `flask add-reference` command:
   ```bash
   cd rest
   export FLASK_APP=ldserver
   flask add-reference <name> <description> <genome build> <samples file> <genotype files>
   ```
   For example, the below command adds a new reference panel named `1000G_GRCh37`. The `samples.txt` file stores list of sample names in `ALL.chr*.bcf` files.
   ```bash
   flask add-reference 1000G "1000 Genomes Project Phase 3" GRCh37 samples.txt ALL.chr*.bcf
   ```
   The genotypes can be stored in VCF, BCF, and SAV formats. For better runtime performance and more compact storage, we highly recommend using SAV format.

##### Add new (sub-)population to the reference

- To define populations in the reference, execute `flask create-population` command:

   ```bash
   cd rest
   export FLASK_APP=ldserver
   flask create-population <genome build> <reference name> <population name> <samples file>
   ```

   For example, the below command defines `AFR` population in the `1000G_GRCh37` reference. The `samples.txt` file stores list of sample names that will be included into `AFR`.

   ```bash
   flask create-population GRCh37 1000G AFR samples.txt
   ```

##### List loaded references

- Run the following command to list all references that are currently loaded into the server:

  ```bash
  flask show-references
  ```

##### List loaded genotype files

- Run the following command to list all genotype files loaded for a specified reference:

  ```bash
  flask show-genotypes <genome build> <reference name>
  ```

  Example:

  ```bash
  flask show-genotypes GRCh37 1000G
  ```

##### List loaded population samples

- Run the following command to list all sample names from a specified population in a specified reference:

  ```bash
  flask show-population <genome build> <reference name> <population name>
  ```

  Example:

  ```bash
  flask show-population GRCh37 1000G EUR
  ```

#### Running the rest app

1. Start Redis cache (provide configuration parameters depending on your needs):

   ```bash
   cd cget/bin/
   ./redis-server
   ```

2. Start LD server API:

   ```bash
   cd rest
   export FLASK_APP=ldserver
   flask run
   ```

   Or with `gunicorn`:

   ```bash
   gunicorn -b 127.0.0.1:[port] -w [n workers] -k gevent --pythonpath rest "ldserver:create_app()"
   ```

   In production you are advised to run with `gunicorn`, as it provides a number of asynchronous workers to handle requests. Using `flask run` as above starts only a single synchronous worker, which can be blocked by long running requests.

### playground app

This app provides a simple developer sandbox for learning the REST API.

#### Running the playground app

```bash
cd rest
export FLASK_APP=playground
flask run [ --port ... ]
```

Or with `gunicorn`:
```bash
gunicorn -b 127.0.0.1:[port] -w [n workers] -k gevent --pythonpath rest "playground:create_app()"
```

### raremetal app

This app serves a more specialized [REST API](docs/raremetal-api.md) for calculating aggregation tests of rare genetic variants. It provides covariance between variants, in addition to score statistics (summary statistics of the association between genetic variants and phenotypes.)

These statistics can be used by the [raremetal.js](https://github.com/statgen/raremetal.js) package to perform the aggregation tests. That package also provides documentation on the methodology.

Note: at times, we may refer to this also as the **raremetal server**.

#### Required data and formats

The server has 2 required types of data:

1. Genotype files (VCF, BCF, or Savvy format)
2. Phenotype files (TAB-delimited, or PED/DAT)

There is one optional type of data, which are mask files that group together variants and assign them to a particular gene or region. This data is optional if you wish to send your definitions of masks directly to the server via the covariance API request.

##### Genotype data

Genotype datasets may be in [VCF, BCF](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.pdf), or [Savvy](https://github.com/statgen/savvy) formats.

VCF files **must** be [bgzipped](http://www.htslib.org/doc/bgzip.html) and [tabixed](http://www.htslib.org/doc/tabix.html). In order to create a tabix-index for your VCF, it **must** be bgzipped first. Gzip is not sufficient. The `bgzip` program is included with `tabix` typically.

These two programs can be installed by downloading [htslib](http://www.htslib.org/download/) and compiling. They can also be installed via Ubuntu package: `sudo apt-get install tabix`, though this package may not be up to date at times depending on your Ubuntu version.

##### Phenotype data

Phenotype datasets are files containing a number of phenotypes (such as BMI, fasting glucose, heart rate, etc.) collected on a set of samples. Each row is a sample/individual, and each column is a phenotype.

There are two common formats for phenotype files - PED is probably the most common, but it is also common to use simple tab-delimited files. We support both formats currently.

PED format is described both on [our site](http://csg.sph.umich.edu/abecasis/merlin/tour/input_files.html) and by [PLINK](http://zzz.bwh.harvard.edu/plink/data.shtml).

For a tab-delimited file, the format is simply one header row denoting the column names, followed by a row for each individual. The first column is assumed to be the sample IDs. For example:

```
    fid      iid  patid  matid     sex  rand_binary   rand_qt
HG00096  HG00096      0      0  female           NA  0.283212
HG00097  HG00097      0      0    male          1.0  0.248001
HG00099  HG00099      0      0    male          0.0  0.691624
HG00100  HG00100      0      0    male          0.0  0.284674
HG00101  HG00101      0      0    male          1.0  0.494104
```

Categorical variables such as `sex` can be given labels as above. `rand_binary` and `rand_qt` are randomly generated phenotypes for this example.

Missing values should be encoded as `NA`.

##### Masks

A mask file maps genetic variants to "groups", which are typically genes (though they could also be arbitrary genomic regions.) Typically mask files are created using variant filters, such as "allele frequency < 0.01" or "protein truncating variants AND loss-of-function variants".

The mask file is tab-delimited, and has one "group" per row. For example:

```
ZBED4   22  50276999  50283712  22:50276999_C/T 22:50277035_C/T 22:50277046_A/C 22:50277051_C/T 22:50277087_T/C 22:50277093_G/A
ALG12   22  50293950  50312062  22:50293950_A/G 22:50293955_G/T 22:50294146_C/T 22:50294162_C/G 22:50294175_G/A 22:50294301_G/A
```

Each row begins with 4 columns:

* The group name (in this case they are genes)
* Chromosome
* Start position of the group
* End position of the group

The remaining values on each row are the variants that are assigned to the group. They should also be tab-delimited.

The mask file must then be bgzipped, and tabix indexed. An easy way to do this is:

```bash
bgzip my_mask_file.tab
tabix -s 2 -b 3 -e 4 my_mask_file.tab.gz
```

##### Summary statistics

Some studies are not able or willing to share their genotypes and/or phenotype data. In such cases, they provide summary statistics instead, which can be used to run various association analyses or meta-analysis.

The LDServer can directly serve these summary statistics (score statistics, covariance matrices) rather than computing them on the fly.

Common programs for computing summary statistics as a frame of reference:

* [rvtests](https://github.com/zhanxw/rvtests#meta-analysis-models)
* [RAREMETALWORKER](https://genome.sph.umich.edu/wiki/RAREMETALWORKER)

Both programs produce two files, one containing score statistics for single variants, the other file containing the covariance of the score statistics. You can add these files to a [Summary statistic record](#summary-statistic-records) under the `score_path` and `cov_path` keys.

#### Configuring the raremetal app

Before running the following commands, you should set the appropriate flask app by doing:

```bash
export FLASK_APP="rest/raremetal"
```

To add datasets to the server, we recommend using [the YAML config approach](#yaml-adding-multiple-datasets-with-one-yaml-config-file). This allows you to specify all of your datasets in a single file, and to maintain them across server restarts.

You can also add datasets via CLI commands for quick testing, but we do not recommend this for production usage.

##### YAML: Adding multiple datasets with one YAML config file

Datasets may be specified in a YAML file, which provides additional features over the CLI interface:

1. Set your own dataset IDs. This allows you to maintain consistency on IDs when reloading the database.
2. Specify columns that are not meant to be used for analysis. These columns will not show up in metadata queries.
3. Specify longer text descriptions for each column.

An example YAML file is located in `data/test.yaml`, and looks like:

```YAML
genotypes:
- id: 1
  name: "1000G"
  description: "1000G chr22 Testing VCF"
  filepath: "data/chr22.test.vcf.gz"
  genome_build: "GRCh37"

phenotypes:
- id: 1
  name: "1000G random phenotypes"
  description: "An example set of randomly generated phenotypes for 1000G"
  genotypes: 1
  filepath: "data/chr22.test.tab"
  delim: "\t"
  columns:
    iid:
      column_type: "TEXT"
      sample_column: true

    sex:
      column_type: "CATEGORICAL"
      for_analysis: false

    rand_binary:
      column_type: "CATEGORICAL"
      description: "A random binary phenotype"

    rand_qt:
      column_type: "FLOAT"
      description: "A random quantitative phenotype"

- id: 2
  name: "1000G random phenotypes II"
  genotypes: 1
  filepath: "data/chr22.more_phenotypes.test.ped"
  description: "Adding a second set of phenotypes for 1000G"
  delim: "\t"
  columns:
    ANOTHER_RAND_QT:
      column_type: "FLOAT"
      description: "Another random quantitative phenotype"

masks:
- id: 1
  name: "AF < 0.01"
  description: "Variants with allele frequency < 1%"
  filepath: "data/mask.epacts.chr22.gencode-exons-AF01.tab.gz"
  genome_build: "GRCh37"
  genotypes: 1
  group_type: "GENE"
  identifier_type: "ENSEMBL"

- id: 2
  name: "AF < 0.05"
  description: "Variants with allele frequency < 5%"
  filepath: "data/mask.epacts.chr22.gencode-exons-AF05.tab.gz"
  genome_build: "GRCh37"
  genotypes: 1
  group_type: "GENE"
  identifier_type: "ENSEMBL"

- id: 3
  name: "All exomic variants"
  description: "All variants falling within exomic regions in TOPMed"
  summary_stats: 1
  genome_build: "GRCh37"
  group_type: "GENE"
  identifier_type: "ENSEMBL"

summary_stats:
- id: 1
  name: "TOPMed 70K Exomes - T2D"
  description: "TOPMed 70K exomes score statistics and covariance matrices for T2D"
  genome_build: "GRCh37"
  score_path: "topmed.metascore.txt.gz"
  cov_path: "topmed.metacov.txt.gz"
```

The file has 4 possible "blocks": `genotypes`, `phenotypes`, `masks`, and `summary-stats`.

###### Genotype records

Each record under `genotypes` looks like:

```YAML
- id: <genotype dataset ID you would like to use>
  name: <short description of genotype dataset, typically a study name>
  description: <long form description of dataset>
  filepath: <path to VCF, BCF, or Savvy file>
  genome_build: <genome build, e.g. GRCh37 or GRCh38>
```

###### Phenotype records

Each record under `phenotypes` looks like:

```YAML
- id: <phenotype dataset ID you would like to use>
  name: <short description of phenotypes>
  description: <long description of phenotypes>
  genotypes: <genotype dataset ID, i.e. which genotypes do these samples line up with?>
  filepath: "data/chr22.test.tab"
  delim: "\t"
  columns:
    <column name as it exists in header or DAT file for PED files>:
      column_type: <type, can be "TEXT" or "CATEGORICAL" or "FLOAT"
      sample_column: <true or false, is this column the sample ID column?>
      for_analysis: <true or false, should this phenotype be used in analysis?>
      description: <long form description of the phenotype>
```

Under `columns`, the following keys are optional:

* `sample_column`: (if no columns specified as the sample ID column, we assume it is the first column in the file.)
* `for_analysis`: all columns are assumed to be used for analysis unless specified otherwise.
* `column_type`: if not specified, the type will be guessed by examining the phenotype file.

Also, specifying columns is entirely optional itself. Column names and types will be deduced if not provided.

For PED files, you do not need to specify anything about the first 5 columns (family ID, individual ID, paternal ID, maternal ID, sex) - they are part of the format and automatically handled. If you specify information for the remaining phenotypes in the file, note that they must match the phenotype as specified in the DAT file.

###### Mask records

Each record under `masks` looks like:

```YAML
- id: <mask ID you would like to use>
  name: <short description of mask>
  description: <long description of mask>
  filepath: "data/mask.epacts.chr22.gencode-exons-AF01.tab.gz"
  genome_build: <genome build, e.g. GRCh37>
  genotypes: <genotype dataset IDs containing the variants specified in this mask>
  summary_stats: <summary stat dataset IDs>
  group_type: <group type, can be "GENE" or "REGION">
  identifier_type: <identifier type, currently only "ENSEMBL" supported>
```

A mask can be linked to any number of genotype datasets, or summary statistic datasets. For example:

```YAML
- id: 1
  name: "All possible protein truncating variants in the genome"
  genotypes: [1, 2, 3]
  summary_stats: [1, 4, 7]
  ...
```

###### Summary statistic records

Each record under `summary_stats` looks like:

```YAML
summary_stats:
- id: 1
  name: "TOPMED-70K-EXOMES"
  description: "TOPMed pre-computed covariance for 70K exomes"
  genome_build: "GRCh37"
  score_path: "topmed.metascore.txt.gz"
  cov_path: "topmed.metacov.txt.gz"
```

This section is for pre-calculated meta-analysis summary statistics, usually generated by a program such as RAREMETALWORKER or rvtests.

For example, in rvtests: https://github.com/zhanxw/rvtests#meta-analysis-models

The score statistic (`score_path`) and covariance matrix (`cov_path`) files must both be tabix indexed. RAREMETALWORKER or rvtest should do this automatically for you (there will be a `.tbi` tabix index file created for each score/cov file.)

If your files are split by chromosome, you can specify both paths as a glob with the `*` character. For example:

```YAML
  score_path: "topmed.chr*.metascore.txt.gz"
  cov_path: "topmed.chr*.metascore.txt.gz"
```

###### Add all records in YAML

The YAML file may have genotype, phenotype, mask, and summary-stats blocks specified in any order, however there may only be 1 block for each type.

Additionally, genotypes will always be processed first, followed by phenotypes, and finally mask files. This is because the latter two are dependent on having at least 1 genotype dataset available.

To add the datasets specified by the YAML file:

```bash
flask add-yaml <path/to/yamlfile>
```

##### CLI: Quickly add data for testing purposes

###### Adding genotype datasets

To add the dataset to the server, use the `add-genotypes` command:

```bash
flask add-genotypes <short label> <long description> <genome build> <VCF/BCF/Savvy file>
```

The parameters:

* `<short label>` is a short description of the dataset, often a study abbreviation.
* `<long description>` a longer description that may include the genotyping platform or sequencing.
* `<genome build>` is the genome build of the positions in the file.
* `<VCF/BCF/Savvy file>` the file containing genotypes for variants over a number of samples. We recommend placing these files (or symbolic links) to these files in the `data/` directory under the application root. You may also provide a glob of files, in the event your genotypes are separated into files by chromosome.

Optional parameters:

* `--samples <file>` provide a file with a list of samples to use, if you do not wish to use all of the samples in the genotype file. One sample per line.

As an example:

```bash
gid=`flask add-genotypes "1000G" "1000G Test VCF" "GRCh37" data/chr*.test.vcf.gz`
```

With the command above, you can capture the genotype dataset ID that was assigned in the database, and use it in later commands.

###### Adding phenotype datasets

To add a phenotype dataset:

```bash
flask add-phenotypes <short label> <long label> <path to ped or tab file> <genotype dataset ID>
```

The parameters:

* `<short label>` is a short description of the dataset, often a study abbreviation.
* `<long description>` a longer description that may describe the set of phenotypes.
* `<path to ped or tab file>` path to either the PED or tab file. If a PED file is supplied, it is assumed there is an accompanying DAT file. If supplying a tab-delimited file, it must end with a `.tab`.
* `<genotype dataset ID>` ID of the genotype dataset to connect this phenotype dataset with.

As an example:

```bash
flask add-phenotypes Test "Test 1" "data/chr22.test.tab" $gid
```

Where `$gid` was set above when adding the genotype dataset for these phenotypes.

Each phenotype dataset added should correspond to a genotype dataset (and hence supplying the genotype dataset ID.) This allows the server to know which phenotypes are available (and valid) for each dataset.

###### Adding masks of genetic variants

To add your mask file to the server:

```bash
flask add-masks <name> <description> <path to mask file> <genome build> <genotype dataset ID> <group type> <identifier type>
```

The parameters:

* `<name>` the name of the mask which will be used in queries. This should ideally be short, for example "AF<0.01 & PTV & LoF".
* `<description>` a longer description of what went into creating the mask. For example: "Allele frequency < 1% and all protein truncating & loss-of-function variants".
* `<path to mask file>` path to the mask file. A tabix-index is assumed to reside next to the file. For example, the mask file is `mask.tab.gz` and this is the file you would provide, there should also be `mask.tab.gz.tbi` as well.
* `<genome build>` genome build the positions are anchored to
* `<genotype dataset ID>` provide the genotype dataset ID this mask corresponds to. Masks are typically generated for a specific set of variants provided in a genotype file.
* `<group type>` can be "GENE", or "REGION"
* `<identifier type>` what is the identifier for each group? Currently only "ENSEMBL" is supported.

For example:

```bash
flask add-masks "AF < 0.01" "Variants with alternate allele freq < 0.01" "data/mask.epacts.chr22.gencode-exons-AF01.tab.gz" "GRCh37" $gid "GENE" "ENSEMBL"
```

Where `$gid` was the genotype dataset ID generated after using the `add-genotypes` command earlier.

###### Listing available datasets

You can see the list of phenotypes, genotypes, and masks that have been given to the server using:

* `flask show-genotypes`
* `flask show-phenotypes`
* `flask show-masks`

#### Running the raremetal app

For quickly starting a server:

```bash
export FLASK_APP=rest/raremetal
flask run
```

You can enable debug mode by also including `export FLASK_DEBUG=1`.

For production, use `gunicorn`:

```bash
gunicorn -b 127.0.0.1:[port] -w [n workers] -k gevent --pythonpath rest "raremetal:create_app()"
```

#### Exposed API endpoints

The full API specification can be found in [a separate document](docs/raremetal-api.md).

[raremetal.js] and [locuszoom.js] have built-in support for these APIs. You may only need to parse the `/aggregation/metadata` endpoint below to provide UI for users selecting which genotypes/phenotypes/masks they would like to use. After that, locuszoom can leverage raremetal.js to make the appropriate calls and parse the response accordingly.

As a brief summary: there are two primary endpoints of interest:

1. GET `/aggregation/metadata` - this endpoint shows available datasets
2. POST `/aggregation/covariance` - computes covariance and score statistics within a region

The metadata endpoint returns JSON that looks like:

```json
{
    "data": [
        {
            "description": "1000G chr22 Testing VCF",
            "genomeBuild": "GRCh37",
            "genotypeDataset": 1,
            "masks": [
                {
                    "id": 1,
                    "name": "AF < 0.01",
                    "description": "Variants with allele frequency < 1%",
                    "groupType": "GENE",
                    "identifierType": "ENSEMBL",
                },
                {
                    "id": 2,
                    "name": "AF < 0.05",
                    "description": "Variants with allele frequency < 5%",
                    "groupType": "GENE",
                    "identifierType": "ENSEMBL",
                }
            ],
            "name": "1000G",
            "phenotypeDatasets": [
                {
                    "description": "An example set of randomly generated phenotypes for 1000G",
                    "name": "1000G random phenotypes",
                    "phenotypeDataset": 1,
                    "phenotypes": [
                        {
                            "description": "A random binary phenotype",
                            "name": "rand_binary"
                        },
                        {
                            "description": "A random quantitative phenotype",
                            "name": "rand_qt"
                        }
                    ]
                },
                {
                    "description": "Adding a second set of phenotypes for 1000G",
                    "name": "1000G random phenotypes II",
                    "phenotypeDataset": 2,
                    "phenotypes": [
                        {
                            "description": "Another random quantitative phenotype",
                            "name": "ANOTHER_RAND_QT"
                        }
                    ]
                }
            ]
        }
    ]
}
```

Each genotype dataset (think VCF) has a number of masks and phenotype datasets associated with it.

The covariance endpoint is a POST request, the query is a JSON object of the following form:

```json
{
	"chrom": "22",
	"start": 50276998,
	"stop": 50357719,
	"genotypeDataset": 1,
	"phenotypeDataset": 1,
	"phenotype": "rand_qt",
	"samples": "ALL",
	"genomeBuild": "GRCh37",
	"masks": [1]
}
```

This tells the covariance endpoint to compute within the region 22:50276998-50357719, using the genotype dataset with ID 1, phenotype dataset with ID 1, and mask with ID 1. The phenotype will be `rand_qt` (in the test/example data, the phenotypes are randomly generated, since we can't release actual data.)

Typically you will always set `samples` to `ALL`. The server supports arbitrary subsets of samples, such as sub-populations, that can be specified with the `flask create-sample-subset` command.

Rather than using a server-side mask and supplying an ID (as in the above example), alternatively you can send your own [mask definitions](https://github.com/statgen/raremetal.js/blob/master/src/docs/portal-api.md#request-1) created in the browser/client to the server.

The response will look like:

```json
{
  "data": {
    "genotypeDataset": 1,
    "description": "1000G chr22 Testing VCF",
    "phenotypeDataset": 1,
    "phenotype": "rand_qt",
    "sigmaSquared": 0.081,
    "nSamples": 2504,
    "variants": [
      {
        "variant": "2:21228642_G/A",
        "altFreq": 0.033,
        "pvalue": 0.000431,
        "score": 0.1
      }
    ],
    "groups": [
      {
        "groupType": "gene",
        "group": "ENSG000001",
        "mask": 1,
        "variants": ["2:21228642_G/A"],
        "covariance": [0.3],
      }
    ]
  }
}
```

The first few entries are simply repeating back the request parameters, namely `genotypeDataset`, `phenotypeDataset`, `phenotype`.

The remaining parameters:

* `sigmaSquared` - this is the variance of the phenotype, in this example `rand_qt` is the phenotype.
* `nSamples` - number of samples that went into the analysis
* `variants` - an array of variant objects, each containing a score statistic, p-value, and alt allele frequency
* `groups` - an array of group objects. Each group object is the result of calculating covariance for each combination of (group, mask). Within each group, covariance is returned as the linearized upper triangle of the covariance matrix.

[raremetal.js]: https://github.com/statgen/raremetal.js
[locuszoom.js]: https://github.com/statgen/locuszoom
