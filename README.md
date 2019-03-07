# LDServer

[![Build Status](https://travis-ci.com/statgen/LDServer.svg?branch=master)](https://travis-ci.com/statgen/LDServer)

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

<!-- TOC depthFrom:1 depthTo:7 withLinks:1 updateOnSave:1 orderedList:0 -->

- [LDServer](#ldserver)
	- [Documentation](#documentation)
	- [Installation](#installation)
	- [Updating](#updating)
	- [Configuring & running the flask apps](#configuring-running-the-flask-apps)
		- [rest app](#rest-app)
			- [Configuring the rest app](#configuring-the-rest-app)
				- [Add new reference](#add-new-reference)
				- [Add new (sub-)population to the reference](#add-new-sub-population-to-the-reference)
				- [List loaded references](#list-loaded-references)
				- [List loaded genotype files](#list-loaded-genotype-files)
				- [List loaded population samples](#list-loaded-population-samples)
			- [Running the rest app](#running-the-rest-app)
		- [playground app](#playground-app)
			- [Running the playground app](#running-the-playground-app)
		- [raremetal app](#raremetal-app)
			- [Configuring the raremetal app](#configuring-the-raremetal-app)
				- [CLI: Adding genotype datasets](#cli-adding-genotype-datasets)
				- [CLI: Adding phenotype datasets](#cli-adding-phenotype-datasets)
				- [CLI: Adding masks of genetic variants](#cli-adding-masks-of-genetic-variants)
				- [CLI: Listing available datasets](#cli-listing-available-datasets)
				- [YAML: Adding multiple datasets with one YAML config file](#yaml-adding-multiple-datasets-with-one-yaml-config-file)
			- [Running the raremetal app](#running-the-raremetal-app)
			- [APIs](#apis)

<!-- /TOC -->

## Installation

Before installing, the following software/packages are required:

* CMake 3+
* gcc 5.4+
* Python 2.7
* A BLAS implementation, such as OpenBLAS or Intel MKL

The following software is optional:

* Redis

Most Linux distributions provide these packages already. For example, in Ubuntu:

```
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
  virtualenv -p python2 venv &&
    source venv/bin/activate &&
    pip install backports.lzma cget &&
    cget install ld
  ```

- *Optional*: Test the LDServer C++ shared library. The '127.0.0.1' and '8888' arguments are the hostname and port for Redis server which will be started during unit tests. If you wish to run Redis server on another port, then change these values correspondingly.

  ```
  cd cget/test/ && ./testAll 127.0.0.1 8888 && cd ../../
  ```

- Install required python packages for the flask apps.

  ```
  # Activate the environment
  source venv/bin/activate

  # Install required packages
  pip install -r rest/requirements.txt
  ```

- *Optional*: Run the flask app test suite.

  ```
  cd rest
  python -m pytest
  ```

- Configure REST API by copying the config file `rest/config/default.py` to `rest/instance/config.py` and then modify values as needed. Modifying `default.py` is not recommended as it will be overwritten when pulling updates from github.

## Updating

Update your files using `git pull`, then recompile the server by doing `cget install --update ld`.

## Configuring & running the flask apps

### rest app

This app serves r2, D', and covariance of genetic variants over a REST API.

#### Configuring the rest app

##### Add new reference

- Execute `flask add-reference` command:
   ```
   cd rest
   export FLASK_APP=rest
   flask add-reference <name> <description> <genome build> <samples file> <genotype files>
   ```
   For example, the below command adds a new reference panel named `1000G_GRCh37`. The `samples.txt` file stores list of sample names in `ALL.chr*.bcf` files.
   ```
   flask add-reference 1000G "1000 Genomes Project Phase 3" GRCh37 samples.txt ALL.chr*.bcf
   ```
   The genotypes can be stored in VCF, BCF, and SAV formats. For better runtime performance and more compact storage, we highly recommend using SAV format.

##### Add new (sub-)population to the reference

- To define populations in the reference, execute `flask create-population` command:

   ```
   cd rest
   export FLASK_APP=rest
   flask create-population <genome build> <reference name> <population name> <samples file>
   ```

   For example, the below command defines `AFR` population in the `1000G_GRCh37` reference. The `samples.txt` file stores list of sample names that will be included into `AFR`.

   ```
   flask create-population GRCh37 1000G AFR samples.txt
   ```

##### List loaded references

- Run the following command to list all references that are currently loaded into the server:

  ```
  flask show-references
  ```

##### List loaded genotype files

- Run the following command to list all genotype files loaded for a specified reference:

  ```
  flask show-genotypes <genome build> <reference name>
  ```

  Example:

  ```
  flask show-genotypes GRCh37 1000G
  ```

##### List loaded population samples

- Run the following command to list all sample names from a specified population in a specified reference:

  ```
  flask show-population <genome build> <reference name> <population name>
  ```

  Example:

  ```
  flask show-population GRCh37 1000G EUR
  ```

#### Running the rest app

1. Start Redis cache (provide configuration parameters depending on your needs):

   ```
   cd cget/bin/
   ./redis-server
   ```

2. Start LD server API:

   ```
   cd rest
   export FLASK_APP=rest/rest
   flask run
   ```

   Or with `gunicorn`:

   ```
   gunicorn -b 127.0.0.1:[port] -w [n workers] -k gevent --pythonpath rest "rest:create_app()"
   ```

   In production you are advised to run with `gunicorn`, as it provides a number of asynchronous workers to handle requests. Using `flask run` as above starts only a single synchronous worker, which can be blocked by long running requests.

### playground app

This app provides a simple developer sandbox for learning the REST API.

#### Running the playground app

```
cd rest
export FLASK_APP=rest/playground
flask run [ --port ... ]
```

Or with `gunicorn`:
```
gunicorn -b 127.0.0.1:[port] -w [n workers] -k gevent --pythonpath rest "playground:create_app()"
```

### raremetal app

This app serves a more specialized [REST API](https://github.com/statgen/raremetal.js/blob/master/docs/portal-api.pdf) for calculating aggregation tests of rare genetic variants. It provides covariance between variants, in addition to score statistics (summary statistics of the association between genetic variants and phenotypes.)

These statistics can be used by the [raremetal.js](https://github.com/statgen/raremetal.js) package to perform the aggregation tests. That package also provides documentation on the methodology.

In addition to genotype files (VCFs, BCFs, Savvy), the raremetal app also requires phenotypes. Phenotype files are typically PED format, or an ad-hoc tab-delimited format.

#### Configuring the raremetal app

Before running the following commands, you should set the appropriate flask app by doing:

```bash
export FLASK_APP="rest/raremetal"
```

There are two methods for inserting datasets:

1. Quickly add datasets with the command-line interface (CLI) commands
2. Provide metadata about all of your datasets at once in a YAML file

**The second option is our recommendation for production**, but to get up and running quickly or for testing, the first option works well.

##### CLI: Adding genotype datasets

Genotype datasets may be in VCF, BCF, or [Savvy](https://github.com/statgen/savvy) formats. VCF files must be [bgzipped](http://www.htslib.org/doc/bgzip.html) and [tabixed](http://www.htslib.org/doc/tabix.html).

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


##### CLI: Adding phenotype datasets

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

##### CLI: Adding masks of genetic variants

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

##### CLI: Listing available datasets

You can see the list of phenotypes, genotypes, and masks that have been given to the server using:

* `flask show-genotypes`
* `flask show-phenotypes`
* `flask show-masks`

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
```

The file has 3 "blocks": `genotypes`, `phenotypes`, and `masks`.

Each record under `genotypes` looks like:

```YAML
- id: <genotype dataset ID you would like to use>
  name: <short description of genotype dataset, typically a study name>
  description: <long form description of dataset>
  filepath: <path to VCF, BCF, or Savvy file>
  genome_build: <genome build, e.g. GRCh37 or GRCh38>
```

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

Each record under `masks` looks like:

```YAML
- id: <mask ID you would like to use>
  name: <short description of mask>
  description: <long description of mask>
  filepath: "data/mask.epacts.chr22.gencode-exons-AF01.tab.gz"
  genome_build: <genome build, e.g. GRCh37>
  genotypes: <genotype dataset ID containing the variants specified in this mask>
  group_type: <group type, can be "GENE" or "REGION">
  identifier_type: <identifier type, currently only "ENSEMBL" supported>
```

The YAML file may have genotype, phenotype, and mask blocks specified in any order, however there may only be 1 block for each type. Additionally, genotypes will always be processed first, followed by phenotypes, and finally mask files. This is because the latter two are dependent on having at least 1 genotype dataset available.

To add the datasets specified by the YAML file:

```bash
flask add-yaml <path/to/yamlfile>
```

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

#### APIs

The [full API specification](https://github.com/statgen/raremetal.js/blob/master/docs/portal-api.pdf) can be found in the [raremetal.js](https://github.com/statgen/raremetal.js) package.

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
