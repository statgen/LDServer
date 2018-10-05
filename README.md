# LDServer
LD Server For Web

## Installation

- Clone repository.
  ```
  git clone https://github.com/statgen/LDServer
  cd LDServer
  ```

- Create and activate virtual environment.
  ```
  python -m virtualenv env
  source env/bin/activate
  ```

- Install and test LDServer C++ shared library. The '127.0.0.1' and '8888' arguments are the hostname and port for Redis server which will be started during unit tests. If you wish to run Redis server on another port, then change these values correspondingly.
  ```
  cget install ld
  cd cget/test/ && ./testAll 127.0.0.1 8888 && cd ../../
  ```

- Install required python packages for REST API.
  ```
  cd rest
  pip install -r requirements.txt
  ```

- Test REST API.
  ```
  python -m pytest
  ```

- Configure REST API by editing `config/default.py` file, or by creating new `instance` folder with `config.py` configuration file.

## Managing references
### Add new reference
- Execute `flask add-reference` command:
   ```
   cd rest
   export FLASK_APP=rest
   flask add-reference <name> <description> <genome build> <samples file> <genotype files>
   ```
   For example, the below command adds a new reference panel named `1000G_GRCh37`. The `samples.txt` file stores list of sample names in `ALL.chr*.bcf` files. 
   ```
   flask add-reference 1000G_GRCh37 "1000 Genomes Project Phase 3" GRCh37 samples.txt ALL.chr*.bcf
   ```
   The genotypes can be stored in VCF, BCF, and SAV formats. For better runtime performance and more compact storage, we highly recommend using SAV format.
### Add new (sub-)population to the reference
- To define populations in the reference, execute `flask create-population` command:
   ```
   cd rest
   export FLASK_APP=rest
   flask create-population <reference name> <population name> <samples file>
   ```
   For example, the below command defines `AFR` population in the `1000G_GRCh37` reference. The `samples.txt` file stores list of sample names that will be included into `AFR`. 
   ```
   flask create-population 1000G_GRCh37 AFR samples.txt
   ```
### List loaded references
- Run the following command to list all references that are currently loaded into the server:
  ```
  flask show-references
  ```  
### List loaded genotype files
- Run the following command to list all genotype files loaded for a specified reference:
  ```
  flask show-genotypes <reference name>
  ```
  Example:
  ```
  flask show-genotypes 1000G_GRCh37
  ```
### List loaded population samples
- Run the following command to list all sample names from a specified population in a specified reference:
  ```
  flask show-population <reference name> <population name>
  ```
  Example:
  ```
  flask show-population 1000G_GRCh37 EUR
  ```
   
## Run REST API
1. Start Redis cache (provide configuration parameters depending on your needs):
   ```
   cd cget/bin/
   ./redis-server
   ```
2. Start LD server API:
   ```
   cd rest
   export FLASK_APP=rest
   flask run
   ```
   with `Gunicorn`
   ```
   gunicorn -b 127.0.0.1:[port] -w [n workers] -k gevent "rest:create_app()"
   ```

## Run Developers Playground
   ```
   cd rest
   export FLASK_APP=playground
   flask run [ --port ... ] 
   ```
   with `Gunicorn`
   ```
   gunicorn -b 127.0.0.1:[port] -w [n workers] -k gevent "playground:create_app()"
   ```
