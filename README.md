# LDServer
LD Server For Web

## Installation

- Create and activate virtual environment:
  ```
  python -m virtualenv env
  source env/bin/activate
  ```

- Install and test LDServer C++ shared library
  ```
  cget install ld
  cd cget/test/ && ./testAll && cd ../../
  ```

- Install required python packages for REST API
  ```
  cd rest
  pip install -r requirements.txt
  ```

- Test REST API
  ```
  python -m pytest
  ```

- Configure REST API by editing `config/default.py` file, or by creating new `instance` folder with `config.py` configuration file.

## Managing references
### Add new reference
- Execute `flask add-reference` command:
   ```
   export FLASK_APP=rest
   flask add-reference <name> <description> <genome build> <samples file> <genotype files>
   ```
   For example, the below command adds a new reference panel named `1000G_GRCh37`. The `samples.txt` file stores list of sample names in `ALL.chr*.bcf` files. 
   ```
   flask add-reference 1000G_GRCh37 "1000 Genomes Project Phase 3" GRCh37 samples.txt ALL.chr*.bcf
   ```
   The genotypes can be stored in VCF, BCF, and SAV formats. For better runtime performance and more compact storage, we highly recommend using SAV format.
### Add new (sub-)population to the reference
- To define populations in the reference, execute `flask add-reference` command:
   ```
   export FLASK_APP=rest
   flask add-reference <reference name> <population name> <samples file>
   ```
   For example, the below command defines `AFR` population in the `1000G_GRCh37` reference. The `samples.txt` file stores list of sample names that will be included into `AFR`. 
   ```
   flask add-reference 1000G_GRCh37 AFR samples.txt
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
```
export FLASK_APP=rest
flask run
```

## Run Developers Playground
```
export FLASK_APP=playground
flask run [ --port ... ] 
```
