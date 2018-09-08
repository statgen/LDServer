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

- Load your LD references from `REFERENCES_JSON` to SQL database:
```
flask load-references
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
