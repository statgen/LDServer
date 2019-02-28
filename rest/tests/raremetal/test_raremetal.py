import re

def test_malformed_chrom(client):
    resp = client.post("/aggregation/covariance", data = {
        "chrom": "22REERERFFAFA",
        "start": 50276998,
        "stop": 50357719,
        "genotypeDataset": 1,
        "phenotypeDataset": 1,
        "phenotype": "rand_qt",
        "samples": "ALL",
        "genomeBuild": "GRCh37",
        "masks": ["AF < 0.01"]
    })

    assert resp.status_code == 400
    assert resp.is_json
    assert re.search("Chromosome.*.not found",resp.json["error"]) is not None

def test_invalid_position(client):
    resp = client.post("/aggregation/covariance", data = {
        "chrom": "22",
        "start": -1,
        "stop": -1,
        "genotypeDataset": 1,
        "phenotypeDataset": 1,
        "phenotype": "rand_qt",
        "samples": "ALL",
        "genomeBuild": "GRCh37",
        "masks": ["AF < 0.01"]
    })

    assert resp.status_code == 400
    assert resp.is_json
    assert re.search("Error while parsing.*start.*", resp.json["error"]) is not None

def test_invalid_mask(client):
    resp = client.post("/aggregation/covariance", data = {
        "chrom": "22",
        "start": 50276998,
        "stop": 50357719,
        "genotypeDataset": 1,
        "phenotypeDataset": 1,
        "phenotype": "rand_qt",
        "samples": "ALL",
        "genomeBuild": "GRCh37",
        "masks": ["bad mask"]
    })

    assert resp.status_code == 400
    assert resp.is_json
    assert re.search("No mask exists by the name", resp.json["error"]) is not None

def test_invalid_phenotype(client):
    resp = client.post("/aggregation/covariance", data = {
        "chrom": "22",
        "start": 50276998,
        "stop": 50357719,
        "genotypeDataset": 1,
        "phenotypeDataset": 1,
        "phenotype": "invalid_phenotype",
        "samples": "ALL",
        "genomeBuild": "GRCh37",
        "masks": ["AF < 0.01"]
    })

    assert resp.status_code == 400
    assert resp.is_json
    assert re.search("Phenotype.*does not exist.*", resp.json["error"]) is not None

def test_covar(client):
    resp = client.post("/aggregation/covariance", data = {
        "chrom": "22",
        "start": 50276998,
        "stop": 50357719,
        "genotypeDataset": 1,
        "phenotypeDataset": 1,
        "phenotype": "rand_qt",
        "samples": "ALL",
        "genomeBuild": "GRCh37",
        "masks": ["AF < 0.01"]
    })

    assert resp.status_code == 200
    assert resp.is_json

    for group in resp.json["data"]["groups"]:
        n_variants = len(group["variants"])
        n_covar = len(group["covariance"])
        assert n_covar == (n_variants * (n_variants + 1) / 2)

    for variant in resp.json["data"]["variants"]:
        assert variant["altFreq"] > 0
        assert variant["pvalue"] > 0
        assert variant["pvalue"] <= 1