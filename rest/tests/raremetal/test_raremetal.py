import re

def is_sha(s):
    return re.search("[A-Fa-f0-9]+",s) is not None

def test_status(client):
    resp = client.get("/status")

    assert resp.status_code == 200
    assert "sha" in resp.json["data"]
    sha = resp.json["data"]["sha"]
    assert sha == "no-git" or is_sha(sha)

def test_trailing_slash(client):
    resp = client.get("/aggregation/metadata/")
    assert resp.status_code == 404

def test_metadata(client):
    resp = client.get("/aggregation/metadata")
    assert resp.status_code == 200

    for genotype_dataset in resp.json["data"]:
        assert "name" in genotype_dataset
        assert "masks" in genotype_dataset
        assert "description" in genotype_dataset
        assert "genomeBuild" in genotype_dataset
        assert "genotypeDataset" in genotype_dataset
        assert "phenotypeDatasets" in genotype_dataset

        for phenotype_dataset in genotype_dataset["phenotypeDatasets"]:
            assert "description" in phenotype_dataset
            assert "phenotypes" in phenotype_dataset
            assert "name" in phenotype_dataset
            assert "phenotypeDataset" in phenotype_dataset

            for phenotype in phenotype_dataset["phenotypes"]:
                assert "description" in phenotype
                assert "name" in phenotype

        for mask in genotype_dataset["masks"]:
            assert "groupType" in mask
            assert "identifierType" in mask
            assert "description" in mask
            assert "name" in mask
            assert "id" in mask

            assert mask["groupType"] == mask["groupType"].upper()
            assert mask["identifierType"] == mask["identifierType"].upper()



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
        "masks": [1]
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
        "masks": [1]
    })

    assert resp.status_code == 400
    assert resp.is_json
    assert re.search("Error while parsing", resp.json["error"]) is not None

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
        "masks": [-9]
    })

    assert resp.status_code == 400
    assert resp.is_json
    assert re.search("No mask exists for ID", resp.json["error"]) is not None

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
        "masks": [1]
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
        "masks": [1]
    })

    assert resp.status_code == 200
    assert resp.is_json

    score_variants = [x["variant"] for x in resp.json["data"]["variants"]]

    for group in resp.json["data"]["groups"]:
        n_variants = len(group["variants"])
        n_covar = len(group["covariance"])
        assert n_covar == (n_variants * (n_variants + 1) / 2)
        assert isinstance(group["mask"], int)
        assert "group" in group
        assert "groupType" in group
        assert group["groupType"] in ("REGION", "GENE")
        assert all([v in score_variants for v in group["variants"]])

    for variant in resp.json["data"]["variants"]:
        assert variant["altFreq"] > 0
        assert variant["pvalue"] > 0
        assert variant["pvalue"] <= 1
        assert "score" in variant

def test_monomorphic(client):
    resp = client.post("/aggregation/covariance", json = {
        "chrom": "22",
        "start": 50276998,
        "stop": 50357719,
        "genotypeDataset": 2, # includes 3 monomorphic variants in chr22.monomorphic_test.vcf.gz
        "phenotypeDataset": 1,
        "phenotype": "rand_qt",
        "samples": "ALL",
        "genomeBuild": "GRCh37",
        "masks": [1]
    })

    assert resp.status_code == 200
    assert resp.is_json

    score_variants = [x["variant"] for x in resp.json["data"]["variants"]]

    for group in resp.json["data"]["groups"]:
        n_variants = len(group["variants"])
        n_covar = len(group["covariance"])
        assert n_covar == (n_variants * (n_variants + 1) / 2)
        assert isinstance(group["mask"], int)
        assert "group" in group
        assert "groupType" in group
        assert group["groupType"] in ("REGION", "GENE")
        assert all([v in score_variants for v in group["variants"]])

    for variant in resp.json["data"]["variants"]:
        assert variant["altFreq"] > 0
        assert variant["pvalue"] > 0
        assert variant["pvalue"] <= 1
        assert "score" in variant