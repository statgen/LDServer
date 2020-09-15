import re
import pytest

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

    for dataset in resp.json["data"]:
        assert "name" in dataset
        assert "masks" in dataset
        assert "description" in dataset
        assert "genomeBuild" in dataset

        assert ("genotypeDataset" in dataset and "phenotypeDatasets" in dataset) or "summaryStatisticDataset" in dataset

        for phenotype_dataset in dataset.get("phenotypeDatasets", []):
            assert "description" in phenotype_dataset
            assert "phenotypes" in phenotype_dataset
            assert "name" in phenotype_dataset
            assert "phenotypeDataset" in phenotype_dataset

            for phenotype in phenotype_dataset["phenotypes"]:
                assert "description" in phenotype
                assert "name" in phenotype

        for mask in dataset["masks"]:
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

    assert resp.status_code == 200
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

    assert "phenotypeDataset" in resp.json["data"]
    assert "genotypeDataset" in resp.json["data"]
    assert "phenotype" in resp.json["data"]

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

def test_summary_stat(client):
    resp = client.post("/aggregation/covariance", data = {
        "chrom": "1",
        "start": 2,
        "stop": 307,
        "summaryStatDataset": 2,
        "genomeBuild": "GRCh37",
        "masks": [3]
    })

    assert resp.status_code == 200
    assert resp.is_json

    score_variants = [x["variant"] for x in resp.json["data"]["variants"]]

    assert "summaryStatDataset" in resp.json["data"]

    for group in resp.json["data"]["groups"]:
        n_variants = len(group["variants"])
        n_covar = len(group["covariance"])
        assert n_variants > 0
        assert n_covar > 0
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

def test_summary_stat_chromglob(client):
    resp = client.post("/aggregation/covariance", data = {
        "chrom": "9",
        "start": 22132,
        "stop": 22142,
        "summaryStatDataset": 3,
        "genomeBuild": "GRCh37",
        "masks": [4]
    })

    assert resp.status_code == 200
    assert resp.is_json
    score_variants = [x["variant"] for x in resp.json["data"]["variants"]]
    assert "summaryStatDataset" in resp.json["data"]

    groups = resp.json["data"]["groups"]
    assert len(groups) > 0
    for group in resp.json["data"]["groups"]:
        n_variants = len(group["variants"])
        n_covar = len(group["covariance"])
        assert n_variants > 0
        assert n_covar > 0
        assert n_covar == (n_variants * (n_variants + 1) / 2)
        assert isinstance(group["mask"], int)
        assert "group" in group
        assert "groupType" in group
        assert group["groupType"] in ("REGION", "GENE")
        assert group["variants"][0].split(":")[0] == '9'
        assert all([v in score_variants for v in group["variants"]])

    assert groups[0]["group"] == "MODX3"
    assert groups[0]["covariance"][0] == pytest.approx(0.0142854)

    variants = resp.json["data"]["variants"]
    assert len(variants) > 0
    for variant in variants:
        assert variant["altFreq"] > 0
        assert variant["pvalue"] > 0
        assert variant["pvalue"] <= 1
        assert "score" in variant

    assert variants[0]["score"] == pytest.approx(-9.93863)
    assert variants[0]["pvalue"] == pytest.approx(0.00881654)

    resp = client.post("/aggregation/covariance", data = {
        "chrom": "1",
        "start": 14895,
        "stop": 14918,
        "summaryStatDataset": 3,
        "genomeBuild": "GRCh37",
        "masks": [4]
    })

    assert resp.status_code == 200
    assert resp.is_json
    score_variants = [x["variant"] for x in resp.json["data"]["variants"]]
    assert "summaryStatDataset" in resp.json["data"]

    groups = resp.json["data"]["groups"]
    assert len(groups) > 0
    for group in resp.json["data"]["groups"]:
        n_variants = len(group["variants"])
        n_covar = len(group["covariance"])
        assert n_variants > 0
        assert n_covar > 0
        assert n_covar == (n_variants * (n_variants + 1) / 2)
        assert isinstance(group["mask"], int)
        assert "group" in group
        assert "groupType" in group
        assert group["groupType"] in ("REGION", "GENE")
        assert group["variants"][0].split(":")[0] == '1'
        assert all([v in score_variants for v in group["variants"]])

    assert groups[0]["group"] == "IDQV5"
    assert groups[0]["covariance"][0] == pytest.approx(0.00140996)

    variants = resp.json["data"]["variants"]
    assert len(variants) > 0
    for variant in variants:
        assert variant["altFreq"] > 0
        assert variant["pvalue"] > 0
        assert variant["pvalue"] <= 1
        assert "score" in variant

    assert variants[0]["score"] == pytest.approx(-1.31314)
    assert variants[0]["pvalue"] == pytest.approx(0.270687)

def test_pheno_bad_float(client):
    resp = client.post("/aggregation/covariance", data = {
        "chrom": "22",
        "start": 50276998,
        "stop": 50357719,
        "genotypeDataset": 1,
        "phenotypeDataset": 3,
        "phenotype": "rand_qt",
        "samples": "ALL",
        "genomeBuild": "GRCh37",
        "masks": [1]
    })

    assert resp.status_code == 500
    assert "An error occurred parsing a phenotype file" in resp.json["error"]

def test_for_analysis_skip_column(client):
    resp = client.post("/aggregation/covariance", data = {
        "chrom": "22",
        "start": 50276998,
        "stop": 50357719,
        "genotypeDataset": 1,
        "phenotypeDataset": 4,
        "phenotype": "rand_qt",
        "samples": "ALL",
        "genomeBuild": "GRCh37",
        "masks": [1]
    })

    assert resp.status_code == 200

def test_user_masks(client):
    resp = client.post("/aggregation/covariance", json = {
        "chrom": "22",
        "start": 50276998,
        "stop": 50357719,
        "genotypeDataset": 1,
        "phenotypeDataset": 1,
        "phenotype": "rand_qt",
        "samples": "ALL",
        "genomeBuild": "GRCh37",
        "maskDefinitions": [
            {
                "id": 10,
                "name": "PTV+LOF<0.01",
                "description": "A mask generated by the user in the browser, including only variants that are PTV or LOF w/ AF < 1%",
                "genome_build": "GRCh37",
                "group_type": "GENE",
                "identifier_type": "ENSEMBL",
                "groups": {
                    "PIM3": ["22:50354416_G/C", "22:50355407_C/T", "22:50356368_C/T", "22:50356386_C/T", "22:50356473_C/T",
                             "22:50356497_G/A", "22:50356731_C/T", "22:50356811_G/T", "22:50356864_G/A", "22:50356875_C/T",
                             "22:50356887_C/T", "22:50356961_C/T", "22:50356965_C/T", "22:50356994_G/A", "22:50357305_C/T",
                             "22:50357350_G/A", "22:50357577_G/A", "22:50357657_A/G", "22:50357667_C/G"],
                }
            }
        ]
    })

    assert resp.status_code == 200
    assert resp.is_json

    score_variants = [x["variant"] for x in resp.json["data"]["variants"]]

    assert "phenotypeDataset" in resp.json["data"]
    assert "genotypeDataset" in resp.json["data"]
    assert "phenotype" in resp.json["data"]

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