import re
import pytest

def is_sha(s):
    return re.search("[A-Fa-f0-9]+",s) is not None

def parse_epacts(v, strict=True):
    """
    Try to parse an EPACTS ID into components.

    Args:
      v (string): variant
      strict (bool): if true, must match an EPACTS ID exactly (chr:pos_ref/alt)
        If false, then ref/alt can be missing, but at least chr:pos must be specified.
        In this case, ref/alt will be None in the returned tuple

    Returns:
      tuple: (chrom, pos, ref, alt)
    """

    split = v.split("_")

    # Split chrom/pos
    # This is the minimum required information. If even this isn't present,
    # it's a bogus ID.
    try:
        chrom, pos = split[0].split(":")
    except:
        raise ValueError("EPACTS ID had no chrom/pos? " + v)

    # Position should be numeric
    try:
        int(pos)
    except:
        raise ValueError("Couldn't recognize position {} for variant {}".format(pos,v))

    # Try to split alleles if they were given
    try:
        ref, alt = split[1].split("/")
    except:
        if strict:
            raise ValueError("No ref/alt alleles found in EPACTS ID " + v)

        ref = None
        alt = None

    return chrom, pos, ref, alt

def parse_colons(v):
    """
    Try to parse a colon separated ID into components.

    Args:
      v (string): variant

    Returns:
      tuple: (chrom, pos, ref, alt)
    """

    chrom, pos, ref, alt = v.split(":")

    # Position should be numeric
    try:
        int(pos)
    except:
        raise ValueError("Couldn't recognize position {} for variant {}".format(pos, v))

    return chrom, pos, ref, alt

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

def test_no_scores_in_region(client):
    resp = client.post("/aggregation/covariance", json = {
        "chrom": "22",
        "start": 1,
        "stop": 2,
        "summaryStatisticDataset": 2,
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
                    "PIM3": ["1:1_A/T"],
                }
            }
        ]
    })

    assert resp.status_code == 400
    assert re.search("No score statistics loaded within genomic region", resp.json["error"]) is not None

def test_tabixpp_file_inaccessible(client):
    resp = client.post("/aggregation/covariance", json = {
        "chrom": "22",
        "start": 1,
        "stop": 2,
        "summaryStatisticDataset": 1,
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
                    "PIM3": ["1:1_A/T"],
                }
            }
        ]
    })

    assert resp.status_code == 400
    assert "[tabix++] error reading file" not in resp.json["error"]

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
        "summaryStatisticDataset": 2,
        "genomeBuild": "GRCh37",
        "masks": [3]
    })

    assert resp.status_code == 200
    assert resp.is_json

    score_variants = [x["variant"] for x in resp.json["data"]["variants"]]

    assert "summaryStatisticDataset" in resp.json["data"]

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

def test_variant_format(client):
    byformat = dict()
    formats = ("EPACTS", "COLONS") # EPACTS must be the first format, add others to the end

    for vfmt in formats:
        resp = client.post("/aggregation/covariance", data = {
            "chrom": "1",
            "start": 2,
            "stop": 307,
            "summaryStatisticDataset": 2,
            "genomeBuild": "GRCh37",
            "masks": [3],
            "variantFormat": vfmt
        })

        assert resp.status_code == 200
        assert resp.is_json

        byformat[vfmt] = resp.json["data"]

    epacts_variants = [parse_epacts(x["variant"]) for x in byformat["EPACTS"]["variants"]]

    for vfmt in formats[1:]:
        score_variants = [x["variant"] for x in byformat[vfmt]["variants"]]
        assert len(score_variants) > 0

        for i, v in enumerate(score_variants):
            if vfmt == "COLONS":
                chrom, pos, ref, alt = parse_colons(v)

            epacts_chrom, epacts_pos, epacts_ref, epacts_alt = epacts_variants[i]

            assert chrom == epacts_chrom
            assert pos == epacts_pos
            assert ref == epacts_ref
            assert alt == epacts_alt

        for g, group in enumerate(byformat[vfmt]["groups"]):
            n_variants = len(group["variants"])
            assert n_variants > 0

            for i, v in enumerate(group["variants"]):
                if vfmt == "COLONS":
                    chrom, pos, ref, alt = parse_colons(v)

                epacts_variant = byformat["EPACTS"]["groups"][g]["variants"][i]
                epacts_chrom, epacts_pos, epacts_ref, epacts_alt = parse_epacts(epacts_variant)

                assert chrom == epacts_chrom
                assert pos == epacts_pos
                assert ref == epacts_ref
                assert alt == epacts_alt

def test_variant_format_invalid(client):
    resp = client.post("/aggregation/covariance", data = {
        "chrom": "1",
        "start": 2,
        "stop": 307,
        "summaryStatisticDataset": 2,
        "genomeBuild": "GRCh37",
        "masks": [3],
        "variantFormat": "NOT_A_FORMAT"
    })

    assert resp.status_code == 400

def test_summary_stat_chromglob(client):
    resp = client.post("/aggregation/covariance", data = {
        "chrom": "9",
        "start": 22132,
        "stop": 22142,
        "summaryStatisticDataset": 3,
        "genomeBuild": "GRCh37",
        "masks": [4]
    })

    assert resp.status_code == 200
    assert resp.is_json
    score_variants = [x["variant"] for x in resp.json["data"]["variants"]]
    assert "summaryStatisticDataset" in resp.json["data"]

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
        "summaryStatisticDataset": 3,
        "genomeBuild": "GRCh37",
        "masks": [4]
    })

    assert resp.status_code == 200
    assert resp.is_json
    score_variants = [x["variant"] for x in resp.json["data"]["variants"]]
    assert "summaryStatisticDataset" in resp.json["data"]

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

    assert resp.status_code == 400
    assert "Error reading phenotype file on line" in resp.json["error"]

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

def test_user_masks_colon_format(client):
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
                    "PIM3": ["22:50354416:G:C", "22:50355407:C:T", "22:50356368:C:T", "22:50356386:C:T", "22:50356473:C:T",
                             "22:50356497:G:A", "22:50356731:C:T", "22:50356811:G:T", "22:50356864:G:A", "22:50356875:C:T",
                             "22:50356887:C:T", "22:50356961:C:T", "22:50356965:C:T", "22:50356994:G:A", "22:50357305:C:T",
                             "22:50357350:G:A", "22:50357577:G:A", "22:50357657:A:G", "22:50357667:C:G"],
                }
            }
        ],
        "variantFormat": "COLONS"
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