import re
import pytest
from copy import deepcopy

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

    assert resp.status_code == 200
    data = resp.json["data"]
    assert len(data["variants"]) == 0
    assert len(data["groups"]) == 0

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

    assert resp.status_code == 200
    data = resp.json["data"]
    assert len(data["variants"]) == 0
    assert len(data["groups"]) == 0
    assert "[tabix++] error reading file" not in resp.json.get("error", "")

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

def test_metastaar(client):
    resp = client.post("/aggregation/covariance", json = {
        "chrom": "1",
        "start": 4957,
        "stop": 5143,
        "summaryStatisticDataset": 4,
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
                    "WYAV7": [
                        "1:4957_C/T", "1:4958_G/A", "1:4959_C/T", "1:4960_G/C", "1:4961_C/G", "1:4962_C/T", "1:4963_G/A",
                        "1:4964_G/C", "1:4965_T/G", "1:4966_G/C", "1:4967_C/A", "1:4968_G/T", "1:4969_T/A", "1:4970_A/G",
                        "1:4971_T/C", "1:4972_C/A", "1:4973_T/C", "1:4974_T/A", "1:4975_C/T", "1:4976_C/T", "1:4977_G/A",
                        "1:4978_T/A", "1:4979_G/A", "1:4980_G/T", "1:4981_T/A", "1:4982_A/G", "1:4983_C/A", "1:4984_G/T",
                        "1:4985_C/G", "1:4986_C/G", "1:4987_G/C", "1:4988_C/A", "1:4989_A/C", "1:4990_A/C", "1:4991_G/T",
                        "1:4992_T/C", "1:4993_T/C", "1:4994_C/T", "1:4995_C/A", "1:4996_G/T", "1:4997_A/G", "1:4998_G/A",
                        "1:4999_C/G", "1:5000_T/C", "1:5001_C/G", "1:5002_G/C", "1:5003_G/T", "1:5004_T/G", "1:5005_C/A",
                        "1:5006_C/A", "1:5007_T/A", "1:5008_C/T", "1:5009_G/C", "1:5010_C/G", "1:5011_G/A", "1:5012_A/G",
                        "1:5013_A/G", "1:5014_A/C", "1:5015_T/G", "1:5016_A/T", "1:5017_C/A", "1:5018_A/C", "1:5019_A/T",
                        "1:5020_T/G", "1:5021_C/G", "1:5022_C/A", "1:5023_G/T", "1:5024_G/A", "1:5025_C/G", "1:5026_A/G",
                        "1:5027_C/A", "1:5028_C/T", "1:5029_T/C", "1:5030_G/A", "1:5031_A/C", "1:5032_C/G", "1:5033_A/G",
                        "1:5034_T/G", "1:5035_T/A", "1:5036_A/C", "1:5037_G/C", "1:5038_A/T", "1:5039_T/C", "1:5040_A/G",
                        "1:5041_G/A", "1:5042_T/G", "1:5043_T/C", "1:5044_C/G", "1:5045_T/G", "1:5046_A/G", "1:5047_A/T",
                        "1:5048_A/C", "1:5049_C/A", "1:5050_T/G", "1:5051_T/C", "1:5052_G/T", "1:5053_A/T", "1:5054_C/A",
                        "1:5055_G/T", "1:5056_T/A", "1:5057_T/G", "1:5058_T/C", "1:5059_G/A", "1:5060_A/T", "1:5061_A/C",
                        "1:5062_C/A", "1:5063_C/A", "1:5064_G/A", "1:5065_T/A", "1:5066_A/C", "1:5067_T/C", "1:5068_T/A",
                        "1:5069_T/A", "1:5070_C/G", "1:5071_C/G", "1:5072_G/A", "1:5073_C/A", "1:5074_C/A", "1:5075_A/G",
                        "1:5076_A/T", "1:5077_T/A", "1:5078_T/G", "1:5079_G/A", "1:5080_A/G", "1:5081_A/G", "1:5082_T/C",
                        "1:5083_C/G", "1:5084_C/T", "1:5085_G/C", "1:5086_A/G", "1:5087_A/G", "1:5088_G/A", "1:5089_T/G",
                        "1:5090_C/A", "1:5091_G/A", "1:5092_A/G", "1:5093_C/T", "1:5094_C/A", "1:5095_G/C", "1:5096_A/G",
                        "1:5097_C/A", "1:5098_A/G", "1:5099_G/T", "1:5100_G/A", "1:5101_G/T", "1:5102_A/C", "1:5103_T/A",
                        "1:5104_C/G", "1:5105_C/A", "1:5106_G/A", "1:5107_T/C", "1:5108_A/C", "1:5109_T/A", "1:5110_A/G",
                        "1:5111_C/T", "1:5112_C/A", "1:5113_C/A", "1:5114_C/G", "1:5115_C/T", "1:5116_G/A", "1:5117_G/A",
                        "1:5118_C/G", "1:5119_G/A", "1:5120_G/A", "1:5121_C/T", "1:5122_G/A", "1:5123_C/T", "1:5124_G/A",
                        "1:5125_T/A", "1:5126_G/C", "1:5127_C/T", "1:5128_A/T", "1:5129_A/C", "1:5130_A/G", "1:5131_T/G",
                        "1:5132_C/G", "1:5133_T/A", "1:5134_A/T", "1:5135_C/T", "1:5136_G/C", "1:5137_G/A", "1:5138_G/T",
                        "1:5139_T/A", "1:5140_C/T", "1:5141_C/A", "1:5142_T/A", "1:5143_A/G"
                    ],
                }
            }
        ]
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

# This test contains a region bound that technically falls within the MetaSTAAR computed bounds,
# but there are no actual variants in the region itself. So the query should return empty.
def test_metastaar_empty_region(client):
    resp = client.post("/aggregation/covariance", json = {
        "chrom": "1",
        "start": 9205,
        "stop": 9210,
        "summaryStatisticDataset": 4,
        "genomeBuild": "GRCh37",
        "maskDefinitions": [
            {
                "id": 10,
                "name": "PTV+LOF<0.01",
                "description": "A mask generated by the user in the browser, including only variants that are PTV or LOF w/ AF < 1%",
                "genome_build": "GRCh37",
                "group_type": "REGION",
                "identifier_type": "COORDINATES",
                "groups": {
                    "bad-region": {
                        "start": 9205,
                        "stop":  9210,
                    }
                }
            }
        ]
    })

    assert resp.status_code == 200
    assert resp.is_json

    data = resp.json["data"]

    assert "summaryStatisticDataset" in data
    assert "groups" in data
    assert "nSamples" in data
    assert len(data["variants"]) == 0
    assert len(data["groups"]) == 0

# This test requests variants in a region that is completely outside what MetaSTAAR tried to compute within this chunk
# In other words, the metadata in the parquet file specifies a range say 5,000 - 10,000 and we ask for 90,000.
def test_metastaar_region_out_of_bounds(client):
    resp = client.post("/aggregation/covariance", json = {
        "chrom": "1",
        "start": 99205,
        "stop": 99210,
        "summaryStatisticDataset": 4,
        "genomeBuild": "GRCh37",
        "maskDefinitions": [
            {
                "id": 10,
                "name": "PTV+LOF<0.01",
                "description": "A mask generated by the user in the browser, including only variants that are PTV or LOF w/ AF < 1%",
                "genome_build": "GRCh37",
                "group_type": "REGION",
                "identifier_type": "COORDINATES",
                "groups": {
                    "bad-region": {
                        "start": 99205,
                        "stop":  99210,
                    }
                }
            }
        ]
    })

    assert resp.status_code == 200
    assert resp.is_json

    data = resp.json["data"]

    assert "summaryStatisticDataset" in data
    assert "groups" in data
    assert "nSamples" in data
    assert len(data["variants"]) == 0
    assert len(data["groups"]) == 0

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

def test_bad_group_or_identifier_types(client):
    json_data = {
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
                "id": 1,
                "name": "Testing regions",
                "description": "Random region picked for testing",
                "genome_build": "GRCh37",
                "group_type": "REGION",
                "identifier_type": "COORDINATES",
                "groups": {
                    "22:50276998-50300000": {
                        "start": 50276998,
                        "stop":  50300000,
                        "filters": [
                            {
                                "field": "maf",
                                "op": "gte",
                                "value": 0.05
                            }
                        ]
                    }
                }
            }
        ]
    }

    test_data1 = deepcopy(json_data)
    test_data1["maskDefinitions"][0]["group_type"] = "region"
    resp1 = client.post("/aggregation/covariance", json=test_data1)

    assert resp1.status_code == 400
    assert resp1.is_json

    test_data2 = deepcopy(json_data)
    test_data2["maskDefinitions"][0]["identifier_type"] = "ensembl"
    resp2 = client.post("/aggregation/covariance", json=test_data2)

    assert resp2.status_code == 400
    assert resp2.is_json

def test_region_covar(client):
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
                "id": 1,
                "name": "Testing regions",
                "description": "Random region picked for testing",
                "genome_build": "GRCh37",
                "group_type": "REGION",
                "identifier_type": "COORDINATES",
                "groups": {
                    "22:50276998-50300000": {
                        "start": 50276998,
                        "stop":  50300000,
                        "filters": [
                            {
                                "field": "maf",
                                "op": "gte",
                                "value": 0.05
                            }
                        ]
                    }
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

def test_region_filter_maf(client):
    resp_filter = client.post("/aggregation/covariance", json = {
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
                "id": 1,
                "name": "Testing regions",
                "description": "Random region picked for testing",
                "genome_build": "GRCh37",
                "group_type": "REGION",
                "identifier_type": "COORDINATES",
                "groups": {
                    "22:50276998-50300000": {
                        "start": 50276998,
                        "stop":  50300000,
                        "filters": [
                            {
                                "field": "maf",
                                "op": "gte",
                                "value": 0.05
                            }
                        ]
                    }
                }
            }
        ]
    })

    resp_nofilter = client.post("/aggregation/covariance", json = {
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
                "id": 1,
                "name": "Testing regions",
                "description": "Random region picked for testing",
                "genome_build": "GRCh37",
                "group_type": "REGION",
                "identifier_type": "COORDINATES",
                "groups": {
                    "22:50276998-50300000": {
                        "start": 50276998,
                        "stop":  50300000
                    }
                }
            }
        ]
    })

    group_filt = resp_filter.json["data"]["groups"][0]
    group_nofilt = resp_nofilter.json["data"]["groups"][0]

    assert len(group_filt["variants"]) < len(group_nofilt["variants"])
    assert len(group_filt["covariance"]) < len(group_nofilt["covariance"])

    def ncovar(n):
        return n * (n + 1) / 2

    assert ncovar(len(group_filt["variants"])) == len(group_filt["covariance"])
    assert ncovar(len(group_nofilt["variants"])) == len(group_nofilt["covariance"])

def test_region_covar_too_big(client):
    resp = client.post("/aggregation/covariance", json = {
        "chrom": "22",
        "start": 50276998,
        "stop": 51300000,
        "genotypeDataset": 1,
        "phenotypeDataset": 1,
        "phenotype": "rand_qt",
        "samples": "ALL",
        "genomeBuild": "GRCh37",
        "maskDefinitions": [
            {
                "id": 1,
                "name": "Testing regions",
                "description": "Random region picked for testing",
                "genome_build": "GRCh37",
                "group_type": "REGION",
                "identifier_type": "COORDINATES",
                "groups": {
                    "22:50276998-51300000": {
                        "start": 50276998,
                        "stop":  51300000,
                        "filters": [
                            {
                                "field": "maf",
                                "op": "gte",
                                "value": 0.05
                            }
                        ]
                    }
                }
            }
        ]
    })

    assert resp.status_code == 400

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

def check_variant_response_ok(resp):
    assert resp.status_code == 200
    assert resp.is_json
    assert len(resp.json["data"]["variants"]) > 0
    for v in resp.json["data"]["variants"]:
        assert "variant" in v
        assert "chrom" in v
        assert "pos" in v
        assert isinstance(v["chrom"], str)
        assert isinstance(v["pos"], int)

def check_variant_response_empty(resp):
    assert resp.status_code == 200
    assert resp.is_json
    assert len(resp.json["data"]["variants"]) == 0
    data = resp.json["data"]
    assert "summaryStatisticDataset" in data or "genotypeDataset" in data

def test_get_variants_from_genotype_file(client):
    resp = client.post("/aggregation/variants", json={
        "chrom": "22",
        "start": 50276998,
        "stop": 50357719,
        "genotypeDataset": 2,
        "genomeBuild": "GRCh37"
    })
    check_variant_response_ok(resp)

def test_get_variants_from_summary_stats(client):
    # Test a rvtest formatted dataset
    resp = client.post("/aggregation/variants", json={
        "chrom": "1",
        "start": 2,
        "stop": 307,
        "summaryStatisticDataset": 2,
        "genomeBuild": "GRCh37"
    })
    check_variant_response_ok(resp)

    # Test a MetaSTAAR dataset
    resp = client.post("/aggregation/variants", json={
        "chrom": "1",
        "start": 4957,
        "stop": 5143,
        "summaryStatisticDataset": 4,
        "genomeBuild": "GRCh37"
    })
    check_variant_response_ok(resp)

    # Test MetaSTAAR potential failure 1
    resp = client.post("/aggregation/variants", json={
        "chrom": "1",
        "start": 9205,
        "stop": 9210,
        "summaryStatisticDataset": 4,
        "genomeBuild": "GRCh37"
    })
    check_variant_response_empty(resp)

    # Test MetaSTAAR potential failure 2
    resp = client.post("/aggregation/variants", json={
        "chrom": "1",
        "start": 99205,
        "stop": 99210,
        "summaryStatisticDataset": 4,
        "genomeBuild": "GRCh37"
    })
    check_variant_response_empty(resp)