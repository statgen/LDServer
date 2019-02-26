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

    assert resp.status_code == 500