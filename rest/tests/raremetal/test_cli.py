from raremetal.model import add_yaml_command, get_genotype_dataset, get_phenotype_dataset, \
                            get_phenotype_column_objects, get_mask_by_id
from ld.pywrapper import VariantGroupType, GroupIdentifierType

def test_add_yaml(app, db):
  with app.app_context():
    db.drop_all()
    db.create_all()

    runner = app.test_cli_runner()
    result = runner.invoke(add_yaml_command, ["../data/test.yaml"])
    gdata = get_genotype_dataset(1)
    assert gdata["name"] == "1000G"
    assert gdata["genome_build"] == "GRCh37"

    pdata = get_phenotype_dataset(1)
    assert pdata["name"] == "1000G random phenotypes"
    assert pdata["description"] == "An example set of randomly generated phenotypes for 1000G"
    assert pdata["sample_column"] == "iid"
    assert pdata["nrows"] == 2504
    assert pdata["ncols"] == 4

    pcols = get_phenotype_column_objects(1)
    for col in pcols:
      assert col["column_type"] in ("TEXT", "CATEGORICAL", "FLOAT", "INTEGER")

    mdata = get_mask_by_id(1)
    assert mdata["name"] == "AF < 0.01"
    assert mdata["genome_build"] == "GRCh37"
    assert mdata["group_type"] == VariantGroupType.GENE
    assert mdata["identifier_type"] == GroupIdentifierType.ENSEMBL

    assert get_phenotype_dataset(2)["sample_column"] == "IID"