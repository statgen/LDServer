from raremetal.model import add_yaml_command, get_genotype_dataset, get_phenotype_dataset, \
                            get_phenotype_column_objects, get_mask_by_id, get_analysis_columns
from core.pywrapper import VariantGroupType, GroupIdentifierType
import traceback

def test_add_yaml(app, db):
  with app.app_context():
    db.drop_all()
    db.create_all()

    runner = app.test_cli_runner()
    result = runner.invoke(add_yaml_command, ["../data/test.yaml"])
    if isinstance(result.exception, Exception):
      traceback.print_tb(result.exc_info[2])
      raise result.exception

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

    analysis_cols = get_analysis_columns(1)
    assert len(analysis_cols) == 2

def test_tab_incorrect_float(app, db):
  with app.app_context():
    db.create_all()
    runner = app.test_cli_runner()
    result = runner.invoke(add_yaml_command, ["../data/test_tab_incorrect_float.yaml"])

    assert isinstance(result.exception, ValueError)
    assert result.exception.message.startswith("Column sex for phenotype dataset")
    assert result.exception.message.startswith("Column sex for phenotype dataset")
    assert result.exception.message.endswith("cannot be coerced to float")

def test_ped_incorrect_float(app, db):
  with app.app_context():
    db.create_all()
    runner = app.test_cli_runner()
    result = runner.invoke(add_yaml_command, ["../data/test_ped_incorrect_float.yaml"])

    assert isinstance(result.exception, ValueError)
    assert result.exception.message.startswith("Column ANCESTRY for PED")
    assert result.exception.message.endswith("cannot be coerced to float")

def test_vcf_missing_tabix(app, db):
  with app.app_context():
    db.create_all()
    runner = app.test_cli_runner()
    result = runner.invoke(add_yaml_command, ["../data/test_no_tabix.yaml"])

    assert isinstance(result.exception, ValueError)
    assert result.exception.message.startswith("Cannot find tabix index for VCF file")

def test_sav_missing_index(app, db):
  with app.app_context():
    db.create_all()
    runner = app.test_cli_runner()
    result = runner.invoke(add_yaml_command, ["../data/test_no_sav_index.yaml"])

    assert isinstance(result.exception, ValueError)
    assert result.exception.message.startswith("Cannot find savvy index")