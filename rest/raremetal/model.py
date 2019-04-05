from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Index, inspect
from sqlalchemy.types import Enum
from flask.cli import with_appcontext
from flask import current_app
from collections import Counter, OrderedDict
from core.pywrapper import ColumnType, ColumnTypeMap, VariantGroupType, GroupIdentifierType, extract_samples
from tabulate import tabulate
from glob import glob
import os
import click
import json
import yaml
import gzip

db = SQLAlchemy()

class Correlation(db.Model):
  __tablename__ = 'correlation'
  id = db.Column(db.Integer, primary_key = True)
  name = db.Column(db.String(), unique = True, nullable = False)
  label = db.Column(db.String(), unique = False, nullable = False)
  descripition = db.Column(db.String(), unique = False, nullable = False)
  type = db.Column(db.String(), unique = False, nullable = False)
  def __repr__(self):
    return '<Correlation %r>' % self.name


class GenotypeDataset(db.Model):
  __tablename__ = 'genotype_datasets'
  id = db.Column(db.Integer, primary_key = True)
  name = db.Column(db.String(), unique = False, nullable = False)
  genome_build = db.Column(db.String(), unique = False, nullable = False)
  description = db.Column(db.String(), unique = False, nullable = False)

  files = db.relationship('File', backref = 'genotype_datasets', lazy = True)
  samples = db.relationship('Sample', backref = 'genotype_datasets', lazy = True)
  phenotypes = db.relationship("PhenotypeDataset", secondary="genotype_phenotype", back_populates="genotypes")
  masks = db.relationship("Mask", back_populates="genotypes")

  __table_args__ = (
    db.UniqueConstraint('name', 'genome_build', name = 'name_build_uc'),
  )

  def __repr__(self):
    return '<GenotypeDataset %r>' % self.name


class File(db.Model):
  __tablename__ = 'file'
  id = db.Column(db.Integer, primary_key = True)
  path = db.Column(db.String(), unique = False, nullable = False)
  genotype_dataset_id = db.Column(db.Integer, db.ForeignKey('genotype_datasets.id'), nullable = False)
  def __repr__(self):
    return '<File %r>' % self.path


class Sample(db.Model):
  __tablename__ = 'sample'
  id = db.Column(db.Integer, primary_key = True)
  subset = db.Column(db.String, unique = False, nullable = False)
  sample = db.Column(db.String, unique = False, nullable = False)
  genotype_dataset_id = db.Column(db.Integer, db.ForeignKey('genotype_datasets.id'), nullable = False)
  def __repr__(self):
    return '<Sample %r %r>' % (self.subset, self.sample)


genotype_phenotype = db.Table("genotype_phenotype",
  db.Column("genotype_dataset_id", db.Integer, db.ForeignKey("genotype_datasets.id"), primary_key = True),
  db.Column("phenotype_dataset_id", db.Integer, db.ForeignKey("phenotype_datasets.id"), primary_key = True)
)


class PhenotypeDataset(db.Model):
  __tablename__ = "phenotype_datasets"
  id = db.Column(db.Integer, primary_key = True)
  name = db.Column(db.String, unique = False, nullable = False)
  description = db.Column(db.String, unique = False, nullable = False)
  filepath = db.Column(db.String, unique = False, nullable = False)
  nrows = db.Column(db.Integer, nullable = False)
  ncols = db.Column(db.Integer, nullable = False)
  sample_column = db.Column(db.String, unique = False, nullable = False)
  delim = db.Column(db.String, unique = False, nullable = False)

  columns = db.relationship("PhenotypeColumn", back_populates = "dataset")
  genotypes = db.relationship("GenotypeDataset", secondary="genotype_phenotype", back_populates="phenotypes")


class PhenotypeColumn(db.Model):
  __tablename__ = "phenotype_columns"
  id = db.Column(db.Integer, primary_key = True)
  phenotype_dataset_id = db.Column(db.Integer, db.ForeignKey('phenotype_datasets.id'))
  column_name = db.Column(db.String, nullable = False)
  column_type = db.Column(Enum(*tuple(ColumnType.values[i].name for i in range(len(ColumnType.values)))), nullable = False, name = "column_type")
  description = db.Column(db.String, nullable = True)
  for_analysis = db.Column(db.Boolean, default = True)

  dataset = db.relationship("PhenotypeDataset", back_populates = "columns")


class Mask(db.Model):
  __tablename__ = "masks"
  id = db.Column(db.Integer, primary_key = True)
  name = db.Column(db.String, unique = True, nullable = False)
  filepath = db.Column(db.String, unique = True, nullable = False)
  description = db.Column(db.String, unique = False, nullable = False)
  genome_build = db.Column(db.String, unique = False, nullable = False)
  group_type = db.Column(Enum(*tuple(VariantGroupType.values[i].name for i in range(len(VariantGroupType.values)))), nullable = False, name = "group_type")
  identifier_type = db.Column(db.String, unique = False, nullable = False)
  genotype_dataset_id = db.Column(db.Integer, db.ForeignKey('genotype_datasets.id'))

  genotypes = db.relationship("GenotypeDataset", back_populates="masks")


Index('genotype_dataset_index_1', GenotypeDataset.genome_build)
Index('file_index', File.genotype_dataset_id, File.path)
Index('sample_index', Sample.genotype_dataset_id, Sample.subset, Sample.sample)

def find_file(relpath):
  if os.path.isfile(relpath):
    return relpath

  data_path = os.path.join(current_app.root_path, "../../", relpath)

  if os.path.isfile(data_path):
    return data_path
  else:
    raise IOError("Could not locate file, tried '{}' and '{}'".format(relpath, data_path))

def get_correlations():
  return [{'name': name, 'label': label, 'description': desc, 'type': type} for name, label, desc, type in db.session.query(Correlation.name, Correlation.label, Correlation.descripition, Correlation.type)]

def has_genome_build(genome_build):
  return db.session.query(GenotypeDataset.genome_build).filter_by(genome_build = genome_build).first() is not None

def get_genome_builds():
  return [x for x, in db.session.query(GenotypeDataset.genome_build).distinct()]

def get_genotype_datasets(genome_build):
  return [{'name': name, 'description': desc } for name, desc, in db.session.query(GenotypeDataset.name, GenotypeDataset.description).filter_by(genome_build = genome_build)]

def get_genotype_dataset_id(genome_build, genotype_dataset_name):
  return db.session.query(GenotypeDataset.id).filter_by(genome_build = genome_build).filter_by(name = genotype_dataset_name).scalar()

def get_genotype_dataset(genotype_dataset_id):
  genotype_dataset = GenotypeDataset.query.filter_by(id = genotype_dataset_id).first()
  if genotype_dataset is not None:
    return {c.key: getattr(genotype_dataset, c.key) for c in inspect(genotype_dataset).mapper.column_attrs}

  return None

def get_phenotype_dataset(phenotype_dataset_id):
  phenotype_dataset = PhenotypeDataset.query.filter_by(id = phenotype_dataset_id).first()
  if phenotype_dataset is not None:
    return {c.key: getattr(phenotype_dataset, c.key) for c in inspect(phenotype_dataset).mapper.column_attrs}

  return None

def get_full_genotype_datasets():
  datasets = []
  for row in db.session.query(GenotypeDataset).all():
    as_dict = {c.key: getattr(row, c.key) for c in inspect(row).mapper.column_attrs}
    as_dict["genotypeDataset"] = as_dict["id"]
    as_dict["genomeBuild"] = as_dict["genome_build"]
    del as_dict["genome_build"]
    del as_dict["id"]
    datasets.append(as_dict)

  return datasets

def has_genotype_dataset(genotype_dataset_id):
  return db.session.query(GenotypeDataset).filter_by(id = genotype_dataset_id).first() is not None

def has_phenotype_dataset(phenotype_dataset_id):
  return db.session.query(PhenotypeDataset.id).filter_by(id = phenotype_dataset_id).first() is not None

def has_phenotype(phenotype_dataset_id, phenotype):
  return db.session.query(PhenotypeColumn).filter_by(phenotype_dataset_id=phenotype_dataset_id, column_name=phenotype).scalar()

def has_mask(mask_id):
  return db.session.query(Mask).filter_by(id = mask_id).scalar() is not None

def get_phenotype_dataset_id(phenotype_dataset_name):
  return db.session.query(PhenotypeDataset.id).filter_by(name = phenotype_dataset_name).scalar()

def has_sample_subset(genotype_dataset_id, subset):
  return db.session.query(Sample.subset).filter_by(genotype_dataset_id = genotype_dataset_id).filter_by(subset = subset).first() is not None

def get_sample_subsets(genotype_dataset_id):
  return [x for x, in db.session.query(Sample.subset).filter_by(genotype_dataset_id = genotype_dataset_id).distinct()]

def has_samples(genotype_dataset_id, subset):
  return db.session.query(Sample.sample).filter_by(genotype_dataset_id = genotype_dataset_id).filter_by(subset = subset).first() is not None

def get_samples_count(genotype_dataset_id, subset):
  return db.session.query(Sample).filter_by(genotype_dataset_id = genotype_dataset_id).filter_by(subset = subset).count()

def get_samples(genotype_dataset_id, subset):
  return [str(x) for x, in db.session.query(Sample.sample).filter_by(genotype_dataset_id = genotype_dataset_id).filter_by(subset = subset)]

def get_files(genotype_dataset_id):
  return [str(x) for x, in db.session.query(File.path).filter_by(genotype_dataset_id = genotype_dataset_id)]

def get_phenotype_file(phenotype_dataset_id):
  p = db.session.query(PhenotypeDataset.filepath).filter_by(id = phenotype_dataset_id).scalar()
  return str(p) if p is not None else None

def get_column_types(phenotype_dataset_id):
  mapping = ColumnTypeMap()
  for row in PhenotypeColumn.query.filter_by(phenotype_dataset_id = phenotype_dataset_id):
    mapping.add(str(row.column_name), ColumnType.names.get(row.column_type))

  return mapping

# TODO: clean this up, should just get object all at once
def get_phenotype_nrows(phenotype_dataset_id):
  return db.session.query(PhenotypeDataset.nrows).filter_by(id = phenotype_dataset_id).scalar()

def get_phenotype_delim(phenotype_dataset_id):
  return db.session.query(PhenotypeDataset.delim).filter_by(id = phenotype_dataset_id).scalar()

def get_phenotype_sample_column(phenotype_dataset_id):
  return db.session.query(PhenotypeDataset.sample_column).filter_by(id = phenotype_dataset_id).scalar()

def get_phenotype_columns(phenotype_dataset_id):
  return [str(x) for x in db.session.query(PhenotypeColumn.column_name).filter_by(id = phenotype_dataset_id).order_by(PhenotypeColumn.id)]

def get_phenotype_column_objects(phenotype_dataset_id):
  results = db.session.query(PhenotypeColumn).filter_by(phenotype_dataset_id=phenotype_dataset_id)
  objs =  [{c.key: getattr(row, c.key) for c in inspect(row).mapper.column_attrs} for row in results]
  return objs

def get_phenotypes_for_genotypes(genotype_dataset_id):
  result = []
  for row in db.session.query(GenotypeDataset).filter_by(id = genotype_dataset_id).scalar().phenotypes:
    as_dict = {c.key: getattr(row, c.key) for c in inspect(row).mapper.column_attrs}
    as_dict["phenotypeDataset"] = as_dict["id"]
    del as_dict["id"]
    del as_dict["nrows"]
    del as_dict["ncols"]
    del as_dict["filepath"]
    del as_dict["sample_column"]
    del as_dict["delim"]

    for col in row.columns:
      if col.for_analysis:
        as_dict.setdefault("phenotypes", []).append({
          "name": col.column_name,
          "description": col.description
        })

    result.append(as_dict)

  return result

def get_masks_for_genotypes(genotype_dataset_id):
  masks = []
  for row in db.session.query(Mask).filter_by(genotype_dataset_id = genotype_dataset_id):
    as_dict = {c.key: getattr(row, c.key) for c in inspect(row).mapper.column_attrs}
    as_dict["groupType"] = as_dict.pop("group_type")
    as_dict["identifierType"] = as_dict.pop("identifier_type")

    del as_dict["genome_build"]
    del as_dict["genotype_dataset_id"]
    del as_dict["filepath"]
    masks.append(as_dict)

  return masks

def get_mask_by_id(mask_id):
  result = db.session.query(Mask).filter_by(id = mask_id).scalar()
  if result is None:
    raise ValueError("No mask exists for ID {}".format(mask_id))

  as_dict = {c.key: getattr(result, c.key) for c in inspect(result).mapper.column_attrs}

  as_dict["group_type"] = VariantGroupType.names.get(result.group_type)
  as_dict["identifier_type"] = GroupIdentifierType.names.get(result.identifier_type)
  return as_dict

def get_mask_by_name(mask_name, genotype_dataset_id):
  result = db.session.query(Mask).filter_by(name = mask_name, genotype_dataset_id = genotype_dataset_id).scalar()
  if result is None:
    raise ValueError("No mask exists by the name {} for genotype dataset {}".format(mask_name, genotype_dataset_id))

  as_dict = {c.key: getattr(result, c.key) for c in inspect(result).mapper.column_attrs}

  as_dict["group_type"] = VariantGroupType.names.get(result.group_type)
  as_dict["identifier_type"] = GroupIdentifierType.names.get(result.identifier_type)
  return as_dict

def load_correlations():
  db.create_all()
  correlations = [str(x) for x, in db.session.query(Correlation.name)]
  if not correlations:
    correlations = [{ 'name': 'r', 'label': 'r', 'description': '', 'type': 'LD' },
            { 'name': 'rsquare', 'label': 'r^2', 'description': '', 'type': 'LD' },
            { 'name': 'cov', 'label': 'covariance', 'description': '', 'type': 'Covariance' }]
    for correlation in correlations:
      db.session.add(Correlation(name = correlation['name'], label = correlation['label'], descripition = correlation['description'], type = correlation['type']))
    db.session.commit()
    return [x['name'] for x in correlations]
  return correlations


def load_genotype_datasets(json_file):
  for t in db.metadata.sorted_tables:
    if t.name != 'correlation':
      t.drop(db.engine)
  db.create_all()
  with open(json_file, 'r') as f:
    for genotype_dataset in json.load(f):
      r = GenotypeDataset(name = genotype_dataset['Name'], description = genotype_dataset['Description'], genome_build = genotype_dataset['Genome build'])
      for path in genotype_dataset['Files']:
        r.files.append(File(path = path))
      for subset, samples in genotype_dataset['Samples'].iteritems():
        for sample in samples:
          r.samples.append(Sample(subset = subset, sample = sample))
      db.session.add(r)
    db.session.commit()


def show_genotype_datasets():
  db.create_all()
  header = ["ID", "Name", "Description", "Genome Build", "Populations", "Files"]
  rows = []
  for genotype_dataset in GenotypeDataset.query.all():
    rows.append([
      genotype_dataset.id,
      genotype_dataset.name,
      genotype_dataset.description,
      genotype_dataset.genome_build,
      ';'.join(list(set([s.subset for s in genotype_dataset.samples]))),
      "\n".join([f.path for f in genotype_dataset.files])
    ])

  print(tabulate(rows, header, 'psql'))

def show_phenotype_datasets():
  db.create_all()
  header = ["ID", "Name", "Description", "Path"]
  rows = []
  for phenotype_dataset in PhenotypeDataset.query.all():
    rows.append([phenotype_dataset.id, phenotype_dataset.name, phenotype_dataset.description, phenotype_dataset.filepath])

  print(tabulate(rows, header, 'psql'))

def show_masks():
  db.create_all()
  header = ["ID", "Name", "Description", "Genome Build", "Group Type", "Identifier Type", "Path"]
  rows = []
  for mask in Mask.query.all():
    rows.append([
      mask.id,
      mask.name,
      mask.description,
      mask.genome_build,
      mask.group_type,
      mask.identifier_type,
      mask.filepath
    ])

  print(tabulate(rows, header, 'psql'))

def add_genotype_dataset(name, description, genome_build, samples_filename, genotype_files, dbid=None):
  db.create_all()
  samples = []
  if samples_filename:
    with open(samples_filename, 'r') as f:
      for sample in f:
        sample = sample.strip()
        if sample:
          samples.append(sample)
  else:
    for f in genotype_files:
      f = find_file(f)
      f_samples = extract_samples(f)

      if len(f_samples) == 0:
        raise ValueError("Extracted 0 samples from {} - is the file missing?".format(f))

      if len(samples) == 0:
        samples = f_samples
      else:
        if set(samples) != set(f_samples):
          raise ValueError("Genotype file {} did not have the same samples as the others")

  args = {"name": name, "description": description, "genome_build": genome_build}
  if dbid is not None:
    args["id"] = dbid

  r = GenotypeDataset(**args)
  for path in genotype_files:
    r.files.append(File(path = path))
  for sample in samples:
    r.samples.append(Sample(subset = 'ALL', sample = sample))
  db.session.add(r)
  db.session.commit()

def guess_type(values):
  guesses = []

  for v in values:
    try:
      f = float(v)
    except ValueError:
      guesses.append(ColumnType.TEXT.name)
      continue

    i = int(f)
    if f == i:
      guesses.append(ColumnType.INTEGER.name)
    else:
      guesses.append(ColumnType.FLOAT.name)

  return Counter(guesses).most_common(1)[0][0]


def _load_phenotype_ped(ped_file, dat_file):
  if ped_file.endswith(".gz"):
    fp_ped = gzip.open(ped_file, "rt")
  else:
    fp_ped = open(ped_file)

  if dat_file.endswith(".gz"):
    fp_dat = gzip.open(dat_file, "rt")
  else:
    fp_dat = open(dat_file)

  nrows = 0
  with fp_ped, fp_dat:
    column_types = OrderedDict()

    # We know some of the column types already in a PED file
    column_types["FID"] = {"column_type": ColumnType.TEXT.name}
    column_types["IID"] = {"column_type": ColumnType.TEXT.name}
    column_types["PID"] = {"column_type": ColumnType.TEXT.name}
    column_types["MID"] = {"column_type": ColumnType.TEXT.name}
    column_types["SEX"] = {"column_type": ColumnType.CATEGORICAL.name}

    for k, v in column_types.items():
      v["column_name"] = k
      v["for_analysis"] = False

    for i, line in enumerate(fp_dat):
      ctype, pheno = line.split()
      if ctype == "A":
        column_types.setdefault(pheno, {})["column_type"] = ColumnType.CATEGORICAL.name
        column_types.setdefault(pheno, {})["column_name"] = pheno
      elif ctype == "T":
        column_types.setdefault(pheno, {})["column_type"] = ColumnType.FLOAT.name
        column_types.setdefault(pheno, {})["column_name"] = pheno
      elif ctype == "M":
        continue
      else:
        ValueError("Unrecognized DAT data type: " + ctype)

    for line in fp_ped:
      if not line.startswith("#"):
        nrows += 1
        break

    nrows += sum(1 for i in fp_ped)

  return column_types, nrows


def _load_phenotype_tab(tab_file):
  if tab_file.endswith(".gz"):
    fp_tab = gzip.open(tab_file, "rt")
  else:
    fp_tab = open(tab_file)

  header = None
  nrows = 0
  with fp_tab:
    # First line is header
    header = next(fp_tab).split("\t")
    header[-1] = header[-1].rstrip()

    column_values = OrderedDict()
    for i, line in enumerate(fp_tab):
      ls = line.split("\t")
      ls[-1] = ls[-1].rstrip()

      for j, v in enumerate(ls):
        column_values.setdefault(header[j], []).append(v)

      nrows += 1

  column_types = OrderedDict()
  for k, v in column_values.items():
    column_types[k] = {"column_name": k, "column_type": guess_type(v)}

  return column_types, nrows


def add_phenotype_dataset(name, description, filepath, genotype_datasets, column_spec=None, delim="\t", pid=None):
  db.create_all()
  fullpath = find_file(filepath)

  if ".ped" in fullpath:
    column_types, nrows = _load_phenotype_ped(fullpath, fullpath.replace(".ped",".dat"))
  elif ".dat" in fullpath:
    column_types, nrows = _load_phenotype_ped(fullpath, fullpath.replace(".dat",".ped"))
  elif ".tab" in fullpath:
    column_types, nrows = _load_phenotype_tab(fullpath)
  else:
    raise ValueError("Must specify PED+DAT or tab-delimited file")

  sample_column = None

  if column_spec is not None:
    for col, col_data in column_spec.items():
      for k, v in col_data.items():
        if k == "sample_column":
          sample_column = col
          column_types[col]["for_analysis"] = False
          continue

        if col not in column_types:
          raise ValueError("Column '{}' was specified in YAML, but it does not exist in phenotype dataset '{}'".format(col, name))

        column_types[col][k] = v

  if ".ped" in filepath or ".dat" in filepath:
    sample_column = "IID"

  args = {
    "name": name,
    "description": description,
    "filepath": filepath,
    "nrows": nrows,
    "ncols": len(column_types),
    "delim": delim,
    "sample_column": sample_column if sample_column is not None else next(iter(column_types)) # assume it is the first column unless told otherwise
  }
  if pid is not None:
    args["id"] = pid

  pheno = PhenotypeDataset(**args)
  for col, ctype in column_types.items():
    pheno.columns.append(PhenotypeColumn(**ctype))

  try:
    genotype_datasets = [int(x) for x in genotype_datasets]
  except TypeError as e:
    if not "'int' object is not iterable" == str(e):
      raise Exception("Unexpected exception when parsing genotype dataset ID for phenotype dataset {}".format(args["id"]))

    genotype_datasets = [int(genotype_datasets)]

  for gd in genotype_datasets:
    genotype_dataset = GenotypeDataset.query.filter_by(id = gd).first()
    if genotype_dataset is not None:
      pheno.genotypes.append(genotype_dataset)
    else:
      raise ValueError("Genotype dataset with ID {} does not exist".format(gd))

  db.session.add(pheno)
  db.session.commit()


def add_masks(name, description, filepath, genome_build, genotype_dataset_id, group_type, identifier_type, mid=None):
  args = locals()
  if mid is not None:
    args["id"] = mid

  if "mid" in args: del args["mid"]

  mask = Mask(**args)
  db.session.add(mask)
  db.session.commit()


def create_subset(genome_build, genotype_dataset_name, subset_name, samples_filename):
  db.create_all()
  samples = []
  with open(samples_filename, 'r') as f:
    for sample in f:
      sample = sample.strip()
      if sample:
        samples.append(sample)
  genotype_dataset = GenotypeDataset.query.filter_by(genome_build = genome_build, name = genotype_dataset_name).first()
  if genotype_dataset is None:
    click.echo('Genotype dataset {} on {} genome build doesn\'t exist'.format(genotype_dataset_name, genome_build))
    return
  for sample in samples:
    genotype_dataset.samples.append(Sample(subset = subset_name, sample = sample))
  db.session.commit()


def show_genotypes(genome_build, genotype_dataset_name):
  db.create_all()
  genotype_dataset = GenotypeDataset.query.filter_by(genome_build = genome_build, name = genotype_dataset_name).first()
  if genotype_dataset is None:
    click.echo('Genotype dataset {} on {} genome build doesn\'t exist'.format(genotype_dataset_name, genome_build))
    return
  print '\t'.join(['GENOME BUILD', 'GENOTYPE DATASET NAME', 'FILE'])
  for file in genotype_dataset.files:
    print '\t'.join([genotype_dataset.genome_build, genotype_dataset.name, file.path])


def show_samples(genome_build, genotype_dataset_name, subset_name):
  db.create_all()
  genotype_dataset = GenotypeDataset.query.filter_by(genome_build = genome_build, name = genotype_dataset_name).first()
  if genotype_dataset is None:
    click.echo('Genotype dataset {} on {} genome build doesn\'t exist'.format(genotype_dataset_name, genome_build))
    return
  samples = Sample.query.with_parent(genotype_dataset).filter_by(subset = subset_name).all()
  if not samples:
    click.echo('Population {} doesn\'t exist in {} genotype_dataset'.format(subset_name, genotype_dataset_name))
    return
  print '\t'.join(['GENOME BUILD', 'GENOTYPE DATASET NAME', 'SAMPLE SUBSET', 'SAMPLE'])
  for sample in samples:
    print '\t'.join([genotype_dataset.genome_build, genotype_dataset.name, sample.subset, sample.sample])


@click.command('load-genotypes')
@click.argument('json', type = click.Path(exists = True))
@with_appcontext
def load_genotypes_command(json):
  """Loads genotype datasets from a JSON file.

  json -- JSON file with genotype datasets to load.
  """
  load_genotype_datasets(json)


@click.command('show-genotypes')
@with_appcontext
def show_genotypes_command():
  """Shows loaded genotype datasets."""
  show_genotype_datasets()


@click.command('show-phenotypes')
@with_appcontext
def show_phenotypes_command():
  """Shows loaded phenotype datasets."""
  show_phenotype_datasets()


@click.command('show-masks')
@with_appcontext
def show_masks_command():
  """Shows loaded masks."""
  show_masks()

def add_from_yaml(filepath):
  if not os.path.isfile(filepath):
    raise IOError("Not a file: " + filepath)

  with open(filepath) as fp:
    data = yaml.safe_load(fp)

  # We need to preprocess the file to make sure their requested database IDs have not already been taken.
  for block_type, block_data in data.iteritems():
    id_lookup = None
    if block_type == "genotypes":
      id_lookup = has_genotype_dataset
    elif block_type == "phenotypes":
      id_lookup = has_phenotype_dataset
    elif block_type == "masks":
      id_lookup = has_mask
    else:
      raise ValueError("Unrecognized block in yaml file: " + block_type)

    requested_ids = set()
    for record in block_data:
      if id_lookup(record.get("id")):
        raise ValueError("Database already has ID {} for this type of record. Please make sure your IDs are unique, if you wish to assign them yourself.".format(record["id"]))

      rec_id = record.get("id")
      if rec_id is not None:
        if rec_id in requested_ids:
          raise ValueError("Cannot have two phenotype datasets with the same ID - {} already seen previously".format(rec_id))
        else:
          requested_ids.add(rec_id)

  # Validate column types for phenotypes.
  for record in data.get("phenotypes", []):
    for col, column_record in record.get("columns", {}).items():
      column_type = column_record.get("type")
      if column_type is not None and column_type not in ColumnType.names:
        raise ValueError("Invalid column type '{}' for column '{}' in phenotype dataset '{}'".format(column_type, column_record["name"], record["name"]))

  # We made it here, the IDs are all acceptable to insert into the database.
  if "genotypes" in data:
    for record in data["genotypes"]:
      if "*" in record["filepath"]:
        genotype_files = glob(record["filepath"])
      else:
        genotype_files = [record["filepath"]]

      add_genotype_dataset(record["name"], record["description"], record["genome_build"], record.get("samples"), genotype_files, record.get("id"))

  if "phenotypes" in data:
    for record in data["phenotypes"]:
      add_phenotype_dataset(record["name"], record["description"], record["filepath"], record["genotypes"], record.get("columns"), record.get("delim"), record.get("id"))

  if "masks" in data:
    for record in data["masks"]:
      add_masks(record["name"], record["description"], record["filepath"], record["genome_build"], record["genotypes"], record["group_type"], record["identifier_type"], record.get("id"))

@click.command("add-yaml")
@click.argument("yaml_path", type=click.Path(exists=True))
@with_appcontext
def add_yaml_command(yaml_path):
  add_from_yaml(yaml_path)

@click.command('add-genotypes')
@click.argument('name')
@click.argument('description')
@click.argument('genome_build')
@click.option('--samples', type = click.Path(exists = True))
@click.argument('genotypes', type = click.Path(exists = True), nargs=-1)
@with_appcontext
def add_genotypes_command(name, description, genome_build, samples, genotypes):
  """Adds new genotype dataset to the database. The specified genome build and genotype dataset short name will uniquely identify new genotype dataset.

  name -- The short name of a genotype dataset.\n
  description -- The description of a genotype dataset.\n
  genome_build -- The genome build version.\n
  samples -- File with a list of sample names. One sample name per line.\n
  genotypes -- Path (glob) to VCF/BCF/SAV file(s).\n
  """
  add_genotype_dataset(name, description, genome_build, samples, genotypes)


@click.command('add-phenotypes')
@click.argument('name')
@click.argument('description')
@click.argument('filepath', type = click.Path(exists = True))
@click.argument('genotype_datasets', nargs=-1)
@with_appcontext
def add_phenotypes_command(name, description, filepath, genotype_datasets):
  """Adds new phenotypes to the database.

  name -- The short name of a phenotype dataset.\n
  description -- The longer form description of the phenotype dataset.\n
  filepath -- Path to PED file or TAB-delimited file containing the phenotypes and samples.\n
  genotype_datasets -- List of genotype datasets (IDs) that are linked to this phenotype dataset.\n
  """
  add_phenotype_dataset(name, description, filepath, genotype_datasets)


@click.command('add-masks')
@click.argument('name')
@click.argument('description')
@click.argument('filepath', type = click.Path(exists = True))
@click.argument('genome_build')
@click.argument('genotype_dataset', type = int)
@click.argument('group_type')
@click.argument('identifier_type')
@with_appcontext
def add_masks_command(name, description, filepath, genome_build, genotype_dataset, group_type, identifier_type):
  """Adds new masks to the database.

  name -- The short name of a mask.\n
  description -- The longer form description of the phenotype dataset.\n
  filepath -- Path to PED file or TAB-delimited file containing the phenotypes and samples.\n
  genotype_dataset -- ID of the genotype dataset that corresponds to this mask.\n
  group_type -- What are the groups? Can be "gene", or "region" for arbitrary region.
  identifier_type -- What is the identifier for each group? If genes, can be "ENSEMBL" or "NCBI".
  """
  add_masks(name, description, filepath, genome_build, genotype_dataset, group_type, identifier_type)


@click.command('create-sample-subset')
@click.argument('genome_build')
@click.argument('genotype_dataset')
@click.argument('sample_subset')
@click.argument('samples', type = click.Path(exists = True))
@with_appcontext
def create_subset_command(genome_build, genotype_dataset, sample_subset, samples):
  """Adds new sample subset to an existing genotype dataset.

  genome_build -- The genome build version.\n
  genotype_dataset -- The short name of a genotype dataset.\n
  sample_subset -- The unique short name of a sample subset.\n
  samples -- File with a list of sample names. One sample name per line.\n
  """
  create_subset(genome_build, genotype_dataset, sample_subset, samples)


@click.command('show-sample-subset')
@click.argument('genome_build')
@click.argument('genotype_dataset')
@click.argument('sample_subset')
@with_appcontext
def show_samples_command(genome_build, genotype_dataset, sample_subset):
  """Shows all samples in a given sample subset.

  genome_build -- The genome build version.\n
  genotype_dataset -- The short name of a genotype dataset.\n
  sample_subset -- The unique short name of a sample subset.\n
  """
  show_samples(genome_build, genotype_dataset, sample_subset)
