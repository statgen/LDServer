from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Index, inspect
from sqlalchemy.types import Enum
from flask.cli import with_appcontext
from collections import Counter, OrderedDict
from ld.pywrapper import ColumnType, ColumnTypeMap, VariantGroupType, GroupIdentifierType, extract_samples
from tabulate import tabulate
import click
import json
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

  columns = db.relationship("PhenotypeColumn", back_populates = "dataset")
  genotypes = db.relationship("GenotypeDataset", secondary="genotype_phenotype", back_populates="phenotypes")


class PhenotypeColumn(db.Model):
  __tablename__ = "phenotype_columns"
  id = db.Column(db.Integer, primary_key = True)
  phenotype_dataset_id = db.Column(db.Integer, db.ForeignKey('phenotype_datasets.id'))
  column_name = db.Column(db.String, nullable = False)
  column_type = db.Column(Enum(*tuple(ColumnType.values[i].name for i in range(len(ColumnType.values)))), nullable = False, name = "column_type")

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

def get_genotype_dataset(genome_build, genotype_dataset_name):
  genotype_dataset = GenotypeDataset.query.filter_by(genome_build = genome_build, name = genotype_dataset_name).first()
  if genotype_dataset is not None:
    return { 'name': genotype_dataset.name, 'description': genotype_dataset.description, 'genome build': genotype_dataset.genome_build }
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

def get_phenotype_nrows(phenotype_dataset_id):
  return db.session.query(PhenotypeDataset.nrows).filter_by(id = phenotype_dataset_id).scalar()

def get_phenotype_columns(phenotype_dataset_id):
  return [str(x) for x in db.session.query(PhenotypeColumn.column_name).filter_by(id = phenotype_dataset_id).order_by(PhenotypeColumn.id)]

def get_phenotypes_for_genotypes(genotype_dataset_id):
  result = []
  for row in db.session.query(GenotypeDataset).filter_by(id = genotype_dataset_id).scalar().phenotypes:
    as_dict = {c.key: getattr(row, c.key) for c in inspect(row).mapper.column_attrs}
    as_dict["phenotypeDataset"] = as_dict["id"]
    del as_dict["id"]
    del as_dict["nrows"]
    del as_dict["ncols"]
    del as_dict["filepath"]
    as_dict["phenotypes"] = [x.column_name for x in row.columns]

    result.append(as_dict)

  return result

def get_masks_for_genotypes(genotype_dataset_id):
  masks = []
  for row in db.session.query(Mask).filter_by(genotype_dataset_id = genotype_dataset_id):
    as_dict = {c.key: getattr(row, c.key) for c in inspect(row).mapper.column_attrs}
    as_dict["id"] = as_dict.pop("name")
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
  header = ["ID", "Name", "Description", "Genome Build", "Populations"]
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

def add_genotype_dataset(name, description, genome_build, samples_filename, genotype_files):
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
      f_samples = extract_samples(f)
      if len(samples) == 0:
        samples = f_samples
      else:
        if set(samples) != set(f_samples):
          raise ValueError("Genotype file {} did not have the same samples as the others")

  r = GenotypeDataset(name = name, description = description, genome_build = genome_build)
  for path in genotype_files:
    r.files.append(File(path = path))
  for sample in samples:
    r.samples.append(Sample(subset = 'ALL', sample = sample))
  db.session.add(r)
  db.session.commit()
  print r.id

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
    column_types["IID"] = ColumnType.TEXT.name
    column_types["FID"] = ColumnType.TEXT.name
    column_types["PID"] = ColumnType.TEXT.name
    column_types["MID"] = ColumnType.TEXT.name
    column_types["SEX"] = ColumnType.CATEGORICAL.name

    for i, line in enumerate(fp_dat):
      ctype, pheno = line.split()
      if ctype == "A":
        column_types[pheno] = ColumnType.CATEGORICAL.name
      elif ctype == "T":
        column_types[pheno] = ColumnType.FLOAT.name
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
    column_types[k] = guess_type(v)

  return column_types, nrows


def add_phenotypes(name, description, filepath, genotype_datasets):
  db.create_all()

  if ".ped" in filepath:
    column_types, nrows = _load_phenotype_ped(filepath, filepath.replace(".ped",".dat"))
  elif ".dat" in filepath:
    column_types, nrows = _load_phenotype_ped(filepath, filepath.replace(".dat",".ped"))
  elif ".tab" in filepath:
    column_types, nrows = _load_phenotype_tab(filepath)
  else:
    raise ValueError("Must specify PED+DAT or tab-delimited file")

  pheno = PhenotypeDataset(name = name, description = description, filepath = filepath, nrows = nrows, ncols = len(column_types))
  for col, ctype in column_types.items():
    pheno.columns.append(PhenotypeColumn(column_name = col, column_type = ctype))

  for gd in genotype_datasets:
    genotype_dataset = GenotypeDataset.query.filter_by(id = gd).first()
    if genotype_dataset is not None:
      pheno.genotypes.append(genotype_dataset)
    else:
      raise ValueError("Genotype dataset with name {} does not exist".format(gd))

  db.session.add(pheno)
  db.session.commit()


def add_masks(name, description, filepath, genome_build, genotype_dataset, group_type, identifier_type):
  mask = Mask(
    name = name,
    description = description,
    filepath = filepath,
    genome_build = genome_build,
    genotype_dataset_id = genotype_dataset,
    group_type = group_type,
    identifier_type = identifier_type
  )
  db.session.add(mask)
  db.session.commit()
  print mask.id


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
  add_phenotypes(name, description, filepath, genotype_datasets)


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