from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Index
from sqlalchemy.types import Enum
from flask.cli import with_appcontext
from collections import Counter, OrderedDict
from ld.pywrapper import ColumnType
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
  __tablename__ = 'genotype_dataset'
  id = db.Column(db.Integer, primary_key = True)
  name = db.Column(db.String(), unique = False, nullable = False)
  genome_build = db.Column(db.String(), unique = False, nullable = False)
  description = db.Column(db.String(), unique = False, nullable = False)
  files = db.relationship('File', backref = 'genotype_dataset', lazy = True)
  samples = db.relationship('Sample', backref = 'genotype_dataset', lazy = True)
  __table_args__ = (
    db.UniqueConstraint('name', 'genome_build', name = 'name_build_uc'),
  )
  def __repr__(self):
    return '<GenotypeDataset %r>' % self.name

class File(db.Model):
  __tablename__ = 'file'
  id = db.Column(db.Integer, primary_key = True)
  path = db.Column(db.String(), unique = False, nullable = False)
  genotype_dataset_id = db.Column(db.Integer, db.ForeignKey('genotype_dataset.id'), nullable = False)
  def __repr__(self):
    return '<File %r>' % self.path

class Sample(db.Model):
  __tablename__ = 'sample'
  id = db.Column(db.Integer, primary_key = True)
  subset = db.Column(db.String, unique = False, nullable = False)
  sample = db.Column(db.String, unique = False, nullable = False)
  genotype_dataset_id = db.Column(db.Integer, db.ForeignKey('genotype_dataset.id'), nullable = False)
  def __repr__(self):
    return '<Sample %r %r>' % (self.subset, self.sample)

class PhenotypeDataset(db.Model):
  __tablename__ = "phenotype_datasets"
  id = db.Column(db.Integer, primary_key = True)
  name = db.Column(db.String, unique = False, nullable = False)
  description = db.Column(db.String, unique = False, nullable = False)
  filepath = db.Column(db.String, unique = False, nullable = False)
  nrows = db.Column(db.Integer, nullable = False)
  ncols = db.Column(db.Integer, nullable = False)

  columns = db.relationship("PhenotypeColumn", back_populates = "dataset")

class PhenotypeColumn(db.Model):
  __tablename__ = "phenotype_columns"
  id = db.Column(db.Integer, primary_key = True)
  phenotype_id = db.Column(db.Integer, db.ForeignKey('phenotype_datasets.id'))
  column_name = db.Column(db.String, nullable = False)
  column_type = db.Column(Enum(*tuple(ColumnType.values[i].name for i in range(len(ColumnType.values)))), nullable = False, name = "column_type")

  dataset = db.relationship("PhenotypeDataset", back_populates = "columns")

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
  print '\t'.join(['NAME', 'DESCRIPTION', 'GENOME BUILD', 'POPULATIONS'])
  for genotype_dataset in GenotypeDataset.query.all():
    print '\t'.join([genotype_dataset.name, genotype_dataset.description, genotype_dataset.genome_build, ';'.join(list(set([s.subset for s in genotype_dataset.samples])))])

def add_genotype_dataset(name, description, genome_build, samples_filename, genotype_files):
  db.create_all()
  samples = []
  with open(samples_filename, 'r') as f:
    for sample in f:
      sample = sample.strip()
      if sample:
        samples.append(sample)
  r = GenotypeDataset(name = name, description = description, genome_build = genome_build)
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

def add_phenotypes(name, description, filepath):
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

  db.session.add(pheno)
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

@click.command('load-genotype-datasets')
@click.argument('json', type = click.Path(exists = True))
@with_appcontext
def load_genotype_datasets_command(json):
  """Loads genotype datasets from a JSON file.

  json -- JSON file with genotype datasets to load.
  """
  load_genotype_datasets(json)

@click.command('show-genotype-datasets')
@with_appcontext
def show_genotype_datasets_command():
  """Shows loaded genotype datasets."""
  show_genotype_datasets()

@click.command('add-genotype-dataset')
@click.argument('name')
@click.argument('description')
@click.argument('genome_build')
@click.argument('samples', type = click.Path(exists = True))
@click.argument('genotypes', type = click.Path(exists = True), nargs=-1)
@with_appcontext
def add_genotype_dataset_command(name, description, genome_build, samples, genotypes):
  """Adds new genotype dataset to the database. The specified genome build and genotype dataset short name will uniquely identify new genotype dataset.

  name -- The short name of a genotype dataset panel.\n
  description -- The description of a genotype dataset panel.\n
  genome_build -- The genome build version.\n
  samples -- File with a list of sample names. One sample name per line.\n
  genotypes -- Path (glob) to VCF/BCF/SAV file(s).\n
  """
  add_genotype_dataset(name, description, genome_build, samples, genotypes)

@click.command('add-phenotypes')
@click.argument('name')
@click.argument('description')
@click.argument('filepath', type = click.Path(exists = True))
@with_appcontext
def add_phenotypes_command(name, description, filepath):
  """Adds new phenotypes to the database.

  name -- The short name of a phenotype dataset.\n
  description -- The longer form description of the phenotype dataset.\n
  filepath -- Path to PED file or TAB-delimited file containing the phenotypes and samples.\n
  """
  add_phenotypes(name, description, filepath)

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

@click.command('show-genotypes')
@click.argument('genome_build')
@click.argument('genotype_dataset')
@with_appcontext
def show_genotypes_command(genome_build, genotype_dataset):
  """Shows raw genotype files.

  genome_build -- The genome build version.\n
  genotype_dataset -- The short name of a genotype dataset.\n
  """
  show_genotypes(genome_build, genotype_dataset)

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
