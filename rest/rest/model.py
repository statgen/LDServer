from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Index
from sqlalchemy.types import Enum
from flask.cli import with_appcontext
from collections import Counter, OrderedDict
import click
import json
import gzip
import enum

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

class Reference(db.Model):
    __tablename__ = 'reference'
    id = db.Column(db.Integer, primary_key = True)
    name = db.Column(db.String(), unique = False, nullable = False)
    genome_build = db.Column(db.String(), unique = False, nullable = False)
    description = db.Column(db.String(), unique = False, nullable = False)
    files = db.relationship('File', backref = 'reference', lazy = True)
    samples = db.relationship('Sample', backref = 'reference', lazy = True)
    __table_args__ = (
        db.UniqueConstraint('name', 'genome_build', name = 'name_build_uc'),
    )
    def __repr__(self):
        return '<Reference %r>' % self.name

class File(db.Model):
    __tablename__ = 'file'
    id = db.Column(db.Integer, primary_key = True)
    path = db.Column(db.String(), unique = False, nullable = False)
    reference_id = db.Column(db.Integer, db.ForeignKey('reference.id'), nullable = False)
    def __repr__(self):
        return '<File %r>' % self.path

class Sample(db.Model):
    __tablename__ = 'sample'
    id = db.Column(db.Integer, primary_key = True)
    subset = db.Column(db.String, unique = False, nullable = False)
    sample = db.Column(db.String, unique = False, nullable = False)
    reference_id = db.Column(db.Integer, db.ForeignKey('reference.id'), nullable = False)
    def __repr__(self):
        return '<Sample %r %r>' % (self.subset, self.sample)

class ColumnType(enum.Enum):
    TEXT = 0
    CATEGORICAL = 1
    INTEGER = 2
    FLOAT = 3

class PhenotypeDataset(db.Model):
    __tablename__ = "phenotype_datasets"
    id = db.Column(db.Integer, primary_key = True)
    name = db.Column(db.String, unique = False, nullable = False)
    description = db.Column(db.String, unique = False, nullable = False)
    filepath = db.Column(db.String, unique = False, nullable = False)

    columns = db.relationship("PhenotypeColumn", back_populates = "dataset")

class PhenotypeColumn(db.Model):
    __tablename__ = "phenotype_columns"
    id = db.Column(db.Integer, primary_key = True)
    phenotype_id = db.Column(db.Integer, db.ForeignKey('phenotype_datasets.id'))
    column_name = db.Column(db.String, nullable = False)
    column_type = db.Column(Enum(ColumnType), nullable = False, name = "column_type")

    dataset = db.relationship("PhenotypeDataset", back_populates = "columns")

Index('reference_index_1', Reference.genome_build)
Index('file_index', File.reference_id, File.path)
Index('sample_index', Sample.reference_id, Sample.subset, Sample.sample)

def get_correlations():
    return [{'name': name, 'label': label, 'description': desc, 'type': type} for name, label, desc, type in db.session.query(Correlation.name, Correlation.label, Correlation.descripition, Correlation.type)]

def has_genome_build(genome_build):
    return db.session.query(Reference.genome_build).filter_by(genome_build = genome_build).first() is not None

def get_genome_builds():
    return [x for x, in db.session.query(Reference.genome_build).distinct()]

def get_references(genome_build):
    return [{'name': name, 'description': desc } for name, desc, in db.session.query(Reference.name, Reference.description).filter_by(genome_build = genome_build)]

def get_reference_id(genome_build, reference_name):
    return db.session.query(Reference.id).filter_by(genome_build = genome_build).filter_by(name = reference_name).scalar()

def get_reference(genome_build, reference_name):
    reference = Reference.query.filter_by(genome_build = genome_build, name = reference_name).first()
    if reference is not None:
        return { 'name': reference.name, 'description': reference.description, 'genome build': reference.genome_build }
    return None

def has_population(reference_id, subset):
    return db.session.query(Sample.subset).filter_by(reference_id = reference_id).filter_by(subset = subset).first() is not None

def get_populations(reference_id):
    return [x for x, in db.session.query(Sample.subset).filter_by(reference_id = reference_id).distinct()]

def has_samples(reference_id, subset):
    return db.session.query(Sample.sample).filter_by(reference_id = reference_id).filter_by(subset = subset).first() is not None

def get_samples_count(reference_id, subset):
    return db.session.query(Sample).filter_by(reference_id = reference_id).filter_by(subset = subset).count()

def get_samples(reference_id, subset):
    return [str(x) for x, in db.session.query(Sample.sample).filter_by(reference_id = reference_id).filter_by(subset = subset)]

def get_files(reference_id):
    return [str(x) for x, in db.session.query(File.path).filter_by(reference_id = reference_id)]

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

def load_references(json_file):
    for t in db.metadata.sorted_tables:
        if t.name != 'correlation':
            t.drop(db.engine)
    db.create_all()
    with open(json_file, 'r') as f:
        for reference in json.load(f):
            r = Reference(name = reference['Name'], description = reference['Description'], genome_build = reference['Genome build'])
            for path in reference['Files']:
                r.files.append(File(path = path))
            for subset, samples in reference['Samples'].iteritems():
                for sample in samples:
                    r.samples.append(Sample(subset = subset, sample = sample))
            db.session.add(r)
        db.session.commit()

def show_references():
    db.create_all()
    print '\t'.join(['NAME', 'DESCRIPTION', 'GENOME BUILD', 'POPULATIONS'])
    for reference in Reference.query.all():
        print '\t'.join([reference.name, reference.description, reference.genome_build, ';'.join(list(set([s.subset for s in reference.samples])))])

def add_reference(name, description, genome_build, samples_filename, genotype_files):
    db.create_all()
    samples = []
    with open(samples_filename, 'r') as f:
        for sample in f:
            sample = sample.strip()
            if sample:
                samples.append(sample)
    r = Reference(name = name, description = description, genome_build = genome_build)
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
            guesses.append(ColumnType.TEXT)
            continue

        i = int(f)
        if f == i:
            guesses.append(ColumnType.INTEGER)
        else:
            guesses.append(ColumnType.FLOAT)

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

    with fp_ped, fp_dat:
        column_types = OrderedDict()

        # We know some of the column types already in a PED file
        column_types["IID"] = ColumnType.TEXT
        column_types["FID"] = ColumnType.TEXT
        column_types["PID"] = ColumnType.TEXT
        column_types["MID"] = ColumnType.TEXT
        column_types["SEX"] = ColumnType.CATEGORICAL

        for i, line in enumerate(fp_dat):
            ctype, pheno = line.split()
            if ctype == "A":
                column_types[pheno] = ColumnType.CATEGORICAL
            elif ctype == "T":
                column_types[pheno] = ColumnType.FLOAT
            elif ctype == "M":
                continue
            else:
                ValueError("Unrecognized DAT data type: " + ctype)

    return column_types

def _load_phenotype_tab(tab_file, lookahead = 1000):
    if tab_file.endswith(".gz"):
        fp_tab = gzip.open(tab_file, "rt")
    else:
        fp_tab = open(tab_file)

    header = None
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

            if i > lookahead:
                break

    column_types = OrderedDict()
    for k, v in column_values.items():
        column_types[k] = guess_type(v)
    return column_types

def add_phenotypes(name, description, filepath):
    db.create_all()

    if ".ped" in filepath:
        column_types = _load_phenotype_ped(filepath, filepath.replace(".ped",".dat"))
    elif ".dat" in filepath:
        column_types = _load_phenotype_ped(filepath, filepath.replace(".dat",".ped"))
    elif ".tab" in filepath:
        column_types = _load_phenotype_tab(filepath)
    else:
        raise ValueError("Must specify PED+DAT or tab-delimited file")

    pheno = PhenotypeDataset(name = name, description = description, filepath = filepath)
    for col, ctype in column_types.items():
        pheno.columns.append(PhenotypeColumn(column_name = col, column_type = ctype))

    db.session.add(pheno)
    db.session.commit()

def create_subset(genome_build, reference_name, subset_name, samples_filename):
    db.create_all()
    samples = []
    with open(samples_filename, 'r') as f:
        for sample in f:
            sample = sample.strip()
            if sample:
                samples.append(sample)
    reference = Reference.query.filter_by(genome_build = genome_build, name = reference_name).first()
    if reference is None:
        click.echo('Reference {} on {} genome build doesn\'t exist'.format(reference_name, genome_build))
        return
    for sample in samples:
        reference.samples.append(Sample(subset = subset_name, sample = sample))
    db.session.commit()

def show_genotypes(genome_build, reference_name):
    db.create_all()
    reference = Reference.query.filter_by(genome_build = genome_build, name = reference_name).first()
    if reference is None:
        click.echo('Reference {} on {} genome build doesn\'t exist'.format(reference_name, genome_build))
        return
    print '\t'.join(['GENOME BUILD', 'REFERENCE NAME', 'FILE'])
    for file in reference.files:
        print '\t'.join([reference.genome_build, reference.name, file.path])

def show_samples(genome_build, reference_name, subset_name):
    db.create_all()
    reference = Reference.query.filter_by(genome_build = genome_build, name = reference_name).first()
    if reference is None:
        click.echo('Reference {} on {} genome build doesn\'t exist'.format(reference_name, genome_build))
        return
    samples = Sample.query.with_parent(reference).filter_by(subset = subset_name).all()
    if not samples:
        click.echo('Population {} doesn\'t exist in {} reference'.format(subset_name, reference_name))
        return
    print '\t'.join(['GENOME BUILD', 'REFERENCE NAME', 'POPULATION', 'SAMPLE'])
    for sample in samples:
        print '\t'.join([reference.genome_build, reference.name, sample.subset, sample.sample])

@click.command('load-references')
@click.argument('json', type = click.Path(exists = True))
@with_appcontext
def load_references_command(json):
    """Loads references from a JSON file.

    json -- JSON file with references to load.
    """
    load_references(json)

@click.command('show-references')
@with_appcontext
def show_references_command():
    """Shows loaded references."""
    show_references()

@click.command('add-reference')
@click.argument('name')
@click.argument('description')
@click.argument('genome_build')
@click.argument('samples', type = click.Path(exists = True))
@click.argument('genotypes', type=click.Path(exists = True), nargs=-1)
@with_appcontext
def add_reference_command(name, description, genome_build, samples, genotypes):
    """Adds new reference to the database. The specified genome build and reference short name will uniquely identify new reference.

    name -- The short name of a reference panel.\n
    description -- The description of a reference panel.\n
    genome_build -- The genome build version.\n
    samples -- File with a list of sample names. One sample name per line.\n
    genotypes -- Path (glob) to VCF/BCF/SAV file(s).\n
    """
    add_reference(name, description, genome_build, samples, genotypes)

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

@click.command('create-population')
@click.argument('genome_build')
@click.argument('reference')
@click.argument('population')
@click.argument('samples', type = click.Path(exists = True))
@with_appcontext
def create_subset_command(genome_build, reference, population, samples):
    """Adds new population to the existing reference.

    genome_build -- The genome build version.\n
    reference -- The short name of a reference panel.\n
    population -- The unique short name of a population.\n
    samples -- File with a list of sample names. One sample name per line.\n
    """
    create_subset(genome_build, reference, population, samples)

@click.command('show-genotypes')
@click.argument('genome_build')
@click.argument('reference')
@with_appcontext
def show_genotypes_command(genome_build, reference):
    """Shows raw genotypes files.

    genome_build -- The genome build version.\n
    reference -- The short name of a reference panel.\n
    """
    show_genotypes(genome_build, reference)

@click.command('show-population')
@click.argument('genome_build')
@click.argument('reference')
@click.argument('population')
@with_appcontext
def show_samples_command(genome_build, reference, population):
    """Shows all samples from a given population.

    genome_build -- The genome build version.\n
    reference -- The short name of a reference panel.\n
    population -- The unique short name of a population.\n
    """
    show_samples(genome_build, reference, population)
