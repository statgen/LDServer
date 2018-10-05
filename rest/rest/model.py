from flask_sqlalchemy import SQLAlchemy
from flask import current_app
from flask.cli import with_appcontext
import click
import json

db = SQLAlchemy()

class Reference(db.Model):
    __tablename__ = 'reference'
    id = db.Column(db.Integer, primary_key = True)
    name = db.Column(db.String(), unique = True, nullable = False)
    description = db.Column(db.String(), unique = False, nullable = False)
    genome_build = db.Column(db.String(), unique = False, nullable = False)
    files = db.relationship('File', backref = 'reference', lazy = True)
    samples = db.relationship('Sample', backref = 'reference', lazy = True)
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

def load_references(json_file):
    db.drop_all()
    db.create_all()
    with open(json_file, 'r') as f:
        for reference in json.load(f):
            r = Reference(name = reference['Name'], description = reference['Description'], genome_build = reference['Genome build'])
            for path in reference['Files']:
                r.files.append(File(path = path))
            for subset, samples in reference['Samples'].iteritems():
                for sample in samples:
                    r.samples.append(Sample(subset = subset, sample = sample));
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

def create_subset(reference_name, subset_name, samples_filename):
    db.create_all()
    samples = []
    with open(samples_filename, 'r') as f:
        for sample in f:
            sample = sample.strip()
            if sample:
                samples.append(sample)
    reference = Reference.query.filter_by(name = reference_name).first()
    if reference is None:
        click.echo('Reference {} doesn\'t exist'.format(reference_name))
        return
    for sample in samples:
        reference.samples.append(Sample(subset = subset_name, sample = sample))
    db.session.commit()

def show_genotypes(reference_name):
    db.create_all()
    reference = Reference.query.filter_by(name = reference_name).first()
    if reference is None:
        click.echo('Reference {} doesn\'t exist'.format(reference_name))
        return
    print '\t'.join(['REFERENCE', 'FILE'])
    for file in reference.files:
        print '\t'.join([reference.name, file.path])

def show_samples(reference_name, subset_name):
    db.create_all()
    reference = Reference.query.filter_by(name = reference_name).first()
    if reference is None:
        click.echo('Reference {} doesn\'t exist'.format(reference_name))
        return
    samples = Sample.query.with_parent(reference).filter_by(subset = subset_name).all()
    if not samples:
        click.echo('Population {} doesn\'t exist in {} reference'.format(subset_name, reference_name))
        return
    print '\t'.join(['REFERENCE', 'POPULATION', 'SAMPLE'])
    for sample in samples:
        print '\t'.join([reference.name, sample.subset, sample.sample])

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
    """Adds new reference to the database.

    name -- The unique short name of a reference panel.\n
    description -- The description of a reference panel.\n
    genome_build -- The genome build version.\n
    samples -- File with a list of sample names. One sample name per line.\n
    genotypes -- Path (glob) to VCF/BCF/SAV file(s).\n
    """
    add_reference(name, description, genome_build, samples, genotypes)

@click.command('create-population')
@click.argument('reference')
@click.argument('population')
@click.argument('samples', type = click.Path(exists = True))
@with_appcontext
def create_subset_command(reference, population, samples):
    """Adds new population to the existing reference.

    reference -- The unique short name of a reference panel.\n
    population -- The unique short name of a population.\n
    samples -- File with a list of sample names. One sample name per line.\n
    """
    create_subset(reference, population, samples)

@click.command('show-genotypes')
@click.argument('reference')
@with_appcontext
def show_genotypes_command(reference):
    """Shows raw genotypes files.

    reference -- The unique short name of a reference panel.\n
    """
    show_genotypes(reference)

@click.command('show-population')
@click.argument('reference')
@click.argument('population')
@with_appcontext
def show_samples_command(reference, population):
    """Shows all samples from a given population.

    reference -- The unique short name of a reference panel.\n
    population -- The unique short name of a population.\n
    """
    show_samples(reference, population)