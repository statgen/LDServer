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

def load_references():
    db.drop_all()
    db.create_all()
    with open(current_app.config['REFERENCES_JSON']) as f:
        for reference in json.load(f):
            r = Reference(name = reference['Name'], description = reference['Description'], genome_build = reference['Genome build'])
            for path in reference['Files']:
                r.files.append(File(path = path))
            for subset, samples in reference['Samples'].iteritems():
                for sample in samples:
                    r.samples.append(Sample(subset = subset, sample = sample));
            db.session.add(r)
            db.session.commit()

@click.command('load-references')
@with_appcontext
def load_references_command():
    load_references()
    click.echo('References loaded.')