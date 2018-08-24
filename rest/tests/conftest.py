import os
from rest import create_app
import pytest

@pytest.fixture
def app():
    app = create_app({
        'TESTING': True,
        'REFERENCES_JSON': os.path.join(os.path.dirname(__file__), 'datasets.json'),
        'SQLALCHEMY_DATABASE_URI': 'sqlite:///' + os.path.join(os.path.dirname(__file__), 'sql.db'),
        'SQLALCHEMY_TRACK_MODIFICATIONS': False,
        'API_MAX_PAGE_SIZE': 1000
    })

    from rest.model import load_references
    with app.app_context():
        load_references()

    yield app

@pytest.fixture
def client(app):
    return app.test_client()

@pytest.fixture
def config(app):
    return app.config

@pytest.fixture
def goldstandard_ld():
    def _goldstandard_ld(filename):
        goldstandard = dict()
        with open(filename, 'r') as f:
            header = f.readline().rstrip().split('\t')
            for line in f:
                record = dict(zip(header, line.rstrip().split('\t')))
                goldstandard[record['POS1'] + '_' + record['POS2']] = float(record['R^2'])
        return goldstandard
    yield _goldstandard_ld