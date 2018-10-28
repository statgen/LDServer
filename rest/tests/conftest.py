import os
from rest import create_app
import pytest

@pytest.fixture
def app():
    app = create_app({
        'TESTING': True,
        'SQLALCHEMY_DATABASE_URI': 'sqlite:///' + os.path.join(os.path.dirname(__file__), 'sql.db'),
        'SQLALCHEMY_TRACK_MODIFICATIONS': False,
        'PROXY_PASS': None,
        'API_MAX_PAGE_SIZE': 1000,
        'SEGMENT_SIZE_BP': 1000,
        'CACHE_ENABLED': False,
        'CACHE_REDIS_HOSTNAME': '127.0.0.1',
        'CACHE_REDIS_PORT': 6379
    })

    app.config['REFERENCES_JSON'] = os.path.join(os.path.dirname(__file__), 'datasets.json')

    from rest.model import load_references
    with app.app_context():
        load_references(app.config['REFERENCES_JSON'])

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