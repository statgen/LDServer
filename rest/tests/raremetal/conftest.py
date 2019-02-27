import os
from raremetal import create_app
import pytest

@pytest.fixture
def app():
    app = create_app({
        'TESTING': True,
        'SQLALCHEMY_DATABASE_URI': 'sqlite:///' + os.path.join(os.path.dirname(__file__), '../../raremetal/sql.db'),
        'SQLALCHEMY_TRACK_MODIFICATIONS': False,
        'PROXY_PASS': None,
        'API_MAX_PAGE_SIZE': 1000,
        'API_MAX_REGION_SIZE': 4000000,
        'SEGMENT_SIZE_BP': 1000,
        'CACHE_ENABLED': False,
        'CACHE_REDIS_HOSTNAME': '127.0.0.1',
        'CACHE_REDIS_PORT': 6379,
        'GZIP_COMPRESSION': True,
        'SENTRY_DSN': None,
        'SENTRY_ENV': None
    })
    print "Using database: " + app.config["SQLALCHEMY_DATABASE_URI"]
    yield app

@pytest.fixture
def client(app):
    return app.test_client()

@pytest.fixture
def config(app):
    return app.config