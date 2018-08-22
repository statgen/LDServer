import pytest

def test_root(client):
    response = client.get('/').get_json()
    assert 'references' in response
    assert 'api_version' in response
    for reference in response['references']:
        assert  'name' in reference
        assert  'description' in reference
        assert  'genome build' in reference
        assert  'populations' in reference

def test_reference(client):
    response = client.get('/1000G_GRCh37').get_json()
    assert  'name' in response
    assert  'description' in response
    assert  'genome build' in response
    assert  'populations' in response
    response = client.get('/SOMETHING_BAD')
    assert 404 == response.status_code

def test_reference_population(client):
    response = client.get('/1000G_GRCh37/EUR').get_json()
    assert 'name' in response
    assert 'size' in response
    response = client.get('/SOMETHING_BAD/EUR')
    assert 404 == response.status_code
    response = client.get('/1000G_GRCh37/SOMETHING_BAD')
    assert 404 == response.status_code

def test_region_ld(client):
    response = client.get('/1000G_GRCh37/ALL/ld/region?chrom=22&start=51241101&stop=51241385').get_json()
    assert 'data' in response
    assert 'next' in response
    assert 'variant1' in response['data']
    assert 'chromosome1' in response['data']
    assert 'position1' in response['data']
    assert 'chromosome2' in response['data']
    assert 'position2' in response['data']
    assert 'rsquare' in response['data']

    response = client.get('/1000G_GRCh37/AFR/ld/region?chrom=22&start=51241101&stop=51241385').get_json()
    assert 'data' in response
    assert 'next' in response
    assert 'variant1' in response['data']
    assert 'chromosome1' in response['data']
    assert 'position1' in response['data']
    assert 'chromosome2' in response['data']
    assert 'position2' in response['data']
    assert 'rsquare' in response['data']

def test_variant_ld(client):
    response = client.get('/1000G_GRCh37/ALL/ld/variant?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385').get_json()
    assert 'data' in response
    assert 'next' in response
    assert 'variant1' in response['data']
    assert 'chromosome1' in response['data']
    assert 'position1' in response['data']
    assert 'chromosome2' in response['data']
    assert 'position2' in response['data']
    assert 'rsquare' in response['data']