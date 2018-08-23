import pytest


def test_root(client):
    response = client.get('/')
    assert response.status_code == 200
    result = response.get_json()
    assert 'references' in result
    assert 'api_version' in result
    for reference in result['references']:
        assert all(x in reference for x in ['name', 'description', 'genome build', 'populations'])


def test_reference(client):
    response = client.get('/1000G_GRCh37')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['name', 'description', 'genome build', 'populations'])
    response = client.get('/SOMETHING_BAD')
    assert response.status_code == 404


def test_reference_population(client):
    response = client.get('/1000G_GRCh37/EUR')
    assert response.status_code == 200
    result = response.get_json()
    assert 'name' in result
    assert 'size' in result
    response = client.get('/SOMETHING_BAD/EUR')
    assert response.status_code == 404
    response = client.get('/1000G_GRCh37/SOMETHING_BAD')
    assert response.status_code == 404


def test_region_ld(client, goldstandard_ld):
    response = client.get('/1000G_GRCh37/ALL/ld/region?chrom=22&start=51241101&stop=51241385')
    assert response.status_code == 200
    result = response.get_json()
    assert 'data' in result
    assert 'next' in result
    assert all(x in result['data'] for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'rsquare'])
    assert result['next'] is None
    assert all(len(result['data']['rsquare']) == len(result['data'][x]) for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'rsquare'])
    goldstandard = goldstandard_ld('../data/region_ld_22_51241101_51241385.hap.ld')
    assert len(goldstandard) == len(result['data']['rsquare'])
    for i in xrange(0, len(goldstandard)):
        key = str(result['data']['position1'][i]) + '_' + str(result['data']['position2'][i])
        assert key in goldstandard
        assert pytest.approx(result['data']['rsquare'][i], 0.00001) == goldstandard[key]

    response = client.get('/1000G_GRCh37/AFR/ld/region?chrom=22&start=51241101&stop=51241385')
    assert response.status_code == 200
    result = response.get_json()
    assert 'data' in result
    assert 'next' in result
    assert all(x in result['data'] for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'rsquare'])
    assert result['next'] is None
    assert all(len(result['data']['rsquare']) == len(result['data'][x]) for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'rsquare'])
    goldstandard = goldstandard_ld('../data/region_ld_22_51241101_51241385.AFR.hap.ld')
    assert len(goldstandard) == len(result['data']['rsquare'])
    for i in xrange(0, len(goldstandard)):
        key = str(result['data']['position1'][i]) + '_' + str(result['data']['position2'][i])
        assert key in goldstandard
        assert pytest.approx(result['data']['rsquare'][i], 0.00001) == goldstandard[key]


def test_variant_ld(client, goldstandard_ld):
    response = client.get('/1000G_GRCh37/ALL/ld/variant?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385')
    assert response.status_code == 200
    result = response.get_json()
    assert 'data' in result
    assert 'next' in result
    assert all(x in result['data'] for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'rsquare'])
    assert result['next'] is None
    assert all(len(result['data']['rsquare']) == len(result['data'][x]) for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'rsquare'])
    goldstandard = goldstandard_ld('../data/variant_ld_22_51241101_vs_51241101_51241385.hap.ld')
    assert len(goldstandard) == len(result['data']['rsquare'])
    for i in xrange(0, len(goldstandard)):
        key = str(result['data']['position1'][i]) + '_' + str(result['data']['position2'][i])
        assert key in goldstandard
        assert pytest.approx(result['data']['rsquare'][i], 0.00001) == goldstandard[key]

