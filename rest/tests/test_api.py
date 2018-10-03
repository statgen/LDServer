import pytest
import json

def test_references(client, config):
    response = client.get('/references')
    assert response.status_code == 200
    result = response.get_json()
    with open(config['REFERENCES_JSON'], 'r') as f:
        references_json = json.load(f)
    assert len(references_json) == len(result)
    for a, b in zip(references_json, result):
        assert a['Name'] == b['name']
        assert a['Description'] == b['description']
        assert a['Genome build'] == b['genome build']
        assert len(a['Samples'].keys()) == len(b['populations'])
        assert all(x in b['populations'] for x in a['Samples'].keys())

def test_reference(client, config):
    response = client.get('/references/1000G_GRCh37')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['name', 'description', 'genome build', 'populations'])
    with open(config['REFERENCES_JSON'], 'r') as f:
        reference_json = [x for x in json.load(f) if x['Name'] == result['name']][0]
    assert reference_json['Name'] == result['name']
    assert reference_json['Description'] == result['description']
    assert reference_json['Genome build'] == result['genome build']
    assert len(reference_json['Samples'].keys()) == len(result['populations'])
    assert all(x in result['populations'] for x in reference_json['Samples'].keys())
    response = client.get('/references/SOMETHING_BAD')
    assert response.status_code == 404

def test_reference_populations(client, config):
    response = client.get('/references/1000G_GRCh37/populations')
    assert response.status_code == 200
    result = response.get_json()
    with open(config['REFERENCES_JSON'], 'r') as f:
        reference_json = [x for x in json.load(f) if x['Name'] == '1000G_GRCh37'][0]
        populations = reference_json['Samples'].keys()
    assert len(result) == len(populations)
    assert all(x in populations for x in result)

def test_reference_population(client, config):
    response = client.get('/references/1000G_GRCh37/populations/EUR')
    assert response.status_code == 200
    result = response.get_json()
    assert 'name' in result
    assert 'size' in result
    with open(config['REFERENCES_JSON'], 'r') as f:
        reference_json = [x for x in json.load(f) if x['Name'] == '1000G_GRCh37'][0]
    assert len(reference_json['Samples']['EUR']) == result['size']
    response = client.get('/references/1000G_GRCh37/populations/SOMETHING_BAD')
    assert response.status_code == 404

def test_region_ld(client, goldstandard_ld):
    # response = client.post('/references/1000G_GRCh37/populations/ALL/regions/compute', data = json.dumps({ 'chrom': '22', 'start': 51241101, 'stop': 51241385, 'correlation': 'rsquare' }), content_type = 'application/json')
    response = client.get('/references/1000G_GRCh37/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare')
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

    # response = client.post('/references/1000G_GRCh37/populations/AFR/regions/compute', data = json.dumps({ 'chrom': '22', 'start': 51241101, 'stop': 51241385, 'correlation': 'rsquare'}), content_type = 'application/json')
    response = client.get('/references/1000G_GRCh37/populations/AFR/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare')
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
    response = client.get('/references/1000G_GRCh37/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare')
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

