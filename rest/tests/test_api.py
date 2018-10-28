import pytest
import json
import StringIO
import gzip

def test_correlations(client, config):
    response = client.get('/correlations')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is None
    data = result['data']
    for correlation in data:
        assert all(x in correlation for x in ['name', 'label', 'description', 'type'])


def test_references(client, config):
    response = client.get('/references')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is None
    data = result['data']
    with open(config['REFERENCES_JSON'], 'r') as f:
        references_json = json.load(f)
    assert len(references_json) == len(data)
    for a, b in zip(references_json, data):
        assert a['Name'] == b['name']
        assert a['Description'] == b['description']
        assert a['Genome build'] == b['genome build']
        assert len(a['Samples'].keys()) == len(b['populations'])
        assert all(x in b['populations'] for x in a['Samples'].keys())


def test_reference(client, config):
    response = client.get('/references/1000G_GRCh37')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is None
    data = result['data']
    assert all(x in data for x in ['name', 'description', 'genome build', 'populations'])
    with open(config['REFERENCES_JSON'], 'r') as f:
        reference_json = [x for x in json.load(f) if x['Name'] == data['name']][0]
    assert reference_json['Name'] == data['name']
    assert reference_json['Description'] == data['description']
    assert reference_json['Genome build'] == data['genome build']
    assert len(reference_json['Samples'].keys()) == len(data['populations'])
    assert all(x in data['populations'] for x in reference_json['Samples'].keys())
    response = client.get('/references/SOMETHING_BAD')
    assert response.status_code == 404


def test_reference_populations(client, config):
    response = client.get('/references/1000G_GRCh37/populations')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is None
    data = result['data']
    with open(config['REFERENCES_JSON'], 'r') as f:
        reference_json = [x for x in json.load(f) if x['Name'] == '1000G_GRCh37'][0]
        populations = reference_json['Samples'].keys()
    assert len(data) == len(populations)
    assert all(x in populations for x in data)


def test_reference_population(client, config):
    response = client.get('/references/1000G_GRCh37/populations/EUR')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is None
    data = result['data']
    assert all(x in data for x in ['name', 'size'])
    with open(config['REFERENCES_JSON'], 'r') as f:
        reference_json = [x for x in json.load(f) if x['Name'] == '1000G_GRCh37'][0]
    assert len(reference_json['Samples']['EUR']) == data['size']
    response = client.get('/references/1000G_GRCh37/populations/SOMETHING_BAD')
    assert response.status_code == 404


def test_region_ld(client, goldstandard_ld):
    # response = client.post('/references/1000G_GRCh37/populations/ALL/regions/compute', data = json.dumps({ 'chrom': '22', 'start': 51241101, 'stop': 51241385, 'correlation': 'rsquare' }), content_type = 'application/json')
    response = client.get('/references/1000G_GRCh37/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is None
    data = result['data']
    assert result['next'] is None
    assert all(x in data for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])
    assert all(len(data['correlation']) == len(data[x]) for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])
    goldstandard = goldstandard_ld('../data/region_ld_22_51241101_51241385.hap.ld')
    assert len(goldstandard) == len(data['correlation'])
    for i in xrange(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        assert pytest.approx(data['correlation'][i], 0.00001) == goldstandard[key]

    # response = client.post('/references/1000G_GRCh37/populations/AFR/regions/compute', data = json.dumps({ 'chrom': '22', 'start': 51241101, 'stop': 51241385, 'correlation': 'rsquare'}), content_type = 'application/json')
    response = client.get('/references/1000G_GRCh37/populations/AFR/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is None
    data = result['data']
    assert result['next'] is None
    assert all(x in data for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])
    assert all(len(data['correlation']) == len(data[x]) for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])
    goldstandard = goldstandard_ld('../data/region_ld_22_51241101_51241385.AFR.hap.ld')
    assert len(goldstandard) == len(data['correlation'])
    for i in xrange(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        assert pytest.approx(data['correlation'][i], 0.00001) == goldstandard[key]


def test_variant_ld(client, goldstandard_ld):
    response = client.get('/references/1000G_GRCh37/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is None
    data = result['data']
    assert result['next'] is None
    assert all(x in data for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])
    assert all(len(data['correlation']) == len(data[x]) for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])
    goldstandard = goldstandard_ld('../data/variant_ld_22_51241101_vs_51241101_51241385.hap.ld')
    assert len(goldstandard) == len(data['correlation'])
    for i in xrange(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        assert pytest.approx(data['correlation'][i], 0.00001) == goldstandard[key]


def test_compression(client):
    response = client.get('/references/1000G_GRCh37/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare', headers=[('Accept-Encoding', 'gzip')])
    assert response.headers['Content-Encoding'] == 'gzip'
    compressed_payload = StringIO.StringIO(response.data)
    uncompressed_payload = gzip.GzipFile(fileobj = compressed_payload, mode = 'rb').read()
    assert len(uncompressed_payload) > 0
    result = json.loads(uncompressed_payload)
    assert all(x in result for x in ['data', 'error', 'next'])


def test_malformed_region_ld(client):
    response = client.get('/references/1000G_GRCh37/populations/ALL/regions?start=100&stop=1000&correlation=rsquare')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0

    response = client.get('/references/1000G_GRCh37/populations/ALL/regions?chrom=chr22&start=-100&stop=1000&correlation=rsquare')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0

    response = client.get('/references/1000G_GRCh37/populations/ALL/regions?chrom=chr22&start=100&stop=-1000&correlation=rsquare')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0

    response = client.get('/references/1000G_GRCh37/populations/ALL/regions?chrom=chr22&start=1000&stop=100&correlation=rsquare')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0

    response = client.get('/references/1000G_GRCh37/populations/ALL/regions?chrom=chr22&start=100&stop=1000&correlation=nonsense')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0

    response = client.get('/references/1000G_GRCh37/populations/ALL/regions?chrom=chr22&start=100&stop=1000&correlation=rsquare&limit=0')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0


def test_malformed_variant_ld(client):
    response = client.get('/references/1000G_GRCh37/populations/ALL/variants?chrom=22&start=51241101&stop=51241385&correlation=rsquare')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0

    response = client.get('/references/1000G_GRCh37/populations/ALL/variants?variant=22:51241101_A/T&start=51241101&stop=51241385&correlation=rsquare')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0

    response = client.get('/references/1000G_GRCh37/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=-1000&stop=51241385&correlation=rsquare')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0

    response = client.get('/references/1000G_GRCh37/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=1000&stop=-100&correlation=rsquare')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0

    response = client.get('/references/1000G_GRCh37/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=1000&stop=100&correlation=rsquare')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0

    response = client.get('/references/1000G_GRCh37/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=1000&stop=10000&correlation=nonsense')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0

    response = client.get('/references/1000G_GRCh37/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=1000&stop=10000&correlation=rsquare&limit=0')
    assert response.status_code == 422
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['data'] is None
    assert result['error'] is not None and len(result['error']) > 0