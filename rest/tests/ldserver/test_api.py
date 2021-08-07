import pytest
import json
import io
import gzip
import msgpack

def test_correlations(client):
    response = client.get('/correlations')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is None
    data = result['data']
    assert len(data) == 4
    for correlation in data:
        assert all(x in correlation for x in ['name', 'label', 'description', 'type'])


def test_genome_builds(client):
    response = client.get('/genome_builds')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])


def test_references(client, config):
    response = client.get('/genome_builds/GRCh37/references')
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

    response = client.get('/genome_builds/SOMETHING/references')
    assert response.status_code == 404
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is not None
    assert result['data'] is None


def test_reference(client, config):
    response = client.get('/genome_builds/GRCh37/references/1000G')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is None
    data = result['data']
    assert all(x in data for x in ['name', 'description', 'genome build'])
    with open(config['REFERENCES_JSON'], 'r') as f:
        reference_json = [x for x in json.load(f) if x['Name'] == data['name']][0]
    assert reference_json['Name'] == data['name']
    assert reference_json['Description'] == data['description']
    assert reference_json['Genome build'] == data['genome build']

    for bad_url in ['genome_builds/GRCh37/references/SOMETHING_BAD', 'genome_builds/SOMETHING_BAD/references/1000G']:
        response = client.get(bad_url)
        assert response.status_code == 404
        result = response.get_json()
        assert all(x in result for x in ['data', 'error'])
        assert result['error'] is not None
        assert result['data'] is None


def test_reference_populations(client, config):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is None
    data = result['data']
    with open(config['REFERENCES_JSON'], 'r') as f:
        reference_json = [x for x in json.load(f) if x['Name'] == '1000G'][0]
        populations = list(reference_json['Samples'].keys())
    assert len(data) == len(populations)
    assert all(x in populations for x in data)

    for bad_url in ['genome_builds/GRCh37/references/SOMETHING_BAD/populations', 'genome_builds/SOMETHING_BAD/references/1000G/populations']:
        response = client.get(bad_url)
        assert response.status_code == 404
        result = response.get_json()
        assert all(x in result for x in ['data', 'error'])
        assert result['error'] is not None
        assert result['data'] is None


def test_reference_population(client, config):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/EUR')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is None
    data = result['data']
    assert all(x in data for x in ['name', 'size'])
    with open(config['REFERENCES_JSON'], 'r') as f:
        reference_json = [x for x in json.load(f) if x['Name'] == '1000G'][0]
    assert len(reference_json['Samples']['EUR']) == data['size']

    for bad_url in ['genome_builds/GRCh37/references/SOMETHING_BAD/populations/EUR',
                    'genome_builds/SOMETHING_BAD/references/1000G/populations/EUR',
                    '/genome_builds/GRCh37/references/1000G/populations/SOMETHING_BAD']:
        response = client.get(bad_url)
        assert response.status_code == 404
        result = response.get_json()
        assert all(x in result for x in ['data', 'error'])
        assert result['error'] is not None
        assert result['data'] is None


def test_chromosomes(client):
    response = client.get('/genome_builds/GRCh37/references/1000G/chromosomes')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is None
    data = result['data']
    assert len(data) > 0

    for bad_url in ['/genome_builds/SOMETHING_BAD/references/1000G/chromosomes',
                    '/genome_builds/GRCh37/references/SOMETHING_BAD/chromosomes']:
        response = client.get(bad_url)
        assert response.status_code == 404
        result = response.get_json()
        assert all(x in result for x in ['data', 'error'])
        assert result['error'] is not None
        assert result['data'] is None

def test_region_ld_classic(client, goldstandard_ld):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare')
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
    for i in range(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        print((key, goldstandard[key], data['correlation'][i]))
        assert pytest.approx(data['correlation'][i], 0.00001) == goldstandard[key]

    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=r')
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
    for i in range(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        print((key, goldstandard[key], data['correlation'][i]**2))
        assert pytest.approx(data['correlation'][i]**2, 0.00001) == goldstandard[key]

    response = client.get('/genome_builds/GRCh37/references/1000G/populations/AFR/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare')
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
    for i in range(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        assert pytest.approx(data['correlation'][i], 0.00001) == goldstandard[key]


def test_region_ld_empty_classic(client):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=10&start=51241101&stop=51241385&correlation=rsquare')
    assert response.status_code == 200
    result = response.get_json()
    assert result['error'] is None
    data = result['data']
    assert all(x in data for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])
    assert all(len(data[x]) == 0 for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])


def test_region_ld_with_paging_classic(client):
    url = '/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare'
    response = client.get(url)
    assert response.status_code == 200
    data = response.get_json()['data']
    single_page_result = [x for x in zip(data['variant1'], data['variant2'], data['correlation'])]
    url += '&limit=2'
    multi_page_result = []
    while True:
        response = client.get(url)
        assert response.status_code == 200
        result = response.get_json()
        data = result['data']
        multi_page_result.extend([x for x in zip(data['variant1'], data['variant2'], data['correlation'])])
        url = result['next']
        if url is None:
            break
    assert len(single_page_result) == len(multi_page_result)
    assert all(x == y for x, y in zip(single_page_result, multi_page_result))


def test_region_ld_compact(client, goldstandard_ld):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    result = response.get_json()
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['variants', 'chromosomes', 'positions', 'offsets', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'positions', 'offsets', 'correlations'])
    n_correlations = sum(len(x) for x in data['correlations'])
    goldstandard = goldstandard_ld('../data/region_ld_22_51241101_51241385.hap.ld')
    assert len(goldstandard) == n_correlations
    for i1, correlations1 in enumerate(data['correlations']):
        offset1 = data['offsets'][i1]
        if len(correlations1) == 0:
            assert offset1 < 0
            continue
        assert (offset1 >= i1) and (offset1 < len(data['variants']))
        for i2, value in enumerate(correlations1, offset1):
            key = str(data['positions'][i1]) + '_' + str(data['positions'][i2])
            assert key in goldstandard
            assert pytest.approx(value, 0.00001) == goldstandard[key]

    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare_approx&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    result = response.get_json()
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['variants', 'chromosomes', 'positions', 'offsets', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'positions', 'offsets', 'correlations'])
    n_correlations = sum(len(x) for x in data['correlations'])
    goldstandard = goldstandard_ld('../data/region_ld_22_51241101_51241385.hap.ld')
    assert len(goldstandard) == n_correlations
    for i1, correlations1 in enumerate(data['correlations']):
        offset1 = data['offsets'][i1]
        if len(correlations1) == 0:
            assert offset1 < 0
            continue
        assert (offset1 >= i1) and (offset1 < len(data['variants']))
        for i2, value in enumerate(correlations1, offset1):
            key = str(data['positions'][i1]) + '_' + str(data['positions'][i2])
            print(key, value, goldstandard[key])
            assert key in goldstandard
            assert pytest.approx(value, 0.00001) == goldstandard[key]

    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=r&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    result = response.get_json()
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['variants', 'chromosomes', 'positions', 'offsets', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'offsets', 'positions', 'correlations'])
    n_correlations = sum(len(x) for x in data['correlations'])
    goldstandard = goldstandard_ld('../data/region_ld_22_51241101_51241385.hap.ld')
    assert len(goldstandard) == n_correlations
    for i1, correlations1 in enumerate(data['correlations']):
        offset1 = data['offsets'][i1]
        if len(correlations1) == 0:
            assert offset1 < 0
            continue
        assert (offset1 >= i1) and (offset1 < len(data['variants']))
        for i2, value in enumerate(correlations1, offset1):
            key = str(data['positions'][i1]) + '_' + str(data['positions'][i2])
            assert key in goldstandard
            assert pytest.approx(value**2, 0.00001) == goldstandard[key]

    response = client.get('/genome_builds/GRCh37/references/1000G/populations/AFR/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    result = response.get_json()
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['variants', 'chromosomes', 'positions', 'offsets', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'offsets', 'positions', 'correlations'])
    n_correlations = sum(len(x) for x in data['correlations'])
    goldstandard = goldstandard_ld('../data/region_ld_22_51241101_51241385.AFR.hap.ld')
    assert len(goldstandard) == n_correlations
    for i1, correlations1 in enumerate(data['correlations']):
        offset1 = data['offsets'][i1]
        if len(correlations1) == 0:
            assert offset1 < 0
            continue
        assert (offset1 >= i1) and (offset1 < len(data['variants']))
        for i2, value in enumerate(correlations1, offset1):
            key = str(data['positions'][i1]) + '_' + str(data['positions'][i2])
            assert key in goldstandard
            assert pytest.approx(value, 0.00001) == goldstandard[key]


def test_region_ld_with_msgpack_compact(client, goldstandard_ld):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare&msgpack=1&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/msgpack'
    result = msgpack.unpackb(response.get_data(), strict_map_key = False)
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['variants', 'chromosomes', 'positions', 'offsets', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'positions', 'offsets', 'correlations'])
    n_correlations = sum(len(x) for x in data['correlations'])
    goldstandard = goldstandard_ld('../data/region_ld_22_51241101_51241385.hap.ld')
    assert len(goldstandard) == n_correlations
    for i1, correlations1 in enumerate(data['correlations']):
        offset1 = data['offsets'][i1]
        if len(correlations1) == 0:
            assert offset1 < 0
            continue
        assert (offset1 >= i1) and (offset1 < len(data['variants']))
        for i2, value in enumerate(correlations1, offset1):
            key = str(data['positions'][i1]) + '_' + str(data['positions'][i2])
            assert key in goldstandard
            print((key, goldstandard[key], value))
            assert pytest.approx(value, 0.00001) == goldstandard[key]

    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=r&msgpack=1&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/msgpack'
    result = msgpack.unpackb(response.get_data(), strict_map_key = False)
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['variants', 'chromosomes', 'positions', 'offsets', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'offsets', 'positions', 'correlations'])
    n_correlations = sum(len(x) for x in data['correlations'])
    goldstandard = goldstandard_ld('../data/region_ld_22_51241101_51241385.hap.ld')
    assert len(goldstandard) == n_correlations
    for i1, correlations1 in enumerate(data['correlations']):
        offset1 = data['offsets'][i1]
        if len(correlations1) == 0:
            assert offset1 < 0
            continue
        assert (offset1 >= i1) and (offset1 < len(data['variants']))
        for i2, value in enumerate(correlations1, offset1):
            key = str(data['positions'][i1]) + '_' + str(data['positions'][i2])
            assert key in goldstandard
            print((key, goldstandard[key], value))
            assert pytest.approx(value**2, 0.00001) == goldstandard[key]

    response = client.get('/genome_builds/GRCh37/references/1000G/populations/AFR/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare&msgpack=1&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/msgpack'
    result = msgpack.unpackb(response.get_data(), strict_map_key = False)
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['variants', 'chromosomes', 'positions', 'offsets', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'offsets', 'positions', 'correlations'])
    n_correlations = sum(len(x) for x in data['correlations'])
    goldstandard = goldstandard_ld('../data/region_ld_22_51241101_51241385.AFR.hap.ld')
    assert len(goldstandard) == n_correlations
    for i1, correlations1 in enumerate(data['correlations']):
        offset1 = data['offsets'][i1]
        if len(correlations1) == 0:
            assert offset1 < 0
            continue
        assert (offset1 >= i1) and (offset1 < len(data['variants']))
        for i2, value in enumerate(correlations1, offset1):
            key = str(data['positions'][i1]) + '_' + str(data['positions'][i2])
            assert key in goldstandard
            print((key, goldstandard[key], value))
            assert pytest.approx(value, 0.00001) == goldstandard[key]


def test_region_ld_empty_compact(client):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=10&start=51241101&stop=51241385&correlation=rsquare&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    result = response.get_json()
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['variants', 'chromosomes', 'positions',  'offsets', 'correlations'])
    assert all(len(data[x]) == 0 for x in ['variants', 'chromosomes', 'positions', 'offsets', 'correlations'])


def test_region_ld_with_paging_compact(client):
    url = '/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare&format=compact'
    response = client.get(url)
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    data = response.get_json()['data']
    single_page_result = {}
    for i1, correlations1 in enumerate(data['correlations']):
        offset1 = data['offsets'][i1]
        if len(correlations1) == 0:
            continue
        for i2, value in enumerate(correlations1, offset1):
            key = str(data['positions'][i1]) + '_' + str(data['positions'][i2])
            assert key not in single_page_result
            single_page_result[key] = value
    url += '&limit=2'
    multi_page_result = {}
    while True:
        response = client.get(url)
        assert response.status_code == 200
        assert response.content_type == 'application/json'
        result = response.get_json()
        data = result['data']
        for i1, correlations1 in enumerate(data['correlations']):
            offset1 = data['offsets'][i1]
            if len(correlations1) == 0:
                continue
            for i2, value in enumerate(correlations1, offset1):
                key = str(data['positions'][i1]) + '_' + str(data['positions'][i2])
                assert key not in multi_page_result
                multi_page_result[key] = value
        url = result['next']
        if url == '':
            break
    assert len(single_page_result) == len(multi_page_result)
    for x, y in single_page_result.items():
        assert x in multi_page_result
        assert y == multi_page_result[x]

def test_variant_ld_classic(client, goldstandard_ld):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare')
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
    for i in range(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        assert pytest.approx(data['correlation'][i], 0.00001) == goldstandard[key]

    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=r')
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
    for i in range(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        assert pytest.approx(data['correlation'][i]**2, 0.00001) == goldstandard[key]

def test_variant_ld_compact(client, goldstandard_ld):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    result = response.get_json()
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['index_variant', 'index_chromosome', 'index_position', 'variants', 'chromosomes', 'positions', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'positions', 'correlations'])
    n_correlations = len(data['correlations'])
    goldstandard = goldstandard_ld('../data/variant_ld_22_51241101_vs_51241101_51241385.hap.ld')
    assert len(goldstandard) == n_correlations
    for i, correlation in enumerate(data['correlations']):
        key = str(data['index_position']) + '_' + str(data['positions'][i])
        assert key in goldstandard
        assert pytest.approx(correlation, 0.00001) == goldstandard[key]

    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=r&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    result = response.get_json()
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['index_variant', 'index_chromosome', 'index_position', 'variants', 'chromosomes', 'positions', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'positions', 'correlations'])
    n_correlations = len(data['correlations'])
    goldstandard = goldstandard_ld('../data/variant_ld_22_51241101_vs_51241101_51241385.hap.ld')
    assert len(goldstandard) == n_correlations
    for i, correlation in enumerate(data['correlations']):
        key = str(data['index_position']) + '_' + str(data['positions'][i])
        assert key in goldstandard
        assert pytest.approx(correlation**2, 0.00001) == goldstandard[key]


def test_variant_ld_with_msgpack(client, goldstandard_ld):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare&msgpack=1&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/msgpack'
    result = msgpack.unpackb(response.get_data(), strict_map_key = False)
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['index_variant', 'index_chromosome', 'index_position', 'variants', 'chromosomes', 'positions', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'positions', 'correlations'])
    n_correlations = len(data['correlations'])
    goldstandard = goldstandard_ld('../data/variant_ld_22_51241101_vs_51241101_51241385.hap.ld')
    assert len(goldstandard) == n_correlations
    for i, correlation in enumerate(data['correlations']):
        key = str(data['index_position']) + '_' + str(data['positions'][i])
        assert key in goldstandard
        assert pytest.approx(correlation, 0.00001) == goldstandard[key]

    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=r&msgpack=1&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/msgpack'
    result = msgpack.unpackb(response.get_data(), strict_map_key = False)
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['index_variant', 'index_chromosome', 'index_position', 'variants', 'chromosomes', 'positions', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'positions', 'correlations'])
    n_correlations = len(data['correlations'])
    goldstandard = goldstandard_ld('../data/variant_ld_22_51241101_vs_51241101_51241385.hap.ld')
    assert len(goldstandard) == n_correlations
    for i, correlation in enumerate(data['correlations']):
        key = str(data['index_position']) + '_' + str(data['positions'][i])
        assert key in goldstandard
        assert pytest.approx(correlation**2, 0.00001) == goldstandard[key]


def test_variant_ld_precision_compact(client, goldstandard_ld):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare&precision=7&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    result = response.get_json()
    assert all(x in result for x in ['data', 'error', 'next'])
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['index_variant', 'index_chromosome', 'index_position', 'variants', 'chromosomes', 'positions', 'correlations'])
    assert all(len(data['variants']) == len(data[x]) for x in ['chromosomes', 'positions', 'correlations'])
    n_correlations = len(data['correlations'])
    goldstandard = goldstandard_ld('../data/variant_ld_22_51241101_vs_51241101_51241385.hap.ld')
    assert len(goldstandard) == n_correlations
    for i, correlation in enumerate(data['correlations']):
        key = str(data['index_position']) + '_' + str(data['positions'][i])
        assert key in goldstandard
        assert pytest.approx(correlation, 0.00001) == round(goldstandard[key], 7)


def test_variant_ld_precision_server_side(client):
    client.application.config["LDSERVER_PRECISION"] = 1
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare&precision=7&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    result = response.get_json()
    data = result['data']
    assert data["correlations"] == [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

    client.application.config["LDSERVER_PRECISION"] = 8
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare&precision=7&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    result = response.get_json()
    data = result['data']
    assert data["correlations"] == [1.0, 1.0, 7.59e-06, 2.8e-07, 2.8e-07, 2.8e-07]


def test_variant_ld_different_chromosomes_compact(client):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=21&start=51241101&stop=51241385&correlation=rsquare&format=compact')
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    result = response.get_json()
    assert result['error'] is ''
    assert result['next'] is ''
    data = result['data']
    assert all(x in data for x in ['index_variant', 'index_chromosome', 'index_position', 'variants', 'chromosomes', 'positions', 'correlations'])
    assert all(len(data[x]) == 0 for x in ['variants', 'chromosomes', 'positions', 'correlations'])
    assert data['index_variant'] is None
    assert data['index_chromosome'] is None
    assert data['index_position'] is None


def test_variant_ld_with_paging_compact(client):
    url = '/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare&format=compact'
    response = client.get(url)
    assert response.status_code == 200
    assert response.content_type == 'application/json'
    data = response.get_json()['data']
    single_page_result = {}
    for i, correlation in enumerate(data['correlations']):
        key = str(data['index_position']) + '_' + str(data['positions'][i])
        single_page_result[key] = correlation

    url += '&limit=2'
    multi_page_result = {}
    while True:
        response = client.get(url)
        assert response.status_code == 200
        assert response.content_type == 'application/json'
        result = response.get_json()
        data = result['data']
        for i, correlation in enumerate(data['correlations']):
            key = str(data['index_position']) + '_' + str(data['positions'][i])
            assert key not in multi_page_result
            multi_page_result[key] = correlation
        url = result['next']
        if url == '':
            break
    assert len(single_page_result) == len(multi_page_result)
    for x, y in single_page_result.items():
        assert x in multi_page_result
        assert y == multi_page_result[x]

def test_variant_ld_different_chromosomes_classic(client):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=21&start=51241101&stop=51241385&correlation=rsquare')
    assert response.status_code == 200
    result = response.get_json()
    assert result['error'] is None
    data = result['data']
    assert all(x in data for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])
    assert all(len(data[x]) == 0 for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])


def test_variant_ld_with_paging_classic(client):
    url = '/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare'
    response = client.get(url)
    assert response.status_code == 200
    data = response.get_json()['data']
    single_page_result = [x for x in zip(data['variant1'], data['variant2'], data['correlation'])]
    url += '&limit=2'
    multi_page_result = []
    while True:
        response = client.get(url)
        assert response.status_code == 200
        result = response.get_json()
        data = result['data']
        multi_page_result.extend([x for x in zip(data['variant1'], data['variant2'], data['correlation'])])
        url = result['next']
        if url is None:
            break
    assert len(single_page_result) == len(multi_page_result)
    assert all(x == y for x, y in zip(single_page_result, multi_page_result))

def test_compression(client):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare', headers=[('Accept-Encoding', 'gzip')])
    assert response.headers['Content-Encoding'] == 'gzip'
    compressed_payload = io.BytesIO(response.data)
    uncompressed_payload = gzip.GzipFile(fileobj = compressed_payload, mode = 'rb').read()
    assert len(uncompressed_payload) > 0
    result = json.loads(uncompressed_payload)
    assert all(x in result for x in ['data', 'error', 'next'])


def test_malformed_region_ld(client):
    for bad_url in ['/genome_builds/SOMETHING_BAD/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare',
                    '/genome_builds/GRCh37/references/SOMETHING_BAD/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/SOMETHING_BAD/regions?chrom=22&start=51241101&stop=51241385&correlation=rsquare']:
        response = client.get(bad_url)
        assert response.status_code == 404
        result = response.get_json()
        assert all(x in result for x in ['data', 'error'])
        assert result['data'] is None
        assert result['error'] is not None and len(result['error']) > 0

    for bad_url in ['/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=22&start=51241101&stop=51241385&variant=10:114758349_C/T&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/regions?start=100&stop=1000&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=chr22&start=-100&stop=1000&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=chr22&start=100&stop=-1000&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=chr22&start=1000&stop=100&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=chr22&start=100&stop=1000&correlation=nonsense',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=chr22&start=100&stop=1000&correlation=rsquare&limit=0']:
        response = client.get(bad_url)
        assert response.status_code == 422
        result = response.get_json()
        assert all(x in result for x in ['data', 'error'])
        assert result['data'] is None
        assert result['error'] is not None and len(result['error']) > 0


def test_malformed_variant_ld(client):
    for bad_url in ['/genome_builds/SOMETHING_BAD/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare',
                    '/genome_builds/GRCh37/references/SOMETHING_BAD/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/SOMETHING_BAD/variants?variant=22:51241101_A/T&chrom=22&start=51241101&stop=51241385&correlation=rsquare']:
        response = client.get(bad_url)
        assert response.status_code == 404
        result = response.get_json()
        assert all(x in result for x in ['data', 'error'])
        assert result['data'] is None
        assert result['error'] is not None and len(result['error']) > 0

    for bad_url in ['/genome_builds/GRCh37/references/1000G/populations/ALL/variants?chrom=22&start=51241101&stop=51241385&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&start=51241101&stop=51241385&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=-1000&stop=51241385&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=1000&stop=-100&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=1000&stop=100&correlation=rsquare',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=1000&stop=10000&correlation=nonsense',
                    '/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=22&start=1000&stop=10000&correlation=rsquare&limit=0']:
        response = client.get(bad_url)
        assert response.status_code == 422
        result = response.get_json()
        assert all(x in result for x in ['data', 'error'])
        assert result['data'] is None
        assert result['error'] is not None and len(result['error']) > 0