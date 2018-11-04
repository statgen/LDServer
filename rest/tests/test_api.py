import pytest
import json
import StringIO
import gzip


def test_correlations(client):
    response = client.get('/correlations')
    assert response.status_code == 200
    result = response.get_json()
    assert all(x in result for x in ['data', 'error'])
    assert result['error'] is None
    data = result['data']
    assert len(data) == 3
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
        populations = reference_json['Samples'].keys()
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


def test_region_ld(client, goldstandard_ld):
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
    for i in xrange(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        print key, goldstandard[key], data['correlation'][i]
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
    for i in xrange(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        print key, goldstandard[key], data['correlation'][i]**2
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
    for i in xrange(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        assert pytest.approx(data['correlation'][i], 0.00001) == goldstandard[key]


def test_region_ld_empty(client):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/regions?chrom=10&start=51241101&stop=51241385&correlation=rsquare')
    assert response.status_code == 200
    result = response.get_json()
    assert result['error'] is None
    data = result['data']
    assert all(x in data for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])
    assert all(len(data[x]) == 0 for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])


def test_region_ld_with_paging(client):
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


def test_variant_ld(client, goldstandard_ld):
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
    for i in xrange(0, len(goldstandard)):
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
    for i in xrange(0, len(goldstandard)):
        key = str(data['position1'][i]) + '_' + str(data['position2'][i])
        assert key in goldstandard
        assert pytest.approx(data['correlation'][i]**2, 0.00001) == goldstandard[key]


def test_variant_ld_different_chromosomes(client):
    response = client.get('/genome_builds/GRCh37/references/1000G/populations/ALL/variants?variant=22:51241101_A/T&chrom=21&start=51241101&stop=51241385&correlation=rsquare')
    assert response.status_code == 200
    result = response.get_json()
    assert result['error'] is None
    data = result['data']
    assert all(x in data for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])
    assert all(len(data[x]) == 0 for x in ['chromosome1', 'position1', 'variant1', 'chromosome2', 'position2', 'variant2', 'correlation'])


def test_variant_ld_with_paging(client):
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
    compressed_payload = StringIO.StringIO(response.data)
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