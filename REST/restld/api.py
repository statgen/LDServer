from flask import Flask, request, jsonify, make_response, abort
from webargs.flaskparser import parser
from webargs import fields
import PyLDServer

# app = Flask(__name__, instance_relative_config = True)

# app.config.from_object('config.default')
# app.config.from_pyfile('config.py', silent = True)
# app.config.from_envvar('RESTLD_CONFIG_FILE', silent = True)


API_VERSION = '1.0'
API_MAX_PAGE_SIZE = 1000


datasets = dict()
for dataset in app.config['DATASETS']:
    datasets[dataset['Name']] = { 'LD': PyLDServer.LDServer(), 'Meta': dataset } 
    datasets[dataset['Name']]['LD'].set_file(dataset['Datafile'])


@app.route('/')
def get_all_datasets():
    response = { 'api_version': API_VERSION, 'references': [] }
    for dataset in app.config['DATASETS']:
        response['references'].append({ 'Name': dataset['Name'], 'Description': dataset['Description'], 'Genome build': dataset['Genome build'], 'Populations': dataset['Samples'].keys() })
    return make_response(jsonify(response), 200)


@app.route('/<reference>', methods = ['GET'])
def get_reference_info(reference):
    dataset = datasets.get(reference, None)
    if dataset:
        response = { 'name': dataset['Meta']['Name'], 'description': dataset['Meta']['Description'], 'genome build': dataset['Meta']['Genome build'], 'populations': dataset['Meta']['Samples'].keys() } 
        return make_response(jsonify(response), 200)
    abort(404)


@app.route('/<reference>/<population>', methods = ['GET'])
def get_population_info(reference, population):
    dataset = datasets.get(reference, None)
    if dataset is None:
        abort(404)
    samples = dataset['Meta']['Samples'].get(population, None)
    if samples is None:
        abort(404)
    response = { 'name': population, 'size': len(samples) }
    return make_response(jsonify(response), 200)


def validate_query(value):
    if 'start' in value and 'stop' in value:
        if value['stop'] <= value['start']:
            return False
    return True


def build_link_next(args, result):
    if result.has_next():
        link_next = request.base_url + '?' + '&'.join(('{}={}'.format(arg, value) for arg, value in request.args.iteritems(True) if arg != 'last'))
        link_next += '&last={}'.format(result.get_last())
        return link_next
    return None


@app.route('/<reference>/<population>/ld/region', methods = ['GET'])
def get_region_ld(reference, population):
    arguments = {
        'chrom': fields.Str(required = True, validate = lambda x: len(x) > 0),
        'start': fields.Int(required = True, validate = lambda x: x >= 0),
        'stop': fields.Int(required = True, validate = lambda x: x > 0),
        'limit': fields.Int(required = False, validate = lambda x: x > 0, missing = API_MAX_PAGE_SIZE),
        'last': fields.Str(required = False, validate = lambda x: len(x) > 0)
    }
    args = parser.parse(arguments, validate = validate_query)
    print args
    dataset = datasets.get(reference, None)
    if dataset is None:
        abort(404)
    if population not in dataset['Meta']['Samples']:
        abort(404)
    response = {}
    if 'last' in args:
        result = PyLDServer.LDQueryResult(args['limit'], str(args['last']))
    else:
        result = PyLDServer.LDQueryResult(args['limit'])
    dataset['LD'].compute_region_ld(str(args['chrom']), args['start'], args['stop'], result, str(population))
    response['data'] = {
        'chromosome1': [str(args['chrom'])] * len(result.data),
        'variant1': [None] * len(result.data),
        'position1': [None] * len(result.data),
        'chromosome2': [str(args['chrom'])] * len(result.data),
        'variant2': [None] * len(result.data),
        'position2': [None] * len(result.data),
        'r': [None] * len(result.data),
        'rsquare': [None] * len(result.data)
    }
    for i, pair in enumerate(result.data):
        response['data']['variant1'][i] = pair.variant1
        response['data']['position1'][i] = pair.position1
        response['data']['variant2'][i] = pair.variant2
        response['data']['position2'][i] = pair.position2
        response['data']['r'][i] = pair.r
        response['data']['rsquare'][i] = pair.rsquare
    response['next'] = build_link_next(args, result)
    return make_response(jsonify(response), 200)


@app.route('/<reference>/<population>/ld/variant', methods = ['GET'])
def get_variant_ld(reference, population):
    arguments = {
        'variant': fields.Str(required = True, validate = lambda x: len(x) > 0),
        'chrom': fields.Str(required = True, validate = lambda x: x > 0),
        'start': fields.Int(required = True, validate = lambda x: x >= 0),
        'stop': fields.Int(required = True, validate = lambda x: x > 0),
        'limit': fields.Int(required = False, validate = lambda x: x > 0, missing = API_MAX_PAGE_SIZE),
        'last': fields.Str(required = False, validate = lambda x: len(x) > 0)
    }
    args = parser.parse(arguments, validate = validate_query)
    dataset = datasets.get(reference, None)
    if dataset is None:
        abort(404)
    if population not in dataset['Meta']['Samples']:
        abort(404)
    response = {}
    if 'last' in args:
        result = PyLDServer.LDQueryResult(args['limit'], str(args['last']))
    else:
        result = PyLDServer.LDQueryResult(args['limit'])
    dataset['LD'].compute_variant_ld(str(args['variant']), str(args['chrom']), args['start'], args['stop'], result, str(population))
    response['data'] = {
        'chromosome1': [str(args['chrom'])] * len(result.data),
        'variant1': [None] * len(result.data),
        'position1': [None] * len(result.data),
        'chromosome2': [str(args['chrom'])] * len(result.data),
        'variant2': [None] * len(result.data),
        'position2': [None] * len(result.data),
        'r': [None] * len(result.data),
        'rsquare': [None] * len(result.data)
    }
    for i, pair in enumerate(result.data):
        response['data']['variant1'][i] = pair.variant1
        response['data']['position1'][i] = pair.position1
        response['data']['variant2'][i] = pair.variant2
        response['data']['position2'][i] = pair.position2
        response['data']['r'][i] = pair.r
        response['data']['rsquare'][i] = pair.rsquare
    response['next'] = build_link_next(args, result)
    return make_response(jsonify(response), 200)


# TODO: LD between arbitrary variants
# @app.route('/<reference>/<population>/ld/cartesian', methods = ['GET'])
# def get_ld(reference, population):
#     arguments = {
#         'variants1': fields.DelimitedList(fields.Str()),
#         'variants2': fields.DelimitedList(fields.Str()),
#         'limit': fields.Int(required = False, validate = lambda x: x > 0, missing = API_MAX_PAGE_SIZE),
#         'last': fields.Function(deserialize = deserialize_query_last)
#     }
#     args = parser.parse(arguments, validate = validate_query)
#     dataset = datasets.get(reference, None)
#     if dataset is None:
#         abort(404)
#     if population not in dataset['Meta']['Samples']:
#         abort(404)
#     response = {}
#     pairs = PyLDServer.Pairs()
#     # compute ld
#     response['data'] = {
#         'chromosome1': [str(args['chrom'])] * len(pairs),
#         'variant1': [None] * len(pairs),
#         'position1': [None] * len(pairs),
#         'chromosome2': [str(args['chrom'])] * len(pairs),
#         'variant2': [None] * len(pairs),
#         'position2': [None] * len(pairs),
#         'r': [None] * len(pairs),
#         'rsquare': [None] * len(pairs)
#     }
#     for i, pair in enumerate(pairs):
#         response['data']['variant1'][i] = pair.name1
#         response['data']['position1'][i] = pair.position1
#         response['data']['variant2'][i] = pair.name2
#         response['data']['position2'][i] = pair.position2
#         response['data']['r'][i] = pair.r
#         response['data']['rsquare'][i] = pair.rsquare
#     response['next'] = None
#     return make_response(jsonify(response), 200)


# if __name__ == '__main__':
#     app.run()
