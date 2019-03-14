from flask import current_app, Blueprint, request, jsonify, make_response, abort
from flask_cors import CORS
from flask_compress import Compress
from webargs.flaskparser import parser
from webargs import fields, ValidationError
from functools import partial
from model import Reference, File, Sample
import model
from core.pywrapper import LDServer, LDQueryResult, StringVec, correlation
import time

API_VERSION = "1.0"

bp = Blueprint('api', __name__)
CORS(bp)

compress = Compress()

@parser.error_handler
def handle_parsing_error(error, request, schema):
    for field, message in error.messages.iteritems():
        message = 'Error while parsing \'{}\' query parameter: {}'.format(field, message[0])
        break
    response = jsonify({'data': None, 'error': message })
    response.status_code = 422
    abort(response)


def validate_query(parsed_fields, all_fields):
    for key, value in request.args.iteritems():
        if key not in all_fields:
            raise ValidationError({key: ['Unknown parameter.']})
    if 'start' in parsed_fields and 'stop' in parsed_fields:
        if parsed_fields['stop'] <= parsed_fields['start']:
            raise ValidationError({'start': ['Start position must be greater than stop position.']})
    return True


@bp.route('/correlations', methods = ['GET'])
def get_correlations():
    response = { 'data': model.get_correlations(), 'error': None }
    return make_response(jsonify(response), 200)


@bp.route('/genome_builds', methods = ['GET'])
def get_genome_build():
    data = model.get_genome_builds()
    response = { 'data': data, 'error': None}
    return make_response(jsonify(response), 200)


@bp.route('/genome_builds/<genome_build>/references', methods = ['GET'])
def get_references(genome_build):
    response = { 'data': None, 'error': None }
    if not model.has_genome_build(genome_build):
        response['error'] = 'Genome build \'{}\' was not found.'.format(genome_build)
    else:
        response['data'] = model.get_references(genome_build)
    return make_response(jsonify(response), 200 if response['data'] is not None else 404)


@bp.route('/genome_builds/<genome_build>/references/<reference_name>', methods = ['GET'])
def get_reference(genome_build, reference_name):
    response = { 'data': None, 'error': None }
    if not model.has_genome_build(genome_build):
        response['error'] = 'Genome build \'{}\' was not found.'.format(genome_build)
    else:
        data = model.get_reference(genome_build, reference_name)
        if data is None:
            response['error'] = 'Reference panel \'{}\' was not found in {} genome build.'.format(reference_name, genome_build)
        else:
            response['data'] = data
    return make_response(jsonify(response), 200 if response['data'] is not None else 404)


@bp.route('/genome_builds/<genome_build>/references/<reference_name>/populations', methods = ['GET'])
def get_populations(genome_build, reference_name):
    response = { 'data': None, 'error': None }
    if not model.has_genome_build(genome_build):
        response['error'] = 'Genome build \'{}\' was not found.'.format(genome_build)
    else:
        reference_id = model.get_reference_id(genome_build, reference_name)
        if not reference_id:
            response['error'] = 'Reference panel \'{}\' was not found in {} genome build.'.format(reference_name, genome_build)
        else:
            response['data'] = model.get_populations(reference_id)
    return make_response(jsonify(response), 200 if response['data'] is not None else 404)


@bp.route('/genome_builds/<genome_build>/references/<reference_name>/populations/<population_name>', methods = ['GET'])
def get_population(genome_build, reference_name, population_name):
    response = { 'data': None, 'error': None }
    if not model.has_genome_build(genome_build):
        response['error'] = 'Genome build \'{}\' was not found.'.format(genome_build)
    else:
        reference_id = model.get_reference_id(genome_build, reference_name)
        if not reference_id:
            response['error'] = 'Reference panel \'{}\' was not found in {} genome build.'.format(reference_name, genome_build)
        elif not model.has_population(reference_id, population_name):
            response['error'] = 'Population \'{}\' was not found in {} reference panel.'.format(population_name, reference_name)
        else:
            response['data'] = { 'name': population_name, 'size': model.get_samples_count(reference_id, population_name) }
    return make_response(jsonify(response), 200 if response['data'] is not None else 404)


@bp.route('/genome_builds/<genome_build>/references/<reference_name>/chromosomes', methods = ['GET'])
def get_chromosomes(genome_build, reference_name):
    response = { 'data': None, 'error': None }
    if not model.has_genome_build(genome_build):
        response['error'] = 'Genome build \'{}\' was not found.'.format(genome_build)
    else:
        reference_id = model.get_reference_id(genome_build, reference_name)
        if not reference_id:
            response['error'] = 'Reference panel \'{}\' was not found in {} genome build.'.format(reference_name, genome_build)
        else:
            ldserver = LDServer(current_app.config['SEGMENT_SIZE_BP'])
            for f in model.get_files(reference_id):
                ldserver.set_file(f)
            response['data'] = [x for x in ldserver.get_chromosomes()]
    return make_response(jsonify(response), 200 if response['data'] is not None else 404)

def correlation_type(correlation_name):
    return { 'r': correlation.ld_r, 'rsquare': correlation.ld_rsquare, 'cov': correlation.cov }[correlation_name]

@bp.route('/genome_builds/<genome_build>/references/<reference_name>/populations/<population_name>/regions', methods = ['GET'])
def get_region_ld(genome_build, reference_name, population_name):
    arguments = {
        'chrom': fields.Str(required = True, validate = lambda x: len(x) > 0, error_messages = {'validator_failed': 'Value must be a non-empty string.'}),
        'start': fields.Int(required = True, validate = lambda x: x >= 0, error_messages = {'validator_failed': 'Value must be greater than or equal to 0.'}),
        'stop': fields.Int(required = True, validate = lambda x: x > 0, error_messages = {'validator_failed': 'Value must be greater than 0.'}),
        'correlation': fields.Str(required = True, validate = lambda x: x in current_app.config['CORRELATIONS'], error_messages = {'validator_failed': 'Value must be one of the following: {}.'.format(', '.join(current_app.config['CORRELATIONS']))}),
        'limit': fields.Int(required = False, validate = lambda x: x > 0, missing = current_app.config['API_MAX_PAGE_SIZE'], error_messages = {'validator_failed': 'Value must be greater than 0.'}),
        'last': fields.Str(required = False, validate = lambda x: len(x) > 0, error_messages = {'validator_failed': 'Value must be a non-empty string.'})
    }
    args = parser.parse(arguments, validate = partial(validate_query, all_fields = ['chrom', 'start', 'stop', 'correlation', 'limit', 'last']))
    if args['limit'] > current_app.config['API_MAX_PAGE_SIZE']:
        args['limit'] = current_app.config['API_MAX_PAGE_SIZE']
    if not model.has_genome_build(genome_build):
        response = { 'data': None, 'error': 'Genome build \'{}\' was not found.'.format(genome_build) }
        return make_response(jsonify(response), 404)
    reference_id = model.get_reference_id(genome_build, reference_name)
    if not reference_id:
        response = { 'data': None, 'error': 'Reference panel \'{}\' was not found in {} genome build.'.format(reference_name, genome_build) }
        return make_response(jsonify(response), 404)
    if not model.has_samples(reference_id, population_name):
        response = { 'data': None, 'error': 'Population \'{}\' was not found in {} reference panel.'.format(population_name, reference_name) }
        return make_response(jsonify(response), 404)
    ldserver = LDServer(current_app.config['SEGMENT_SIZE_BP'])
    for f in model.get_files(reference_id):
        ldserver.set_file(f)
    if 'last' in args:
        result = LDQueryResult(args['limit'], str(args['last']))
    else:
        result = LDQueryResult(args['limit'])
    if population_name != 'ALL':
        s = StringVec()
        s.extend(model.get_samples(reference_id, population_name))
        ldserver.set_samples(str(population_name), s)
    if current_app.config['CACHE_ENABLED']:
        ldserver.enable_cache(reference_id, current_app.config['CACHE_REDIS_HOSTNAME'], current_app.config['CACHE_REDIS_PORT'])
    #start = time.time()
    ldserver.compute_region_ld(str(args['chrom']), args['start'], args['stop'], correlation_type(args['correlation']), result, str(population_name))
    #print "Computed results in {} seconds.".format("%0.4f" % (time.time() - start))
    #start = time.time()
    if current_app.config['PROXY_PASS']:
        base_url = '/'.join(x.strip('/') for x in [current_app.config['PROXY_PASS'], request.path])
    else:
        base_url = request.base_url
    base_url += '?' + '&'.join(('{}={}'.format(arg, value) for arg, value in request.args.iteritems(True) if arg != 'last'))
    #print "Jsonified result in {} seconds.".format("%0.4f" % (time.time() - start))
    #start = time.time()
    r = make_response(result.get_json(str(base_url)), 200)
    r.mimetype = 'application/json'
    #print "Response created in {} seconds.".format("%0.4f" % (time.time() - start))
    return r


@bp.route('/genome_builds/<genome_build>/references/<reference_name>/populations/<population_name>/variants', methods = ['GET'])
def get_variant_ld(genome_build, reference_name, population_name):
    arguments = {
        'variant': fields.Str(required = True, validate = lambda x: len(x) > 0, error_messages = {'validator_failed': 'Value must be a non-empty string.'}),
        'chrom': fields.Str(required = True, validate = lambda x: len(x) > 0, error_messages = {'validator_failed': 'Value must be a non-empty string.'}),
        'start': fields.Int(required = True, validate = lambda x: x >= 0, error_messages = {'validator_failed': 'Value must be greater than or equal to 0.'}),
        'stop': fields.Int(required = True, validate = lambda x: x > 0, error_messages = {'validator_failed': 'Value must be greater than 0.'}),
        'correlation': fields.Str(required = True, validate = lambda x: x in current_app.config['CORRELATIONS'], error_messages = {'validator_failed': 'Value must be one of the following: {}.'.format(', '.join(current_app.config['CORRELATIONS']))}),
        'limit': fields.Int(required = False, validate = lambda x: x > 0, missing = current_app.config['API_MAX_PAGE_SIZE'], error_messages = {'validator_failed': 'Value must be greater than 0.'}),
        'last': fields.Str(required = False, validate = lambda x: len(x) > 0, error_messages = {'validator_failed': 'Value must be a non-empty string.'})
    }
    args = parser.parse(arguments, validate = partial(validate_query, all_fields = ['variant', 'chrom', 'start', 'stop', 'correlation', 'limit', 'last']))
    if args['limit'] > current_app.config['API_MAX_PAGE_SIZE']:
        args['limit'] = current_app.config['API_MAX_PAGE_SIZE']
    if not model.has_genome_build(genome_build):
        response = { 'data': None, 'error': 'Genome build \'{}\' was not found.'.format(genome_build) }
        return make_response(jsonify(response), 404)
    reference_id = model.get_reference_id(genome_build, reference_name)
    if not reference_id:
        response = { 'data': None, 'error': 'Reference panel \'{}\' was not found in {} genome build.'.format(reference_name, genome_build) }
        return make_response(jsonify(response), 404)
    if not model.has_samples(reference_id, population_name):
        response = { 'data': None, 'error': 'Population \'{}\' was not found in {} reference panel.'.format(population_name, reference_name) }
        return make_response(jsonify(response), 404)
    ldserver = LDServer(current_app.config['SEGMENT_SIZE_BP'])
    for f in model.get_files(reference_id):
        ldserver.set_file(f)
    if 'last' in args:
        result = LDQueryResult(args['limit'], str(args['last']))
    else:
        result = LDQueryResult(args['limit'])
    if population_name != 'ALL':
        s = StringVec()
        s.extend(model.get_samples(reference_id, population_name))
        ldserver.set_samples(str(population_name), s)
    if current_app.config['CACHE_ENABLED']:
        ldserver.enable_cache(reference_id, current_app.config['CACHE_REDIS_HOSTNAME'], current_app.config['CACHE_REDIS_PORT'])
    start = time.time()
    ldserver.compute_variant_ld(str(args['variant']), str(args['chrom']), args['start'], args['stop'], correlation_type(args['correlation']), result, str(population_name))
    print "Computed results in {} seconds.".format("%0.4f" % (time.time() - start))
    if current_app.config['PROXY_PASS']:
        base_url = '/'.join(x.strip('/') for x in [current_app.config['PROXY_PASS'], request.path])
    else:
        base_url = request.base_url
    base_url += '?' + '&'.join(('{}={}'.format(arg, value) for arg, value in request.args.iteritems(True) if arg != 'last'))
    r = make_response(result.get_json(str(base_url)), 200)
    r.mimetype = 'application/json'
    return r

# TODO: LD between arbitrary variants
# @bp.route('/<reference>/<population>/ld/cartesian', methods = ['GET'])
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
