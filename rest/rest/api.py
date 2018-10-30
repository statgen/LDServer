from flask import current_app, Blueprint, request, jsonify, make_response, abort
from flask_cors import CORS
from flask_compress import Compress
from webargs.flaskparser import parser
from webargs import fields, ValidationError
from model import Reference, File, Sample
import model
from ld.pywrapper import LDServer, LDQueryResult, StringVec, correlation
import time

API_VERSION = "1.0"

bp = Blueprint('api', __name__)
CORS(bp)

compress = Compress()

@parser.error_handler
def handle_parsing_error(error, request):
    for field, message in error.messages.iteritems():
        message = 'Error while parsing \'{}\' query parameter: {}'.format(field, message[0])
        print message
        break
    response = jsonify({'data': None, 'error': message })
    response.status_code = 422
    abort(response)


def validate_query(value):
    if 'start' in value and 'stop' in value:
        if value['stop'] <= value['start']:
            raise ValidationError({'start': ['Start position must be greater than stop position.']})
    return True


def build_link_next(args, result):
    if result.has_next():
        link_next = request.base_url + '?' + '&'.join(('{}={}'.format(arg, value) for arg, value in request.args.iteritems(True) if arg != 'last'))
        link_next += '&last={}'.format(result.get_last())
        return link_next
    return None


@bp.route('/correlations', methods = ['GET'])
def get_correlations():
    data = [
        { 'name': 'r', 'label': 'r', 'description': '', 'type': 'LD' },
        { 'name': 'rsquare', 'label': 'r^2', 'description': '', 'type': 'LD' },
        { 'name': 'covariance', 'label': 'cov', 'description': '', 'type': 'Covariance' }
    ]
    response = { 'data': data, 'error': None }
    return make_response(jsonify(response), 200)


@bp.route('/genome_builds', methods = ['GET'])
def get_genome_build():
    data = model.get_genome_builds()
    response = { 'data': data, 'error': None}
    return make_response(jsonify(response), 200)


@bp.route('/genome_builds/<genome_build>/references', methods = ['GET'])
def get_references(genome_build):
    data = model.get_references(genome_build)
    response = { 'data': data, 'error': None }
    return make_response(jsonify(response), 200)


@bp.route('/genome_builds/<genome_build>/references/<reference_name>', methods = ['GET'])
def get_reference(genome_build, reference_name):
    data = model.get_reference(genome_build, reference_name)
    response = { 'data': data, 'error': None }
    return make_response(jsonify(response), 200)
    #reference = Reference.query.filter_by(name = reference_name).first()
    #if reference:
    #    data = { 'name': reference.name, 'description': reference.description, 'genome build': reference.genome_build, 'populations': list(set([s.subset for s in reference.samples])) }
    #    response = { 'data': data, 'error': None }
    #    return make_response(jsonify(response), 200)
    #bort(404)


@bp.route('/genome_builds/<genome_build>/references/<reference_name>/populations', methods = ['GET'])
def get_populations(genome_build, reference_name):
    data = model.get_populations(genome_build, reference_name)
    response = { 'data': data, 'error': None }
    return make_response(jsonify(response), 200)
    # reference = Reference.query.filter_by(name = reference_name).first()
    # if reference:
    #     data = list(set([s.subset for s in reference.samples]))
    #     response = { 'data': data, 'error': None }
    #     return make_response(jsonify(response), 200)
    # abort(404)


@bp.route('/genome_builds/<genome_build>/references/<reference_name>/populations/<population_name>', methods = ['GET'])
def get_population(genome_build, reference_name, population_name):
    data = model.get_population(genome_build, reference_name, population_name)
    response = { 'data': data, 'error': None }
    return make_response(jsonify(response), 200)
    # reference = Reference.query.filter_by(name = reference_name).first()
    # if reference is not None:
    #     samples = Sample.query.with_parent(reference).filter_by(subset = population_name).all()
    #     if samples:
    #         data = { 'name': population_name, 'size': len(samples) }
    #         response = { 'data': data, 'error': None }
    #         return make_response(jsonify(response), 200)
    # abort(404)


@bp.route('/genome_builds/<genome_build>/references/<reference_name>/chromosomes', methods = ['GET'])
def get_chromosomes(genome_build, reference_name):
    reference = Reference.query.filter_by(genome_build = genome_build, name = reference_name).first()
    if reference is None:
        abort(404)
    files = File.query.with_parent(reference).all()
    if not files:
        abort(404)
    ldserver = LDServer(current_app.config['SEGMENT_SIZE_BP'])
    for file in files:
        ldserver.set_file(str(file.path))
    response = { 'data': [x for x in ldserver.get_chromosomes()], 'error': None }
    return make_response(jsonify(response), 200)


@bp.route('/genome_builds/<genome_build>/references/<reference_name>/populations/<population_name>/regions', methods = ['GET'])
def get_region_ld(genome_build, reference_name, population_name):
    arguments = {
        'chrom': fields.Str(required = True, validate = lambda x: len(x) > 0, error_messages = {'validator_failed': 'Value must be a non-empty string.'}),
        'start': fields.Int(required = True, validate = lambda x: x >= 0, error_messages = {'validator_failed': 'Value must be greater than or equal to 0.'}),
        'stop': fields.Int(required = True, validate = lambda x: x > 0, error_messages = {'validator_failed': 'Value must be greater than 0.'}),
        'correlation': fields.Str(required = True, validate = lambda x: x in ['r', 'rsquare'], error_messages = {'validator_failed': 'Value must be equal to \'r\' or \'rsquare\'.'}),
        'limit': fields.Int(required = False, validate = lambda x: x > 0, missing = current_app.config['API_MAX_PAGE_SIZE'], error_messages = {'validator_failed': 'Value must be greater than 0.'}),
        'last': fields.Str(required = False, validate = lambda x: len(x) > 0, error_messages = {'validator_failed': 'Value must be a non-empty string.'})
    }
    args = parser.parse(arguments, validate = validate_query)
    if args['limit'] > current_app.config['API_MAX_PAGE_SIZE']:
        args['limit'] = current_app.config['API_MAX_PAGE_SIZE']
    if args['correlation'] == 'r':
        correlation_type = correlation.ld_r
    elif args['correlation'] == 'rsquare':
        correlation_type = correlation.ld_rsquare
    else:
        abort(404)
    reference = Reference.query.filter_by(genome_build = genome_build, name = reference_name).first()
    if reference is None:
        abort(404)
    samples = Sample.query.with_parent(reference).filter_by(subset = population_name).all()
    if not samples:
        abort(404)
    files = File.query.with_parent(reference).all()
    if not files:
        abort(404)
    #start = time.time()
    ldserver = LDServer(current_app.config['SEGMENT_SIZE_BP'])
    #end = time.time()
    #print "Created LD server in {} seconds.".format("%0.4f" % (end - start))
    #start = time.time()
    for file in files:
        ldserver.set_file(str(file.path))
    #end = time.time()
    #print "Files initialized in {} seconds.".format("%0.4f" % (end - start))
    if 'last' in args:
        result = LDQueryResult(args['limit'], str(args['last']))
    else:
        result = LDQueryResult(args['limit'])
    if population_name != 'ALL':
        s = StringVec()
        s.extend(str(sample.sample) for sample in samples)
        ldserver.set_samples(str(population_name), s)
    if current_app.config['CACHE_ENABLED']:
        ldserver.enable_cache(file.reference_id, current_app.config['CACHE_REDIS_HOSTNAME'], current_app.config['CACHE_REDIS_PORT'])
    #start = time.time()
    ldserver.compute_region_ld(str(args['chrom']), args['start'], args['stop'], correlation_type, result, str(population_name))
    #print "Computed results in {} seconds.".format("%0.4f" % (time.time() - start))
    #start = time.time()
    if current_app.config['PROXY_PASS']:
        base_url = '/'.join(x.strip('/') for x in [current_app.config['PROXY_PASS'], request.path])
    else:
        base_url = request.base_url
    base_url += '?' + '&'.join(('{}={}'.format(arg, value) for arg, value in request.args.iteritems(True) if arg != 'last'))
    j = result.get_json(str(base_url))
    #print "Jsonified result in {} seconds.".format("%0.4f" % (time.time() - start))
    #start = time.time()
    r = make_response(j, 200)
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
        'correlation': fields.Str(required = True, validate = lambda x: x in ['r', 'rsquare'], error_messages = {'validator_failed': 'Value must be equal to \'r\' or \'rsquare\'.'}),
        'limit': fields.Int(required = False, validate = lambda x: x > 0, missing = current_app.config['API_MAX_PAGE_SIZE'], error_messages = {'validator_failed': 'Value must be greater than 0.'}),
        'last': fields.Str(required = False, validate = lambda x: len(x) > 0, error_messages = {'validator_failed': 'Value must be a non-empty string.'})
    }
    args = parser.parse(arguments, validate = validate_query)
    if args['limit'] > current_app.config['API_MAX_PAGE_SIZE']:
        args['limit'] = current_app.config['API_MAX_PAGE_SIZE']
    if args['correlation'] == 'r':
        correlation_type = correlation.ld_r
    elif args['correlation'] == 'rsquare':
        correlation_type = correlation.ld_rsquare
    else:
        abort(404)
    reference = Reference.query.filter_by(genome_build = genome_build, name = reference_name).first()
    if reference is None:
        abort(404)
    samples = Sample.query.with_parent(reference).filter_by(subset = population_name).all()
    if not samples:
        abort(404)
    files = File.query.with_parent(reference).all()
    if not files:
        abort(404)
    ldserver = LDServer(current_app.config['SEGMENT_SIZE_BP'])
    for file in files:
        ldserver.set_file(str(file.path))
    if 'last' in args:
        result = LDQueryResult(args['limit'], str(args['last']))
    else:
        result = LDQueryResult(args['limit'])
    if population_name != 'ALL':
        s = StringVec()
        s.extend(str(sample.sample) for sample in samples)
        ldserver.set_samples(str(population_name), s)
    if current_app.config['CACHE_ENABLED']:
        ldserver.enable_cache(file.reference_id, current_app.config['CACHE_REDIS_HOSTNAME'], current_app.config['CACHE_REDIS_PORT'])
    start = time.time()
    ldserver.compute_variant_ld(str(args['variant']), str(args['chrom']), args['start'], args['stop'], correlation_type, result, str(population_name))
    print "Computed results in {} seconds.".format("%0.4f" % (time.time() - start))
    start = time.time()
    if current_app.config['PROXY_PASS']:
        base_url = '/'.join(x.strip('/') for x in [current_app.config['PROXY_PASS'], request.path])
    else:
        base_url = request.base_url
    base_url += '?' + '&'.join(('{}={}'.format(arg, value) for arg, value in request.args.iteritems(True) if arg != 'last'))
    j = result.get_json(str(base_url))
    print "Jsonified result in {} seconds.".format("%0.4f" % (time.time() - start))
    start = time.time()
    r = make_response(j, 200)
    r.mimetype = 'application/json'
    print "Response created in {} seconds.".format("%0.4f" % (time.time() - start))
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
