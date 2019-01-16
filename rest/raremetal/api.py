from flask import current_app, Blueprint, request, jsonify, make_response, abort
from flask_cors import CORS
from flask_compress import Compress
from webargs.flaskparser import parser
from webargs import fields, ValidationError
from functools import partial
from .model import GenotypeDataset, File, Sample
import model
from ld.pywrapper import LDServer, LDQueryResult, StringVec, correlation
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


@bp.route('/constants/correlations', methods = ['GET'])
def get_correlations():
  response = { 'data': model.get_correlations(), 'error': None }
  return make_response(jsonify(response), 200)


def correlation_type(correlation_name):
  return { 'r': correlation.ld_r, 'rsquare': correlation.ld_rsquare, 'cov': correlation.cov }[correlation_name]


@bp.route("/aggregation/datasets", methods=["GET"])
def get_metadata():
  """
  Endpoint to describe all available datasets on which covariance may be computed.
  :return:
  """

  pass


@bp.route('/aggregation/covariance', methods = ['GET'])
def get_covariance(genome_build, genotype_dataset_name, sample_subset_name, phenotype):
  """
  This endpoint returns covariance and score statistics within a given region.
  :param genome_build:
  :param genotype_dataset_name:
  :param sample_subset_name:
  :param phenotype:
  :return:
  """

  args_defined = {
    'chrom': fields.Str(required = True, validate = lambda x: len(x) > 0, error_messages = {'validator_failed': 'Value must be a non-empty string.'}),
    'start': fields.Int(required = True, validate = lambda x: x >= 0, error_messages = {'validator_failed': 'Value must be greater than or equal to 0.'}),
    'stop': fields.Int(required = True, validate = lambda x: x > 0, error_messages = {'validator_failed': 'Value must be greater than 0.'}),
    'correlation': fields.Str(
      required = True,
      validate = lambda x: x in current_app.config['CORRELATIONS'],
      error_messages = {'validator_failed': 'Value must be one of the following: {}.'.format(', '.join(current_app.config['CORRELATIONS']))}
    ),
    'limit': fields.Int(required = False, validate = lambda x: x > 0, missing = current_app.config['API_MAX_PAGE_SIZE'], error_messages = {'validator_failed': 'Value must be greater than 0.'}),
    'last': fields.Str(required = False, validate = lambda x: len(x) > 0, error_messages = {'validator_failed': 'Value must be a non-empty string.'})
  }
  args = parser.parse(args_defined, validate = partial(validate_query, all_fields = ['chrom', 'start', 'stop', 'correlation', 'limit', 'last']))

  if args['limit'] > current_app.config['API_MAX_PAGE_SIZE']:
    args['limit'] = current_app.config['API_MAX_PAGE_SIZE']

  if not model.has_genome_build(genome_build):
    response = { 'data': None, 'error': 'Genome build \'{}\' was not found.'.format(genome_build) }
    return make_response(jsonify(response), 404)

  genotype_dataset_id = model.get_genotype_dataset_id(genome_build, genotype_dataset_name)
  if not genotype_dataset_id:
    response = { 'data': None, 'error': 'No genotype dataset \'{}\' available for genome build {}.'.format(genotype_dataset_name, genome_build) }
    return make_response(jsonify(response), 404)
  if not model.has_samples(genotype_dataset_id, sample_subset_name):
    response = { 'data': None, 'error': 'Population \'{}\' was not found in {} genotype_dataset panel.'.format(sample_subset_name, genotype_dataset_name) }
    return make_response(jsonify(response), 404)

  ldserver = LDServer(current_app.config['SEGMENT_SIZE_BP'])
  for f in model.get_files(genotype_dataset_id):
    ldserver.set_file(f)

  if 'last' in args:
    result = LDQueryResult(args['limit'], str(args['last']))
  else:
    result = LDQueryResult(args['limit'])

  if sample_subset_name != 'ALL':
    s = StringVec()
    s.extend(model.get_samples(genotype_dataset_id, sample_subset_name))
    ldserver.set_samples(str(sample_subset_name), s)

  if current_app.config['CACHE_ENABLED']:
    ldserver.enable_cache(genotype_dataset_id, current_app.config['CACHE_REDIS_HOSTNAME'], current_app.config['CACHE_REDIS_PORT'])

  #start = time.time()
  ldserver.compute_region_ld(str(args['chrom']), args['start'], args['stop'], correlation_type(args['correlation']), result, str(sample_subset_name))
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