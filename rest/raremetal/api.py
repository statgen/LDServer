from flask import current_app, Blueprint, request, jsonify, make_response, abort
from flask_cors import CORS
from flask_compress import Compress
from webargs.flaskparser import parser
from webargs import fields, ValidationError
from functools import partial
import model
import os
from ld.pywrapper import LDServer, LDQueryResult, StringVec, correlation, ScoreServer, ScoreStatQueryResult, \
                         ColumnType, ColumnTypeMap, make_shared_segment_vector
from pysam import TabixFile

API_VERSION = "1.0"
MAX_UINT32 = 4294967295

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

def return_error(error_message, http_code):
  response = {"data": None, "error": error_message}
  return make_response(jsonify(response), http_code)

@bp.route('/aggregation/covariance', methods = ['POST'])
def get_covariance():
  """
  This endpoint returns covariance and score statistics within a given region.
  """

  args_defined = {
    'chrom': fields.Str(required=True, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'start': fields.Int(required=True, validate=lambda x: x >= 0, error_messages={'validator_failed': 'Value must be greater than or equal to 0.'}),
    'stop': fields.Int(required=True, validate=lambda x: x > 0, error_messages={'validator_failed': 'Value must be greater than 0.'}),
    'genotype_dataset': fields.Int(required=True, validate=lambda x: x > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'phenotype_dataset': fields.Int(required=True, validate=lambda x: x > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'phenotype': fields.Str(required=True, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'samples': fields.Str(required=True, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'masks': fields.DelimitedList(fields.Str(), required=True, validate=lambda x: len(x) > 0, error_messages={'validator_failed': "Must provide at least 1 mask ID"}),
    'genome_build': fields.Str(required=True, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'limit': fields.Int(required=False, validate=lambda x: x > 0, missing=current_app.config['API_MAX_PAGE_SIZE'], error_messages={'validator_failed': 'Value must be greater than 0.'}),
    'last': fields.Str(required=False, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'})
  }

  args = parser.parse(
    args_defined,
    validate = partial(
      validate_query,
      all_fields = ['chrom', 'start', 'stop', 'limit', 'last']
    )
  )

  chrom = str(args["chrom"])
  start = args["start"]
  stop = args["stop"]
  build = args["genome_build"]
  genotype_dataset_id = args["genotype_dataset"]
  phenotype_dataset_id = args["phenotype_dataset"]
  sample_subset = str(args["samples"])
  phenotype = str(args["phenotype"])
  masks = [str(x) for x in args["masks"]]
  limit = args.get("limit")
  last = args.get("last")

  if limit > current_app.config['API_MAX_PAGE_SIZE']:
    limit = current_app.config['API_MAX_PAGE_SIZE']

  if not model.has_genome_build(build):
    response = { 'data': None, 'error': 'Genome build \'{}\' was not found.'.format(build) }
    return make_response(jsonify(response), 404)

  if not model.has_genotype_dataset(build, genotype_dataset_id):
    response = { 'data': None, 'error': 'No genotype dataset \'{}\' available for genome build {}.'.format(genotype_dataset_id, build) }
    return make_response(jsonify(response), 404)

  if not model.has_phenotype_dataset(phenotype_dataset_id):
    response = { 'data': None, 'error': 'No phenotype dataset \'{}\' available for genome build {}.'.format(phenotype_dataset_id, build) }
    return make_response(jsonify(response), 404)

  if not model.has_samples(genotype_dataset_id, sample_subset):
    response = { 'data': None, 'error': 'Sample subset \'{}\' was not found in genotype dataset {}.'.format(sample_subset, genotype_dataset_id) }
    return make_response(jsonify(response), 404)

  ldserver = LDServer(current_app.config['SEGMENT_SIZE_BP'])
  score_server = ScoreServer(current_app.config["SEGMENT_SIZE_BP"])
  for f in model.get_files(genotype_dataset_id):
    ldserver.set_file(f)
    score_server.set_genotypes_file(f, genotype_dataset_id)

  if 'last' in args:
    if "," not in last:
      raise Exception("This shouldn't happen")

    last_ld, last_score = last.split(",")
    result = LDQueryResult(limit, last_ld)
    score_result = ScoreStatQueryResult(limit, last_score)
  else:
    result = LDQueryResult(limit)
    score_result = ScoreStatQueryResult(limit)

  if sample_subset != 'ALL':
    s = StringVec()
    s.extend(model.get_samples(genotype_dataset_id, sample_subset))
    ldserver.set_samples(str(sample_subset), s)
    score_server.set_samples(str(sample_subset), s)

  # Load phenotypes.
  # This should be done after genotypes/samples have been loaded, because the phenotype object will need to match
  # those samples.
  phenotype_file = model.get_phenotype_file(phenotype_dataset_id)
  if phenotype_file.endswith(".ped"):
    # PED files specify everything we need in the DAT file, so we don't need column types / nrows like we do below
    # for tab files.
    score_server.load_phenotypes_ped(phenotype_file, phenotype_dataset_id)
  elif phenotype_file.endswith(".tab"):
    # Get the other information necessary for parsing the tab file.
    column_types = model.get_column_types(phenotype_dataset_id)
    nrows = model.get_phenotype_nrows(phenotype_dataset_id)
    score_server.load_phenotypes_tab(phenotype_file, column_types, nrows, phenotype_dataset_id)
  else:
    return_error("File format for phenotype file '{}' not supported".format(phenotype_file), 500)

  # Set the phenotype to calculate score stats / p-values for.
  score_server.set_phenotype(phenotype)

  # if current_app.config['CACHE_ENABLED']:
  #   ldserver.enable_cache(genotype_dataset_id, current_app.config['CACHE_REDIS_HOSTNAME'], current_app.config['CACHE_REDIS_PORT'])
  #   score_server.enable_cache(current_app.config['CACHE_REDIS_HOSTNAME'], current_app.config['CACHE_REDIS_PORT'])

  # Determine regions in which to compute LD/scores.
  # This is determined by the mask file, and the overall window requested.
  for mask_id in masks:
    # id = db.Column(db.Integer, primary_key = True)
    # name = db.Column(db.String, unique = True, nullable = False)
    # filepath = db.Column(db.String, unique = True, nullable = False)
    # description = db.Column(db.String, unique = False, nullable = False)
    # genome_build = db.Column(db.String, unique = False, nullable = False)
    # group_type = db.Column(db.String, unique = False, nullable = False)
    # identifier_type = db.Column(db.String, unique = False, nullable = False)
    # genotype_dataset_id = db.Column(db.Integer, db.ForeignKey('genotype_datasets.id'))

    mask = model.get_mask(mask_id)

    if not os.path.isfile(mask.filepath):
      return_error("Could not find mask file on server for ID {}".format(mask_id),500)

    tb = TabixFile(mask.filepath)
    window = "{}:{}-{}".format(chrom, start, stop)
    for row in tb.fetch(window):
      row = row.split("\t")
      row_group, row_chrom, row_start, row_end = row[0:4]
      row_variants = row[4:]

      # Tell the score server and LD server to compute results.
      ld_result = LDQueryResult(MAX_UINT32)
      score_result = ScoreStatQueryResult(MAX_UINT32)
      shared_segments = make_shared_segment_vector()
      ldserver.compute_region_ld(chrom, row_start, row_end, correlation.cov, ld_result, sample_subset, shared_segments)
      score_server.compute_scores(chrom, row_start, row_end, score_result, sample_subset, shared_segments)

      # We need to sort the results to make sure they come back in order, and matched up.
      ld_result.sort_by_variant()
      score_result.sort_by_variant()

  if current_app.config['PROXY_PASS']:
    base_url = '/'.join(x.strip('/') for x in [current_app.config['PROXY_PASS'], request.path])
  else:
    base_url = request.base_url

  base_url += '?' + '&'.join(('{}={}'.format(arg, value) for arg, value in request.args.iteritems(True) if arg != 'last'))
  r = make_response(result.get_json(str(base_url)), 200)
  r.mimetype = 'application/json'
  return r