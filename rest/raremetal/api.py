from flask import current_app, Blueprint, request, jsonify, make_response, abort
from flask_cors import CORS
from flask_compress import Compress
from webargs.flaskparser import parser
from webargs import fields, ValidationError
from functools import partial
import model
import os
from ld.pywrapper import StringVec, ColumnType, ColumnTypeMap, Mask, MaskVec, VariantGroupType, \
                         ScoreCovarianceRunner, ScoreCovarianceConfig, GroupIdentifierType

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


@bp.route("/aggregation/metadata", methods=["GET"])
def get_metadata():
  """
  Endpoint to describe all available datasets on which covariance may be computed.
  :return:
  """

  outer = {"data": model.get_full_genotype_datasets()}
  for gd in outer["data"]:
    gd["masks"] = model.get_masks_for_genotypes(gd["genotypeDataset"])
    gd["phenotypeDatasets"] = model.get_phenotypes_for_genotypes(gd["genotypeDataset"])

  resp = make_response(jsonify(outer), 200)
  resp.mimetype = 'application/json'

  return resp

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
    'genotypeDataset': fields.Int(required=True, validate=lambda x: x > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'phenotypeDataset': fields.Int(required=True, validate=lambda x: x > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'phenotype': fields.Str(required=True, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'samples': fields.Str(required=True, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'masks': fields.DelimitedList(fields.Str(), required=True, validate=lambda x: len(x) > 0, error_messages={'validator_failed': "Must provide at least 1 mask ID"}),
    'genomeBuild': fields.Str(required=True, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
  }

  args = parser.parse(
    args_defined,
    validate = partial(
      validate_query,
      all_fields = ['chrom', 'start', 'stop']
    )
  )

  config = ScoreCovarianceConfig()
  config.segment_size = current_app.config["SEGMENT_SIZE_BP"]

  chrom = str(args["chrom"])
  start = args["start"]
  stop = args["stop"]
  build = args["genomeBuild"]
  genotype_dataset_id = args["genotypeDataset"]
  phenotype_dataset_id = args["phenotypeDataset"]
  sample_subset = str(args["samples"])
  phenotype = str(args["phenotype"])
  masks = [str(x) for x in args["masks"]]

  config.chrom = chrom
  config.start = start
  config.stop = stop
  config.sample_subset = str(args["samples"])

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

  genotype_files = StringVec()
  genotype_files.extend(model.get_files(genotype_dataset_id))
  config.genotype_files = genotype_files
  config.genotype_dataset_id = genotype_dataset_id

  if sample_subset != 'ALL':
    s = StringVec()
    s.extend(model.get_samples(genotype_dataset_id, sample_subset))
    config.samples = s

  config.phenotype_file = model.get_phenotype_file(phenotype_dataset_id)
  config.column_types = model.get_column_types(phenotype_dataset_id)
  config.nrows = model.get_phenotype_nrows(phenotype_dataset_id)
  config.columns = model.get_phenotype_columns(phenotype_dataset_id)

  config.phenotype_dataset_id = phenotype_dataset_id
  config.phenotype = phenotype

  if current_app.config["CACHE_ENABLED"]:
    config.redis_hostname = current_app.config["CACHE_REDIS_HOSTNAME"]
    config.redis_port = current_app.config["CACHE_REDIS_PORT"]

  # Determine regions in which to compute LD/scores.
  # This is determined by the mask file, and the overall window requested.
  # TODO: check mask is valid for genome build
  # TODO: check mask is valid for genotype dataset ID
  mask_vec = MaskVec()
  for mask_name in masks:
    mask = model.get_mask_by_name(mask_name, genotype_dataset_id)

    if not os.path.isfile(mask["filepath"]):
      return_error(
        "Could not find mask file on server for mask {} with genotype dataset ID {}".format(mask_name, genotype_dataset_id),
        500
      )

    tb = Mask(str(mask["filepath"]), str(mask["name"]), mask["group_type"], mask["identifier_type"], chrom, start, stop)
    mask_vec.append(tb)

  config.masks = mask_vec

  runner = ScoreCovarianceRunner(config)
  runner.run()
  json = runner.getJSON()

  resp = make_response(json, 200)
  resp.mimetype = 'application/json'

  return resp