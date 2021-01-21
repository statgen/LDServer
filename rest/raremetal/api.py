from flask import current_app, Blueprint, request, jsonify, make_response, abort
from flask_cors import CORS
from flask_compress import Compress
from .errors import FlaskException
from webargs.flaskparser import parser
from marshmallow import Schema
from marshmallow.validate import OneOf
from webargs import fields, ValidationError
from functools import partial
from raven.versioning import fetch_git_sha
from . import model
import os
from core.pywrapper import (
  StringVec, ColumnType, ColumnTypeMap, Mask, MaskVec, VariantGroupType, LDServerGenericException,
  ScoreCovarianceRunner, ScoreCovarianceConfig, GroupIdentifierType, VariantGroup, VariantGroupVector, VariantFormat,
  VariantFilter
)

API_VERSION = "1.0"
MAX_UINT32 = 4294967295

bp = Blueprint('api', __name__)
CORS(bp)

compress = Compress()

def makeStringVec(arr=[]):
  s = StringVec()
  s.extend(arr)
  return s

@parser.error_handler
def handle_parsing_error(error, request, schema, error_status_code=None, error_headers=None):
  for location, field_dict in error.messages.items():
    for field, message in field_dict.items():
      message = 'Error while parsing \'{}\' query parameter: {}'.format(field, message[0])
      break

  raise FlaskException(message, 400)


def validate_query(parsed_fields, all_fields):
  for key, value in request.args.items():
    if key not in all_fields:
      raise ValidationError({key: ['Unknown parameter.']})

  if 'start' in parsed_fields and 'stop' in parsed_fields:
    if parsed_fields['stop'] <= parsed_fields['start']:
      raise ValidationError({'start': ['Start position must be greater than stop position.']})

  return True


@bp.route("/status", methods=["GET"])
def get_status():
  try:
    sha = fetch_git_sha(os.path.join(current_app.root_path, "../../"))
  except:
    sha = "no-git"

  json = {
    "data": {
      "sha": sha
    }
  }
  resp = make_response(jsonify(json), 200)
  resp.mimetype = "application/json"
  return resp


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

  sum_stat_datasets = model.get_full_summary_stat_datasets()
  for sd in sum_stat_datasets:
    sd["masks"] = model.get_masks_for_summary_stats(sd["summaryStatisticDataset"])

  outer["data"].extend(sum_stat_datasets)

  resp = make_response(jsonify(outer), 200)
  resp.mimetype = 'application/json'

  return resp

def return_error(error_message, http_code):
  response = {"data": None, "error": error_message}
  return make_response(jsonify(response), http_code)

ALLELES = ("A","C","G","T","N","U")

def validate_variant(value):
  try:
    if "_" in value:
      s1 = value.split(":")
      chrom = s1[0]
      pos, alleles = s1[1].split("_")
      pos = int(pos)
      ref, alt = alleles.split("/")
    else:
      chrom, pos, ref, alt = value.split(":")
      pos = int(pos)
  except:
    raise ValidationError("Invalid variant {}, should be of format: CHROM:POS_REF/ALT or CHROM:POS:REF:ALT".format(value))

  if chrom.startswith("chr"):
    raise ValidationError("Variant chromosome should not contain 'chr': " + value)

  if not all((x in ALLELES for x in ref)):
    raise ValidationError("Variant {} had invalid alleles".format(value))

  if not all((x in ALLELES for x in alt)):
    raise ValidationError("Variant {} had invalid alleles".format(value))

  return value

class GroupDefinitionField(fields.Field):
  keys = ("start", "stop", "filters")
  ops = ("gte", "lte", "eq")

  def _serialize(self, value, attr, obj, **kwargs):
    validated = str(self._validated(value)) if value is not None else None
    return super(fields.Field, self)._serialize(validated, attr, obj)

  def _deserialize(self, value, attr, data, **kwargs):
    return self._validated(value)

  def _validated(self, value):
    if not (isinstance(value, list) or isinstance(value, dict)):
      raise ValidationError("Each group must be either a list of variants, or dict (key/value) pairs specifying region parameters")

    if isinstance(value, dict):
      if "start" not in value:
        raise ValidationError("Must provide start position for region")

      if "stop" not in value:
        raise ValidationError("Must provide stop position for region")

      if "filter" in value:
        for f in value["filters"]:
          if f["op"] not in GroupDefinitionField.ops:
            raise ValidationError(f"Invalid op, must be one of: {', '.join(GroupDefinitionField.ops)}")

          if f["field"] == "maf":
            if f["value"] < 0 or f["value"] > 1:
              raise ValidationError("Filter for 'maf' must specify value between 0 and 1 inclusive")

    if isinstance(value, list):
      if len(value) == 0:
        raise ValidationError("List of variants in group definition must have at least 1 variant")

      for v in value:
        validate_variant(v)

    return value

class MaskSchema(Schema):
  id = fields.Int(required=True)
  name = fields.Str(required=True)
  description = fields.Str(required=True)
  genome_build = fields.Str(required=True)
  group_type = fields.Str(required=True)
  identifier_type = fields.Str(required=True)
  groups = fields.Dict(keys=fields.Str(), values=GroupDefinitionField)

@bp.route('/aggregation/covariance', methods = ['POST'])
def get_covariance():
  """
  This endpoint returns covariance and score statistics within a given region.
  """

  if request.content_type not in ("application/x-www-form-urlencoded", "application/json"):
    raise FlaskException("Content-Type must be application/json or application/x-www-form-urlencoded", 415)

  args_defined = {
    'chrom': fields.Str(required=True, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'start': fields.Int(required=True, validate=lambda x: x >= 0, error_messages={'validator_failed': 'Value must be greater than or equal to 0.'}),
    'stop': fields.Int(required=True, validate=lambda x: x > 0, error_messages={'validator_failed': 'Value must be greater than 0.'}),
    'genotypeDataset': fields.Int(required=False, validate=lambda x: x > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'phenotypeDataset': fields.Int(required=False, validate=lambda x: x > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'summaryStatisticDataset': fields.Int(required=False, validate=lambda x: x > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'phenotype': fields.Str(required=False, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'samples': fields.Str(required=False, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'masks': fields.List(fields.Int(), validate=lambda x: len(x) > 0, error_messages={'validator_failed': "Must provide at least 1 mask ID"}),
    'maskDefinitions': fields.Nested(MaskSchema, many=True),
    'genomeBuild': fields.Str(required=True, validate=lambda x: len(x) > 0, error_messages={'validator_failed': 'Value must be a non-empty string.'}),
    'variantFormat': fields.Str(required=False, default="EPACTS", validate=OneOf(["EPACTS", "COLONS"]))
  }

  args = parser.parse(
    args_defined,
    validate = partial(
      validate_query,
      all_fields = ['chrom', 'start', 'stop']
    ),
    location = "json_or_form"
  )

  if not (bool(args.get("masks")) ^ bool(args.get("maskDefinitions"))):
    raise FlaskException("Must provide either 'masks' or 'maskDefinitions' in request, and not both.", 400)

  calc_mode = (args.get("genotypeDataset") is not None) and (args.get("phenotypeDataset") is not None)
  precalc_mode = args.get("summaryStatisticDataset") is not None

  if calc_mode and precalc_mode:
    raise FlaskException("Must give either genotypeDataset and phenotypeDataset, or summaryStatisticDataset by itself", 400)
  if not calc_mode and not precalc_mode:
    raise FlaskException("No genotypeDataset, phenotypeDataset, or summaryStatisticDataset was provided in the request", 400)

  config = ScoreCovarianceConfig()
  config.segment_size = current_app.config["SEGMENT_SIZE_BP"]

  chrom = str(args["chrom"])
  start = args["start"]
  stop = args["stop"]
  build = args["genomeBuild"]
  genotype_dataset_id = args.get("genotypeDataset")
  phenotype_dataset_id = args.get("phenotypeDataset")
  summary_stat_dataset_id = args.get("summaryStatisticDataset")
  sample_subset = str(args.get("samples"))
  phenotype = str(args.get("phenotype"))
  masks = args.get("masks")
  mask_definitions = args.get("maskDefinitions")
  variant_format = args.get("variantFormat", "EPACTS")

  config.chrom = chrom
  config.start = start
  config.stop = stop
  config.sample_subset = str(args.get("samples"))

  if variant_format == "EPACTS":
    config.variant_format = VariantFormat.EPACTS
  elif variant_format == "COLONS":
    config.variant_format = VariantFormat.COLONS
  else:
    raise FlaskException("Invalid variant format given by variantFormat: {}".format(variant_format), 400)

  if (stop - start) > current_app.config["API_MAX_REGION_SIZE"]:
    raise FlaskException("Region requested for analysis exceeds maximum width of {}".format(current_app.config["API_MAX_REGION_SIZE"]), 400)

  if genotype_dataset_id or phenotype_dataset_id:
    if not (genotype_dataset_id and phenotype_dataset_id and phenotype):
      raise FlaskException("Must specify genotype dataset ID, phenotype dataset ID, and phenotype together")

  if genotype_dataset_id and phenotype_dataset_id:
    if not model.has_genome_build(build):
      raise FlaskException('Genome build \'{}\' was not found.'.format(build), 400)

    if not model.has_genotype_dataset(genotype_dataset_id):
      raise FlaskException('No genotype dataset \'{}\' available for genome build {}.'.format(genotype_dataset_id, build), 400)

    if not model.has_phenotype_dataset(phenotype_dataset_id):
      raise FlaskException('No phenotype dataset \'{}\' available for genome build {}.'.format(phenotype_dataset_id, build), 400)

    if not model.has_samples(genotype_dataset_id, sample_subset):
      raise FlaskException('Sample subset \'{}\' was not found in genotype dataset {}.'.format(sample_subset, genotype_dataset_id), 400)

    if not model.has_phenotype(phenotype_dataset_id, phenotype):
      raise FlaskException("Phenotype '{}' does not exist in phenotype dataset {}".format(phenotype, phenotype_dataset_id), 400)

    genotype_files = StringVec()
    genotype_files.extend([model.find_file(x) for x in model.get_genotype_files(genotype_dataset_id)])
    config.genotype_files = genotype_files
    config.genotype_dataset_id = genotype_dataset_id

    if sample_subset != 'ALL':
      s = StringVec()
      s.extend(model.get_samples(genotype_dataset_id, sample_subset))
      config.samples = s

    config.phenotype_file = model.find_file(model.get_phenotype_file(phenotype_dataset_id))
    config.phenotype_column_types = model.get_column_types(phenotype_dataset_id)
    config.phenotype_nrows = model.get_phenotype_nrows(phenotype_dataset_id)
    config.phenotype_sample_column = str(model.get_phenotype_sample_column(phenotype_dataset_id))
    config.phenotype_delim = str(model.get_phenotype_delim(phenotype_dataset_id))
    config.phenotype_dataset_id = phenotype_dataset_id
    config.phenotype = phenotype
    analysis_cols = model.get_analysis_columns(phenotype_dataset_id)
    config.phenotype_analysis_columns = makeStringVec(analysis_cols)

  elif summary_stat_dataset_id:
    config.summary_stat_dataset_id = summary_stat_dataset_id

    score_files = model.get_score_files(summary_stat_dataset_id)
    cov_files = model.get_cov_files(summary_stat_dataset_id)

    score_files = [model.find_file(f) for f in score_files]
    cov_files = [model.find_file(f) for f in cov_files]

    score_vec = StringVec()
    cov_vec = StringVec()
    score_vec.extend(score_files)
    cov_vec.extend(cov_files)

    config.summary_stat_score_files = score_vec
    config.summary_stat_cov_files = cov_vec

  if current_app.config["CACHE_ENABLED"]:
    config.redis_hostname = current_app.config["CACHE_REDIS_HOSTNAME"]
    config.redis_port = current_app.config["CACHE_REDIS_PORT"]

  # Determine regions in which to compute LD/scores.
  # This is determined by the mask file, and the overall window requested.
  mask_vec = MaskVec()
  if masks:
    for mask_id in masks:
      mask = model.get_mask_by_id(mask_id)
      mask_path = model.find_file(mask["filepath"])

      if not os.path.isfile(mask_path):
        raise FlaskException("Could not find mask file on server for mask ID {}".format(mask_id), 400)

      if mask["genome_build"] != build:
        raise FlaskException("Mask ID {} is invalid for genome build {}".format(mask_id, build), 400)

      if genotype_dataset_id and (genotype_dataset_id not in [g.id for g in mask["genotypes"]]):
        raise FlaskException("Mask ID {} is invalid for genotype dataset ID {}".format(mask_id, genotype_dataset_id), 400)

      if summary_stat_dataset_id and (summary_stat_dataset_id not in [s.id for s in mask["sumstats"]]):
        raise FlaskException("Mask ID {} is invalid for summary statistic dataset ID {}".format(mask_id, summary_stat_dataset_id), 400)

      tb = Mask(str(mask_path), mask_id, mask["group_type"], mask["identifier_type"], chrom, start, stop)
      mask_vec.append(tb)

  elif mask_definitions:
    for mask in mask_definitions:
      vg_vec = VariantGroupVector()
      for group_name, group_def in mask["groups"].items():
        vg = VariantGroup()
        vg.name = str(group_name)

        if isinstance(group_def, list):
          # This type of group definition is a list of variants.
          for v in group_def:
            if variant_format == "COLONS":
              # Internally we still use EPACTS format but translate on read/write. Later we will move to using a
              # native variant object type internally.
              v_chrom, pos, ref, alt = v.split(":")
              v = "{}:{}_{}/{}".format(v_chrom, pos, ref, alt)

            vg.add_variant(str(v))
        else:
          # This type of group definition is a dictionary, with keys for start/stop position, and optional filters.
          vg.chrom = chrom
          vg.start = group_def["start"]
          vg.stop = group_def["stop"]

          if "filters" in group_def:
            for f in group_def["filters"]:
              vf = VariantFilter()
              vf.op = f["op"]
              vf.field = f["field"]
              vf.set_value(f["value"])

              vg.filters.append(vf)

        # vg.stop and vg.start are either specified by the query (in the case of regions)
        # or computed by VariantGroup.add_variant (list of variants)
        if (vg.stop - vg.start) > current_app.config["API_MAX_COV_REGION_SIZE"]:
          raise FlaskException("Region requested for analysis exceeds maximum width of {}".format(current_app.config["API_MAX_COV_REGION_SIZE"]), 400)

        vg_vec.append(vg)

      tb = Mask(
        mask["id"],
        VariantGroupType.names.get(mask["group_type"]),
        GroupIdentifierType.names.get(mask["identifier_type"]),
        vg_vec
      )

      mask_vec.append(tb)

  config.masks = mask_vec

  runner = ScoreCovarianceRunner(config)
  runner.run()

  json = runner.getJSON()

  resp = make_response(json, 200)
  resp.mimetype = 'application/json'

  return resp
