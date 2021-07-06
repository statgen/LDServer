#include <boost/python.hpp>
//#include <boost/python/def.hpp>
//#include <boost/python/return_by_value/vector_indexing_suite.hpp>
//#include <boost/python/exception_translator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <Python.h>
#include "LDServer.h"
#include "ScoreServer.h"
#include "Phenotypes.h"
#include "Mask.h"
#include "ScoreCovarianceRunner.h"
#include "Raw.h"
using namespace boost::python;

//template <class E, class... Policies, class... Args>
//class_<E, Policies...> exception_(Args&&... args) {
//  class_<E, Policies...> cls(std::forward<Args>(args)...);
//
//  string scope_name = extract<string>(scope().attr("__name__"));
//  string class_name = extract<string>(cls.attr("__name__"));
//  string exc_name = scope_name + "." + class_name;
//  PyObject* type = PyErr_NewException(exc_name.c_str(), PyExc_Exception, nullptr);
//  if (!type) throw_error_already_set();
//  scope().attr(exc_name.c_str()) = handle<>(borrowed(type));
//
//  register_exception_translator<E>([ptr=cls.ptr(), type](E const& e) {
//    PyErr_SetObject(type, object(e).ptr());
//  });
//
//  cls.def("__str__", &E::what);
//  return cls;
//}

// General purpose boost python translator for exceptions deriving from std::exception.
// Example:
//   class UserDefinedException : public std::exception { using std::exception::exception; };
//   exception_<UserDefinedException>("UserDefinedException", init<std::string>());
template <class E, class... Policies, class... Args>
class_<E, Policies...> exception_(Args&&... args) {
  class_<E, Policies...> cls(std::forward<Args>(args)...);

  register_exception_translator<E>([cls](E const& e) {
    PyErr_SetObject(PyExc_Exception, object(e).ptr());
  });

  cls.def("__str__", &E::what);
  return cls;
}

// This means generate a thin set of wrappers for versions of this function with or without default arguments
// 2nd to last argument is minimum number of arguments the function should accept
// last argument is the maximum number of arguments the function should accept
// compute_region_ld()'s last argument has a default argument of nullptr
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(compute_region_ld_overloads, compute_region_ld, 5, 8)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(phenotype_load_file_overloads, Phenotypes::load_file, 5, 6)

BOOST_PYTHON_MODULE(pywrapper) {

    exception_<LDServerGenericException>("LDServerGenericException", init<std::string>())
      .def("get_secret", &LDServerGenericException::get_secret);

    boost::python::enum_<correlation>("correlation")
            .value("ld_r", LD_R)
            .value("ld_rsquare", LD_RSQUARE)
            .value("cov", COV)
            .value("ld_rsquare_approx", LD_RSQUARE_APPROX)
            ;

    boost::python::enum_<ColumnType>("ColumnType")
            .value("TEXT", TEXT)
            .value("CATEGORICAL", CATEGORICAL)
            .value("INTEGER", INTEGER)
            .value("FLOAT", FLOAT)
            ;

    boost::python::enum_<VariantGroupType>("VariantGroupType")
            .value("GENE", GENE)
            .value("REGION", REGION)
            ;

    boost::python::enum_<GroupIdentifierType>("GroupIdentifierType")
            .value("ENSEMBL", ENSEMBL)
            .value("COORDINATES", COORDINATES)
            ;

    void (VariantFilter::*variant_filter_set_value_str)(const std::string&) = &VariantFilter::set_value;
    void (VariantFilter::*variant_filter_set_value_double)(const double&) = &VariantFilter::set_value;

    boost::python::class_<VariantFilter, shared_ptr<VariantFilter>, boost::noncopyable>("VariantFilter")
            .def_readwrite("op", &VariantFilter::op)
            .def_readwrite("field", &VariantFilter::field)
            .def_readonly("value_string", &VariantFilter::value_string)
            .def_readonly("value_double", &VariantFilter::value_double)
            .def("set_value", variant_filter_set_value_str)
            .def("set_value", variant_filter_set_value_double)
            ;

    boost::python::class_<VariantGroup, shared_ptr<VariantGroup>, boost::noncopyable>("VariantGroup")
            .def("get_variants", &VariantGroup::get_variants)
            .def("add_variant", &VariantGroup::add_variant)
            .def_readwrite("name", &VariantGroup::name)
            .def_readwrite("chrom", &VariantGroup::chrom)
            .def_readwrite("start", &VariantGroup::start)
            .def_readwrite("stop", &VariantGroup::stop)
            .def_readwrite("filters", &VariantGroup::filters);

    boost::python::class_<LDQueryResultMatrix, boost::noncopyable>("LDQueryResultMatrix")
        .def_readonly("variants", &LDQueryResultMatrix::variants)
        .def_readonly("chromosomes", &LDQueryResultMatrix::chromosomes)
        .def_readonly("positions", &LDQueryResultMatrix::positions)
        .def_readonly("offsets", &LDQueryResultMatrix::offsets)
        .def_readonly("correlations", &LDQueryResultMatrix::correlations)
        ;

    boost::python::class_<LDQueryResultVector, boost::noncopyable>("LDQueryResultVector")
        .def_readonly("index_variant", &LDQueryResultVector::index_variant)
        .def_readonly("index_chromosome", &LDQueryResultVector::index_chromosome)
        .def_readonly("index_position", &LDQueryResultVector::index_position)
        .def_readonly("variants", &LDQueryResultVector::variants)
        .def_readonly("chromosomes", &LDQueryResultVector::chromosomes)
        .def_readonly("positions", &LDQueryResultVector::positions)
        .def_readonly("correlations", &LDQueryResultVector::correlations)
        ;

    boost::python::class_<ScoreResult>("ScoreResult")
            .def_readonly("variant", &ScoreResult::variant)
            .def_readonly("score_stat", &ScoreResult::score_stat)
            .def_readonly("pvalue", &ScoreResult::pvalue)
            .def_readonly("alt_freq", &ScoreResult::alt_freq)
            ;

    boost::python::class_<std::vector<VariantGroup>>("VariantGroupVector")
            .def(boost::python::vector_indexing_suite<std::vector<VariantGroup>>());

    boost::python::class_<std::vector<VariantFilter>>("VariantFilterVector")
            .def(boost::python::vector_indexing_suite<std::vector<VariantFilter>>());

    boost::python::class_<std::vector<ScoreResult>>("ScoreResultVector")
            .def(boost::python::vector_indexing_suite<std::vector<ScoreResult>>())
            ;

    boost::python::class_<std::vector<std::string> >("StringVec")
            .def(boost::python::vector_indexing_suite<std::vector<std::string>>())
            ;

    boost::python::class_<ColumnTypeMap, std::shared_ptr<ColumnTypeMap>>("ColumnTypeMap")
            .def("add", &ColumnTypeMap::add)
            ;

    boost::python::class_<LDQueryResult, boost::noncopyable>("LDQueryResult", boost::python::init<boost::uint32_t>())
            .def(boost::python::init<boost::uint32_t, const string&>())
            .def_readonly("limit", &LDQueryResult::limit)
            .def_readonly("data", &LDQueryResult::data)
            .def_readonly("next", &LDQueryResult::next)
            .def_readonly("error", &LDQueryResult::error)
            .def("is_last", &LDQueryResult::is_last)
            .def("has_next", &LDQueryResult::has_next)
            .def("set_next", &LDQueryResult::set_next)
            .def("get_json", &LDQueryResult::get_json)
            .def("get_messagepack_py", &LDQueryResult::get_messagepack_py)
            ;

    boost::python::class_<SingleVariantLDQueryResult, boost::noncopyable>("SingleVariantLDQueryResult", boost::python::init<boost::uint32_t>())
        .def(boost::python::init<boost::uint32_t, const string&>())
        .def_readonly("limit", &SingleVariantLDQueryResult::limit)
        .def_readonly("data", &SingleVariantLDQueryResult::data)
        .def_readonly("next", &SingleVariantLDQueryResult::next)
        .def_readonly("error", &SingleVariantLDQueryResult::error)
        .def("is_last", &SingleVariantLDQueryResult::is_last)
        .def("has_next", &SingleVariantLDQueryResult::has_next)
        .def("set_next", &SingleVariantLDQueryResult::set_next)
        .def("get_json", &SingleVariantLDQueryResult::get_json)
        .def("get_messagepack_py", &SingleVariantLDQueryResult::get_messagepack_py)
        ;

    boost::python::class_<ScoreStatQueryResult, boost::noncopyable>("ScoreStatQueryResult", boost::python::init<boost::uint32_t>())
            .def(boost::python::init<boost::uint32_t, const string&>())
            .def_readonly("limit", &ScoreStatQueryResult::limit)
            .def_readonly("data", &ScoreStatQueryResult::data)
            .def("has_next", &ScoreStatQueryResult::has_next)
            .def("is_last", &ScoreStatQueryResult::is_last)
            .def("get_last", &ScoreStatQueryResult::get_last)
            .def("get_json", &ScoreStatQueryResult::get_json)
            ;

    boost::python::class_<LDServer, boost::noncopyable>("LDServer", boost::python::init<boost::uint32_t>())
            .def("set_file", &LDServer::set_file)
            .def("set_samples", &LDServer::set_samples)
            .def("enable_cache", &LDServer::enable_cache)
            .def("disable_cache", &LDServer::disable_cache)
            .def("get_chromosomes", &LDServer::get_chromosomes)
            .def("compute_region_ld", &LDServer::compute_region_ld, compute_region_ld_overloads())
            .def("compute_variant_ld", &LDServer::compute_variant_ld)
            ;

    boost::python::class_<ScoreServer, boost::noncopyable>("ScoreServer", boost::python::init<boost::uint32_t>())
            .def("set_genotypes_file", &ScoreServer::set_genotypes_file)
            .def("load_phenotypes_file", &ScoreServer::load_phenotypes_file)
            .def("set_phenotype", &ScoreServer::set_phenotype)
            .def("set_samples", &ScoreServer::set_samples)
            .def("enable_cache", &ScoreServer::enable_cache)
            .def("disable_cache", &ScoreServer::disable_cache)
            .def("get_chromosomes", &ScoreServer::get_chromosomes)
            .def("compute_scores", &ScoreServer::compute_scores)
            ;

    boost::python::class_<Mask, shared_ptr<Mask>, boost::noncopyable>("Mask", boost::python::init<const std::string&, const uint64_t, VariantGroupType, GroupIdentifierType>())
            .def(boost::python::init<const std::string&, const uint64_t, VariantGroupType, GroupIdentifierType, const std::string&, uint64_t, uint64_t>())
            .def(boost::python::init<const uint64_t, VariantGroupType, GroupIdentifierType, const std::vector<VariantGroup>&>())
            .def("get_variant_set", &Mask::get_variant_set)
            .def("get_group_names", &Mask::get_group_names)
            .def("get_group", &Mask::get_group)
            .def("get_id", &Mask::get_id)
            .def("get_group_type", &Mask::get_group_type)
            .def("set_group_type", &Mask::set_group_type)
            .def("get_identifier_type", &Mask::get_identifier_type)
            .def("set_identifier_type", &Mask::set_identifier_type)
            ;

    boost::python::class_<std::vector<Mask> >("MaskVec")
            .def(boost::python::vector_indexing_suite<std::vector<Mask>>())
            ;

    boost::python::class_<shared_ptr<vector<shared_ptr<Segment>>>>("SharedSegmentVector");
    boost::python::def("make_shared_segment_vector", &make_shared_segment_vector);
    boost::python::class_<std::set<std::string>, shared_ptr<std::set<std::string>>>("StringSet");

    boost::python::enum_<VariantFormat>("VariantFormat")
      .value("EPACTS", VariantFormat::EPACTS)
      .value("COLONS", VariantFormat::COLONS);

    boost::python::class_<ScoreCovarianceConfig, shared_ptr<ScoreCovarianceConfig>, boost::noncopyable>("ScoreCovarianceConfig")
            .def("pprint", &ScoreCovarianceConfig::pprint)
            .def_readwrite("chrom", &ScoreCovarianceConfig::chrom)
            .def_readwrite("start", &ScoreCovarianceConfig::start)
            .def_readwrite("stop", &ScoreCovarianceConfig::stop)
            .def_readwrite("genotype_files", &ScoreCovarianceConfig::genotype_files)
            .def_readwrite("genotype_dataset_id", &ScoreCovarianceConfig::genotype_dataset_id)
            .def_readwrite("phenotype_file", &ScoreCovarianceConfig::phenotype_file)
            .def_readwrite("phenotype_dataset_id", &ScoreCovarianceConfig::phenotype_dataset_id)
            .def_readwrite("phenotype", &ScoreCovarianceConfig::phenotype)
            .def_readwrite("phenotype_column_types", &ScoreCovarianceConfig::phenotype_column_types)
            .def_readwrite("phenotype_nrows", &ScoreCovarianceConfig::phenotype_nrows)
            .def_readwrite("phenotype_delim", &ScoreCovarianceConfig::phenotype_delim)
            .def_readwrite("phenotype_sample_column", &ScoreCovarianceConfig::phenotype_sample_column)
            .def_readwrite("phenotype_analysis_columns", &ScoreCovarianceConfig::phenotype_analysis_columns)
            .def_readwrite("masks", &ScoreCovarianceConfig::masks)
            .def_readwrite("sample_subset", &ScoreCovarianceConfig::sample_subset)
            .def_readwrite("samples", &ScoreCovarianceConfig::samples)
            .def_readwrite("summary_stat_dataset_id", &ScoreCovarianceConfig::summary_stat_dataset_id)
            .def_readwrite("summary_stat_score_files", &ScoreCovarianceConfig::summary_stat_score_files)
            .def_readwrite("summary_stat_cov_files", &ScoreCovarianceConfig::summary_stat_cov_files)
            .def_readwrite("summary_stat_format", &ScoreCovarianceConfig::summary_stat_format)
            .def_readwrite("segment_size", &ScoreCovarianceConfig::segment_size)
            .def_readwrite("redis_hostname", &ScoreCovarianceConfig::redis_hostname)
            .def_readwrite("redis_port", &ScoreCovarianceConfig::redis_port)
            .def_readwrite("variant_format", &ScoreCovarianceConfig::variant_format)
            ;

    boost::python::class_<ScoreCovarianceRunner, shared_ptr<ScoreCovarianceRunner>, boost::noncopyable>("ScoreCovarianceRunner", boost::python::init<std::shared_ptr<ScoreCovarianceConfig>>())
            .def("run", &ScoreCovarianceRunner::run)
            .def("getJSON", &ScoreCovarianceRunner::getJSON)
            .def("getPrettyJSON", &ScoreCovarianceRunner::getPrettyJSON)
            ;

    boost::python::def("make_score_covariance_config", &make_score_covariance_config);

    boost::python::def("extract_samples", &extract_samples);
}