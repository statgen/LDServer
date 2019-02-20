#include <boost/python.hpp>
//#include <boost/python/def.hpp>
//#include <boost/python/return_by_value/vector_indexing_suite.hpp>
//#include <boost/python/exception_translator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <python2.7/Python.h>
#include "LDServer.h"
#include "ScoreServer.h"
#include "Phenotypes.h"
#include "Mask.h"
#include "ScoreCovarianceRunner.h"

BOOST_PYTHON_MODULE(pywrapper) {

    boost::python::enum_<correlation>("correlation")
            .value("ld_r", LD_R)
            .value("ld_rsquare", LD_RSQUARE)
            .value("cov", COV)
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
            ;

    boost::python::class_<VariantsPair>("VariantsPair", boost::python::init<const char*, const char*, unsigned long int, const char*, const char*, unsigned long int , double>())
            .def_readonly("variant1", &VariantsPair::variant1)
            .def_readonly("chromosome1", &VariantsPair::chromosome1)
            .def_readonly("position1", &VariantsPair::position1)
            .def_readonly("variant2", &VariantsPair::variant2)
            .def_readonly("chromosome2", &VariantsPair::chromosome2)
            .def_readonly("position2", &VariantsPair::position2)
            .def_readonly("value", &VariantsPair::value)
            ;

    boost::python::class_<ScoreResult>("ScoreResult")
            .def_readonly("variant", &ScoreResult::variant)
            .def_readonly("score_stat", &ScoreResult::score_stat)
            .def_readonly("pvalue", &ScoreResult::pvalue)
            .def_readonly("alt_freq", &ScoreResult::alt_freq)
            ;

    boost::python::class_<std::vector<ScoreResult>>("ScoreResultVector")
            .def(boost::python::vector_indexing_suite<std::vector<ScoreResult>>())
            ;

    boost::python::class_<std::vector<VariantsPair> >("VariantsPairLDVec")
            .def(boost::python::vector_indexing_suite<std::vector<VariantsPair>>())
            ;

    boost::python::class_<std::vector<std::string> >("StringVec")
            .def(boost::python::vector_indexing_suite<std::vector<std::string>>())
            ;

    boost::python::class_<std::map<std::string, ColumnType>>("ColumnTypeMap")
            .def(boost::python::map_indexing_suite<std::map<std::string, ColumnType>>())
            ;

    boost::python::class_<LDQueryResult, boost::noncopyable>("LDQueryResult", boost::python::init<boost::uint32_t>())
            .def(boost::python::init<boost::uint32_t, const string&>())
            .def_readonly("limit", &LDQueryResult::limit)
            .def_readonly("data", &LDQueryResult::data)
            .def("has_next", &LDQueryResult::has_next)
            .def("is_last", &LDQueryResult::is_last)
            .def("get_last", &LDQueryResult::get_last)
            .def("get_json", &LDQueryResult::get_json)
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
            .def("compute_region_ld", &LDServer::compute_region_ld)
            .def("compute_variant_ld", &LDServer::compute_variant_ld)
            ;

    boost::python::class_<ScoreServer, boost::noncopyable>("ScoreServer", boost::python::init<boost::uint32_t>())
            .def("set_genotypes_file", &ScoreServer::set_genotypes_file)
            .def("load_phenotypes_tab", &ScoreServer::load_phenotypes_tab)
            .def("load_phenotypes_ped", &ScoreServer::load_phenotypes_ped)
            .def("set_phenotype", &ScoreServer::set_phenotype)
            .def("set_samples", &ScoreServer::set_samples)
            .def("enable_cache", &ScoreServer::enable_cache)
            .def("disable_cache", &ScoreServer::disable_cache)
            .def("get_chromosomes", &ScoreServer::get_chromosomes)
            .def("compute_scores", &ScoreServer::compute_scores)
            ;

    boost::python::class_<Mask, shared_ptr<Mask>, boost::noncopyable>("Mask", boost::python::init<const std::string&, const std::string&, VariantGroupType, GroupIdentifierType>())
            .def(boost::python::init<const std::string&, const std::string&, VariantGroupType, GroupIdentifierType, const std::string&, uint64_t, uint64_t>())
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
            .def_readwrite("column_types", &ScoreCovarianceConfig::column_types)
            .def_readwrite("nrows", &ScoreCovarianceConfig::nrows)
            .def_readwrite("masks", &ScoreCovarianceConfig::masks)
            .def_readwrite("sample_subset", &ScoreCovarianceConfig::sample_subset)
            .def_readwrite("samples", &ScoreCovarianceConfig::samples)
            .def_readwrite("segment_size", &ScoreCovarianceConfig::segment_size)
            .def_readwrite("redis_hostname", &ScoreCovarianceConfig::redis_hostname)
            .def_readwrite("redis_port", &ScoreCovarianceConfig::redis_port)
            ;

    boost::python::class_<ScoreCovarianceRunner, shared_ptr<ScoreCovarianceRunner>, boost::noncopyable>("ScoreCovarianceRunner", boost::python::init<std::shared_ptr<ScoreCovarianceConfig>>())
            .def("run", &ScoreCovarianceRunner::run)
            .def("getJSON", &ScoreCovarianceRunner::getJSON)
            .def("getPrettyJSON", &ScoreCovarianceRunner::getPrettyJSON)
            ;

    boost::python::def("make_score_covariance_config", &make_score_covariance_config);
}