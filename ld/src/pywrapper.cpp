#include <boost/python.hpp>
//#include <boost/python/def.hpp>
//#include <boost/python/return_by_value/vector_indexing_suite.hpp>
//#include <boost/python/exception_translator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <python2.7/Python.h>
#include "LDServer.h"
#include "Phenotypes.h"

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

    boost::python::class_<VariantsPair>("VariantsPair", boost::python::init<const char*, const char*, unsigned long int, const char*, const char*, unsigned long int , double>())
            .def_readonly("variant1", &VariantsPair::variant1)
            .def_readonly("chromosome1", &VariantsPair::chromosome1)
            .def_readonly("position1", &VariantsPair::position1)
            .def_readonly("variant2", &VariantsPair::variant2)
            .def_readonly("chromosome2", &VariantsPair::chromosome2)
            .def_readonly("position2", &VariantsPair::position2)
            .def_readonly("value", &VariantsPair::value)
            ;

    boost::python::class_<std::vector<VariantsPair> >("VariantsPairLDVec")
            .def(boost::python::vector_indexing_suite<std::vector<VariantsPair>>())
            ;

    boost::python::class_<std::vector<std::string> >("StringVec")
            .def(boost::python::vector_indexing_suite<std::vector<std::string>>())
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

    boost::python::class_<LDServer, boost::noncopyable>("LDServer", boost::python::init<boost::uint32_t>())
            .def("set_file", &LDServer::set_file)
            .def("set_samples", &LDServer::set_samples)
            .def("enable_cache", &LDServer::enable_cache)
            .def("disable_cache", &LDServer::disable_cache)
            .def("get_chromosomes", &LDServer::get_chromosomes)
            .def("compute_region_ld", &LDServer::compute_region_ld)
            .def("compute_variant_ld", &LDServer::compute_variant_ld)
            ;
}