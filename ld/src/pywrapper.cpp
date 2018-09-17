#include <boost/python.hpp>
//#include <boost/python/def.hpp>
//#include <boost/python/return_by_value/vector_indexing_suite.hpp>
//#include <boost/python/exception_translator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <python2.7/Python.h>
#include "LDServer.h"

void (LDServer::*set_file)(const std::string& file) = &LDServer::set_file;
//void (LDServer::*set_samples)(const std::string& name, const std::vector<std::string>& samples) = &LDServer::set_samples;
bool (LDServer::*compute_region_ld)(const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result, const std::string& subset_name) const = &LDServer::compute_region_ld;
bool (LDServer::*compute_variant_ld)(const std::string& index_variant, const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result, const std::string& subset_name) const = &LDServer::compute_variant_ld;

BOOST_PYTHON_MODULE(pywrapper) {

    boost::python::class_<VariantsPairLD>("VariantsPairLD", boost::python::init<const char*, const char*, unsigned long int, const char*, const char*, unsigned long int , double, double>())
            .def_readonly("variant1", &VariantsPairLD::variant1)
            .def_readonly("chromosome1", &VariantsPairLD::chromosome1)
            .def_readonly("position1", &VariantsPairLD::position1)
            .def_readonly("variant2", &VariantsPairLD::variant2)
            .def_readonly("chromosome2", &VariantsPairLD::chromosome2)
            .def_readonly("position2", &VariantsPairLD::position2)
            .def_readonly("r", &VariantsPairLD::r)
            .def_readonly("rsquare", &VariantsPairLD::rsquare)
            ;

    boost::python::class_<std::vector<VariantsPairLD> >("VariantsPairLDVec")
            .def(boost::python::vector_indexing_suite<std::vector<VariantsPairLD>>())
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
            .def("set_file", set_file)
            .def("set_samples", &LDServer::set_samples)
            .def("enable_cache", &LDServer::enable_cache)
            .def("disable_cache", &LDServer::disable_cache)
            .def("compute_region_ld", compute_region_ld)
            .def("compute_variant_ld", compute_variant_ld)
            ;
}