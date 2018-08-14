#include <boost/python.hpp>
//#include <boost/python/def.hpp>
//#include <boost/python/return_by_value/vector_indexing_suite.hpp>
//#include <boost/python/exception_translator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <python2.7/Python.h>
#include "LDServer.h"

void (LDServer::*set_file)(const std::string& file) = &LDServer::set_file;
void (LDServer::*compute_region_ld)(const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result, const std::string& subset_name) = &LDServer::compute_region_ld;
void (LDServer::*compute_variant_ld)(const std::string& index_variant, const std::string& region_chromosome, std::uint64_t region_start_bp, std::uint64_t region_stop_bp, struct LDQueryResult& result, const std::string& subset_name) = &LDServer::compute_variant_ld;

BOOST_PYTHON_MODULE(PyLDServer) {

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

    boost::python::class_<std::vector<VariantsPairLD> >("VariantsPairLD")
            .def(boost::python::vector_indexing_suite<std::vector<VariantsPairLD>>())
            ;

    boost::python::class_<LDQueryResult, boost::noncopyable>("LDQueryResult", boost::python::init<boost::uint32_t>())
            .def(boost::python::init<const string&>())
            .def_readonly("limit", &LDQueryResult::limit)
            .def_readonly("data", &LDQueryResult::data)
            .def("has_next", &LDQueryResult::has_next)
            .def("get_last", &LDQueryResult::get_last)
            ;

    boost::python::class_<LDServer, boost::noncopyable>("LDServer")
            .def("set_file", set_file)
            .def("compute_region_ld", compute_region_ld)
            .def("compute_variant_ld", compute_variant_ld)
            ;
}