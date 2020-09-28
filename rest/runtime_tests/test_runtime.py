import argparse
import requests
import time
import sys
import random

argparser = argparse.ArgumentParser(description = 'Tool for runtime experiments.')
argparser_subparsers = argparser.add_subparsers(help = '', dest = 'command')

argparser.add_argument('-n', '--hostname', metavar = 'name', required = True, type = str, dest = 'hostname', help = 'LD Server hostname.')
argparser.add_argument('-p', '--port', metavar = 'number', required = False, type = int, dest = 'port', help = 'LD Server port number.')
argparser.add_argument('-b', '--genome-build', metavar = 'name', required = True, type = str, dest ='genome_build', help = 'Genome build name.')
argparser.add_argument('-f', '--reference', metavar = 'name', required = True, type = str, dest = 'reference', help = 'LD reference panel name.')

argparser_region = argparser_subparsers.add_parser('region', help = 'Region queries.')
argparser_region.add_argument('-g', '--genes', metavar = 'file', required = True, type = str, dest = 'genes_file', help = 'File with gene coordinates. Must have two columns (without header): start position (bp), stop position (bp).')
argparser_region.add_argument('-s', '--page-size', metavar = 'number', required = True, type = int, dest = 'page_size', default = 10000, help = 'Maximal page size.')
argparser_region.add_argument('-l', '--region-length', metavar = 'number', required = True, type = int, dest = 'region_length', default = 10000, help = 'Region length in base-pairs.')
argparser_region.add_argument('-c', '--queries-count', metavar = 'number', required = True, type = int, dest = 'c', default = 10000, help = 'Number of queries to generate.')

argparser_variant = argparser_subparsers.add_parser('variant', help = 'Variant queries.')
argparser_variant.add_argument('-v', '--variants', metavar = 'file', required = True, type = str, dest = 'variants_file', help = 'File with variants. Must have four columns (without header): chromosome, position (bp), REF allele, ALT allele.')
argparser_variant.add_argument('-s', '--page-size', metavar = 'number', required = True, type = int, dest = 'page_size', default = 10000, help = 'Maximal page size.')
argparser_variant.add_argument('-l', '--region-length', metavar = 'number', required = True, type = int, dest = 'region_length', default = 10000, help = 'Region length in base-pairs.')
argparser_variant.add_argument('-c', '--queries-count', metavar = 'number', required = True, type = int, dest = 'c', default = 10000, help = 'Number of queries to generate.')

argparser_query = argparser_subparsers.add_parser('query', help = 'List of queries to re-run.')
argparser_query.add_argument('-q', '--queries', metavar = 'file', required = True, type = str, dest = 'queries_file', help = 'Output file form previous run.')

if __name__ == '__main__':
    args = argparser.parse_args()
    queries = []
    if args.command == 'region':
        window_bp = 100000
        genes = []
        with open(args.genes_file, 'r') as f:
            for line in f:
                if line:
                    chrom, start_bp, stop_bp = line.rstrip().split()
                    start_bp, stop_bp = list(map(int, (start_bp, stop_bp)))
                    genes.append((chrom, start_bp, stop_bp))
        for i in range(0, args.c):
            chrom, start_bp, stop_bp = random.choice(genes)
            start_bp = random.randrange( start_bp - window_bp if start_bp > window_bp else start_bp, stop_bp + window_bp, 1)
            stop_bp = start_bp + args.region_length
            query = 'genome_builds/{}/references/{}/populations/ALL/regions?correlation=rsquare&chrom={}&start={}&stop={}&limit={}'.format(args.genome_build, args.reference, chrom, start_bp, stop_bp, args.page_size)
            queries.append((query, args.region_length, args.page_size))
    elif args.command == 'variant':
        variants = []
        with open(args.variants_file, 'r') as f:
            for line in f:
                chrom, position, ref, alt = line.rstrip().split()
                variants.append((chrom, int(position), ref, alt))
        for i in range(0, args.c):
            chrom, position, ref, alt = random.choice(variants)
            start_bp = random.randrange(position - args.region_length if position > args.region_length else 0 , position, 1)
            stop_bp = start_bp + args.region_length
            query = 'genome_builds/{}/references/{}/populations/ALL/variants?correlation=rsquare&variant={}:{}_{}/{}&chrom={}&start={}&stop={}&limit={}'.format(args.genome_build, args.reference, chrom, position, ref, alt, chrom, start_bp, stop_bp, args.page_size)
            queries.append((query, args.region_length, args.page_size))
    elif args.command == 'query':
        with open(args.queries_file, 'r') as f:
            header = f.readline().rstrip().split()
            for line in f:
                fields = dict(list(zip(header, line.rstrip().split())))
                queries.append((fields['QUERY'], fields['REGION_LENGTH'], fields['PAGE_SIZE']))
    print('QUERY\tREGION_LENGTH\tPAGE_SIZE\tN_VARIANTS\tN_RESULTS\tN_PAGES\tTOTAL_SECONDS\tRESPONSE_SECONDS\tUNCOMPRESSED_MB\tCOMPRESSED_MB')
    for query, region_length, page_size in queries:
        url = 'http://{}{}/{}'.format(args.hostname, ':' + str(args.port) if args.port else '', query)
        total_time = 0
        total_results = 0
        total_pages = 0
        while url is not None:
            start = time.time()
            response = requests.get(url)
            end = time.time()
            if response.status_code != 200:
                sys.exit('Request failed with code {}.\nQuery: {}'.format(response.status_code, query))
            data = response.json()['data']
            total_time += (end - start)
            assert len(data['variants']) == len(data['chromosomes'])
            assert len(data['variants']) == len(data['positions'])
            total_results += sum([len(x) for _, x in data['correlation'].items()])
            total_pages += 1
            url = response.json()['next']
        print('{}\t{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}'.format(query, region_length, page_size, len(data['variants']), total_results, total_pages, total_time, response.elapsed.total_seconds(), len(response.content) / (1024.0 * 1024.0), int(response.headers['Content-Length']) / (1024.0 * 1024.0)))
