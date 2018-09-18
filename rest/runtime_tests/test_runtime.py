import argparse
import requests
import time
import sys
import random

argparser = argparse.ArgumentParser(description = 'Tool for runtime experiments.')
argparser_subparsers = argparser.add_subparsers(help = '', dest = 'command')

argparser.add_argument('-n', '--hostname', metavar = 'name', required = True, type = str, dest = 'hostname', help = 'LD Server hostname.')
argparser.add_argument('-p', '--port', metavar = 'number', required = False, type = int, dest = 'port', help = 'LD Server port number.')

argparser_region = argparser_subparsers.add_parser('region', help = 'Region queries.')
argparser_region.add_argument('-r', '--repeat', metavar = 'file', required = False, type = str, dest = 'repeat_file', help = 'Output file from previous run. If specified, will rerun all queries.')
argparser_region.add_argument('-g', '--genes', metavar = 'file', required = False, type = str, dest = 'genes_file', help = 'File with gene coordinates. Must have two columns (without header): start position (bp), stop position (bp).')
argparser_region.add_argument('-s', '--page-size', metavar = 'number', required = True, type = int, dest = 'page_size', default = 10000, help = 'Maximal page size.')
argparser_region.add_argument('-l', '--region-length', metavar = 'number', required = True, type = int, dest = 'region_length', default = 10000, help = 'Region length in base-pairs.')
argparser_region.add_argument('-c', '--queries-count', metavar = 'number', required = True, type = int, dest = 'c', default = 10000, help = 'Number of queries to generate.')

argparser_variant = argparser_subparsers.add_parser('variant', help = 'Variant queries.')
argparser_variant.add_argument('-r', '--repeat', metavar = 'file', required = False, type = str, dest = 'repeat_file', help = 'Output file from previous run. If specified, will rerun all queries.')
argparser_variant.add_argument('-v', '--variants', metavar = 'file', required = False, type = str, dest = 'variants_file', help = 'File with variants. Must have four columns (without header): chromosome, position (bp), REF allele, ALT allele.')
argparser_variant.add_argument('-s', '--page-size', metavar = 'number', required = True, type = int, dest = 'page_size', default = 10000, help = 'Maximal page size.')
argparser_variant.add_argument('-l', '--region-length', metavar = 'number', required = True, type = int, dest = 'region_length', default = 10000, help = 'Region length in base-pairs.')
argparser_variant.add_argument('-c', '--queries-count', metavar = 'number', required = True, type = int, dest = 'c', default = 10000, help = 'Number of queries to generate.')

if __name__ == '__main__':
    args = argparser.parse_args()

    if args.command == 'region':
        window_bp = 100000
        if args.repeat_file and args.genes_file:
            sys.exit('Options -r/--repeat and -g/-genes can\'t be used together')
        if args.repeat_file:
            queries = []
            with open(args.repeat_file, 'r') as f:
                header = f.readline().rstrip().split()
                for line in f:
                    fields = dict(zip(header, line.rstrip().split()))
                    queries.append((fields['QUERY'], fields['REGION_LENGTH'], fields['PAGE_SIZE']))
            args.c = len(queries)
        elif args.genes_file:
            genes = []
            with open(args.genes_file, 'r') as f:
                for line in f:
                    if line:
                        start_bp, stop_bp = line.rstrip().split()
                        genes.append((int(start_bp), int(stop_bp)))
        print 'QUERY\tREGION_LENGTH\tPAGE_SIZE\tN_VARIANTS\tN_RESULTS\tN_PAGES\tSECONDS'
        for i in range(0, args.c):
            if args.repeat_file:
                query, args.region_length, args.page_size = queries[i]
            elif args.genes_file:
                start_bp, stop_bp = random.choice(genes)
                start_bp = random.randrange( start_bp - window_bp if start_bp > window_bp else start_bp, stop_bp + window_bp, 1)
                stop_bp = start_bp + args.region_length
                query = '1000G_GRCh37/ALL/ld/region?chrom=22&start={}&stop={}&limit={}'.format(start_bp, stop_bp, args.page_size)
            url = 'http://{}{}/{}'.format(args.hostname, ':' + str(args.port) if args.port else '', query)
            total_time = 0
            all_variants = set()
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
                all_variants.update(data['variant1'])
                all_variants.update(data['variant2'])
                total_results += len(data['rsquare'])
                total_pages += 1
                url = response.json()['next']
            print '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(query, args.region_length, args.page_size, len(all_variants), total_results, total_pages, '%0.4f' % total_time)
    elif args.command == 'variant':
        if args.repeat_file and args.variants_file:
            sys.exit('Options -r/--repeat and -v/--variants can\'t be used together')
        if args.repeat_file:
            queries = []
            with open(args.repeat_file, 'r') as f:
                header = f.readline().rstrip().split()
                for line in f:
                    fields = dict(zip(header, line.rstrip().split()))
                    queries.append((fields['QUERY'], fields['REGION_LENGTH'], fields['PAGE_SIZE']))
            args.c = len(queries)
        elif args.variants_file:
            variants = []
            with open(args.variants_file, 'r') as f:
                for line in f:
                    chrom, position, ref, alt = line.rstrip().split()
                    variants.append((chrom, int(position), ref, alt))
        print 'QUERY\tREGION_LENGTH\tPAGE_SIZE\tN_VARIANTS\tN_RESULTS\tN_PAGES\tSECONDS'
        for i in range(0, args.c):
            if args.repeat_file:
                query, args.region_length, args.page_size = queries[i]
            elif args.variants_file:
                chrom, position, ref, alt = random.choice(variants)
                start_bp = random.randrange(position - args.region_length if position > args.region_length else 0 , position, 1)
                stop_bp = start_bp + args.region_length
                query = '1000G_GRCh37/ALL/ld/variant?variant={}:{}_{}/{}&chrom={}&start={}&stop={}&limit={}'.format(chrom, position, ref, alt, chrom, start_bp, stop_bp, args.page_size)
            url = 'http://{}{}/{}'.format(args.hostname, ':' + str(args.port) if args.port else '', query)
            total_time = 0
            all_variants = set()
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
                all_variants.update(data['variant1'])
                all_variants.update(data['variant2'])
                total_results += len(data['rsquare'])
                total_pages += 1
                url = response.json()['next']
            print '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(query, args.region_length, args.page_size, len(all_variants), total_results, total_pages, '%0.4f' % total_time)



