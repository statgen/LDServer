import argparse
import requests
import time
import sys
import random

argparser = argparse.ArgumentParser(description = '')
argparser.add_argument('-n', '--hostname', metavar = 'name', required = True, type = str, dest = 'hostname', help = 'LD Server hostname.')
argparser.add_argument('-p', '--port', metavar = 'number', required = False, type = int, dest = 'port', help = 'LD Server port number.')
argparser.add_argument('-g', '--genes', metavar = 'file', required = True, type = str, dest = 'genes_file', help = 'File with gene coordinates. Must have two columns (without header): start position (bp), stop position (bp).')
argparser.add_argument('-s', '--page-size', metavar = 'number', required = True, type = int, dest = 'page_size', default = 10000, help = 'Maximal page size.')
argparser.add_argument('-r', '--region-length', metavar = 'number', required = True, type = int, dest = 'region_bp', default = 10000, help = 'Region length in base-pairs.')
argparser.add_argument('-k', '--k-queries', metavar = 'number', required = True, type = int, dest = 'k', default = 10000, help = 'Number of queries to generate.')

if __name__ == '__main__':
    args = argparser.parse_args()

    window_bp = 100000

    genes = []
    with open(args.genes_file, 'r') as f:
        for line in f:
            if line:
                start_bp, stop_bp = line.rstrip().split()
                genes.append((int(start_bp), int(stop_bp), int(stop_bp) - int(start_bp) + 1))

    print 'REGION_BP\tPAGE_SIZE\tN_VARIANTS\tN_RESULTS\tN_PAGES\tSECONDS'

    for i in range(0, args.k):
        start_bp, stop_bp, length_bp = random.choice(genes)
        start_bp = random.randrange( start_bp - window_bp if start_bp > window_bp else start_bp, stop_bp + window_bp, 1)
        stop_bp = start_bp + args.region_bp
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
                sys.exit('Request failed with code {}'.format(response.status_code))
            data = response.json()['data']
            total_time += (end - start)
            all_variants.update(data['variant1'])
            all_variants.update(data['variant2'])
            total_results += len(data['rsquare'])
            total_pages += 1
            url = response.json()['next']

        print '{}\t{}\t{}\t{}\t{}\t{}'.format(stop_bp - start_bp, args.page_size, len(all_variants), total_results, total_pages, '%0.4f' % total_time)

