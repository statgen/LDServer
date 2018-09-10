import argparse
import requests
import time
import sys

argparser = argparse.ArgumentParser(description = '')
argparser.add_argument('-n', '--hostname', metavar = 'name', required = True, type = str, dest = 'hostname', help = 'LD Server hostname.')
argparser.add_argument('-p', '--port', metavar = 'number', required = False, type = int, dest = 'port', help = 'LD Server port number.')


if __name__ == '__main__':
    args = argparser.parse_args()

    query = '1000G_GRCh37/ALL/ld/region?chrom=22&start=51241101&stop=51241385'

    url = 'http://{}{}/{}'.format(args.hostname, ':' + str(args.port) if args.port else '', query)

    start = time.time()
    response = requests.get(url)
    print response
    end = time.time()

    if response.status_code != 200:
        sys.exit('Request failed with code {}'.format(response.status_code))

    data = response.json()['data']
    n_variants = len(set(data['variant1']).union(data['variant2']))
    n_results = len(data['rsquare'])

    print "Received {} LD values for {} variants in {} seconds.".format(n_results, n_variants, "%0.4f" % (end - start))