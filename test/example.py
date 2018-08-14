import PyLDServer

if __name__ == '__main__':
    ld = PyLDServer.LDServer()
    ld.set_file('test/chr22.test.sav')
    result = PyLDServer.LDQueryResult(4)
    print 'NAME_1 NAME_2 R R^2'
    while True:
        ld.compute_region_ld('22', 51241101, 51241385, result, 'ALL')
        for pair in result.data:
            print pair.variant1, pair.variant2, pair.r, pair.rsquare
        if not result.has_next():
            break
