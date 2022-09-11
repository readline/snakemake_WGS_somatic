#!/usr/bin/env python
import sys
import pandas as pd

def somsamplesheet(sspath):
    df = pd.read_csv(sspath, sep='\t')
    pair = {}
    single = {}
    for i in df.index:
        sample = df.loc[i,'Sample']
        control = df.loc[i,'Control']
        somid = '%s.vs.%s'%(sample, control)
        if control == 'PoN':
            single[somid] = [sample,control]
        else:
            pair[somid] = [sample, control]
    return pair, single

if __name__ == '__main__':
    print(somsamplesheet(sys.argv[1]))
