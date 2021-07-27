from os import remove
from os.path import join
from pandas import read_csv

#-------------------------------------------------------------------------------
# INPUT

#output directory
dirout = join('..', 'out')

#target file
fntrials = join(dirout, 'trials.csv')

#-------------------------------------------------------------------------------
# MAIN

#read the table and get rid of the index column
df = read_csv(fntrials)
df.drop('trial', axis=1, inplace=True)

#write columns to binary files
for col, dtype in zip(df.columns, df.dtypes):
    if col == 'classification':
        df[col].values.astype('int8').tofile(join(dirout, col))
    else:
        df[col].values.astype('float32').tofile(join(dirout, col))

#delete the big csv
#remove(fntrials)
#print('removed file: %s' % fntrials)

print('finished')
