import numpy as np
import pandas as pd
import os
import glob
import random
#
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.feature_extraction.text import TfidfVectorizer
from scipy.special import gammaln
#
DIR = r'data_folder/wordcounts'
allfiles = glob.glob(os.path.join(DIR,"*.CSV"))
p=.5
# sample files for train
gen_sample = np.array(sorted(random.sample(xrange(len(allfiles)), int(p * len(allfiles)))))
rand_sample = [ allfiles[i] for i in gen_sample ]
#
# take rest for test
rand_sample2 = []
for i in xrange(len(allfiles)):
    if i not in gen_sample:
        rand_sample2.append(allfiles[i])
#
# train data

np_array_list = []
for file_ in rand_sample:
    df = pd.read_csv(file_,index_col=None, header=0)
    df['source'] = file_
    np_array_list.append(df.as_matrix())
#
# test data
np_array_list_test = []
for file_ in rand_sample2:
    df = pd.read_csv(file_, index_col = None, header = 0)
    df['source'] = file_
    np_array_list_test.append(df.as_matrix())
    
#
# train data frame
comb_np_array = np.vstack(np_array_list)
train_frame = pd.DataFrame(comb_np_array)
train_frame.columns = ['words','count','source']
subless = (train_frame['words'].str.len() > 2)
submore = (train_frame['words'].str.len() < 20)
train_frame = train_frame.loc[subless]
train_frame = train_frame.loc[submore]
train_frame = train_frame.fillna(value = 0)
train_frame = train_frame.pivot(index = 'source',columns = 'words', values = 'count')
train_frame = train_frame.fillna(value = 0)
train_frame = train_frame.loc[:, (train_frame.sum(axis = 0) > 5)]
#

# test data frame
comb_np_array_test = np.vstack(np_array_list_test)
test_frame = pd.DataFrame(comb_np_array_test)
test_frame.columns = ['words','count','source']
test_frame = test_frame.fillna(value=0)
test_frame = test_frame.pivot(index = 'source', columns = 'words', values = 'count')
test_frame = test_frame.fillna(value = 0)

