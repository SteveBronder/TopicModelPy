import numpy as np
import random
from mpmath import mp

DTM = np.matrix('2,0,4,1;1,1,2,3')
print DTM


def initialize(DTM, K):
    '''
    :param DTM: a document term matrix
    :param K: the number of topics
    :return: NMZ, NM, NZ, NZW, M, J, topic assignments
    '''
    topics = []
    M = DTM.shape[0]
    J = DTM.shape[1]
    NMZ = np.zeros(shape=(M, K))
    NM = np.zeros(shape=(M, 1))
    NZ = np.zeros(shape=(K, 1))
    NZW = np.zeros(shape=(K, J))

    for m in range(0, M):
        for n in range(0, J):
            if (DTM[m, n] != 0):
                for word in range(0, DTM[m, n]):
                    # print word

                    k = np.random.randint(K)
                    # print k
                    # print DTM[m,n]
                    # print "DTM",m,",",k,":",DTM[m,k]
                    NMZ[m][k] += 1
                    NM[m] += 1
                    NZ[k] += 1
                    NZW[k][n] += 1
                    topics.append(k)

    topics = np.vstack(topics)

    return NMZ, NM, NZ, NZW, M, J, topics


def fullconditional(NMZ, NM, NZ, NZW, m, n, beta, alpha):
    '''
    :param NMZ: number of times doc M and topic Z interact
    :param NM: marginalizing NMZ over topics
    :param NZ: marginalizing NZW over words
    :param NZW: number of times topic Z and word W interact
    :param M: number of docs
    :param K: number of topics
    :param J: unique number of words
    :param m: a document in M
    :param n: a word in Nm
    :param beta: priors
    :param alpha: priors
    :return: the conditional distribution for the gibbs sampler
    '''
    # this is defined by equation 79 of PETA
    left = (NZW[:, n] + beta) / (NZ + beta * J).transpose()  # NM[m] below
    right = (NMZ[m, :] + alpha) / (NM[m] + alpha * K)  # there is -1 in eq 79
    pz = abs(left * right)[0]
    pz /= np.sum(pz)
    return pz


def draw(p):
    '''
    sample from a multinomial
    '''
    r = random.random()
    for i in range(len(p)):
        r = r - p[i]
        if r < 0:
            return i
    return len(p) - 1


def loglik(NMZ, NM, NZ, NZW, m, n, beta, alpha):
    '''
    compute the log likelihood of the model
    '''
    lik = 0
    for k in range(0, K):
        lik += logdir(NZW[Z:, ] + beta)
        lik -= logdir(beta, J)
    for m in range(0, M):
        lik += logdir(NMZ[m:, ] + alpha)
        lik -= logdir(alpha, K)

    return lik


def logdir(alpha, K):
    return K * gammaln(alpha) - gammaln(K * alpha)


def gibbs(NMZ, NM, NZ, NZW, M, K, J, maxiter):
    '''
    :param NMZ: number of times doc M and topic Z interact
    :param NM: marginalizing NMZ over topics
    :param NZ: marginalizing NZW over words
    :param NZW: number of times topic Z and word W interact
    :param M: number of docs
    :param K: number of topics
    :param J: unique number of words
    :param maxiter: how many times to run gibbs sampler
    :return: phi parameters
    '''
    for iter in range(0, maxiter):
        for m in range(0, M):
            for n in range(0, J):
                if (DTM[m, n] != 0):
                    for word in range(0, DTM[m, n]):
                        idx = word * n
                        k = topics[idx][0]
                        NMZ[m][k] -= 1
                        NM[m] -= 1
                        NZ[k] -= 1
                        NZW[k][n] -= 1

                        p_z = fullconditional(NMZ, NM, NZ, NZW, m, n, 1, 1)
                        print "p_z:\n"
                        print p_z
                        k = draw(p_z)
                        NMZ[m][k] += 1
                        NM[m] += 1
                        NZ[k] += 1
                        NZW[k][n] += 1
                        topics[idx][0] = k
        num = abs(NZW + 1)
        num /= np.sum(num, axis=1)[:, np.newaxis]
        print "this is phi:\n"
        print num
    return num

# test code

K = 3

# test the Gibbs sampler
NMZ, NM, NZ, NZW, M, J, topics = initialize(DTM, K)
dryrun = gibbs(NMZ, NM, NZ, NZW, M, K, J, 100)


def entropy(NMZ, NM, NZ, NZW, M, K, J, alpha, phi):
    '''
    compute perplexity as a function of entropy of the model
    '''
    AK = K * alpha
    N = 0
    ent = 0
    for m, d in enumerate(DTM):
        #print "m:", m
        #print "d", d
        theta = NMZ[m, :] / (M + AK)
        #print theta
        ent -= mp.log(np.inner(dryrun[:,m],theta))
        #print "ent:", ent
        N += M
    return mp.exp(ent/N)

test = entropy(NMZ, NM, NZ, NZW, M, K, J, 1, dryrun)
