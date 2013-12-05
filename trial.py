from mrjob.job import MRJob
from test_parallel import run_matrix_tests
from test_parallel import random_bernoulli  


class TrialJob(MRJob):
    
    def __init__(self, *args, **kwargs):
        super(TrialJob, self).__init__(*args, **kwargs)

    def mapper(self, key, line):
        bigN, n, k, reps, mat_type = line.split(' ')
        if mat_type == "b":
            result = run_matrix_tests(random_bernoulli(int(n), int(bigN)), int(k), int(reps))
        yield (bigN, n, k, reps, mat_type), result

    def reducer(self, key, values):
        yield key, list(values)

if __name__ == '__main__':
    TrialJob().run()
    
    
