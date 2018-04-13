import sys

def infname_to_genome_pair_to_jaccard(infname):
	with open(infname) as infile:
		lines = infile.readlines()
		genome_pair_to_jaccard = {}
		for line in lines:
			toks = line.strip().split()
			genome1 = toks[0]
			genome2 = toks[1]
			jaccard = toks[4]
			a, b = tuple(toks[4].split('/'))
			jaccard = float(a) / float(b)
			if (genome1 > genome2): genome1, genome2 = (genome2, genome1)
			genome_pair_to_jaccard[(genome1, genome2)] = jaccard
	return genome_pair_to_jaccard

expected = infname_to_genome_pair_to_jaccard(sys.argv[1]) # expected
observed = infname_to_genome_pair_to_jaccard(sys.argv[2]) # observed
exp_diff_list = []
for genome_pair in expected:
	jexp = expected[genome_pair]
	jobs = observed[genome_pair]
	diff = jobs - jexp
	exp_diff_list.append((jexp, diff))
exp_diff_list = sorted(exp_diff_list)
for exp, diff in exp_diff_list:
	print('{}\t{}'.format(exp, diff))
	
