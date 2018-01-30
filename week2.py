import argparse
import random
import itertools
from collections import deque, Counter
import subprocess
from lib import choices
import math

RANDOM_ABOVE = 1000
REGION_SIZE = 100000

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return list(itertools.zip_longest(*args, fillvalue=fillvalue))

def upuntil(predicate, iterable):
    for x in iterable:
        if predicate(x):
            yield x
        else:
            yield x
            break

def parse_fasta(fasta_name):
    lines = iter(open(fasta_name).readlines())
    # Get header sans > symbol
    new_header = next(lines)[1:].strip()
    for line in lines:
        header = new_header
        # Take lines until next header
        seq = line.strip()
        items = deque((x.strip() for x in upuntil(lambda x: x[0] != ">", lines)))
        # Get header sans > symbol
        new_header = items.pop()[1:].strip()
        seq = seq.upper() + ''.join(items).upper()
        yield header, seq

def reverse_complement(string):
    return string[::-1]

# Also need to take reverse complement occasionally
def sample_fasta(fn, contig, region, outfile):
    # samtools faidx chr22_157501105.fasta "gi|157501105|gb|CM000483.1|":1000000-1001000
    subprocess.call('samtools faidx %s %s:%s-%s > %s' % (fn, contig, region[0], region[1], outfile), shell=True)

def get_empirical_freqs(seq):
    cnt = dict(Counter(seq))
    freqs = { k: (v/len(seq)) for k,v in cnt.items() }
    return freqs


def generate_new_samples(fns):
    # sample 1
    fn = fns[0][0]
    contig = '"gi|157501105|gb|CM000483.1|"'
    outfile = fns[0][1]


    # get count so that we can calculate max region start
    seqs = list(parse_fasta(fn))
    max_length = len(seqs[0][1]) - REGION_SIZE - RANDOM_ABOVE

    random_start = random.randint(RANDOM_ABOVE, max_length)
    region = (random_start, random_start + REGION_SIZE)

    # sample fasta
    sample_fasta(fn, contig, region, outfile)

    # sample 2
    fn = fns[1][0]
    contig = '"gi|206583718|gb|CM000512.1|"'
    outfile = fns[1][1]

    # sample fasta
    sample_fasta(fn, contig, region, outfile)

def qscore_to_prob(Q):
    return math.pow(10,(Q/(-10)))

def parse_profile(profile):
    lines = dict((k, list(g)) for k,g in itertools.groupby([line.strip().split('\t') for line in open(profile).readlines()], lambda x: x[0]))

    def remove_key(x) :
        removed_key = list(map(lambda z: [int(w) for w in z], [y[1:] for y in x]))
        # group by once again, this time returning a tuple
        return dict((k, tuple([y[1:] for y in g])) for k,g in itertools.groupby(removed_key, lambda x: x[0]))

    # map items in list
    maps = {k: remove_key(v) for k, v in lines.items()}
    return maps

def draw_char(char, index, probs, freqs):
    qscore = choices(probs[char][index][0], cum_weights=probs[char][index][1])[0]
    prob_of_error = qscore_to_prob(qscore)
    is_error = choices([True, False], cum_weights=(prob_of_error, 1.00))[0]
    # If there is an error, draw a different character from empirical frequencies
    if(is_error):
        # draw new character from empirical freqs
        #import pdb;pdb.set_trace()
        return choices(list(freqs.keys()), weights=freqs.values())[0]
    else:
        return char

def simulate_short_reads(fn, bp_length, coverage, probs):
    # store strings from supplied fasta files and their reverse complements
    # use profile provided in ART to simulate errors

    # get the sequence from the file
    seq = list(parse_fasta(fn))[0][1]
    freqs = get_empirical_freqs(seq)

    # get reverse complement
    reverse_seq = reverse_complement(seq)

    # uniformly sample bp_length substrings (len(seq)/bp_length) * coverage times
    start_points = len(seq) - bp_length
    num_to_sample = math.ceil(len(seq)/bp_length) * coverage
    indexes = random.sample(range(start_points), k=num_to_sample)

    # Pull evenly from both reverse and forward strands
    sampled_seqs = [seq[i:i+bp_length] for i in indexes[:math.floor(num_to_sample/2)]]
    sampled_seqs += [reverse_seq[i:i+bp_length] for i in indexes[math.floor(num_to_sample/2):]]

    sampled_seqs_with_err = [''.join([draw_char(x, i, probs, freqs) for i, x in enumerate(seq)]) for seq in sampled_seqs]
    return sampled_seqs_with_err

def write_to_fastq(fn, header, reads):
    f = open(fn, 'w')
    for i, read in enumerate(reads):
        f.write('@' + header + '-' + str(i) + '\n')
        f.write(read + '\n')
        f.write('+\n')
        f.write('I'*len(read) + '\n')

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="filename")
    parser.add_argument("-b", "--basepairs", type=int, help="base pair length")
    parser.add_argument("-c", "--coverage", type=int, help="coverage")
    parser.add_argument("-o", "--output", type=str, help="output")
    args = parser.parse_args()

    #FN = [("chr22_157501105.fasta", "chr22_157501105.sampled.fasta"), ("chr22_206583718.fasta", "chr22_206583718.sampled.fasta")]
    #coverage = 10
    #bp_length = 100

    fn = args.input
    bp_length = args.basepairs
    coverage = args.coverage
    output = args.output

    #generate_new_samples(FN)

    probs = parse_profile("HiSeq2kL100R1.txt")
    reads = simulate_short_reads(fn, bp_length, coverage, probs)

    #write reads to fastq
    header = list(parse_fasta(fn))[0][0]
    write_to_fastq(output, header, reads)


if __name__ == '__main__':
    main()

