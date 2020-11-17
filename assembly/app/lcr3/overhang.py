'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=superfluous-parens
import os.path
import sys

from synbiochem.utils import seq_utils
from assembly.app.lcr3 import vienna_utils


def get_seqs(num_overhangs=16, target_melt_temp=70.0,
             mfe_threshold=0.0, max_repeat_nuc=3, evalue=1e-3):
    '''Get overhangs.'''
    overhangs = {}

    directory = os.path.dirname(os.path.realpath(__file__))
    filename = os.path.join(directory,
                            '%i_%.3f_%.3f_%i_%.3e.txt' %
                            (num_overhangs,
                             target_melt_temp,
                             mfe_threshold,
                             max_repeat_nuc,
                             evalue))

    if os.path.exists(filename):
        with open(filename) as fle:
            overhangs = {idx: line.strip() for idx, line in enumerate(fle)}

    while len(overhangs) < num_overhangs:
        overhang = seq_utils.get_rand_seq_by_melt_temp(target_melt_temp,
                                                       max_repeat_nuc)[0]

        if overhangs and not _is_valid(overhangs, overhang, evalue,
                                       target_melt_temp, mfe_threshold):
            continue

        overhangs[str(len(overhangs))] = overhang

    overhangs = list(overhangs.values())

    with open(filename, 'w') as fle:
        for overhang in overhangs:
            fle.write(overhang + '\n')

    return overhangs[:num_overhangs]


def _get_seqs():
    '''Get seqs.'''


def _is_valid(overhangs, overhang, evalue, anneal_temp, mfe_threshold):
    '''Check validity of queries.'''
    return not _is_blast_similar(overhangs, overhang, evalue) \
        and not _has_secondary_structure(overhangs, overhang, anneal_temp,
                                         mfe_threshold)


def _is_blast_similar(overhangs, overhang, evalue):
    '''Check if overhang is BLAST similar to other overhangs.'''
    for result in seq_utils.do_blast(overhangs,
                                     {'query': overhang},
                                     evalue=evalue,
                                     word_size=4):
        if result.alignments:
            return True

    return False


def _has_secondary_structure(overhangs, overhang, anneal_temp,
                             mfe_threshold):
    '''Check if overhang has secondary structure wrt other overhangs.'''
    for curr_overhang in overhangs.values():
        seqs = [curr_overhang + overhang, overhang + curr_overhang]

        mfes = [vienna_utils.run('mfe', [seq], temp=anneal_temp)
                for seq in seqs]

        for mfe in mfes:
            if mfe[0][0] < mfe_threshold:
                return True

    return False


def main(args):
    '''main method.'''
    get_seqs(int(args[0]), float(args[1]), float(args[2]), int(args[3]),
             float(args[4]))


if __name__ == '__main__':
    main(sys.argv[1:])
