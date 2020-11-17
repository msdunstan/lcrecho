'''
AssemblyGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes
# pylint: disable=wrong-import-order
from collections import defaultdict
from operator import itemgetter
import os
import sys

from Bio.Seq import Seq
from synbiochem.utils import dna_utils, ice_utils, seq_utils

from assembly.app.lcr3 import overhang
import pandas as pd


_H_PRIMER = 'ATAGTTCCCTTCACGATAGCCG'
_L_PRIMER = 'TGCTGGATACGACGCTCTACTC'

_HBB_PRIMER_FORW = 'GGATCCAAACTCGAGTAAGG'
_HBB_PRIMER_REV = 'CTTCTTAAAAGATCTTTTGAATTC'
_LBB_PRIMER_FORW = 'GGATCCAAACTCGAGTAAGG'
_LBB_PRIMER_REV = 'CTTCTTAAAAGATCTTTTGAATTC'


class Lcr3Designer():
    '''Class to design LCR v3 assemblies.'''

    def __init__(self, filename, ice_params, primer_melt_temp=60.0,
                 lcr_melt_temp=70.0, restr_enz=None):
        self.__filename = filename
        self.__primer_melt_temp = primer_melt_temp
        self.__lcr_melt_temp = lcr_melt_temp

        self.__restr_enz = ['MlyI'] if restr_enz is None else restr_enz

        self.__seqs = {}

        self.__ice_client_fact = ice_utils.ICEClientFactory()
        self.__ice_client = \
            self.__ice_client_fact.get_ice_client(ice_params['url'],
                                                  ice_params['username'],
                                                  ice_params['password'])

        self.__primers = defaultdict(dict)
        self.__primers[True]['Hbb'] = _HBB_PRIMER_FORW
        self.__primers[True]['Lbb'] = _LBB_PRIMER_FORW
        self.__primers[False]['Hbb'] = _HBB_PRIMER_REV
        self.__primers[False]['Lbb'] = _LBB_PRIMER_REV

        self.__domino_parts = defaultdict(dict)

        for forw in [True, False]:
            self.__domino_parts[forw]['Hbb'] = \
                self.__get_subseq('SBC010499', lcr_melt_temp, forw, False)
            self.__domino_parts[forw]['Lbb'] = \
                self.__get_subseq('SBC010500', lcr_melt_temp, forw, False)

        self.__overhangs = overhang.get_seqs()
        self.__overhang_idx = 0

        self.__part_overhang = {}

        self.__design_parts = self.__get_design_parts()
        self.__part_primers = self.__get_part_primers()
        self.__pair_dominoes = self.__get_pair_dominoes()

    def close(self):
        '''Close.'''
        self.__ice_client_fact.close()

    def get_design_parts(self):
        '''Get design parts.'''
        return self.__design_parts

    def get_part_primers(self):
        '''Get part primers.'''
        return self.__part_primers

    def get_pair_dominoes(self):
        '''Get pair dominoes.'''
        return self.__pair_dominoes

    def to_csv(self, out_dir):
        '''Write data to csv.'''
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        design_parts_df = pd.DataFrame(self.get_design_parts().items(),
                                       columns=['design', 'parts'])

        design_parts_df.to_csv(os.path.join(out_dir, 'design_parts.csv'),
                               index=False)

        part_primers_df = \
            pd.DataFrame({_get_primer(part, idx, prmr)
                          for part, prmrs in self.get_part_primers().items()
                          for idx, prmr in enumerate(prmrs)},
                         columns=['Name', 'Primer'])

        part_primers_df.to_csv(os.path.join(out_dir, 'part_primers.csv'),
                               index=False)

        pair_dominoes_df = \
            pd.DataFrame(sorted(list(
                {(_get_domino_id(*pair), dmn)
                 for pair, dmn in self.get_pair_dominoes().items()})),
                columns=['id', 'domino'])

        pair_dominoes_df.to_csv(os.path.join(out_dir, 'pair_dominoes.csv'),
                                index=False)

    def __get_design_parts(self):
        '''Get design parts.'''
        design_parts = defaultdict(list)

        with open(self.__filename) as fle:
            designs = [tuple(line.strip().split(',')) for line in fle]

            for design in designs:
                design_parts[design].append(design[0] + 'bb')

                for idx, _id in enumerate(design):
                    if _id not in ['H', 'L', '']:
                        environment = design[idx - 1:idx + 2]

                        part = ('H' if idx > 1 and
                                environment[0] == 'H' else '',
                                environment[1],
                                'L' if len(environment) > 2 and
                                environment[2] == 'L' else '')

                        design_parts[design].append(part)

        return design_parts

    def __get_part_primers(self):
        '''Get part primers.'''
        parts = sorted(list({part for parts in self.__design_parts.values()
                             for part in parts}), key=itemgetter(1, 0, 2))

        return {part: self.__get_primers_for_part(part) for part in parts}

    def __get_pair_dominoes(self):
        '''Get dominoes.'''
        pair_dominoes = {}

        for design_part in self.__design_parts.values():
            for idx, part in enumerate(design_part[:-1]):
                pair = (part, design_part[idx + 1])

                if pair not in pair_dominoes:
                    domino = self.__get_domino(pair)
                    pair_dominoes[pair] = domino

        return pair_dominoes

    def __get_primers_for_part(self, part):
        '''Get primers.'''
        return (self.__get_primer(part, True),
                self.__get_primer(part, False))

    def __get_primer(self, part, forward):
        '''Get primer from ICE id.'''
        primer = None

        if part not in self.__primers[forward]:
            if forward and part[0] == 'H':
                return self.__get_part_overhang(part[:2]) + _H_PRIMER
            if not forward and part[2] == 'L':
                return self.__get_part_overhang(part[1:]) + _L_PRIMER

            # else:
            primer = \
                self.__get_subseq(part[1], self.__primer_melt_temp,
                                  forward, True, self.__restr_enz)

            self.__primers[forward][part] = primer

        return self.__primers[forward][part]

    def __get_subseq(self, part_id, mlt_temp, forward, primer, restr_enz=None):
        '''Get subsequence by melting temperature.'''
        seq = self.__get_seq(part_id)

        if restr_enz:
            dna = dna_utils.DNA(seq=seq, name=part_id, desc=part_id)
            restrict_dnas = dna_utils.apply_restricts(dna, restr_enz)

            # This is a bit fudgy...
            # Essentially, return the longest fragment remaining after
            # digestion.
            # Assumes prefix and suffix are short sequences that are cleaved
            # off.
            restrict_dnas.sort(key=lambda x: len(x['seq']), reverse=True)
            seq = restrict_dnas[0]['seq']

        reagent_concs = {seq_utils.MG: 0.0} if primer else None

        subseq = seq_utils.get_seq_by_melt_temp(seq, mlt_temp,
                                                forward=forward,
                                                reagent_concs=reagent_concs)[0]

        if primer and not forward:
            return str(Seq(subseq).reverse_complement())

        return subseq

    def __get_seq(self, ice_id):
        '''Get seq.'''
        if ice_id not in self.__seqs:
            self.__seqs[ice_id] = \
                self.__ice_client.get_ice_entry(ice_id).get_seq()

        return self.__seqs[ice_id]

    def __get_part_overhang(self, part):
        '''Get part overhang.'''
        if part not in self.__part_overhang:
            self.__part_overhang[part] = self.__get_next_overhang()

        return self.__part_overhang[part]

    def __get_next_overhang(self):
        '''Get next overhang.'''
        ovrhng = self.__overhangs[self.__overhang_idx]
        self.__overhang_idx += 1
        return ovrhng.lower()

    def __get_domino(self, pair):
        '''Get domino.'''
        part1 = self.__get_domino_part(pair[0], True)
        part2 = self.__get_domino_part(pair[1], False)

        return part1 + part2

    def __get_domino_part(self, part, left):
        '''Get domino part from ICE id.'''
        primer = None

        if part not in self.__domino_parts[left]:
            if not left and part[0] == 'H':
                return self.__get_part_overhang(part[:2])
            if left and part[2] == 'L':
                return self.__get_part_overhang(part[1:])

            # else:
            primer = \
                self.__get_subseq(part[1], self.__lcr_melt_temp,
                                  not left, False, self.__restr_enz)

            self.__domino_parts[left][part] = primer

        return self.__domino_parts[left][part]


def _get_primer(part, reverse, prmr):
    '''Get primer.'''
    return (_get_part_id(part, not reverse) + ('-R' if reverse else '-F'),
            '/5Phos/' + prmr)


def _get_domino_id(left_part, right_part):
    '''Get domino id.'''
    return '_'.join([_get_part_id(left_part, False),
                     _get_part_id(right_part, True)])


def _get_part_id(part, start):
    '''Return part id.'''
    if isinstance(part, str):
        return part

    if start:
        return ''.join(part[:2])

    return ''.join(part[1:])


def main(args):
    '''main method.'''
    designer = Lcr3Designer(args[0],
                            {'url': args[1],
                             'username': args[2],
                             'password': args[3]})

    designer.to_csv(args[4])

    design_parts = designer.get_design_parts()
    part_primers = designer.get_part_primers()
    pair_dominoes = designer.get_pair_dominoes()
    designer.close()

    with open(args[5], 'w') as fle:
        fle.write('Designs and parts\n')

        for design, prts in design_parts.items():
            fle.write('%s\t%s\n' % (design, prts))

        fle.write('\nParts and primers\n')

        for part, primers in part_primers.items():
            fle.write('%s\t%s\n' % (part, primers))

        fle.write('\nPart pairs and dominoes\n')

        for pair, domino in pair_dominoes.items():
            fle.write('%s\t%s\n' % (pair, domino))

        fle.write('\nDominoes\n')

        for domino in sorted(list(set(pair_dominoes.values()))):
            fle.write('%s\n' % domino)


if __name__ == '__main__':
    main(sys.argv[1:])
