'''
AssemblyGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=too-few-public-methods
from assembly.graph_writer import GraphWriter


class PartDigestWriter(GraphWriter):
    '''Class for generating Part digest worklist graphs.'''

    def __init__(self, part_ids, pcr_numbers, output_name='part_dig'):
        self.__part_ids = part_ids
        self.__pcr_numbers = pcr_numbers
        GraphWriter.__init__(self, output_name)

    def _initialise(self):
        for part_id in self.__part_ids:
            for idx in range(self.__pcr_numbers[part_id]):
                part = self._add_vertex('%s_dig_%i' % (part_id, idx + 1),
                                        {'is_reagent': False})
                primer_mix = self._add_vertex('mm_dig',
                                              {'is_reagent': True})
                self._add_edge(primer_mix, part, {'Volume': 75.0})


class PartPoolWriter(GraphWriter):
    '''Class for generating pools of digested Parts worklist graphs.'''

    def __init__(self, pcr_numbers, output_name='pool_pcr'):
        self.__pcr_numbers = pcr_numbers
        GraphWriter.__init__(self, output_name)

    def _initialise(self):
        for part_id, count in self.__pcr_numbers.items():
            pool = self._add_vertex('%s_dig' % (part_id),
                                    {'is_reagent': True})

            for idx in range(count):
                part = self._add_vertex('%s_dig_%i' % (part_id, idx + 1),
                                        {'is_reagent': False})

                self._add_edge(part, pool, {'Volume': 1000.0})
