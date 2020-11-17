'''
AssemblyGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=too-few-public-methods
from assembly.graph_writer import GraphWriter


class PartQcWriter(GraphWriter):
    '''Class for generating Part qc worklist graphs.'''

    def __init__(self, part_ids, pcr_numbers, output_name='part_qc'):
        self.__part_ids = part_ids
        self.__pcr_numbers = pcr_numbers
        GraphWriter.__init__(self, output_name)

    def _initialise(self):
        ladder = self._add_vertex('ladder',
                                  {'is_reagent': True, 'well_fixed': 'H12'})

        bffer = self._add_vertex('buffer', {'is_reagent': True})

        ladder_product = self._add_vertex('ladder_product',
                                          {'is_reagent': False,
                                           'well_fixed': 'H12'})
        self._add_edge(ladder, ladder_product, {'Volume': 2.0})
        self._add_edge(bffer, ladder_product, {'Volume': 22.0})

        for part_id in self.__part_ids:
            for idx in range(self.__pcr_numbers[part_id]):
                part = self._add_vertex('%s_dig_%i' % (part_id, idx + 1),
                                        {'is_reagent': False})
                product = \
                    self._add_vertex('%s_product_%i' % (part_id, idx + 1),
                                     {'is_reagent': False})

                self._add_edge(part, product, {'Volume': 1.0})
                self._add_edge(bffer, product, {'Volume': 23.0})
