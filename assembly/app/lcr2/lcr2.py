'''
AssemblyGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-locals
# pylint: disable=wrong-import-order
from assembly.graph_writer import GraphWriter


class DominoPoolWriter(GraphWriter):
    '''Class for generating LCR worklist graphs.'''

    def __init__(self, plasmid_parts, output_name='dom'):
        self.__plasmid_parts = plasmid_parts
        GraphWriter.__init__(self, output_name)

    def _initialise(self):
        water = self._add_vertex('water', {'is_reagent': True})
        domino_vol = 0.15

        for plasmid_id, parts_map in self.__plasmid_parts.items():
            dom_pool_water_vol = 20

            domino_pool = self._add_vertex(plasmid_id + '_dominoes',
                                           {'is_reagent': False})

            for ice_id, part_ice in parts_map.items():
                if part_ice.get_parameter('Type') == 'DOMINO':
                    domino = self._add_vertex(ice_id, {'is_reagent': False})

                    self._add_edge(domino, domino_pool, {'Volume': domino_vol})

                    dom_pool_water_vol -= domino_vol

            self._add_edge(water, domino_pool, {'Volume': dom_pool_water_vol})


class LcrWriter(GraphWriter):
    '''Class for generating LCR worklist graphs.'''

    def __init__(self, plasmid_parts, part_vol, output_name='lcr'):
        self.__plasmid_parts = plasmid_parts
        self.__part_vol = part_vol
        GraphWriter.__init__(self, output_name)

    def _initialise(self):
        water = self._add_vertex('water', {'is_reagent': True})
        ampligase = self._add_vertex('ampligase', {'is_reagent': True,
                                                   'well_fixed': 'G12'})
        mm = self._add_vertex('mm_lcr', {'is_reagent': True})

        mm_vol = 7.0
        ampligase_vol = 1.5
        domino_pool_vol = 1.0

        for plasmid_id, parts_map in self.__plasmid_parts.items():
            plasmid_water_vol = 25

            domino_pool = self._add_vertex(plasmid_id + '_dominoes',
                                           {'is_reagent': False})

            # Make lcr plate:
            plasmid = self._add_vertex(plasmid_id, {'is_reagent': False})

            self._add_edge(mm, plasmid, {'Volume': mm_vol})
            self._add_edge(ampligase, plasmid, {'Volume': ampligase_vol})
            self._add_edge(domino_pool, plasmid, {'Volume': domino_pool_vol})

            plasmid_water_vol -= (mm_vol + ampligase_vol + domino_pool_vol)

            for ice_id, part_ice in parts_map.items():
                if part_ice.get_parameter('Type') != 'DOMINO':
                    part = self._add_vertex(ice_id + '_dig',
                                            {'is_reagent': False})
                    self._add_edge(part, plasmid, {'Volume': self.__part_vol})
                    plasmid_water_vol -= self.__part_vol

            self._add_edge(water, plasmid, {'Volume': plasmid_water_vol})
