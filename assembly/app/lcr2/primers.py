'''
AssemblyGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=wrong-import-order
from collections import defaultdict
import sys

from Bio.Seq import Seq
from synbiochem.utils import dna_utils, seq_utils

from assembly import plate
from assembly.app.lcr2 import utils
import pandas as pd


class PrimerDesigner():
    '''Class to design primers.'''

    def __init__(self, ice_details):
        self.__ice_helper = utils.ICEHelper(ice_details['url'],
                                            ice_details['username'],
                                            ice_details['password'])

    def close(self):
        '''Close.'''
        self.__ice_helper.close()

    def get_primers_from_plasmids(self, pathway_plasmid_ids, plates,
                                  restr_enz, tm, mg_conc=0.0):
        '''Get primers for Parts from pathway plasmid ids.'''
        plasmid_parts = \
            self.__ice_helper.get_plasmid_parts(pathway_plasmid_ids,
                                                type_filter='PART')

        parts = {}

        for dct in plasmid_parts.values():
            parts.update(dct)

        return self.get_primers_from_parts(parts, plates, restr_enz, tm,
                                           mg_conc)

    def get_primers_from_parts(self, parts, plates, restr_enz, tm,
                               mg_conc=0.0):
        '''Get primers from parts.'''
        primers = _get_primers(parts, restr_enz, tm, mg_conc)
        return self.__get_plates(primers, plates)

    def __get_plates(self, primers, plates):
        '''Map primers to plates.'''
        primer_plates = defaultdict(list)
        new_plates = [plate.Plate('primer_1')]

        for _, values in primers.items():
            primer_id, plate_id, row_col = \
                self.__get_location(values[0], plates, new_plates)

            row, col = plate.get_indices(row_col)

            nonphospho = primer_plates[plate_id + '_primer_nonphospho']
            phospho = primer_plates[plate_id + '_primer_phospho']
            nonphospho.append(
                [row_col, primer_id + '_F', values[1][0], row, col])
            nonphospho.append(
                [row_col, primer_id + '_R', values[2][0], row, col])
            phospho.append(
                [row_col, primer_id + '_F', '/5Phos/' + values[1][0],
                 row, col])
            phospho.append(
                [row_col, primer_id + '_R', '/5Phos/' + values[2][0],
                 row, col])

        columns = ['Well Position', 'Sequence Name',
                   'Sequence', 'row', 'column']

        for plate_id in primer_plates:
            df = pd.DataFrame(primer_plates[plate_id], columns=columns)
            df.sort_values(['row', 'column'], inplace=True)
            df.drop(['row', 'column'], axis=1, inplace=True)
            primer_plates[plate_id] = df

        return primer_plates

    def __get_location(self, part, existing_plates, new_plates):
        '''Get location.'''
        parent_plasmid = self.__get_parent_plasmid(part)
        parent_plasmid_id = parent_plasmid.get_ice_id()

        part_locs = plate.find(existing_plates, parent_plasmid_id)

        if part_locs:
            # TODO: check whether part_locs can be accessed with keys by index
            plate_id = part_locs.keys()[0]
            row_col = part_locs.values()[0][0]
        else:
            # No location found: return default plate entry:
            plt = new_plates[-1]
            plate_id = plt.get_name()
            row_col = plt.add({'id': parent_plasmid_id})

        return parent_plasmid_id, plate_id, row_col

    def __get_parent_plasmid(self, part_ice):
        '''Get "parent" Plasmid from Part.'''
        part_metadata = part_ice.get_metadata()

        for parent in part_metadata['parents']:
            if parent['visible'] == 'OK':
                parent = self.__ice_helper.get_ice_entry(parent['id'])
                linked_part_ids = \
                    [linked_part['id']
                     for linked_part in parent.get_metadata()['linkedParts']]

                if len(linked_part_ids) == 2 and \
                        part_metadata['id'] in linked_part_ids:
                    return parent

        return None


def _get_primers(parts, restr_enz, tm, mg_conc):
    '''Design primers.'''
    primers = {}
    reag_conc = {seq_utils.MG: mg_conc}

    for part_id, part in parts.items():
        digest = _apply_restricts(part.get_dna(), restr_enz)['seq']

        primers[part_id] = [
            part,
            seq_utils.get_seq_by_melt_temp(digest, tm,
                                           reagent_concs=reag_conc),
            seq_utils.get_seq_by_melt_temp(Seq(digest).reverse_complement(),
                                           tm,
                                           reagent_concs=reag_conc)
        ]

    return primers


def _apply_restricts(dna, restr_enz):
    '''Apply restruction enzyme.'''
    if not restr_enz:
        return dna

    restrict_dnas = dna_utils.apply_restricts(dna, restr_enz)

    # This is a bit fudgy...
    # Essentially, return the longest fragment remaining after digestion.
    # Assumes prefix and suffix are short sequences that are cleaved off.
    restrict_dnas.sort(key=lambda x: len(x['seq']), reverse=True)
    return restrict_dnas[0]


def main(args):
    '''main method.'''
    ice_details = dict(zip(['url', 'username', 'password'], args[:3]))
    designer = PrimerDesigner(ice_details)

    plates = [plate.from_table(pd.read_csv(filename), filename)
              for filename in args[5].split(',')] \
        if args[5] != 'None' else {}

    primer_plates = \
        designer.get_primers_from_plasmids(args[6:], plates,
                                           args[3].split(','),
                                           float(args[4]))

    for plate_id, plt in primer_plates.items():
        plt.to_csv(plate_id + '.csv', index=False, encoding='utf-8')

    designer.close()


if __name__ == '__main__':
    main(sys.argv[1:])
