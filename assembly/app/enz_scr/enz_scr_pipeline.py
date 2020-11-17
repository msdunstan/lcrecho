'''
AssemblyGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=wrong-import-order
from collections import Counter
import os
import sys
from time import gmtime, strftime

from assembly import pipeline, plate, worklist
from assembly.graph_writer import GraphWriter
import pandas as pd


class EnzymeScreenWriter(GraphWriter):
    '''Class for generating enzyme screen worklist graphs.'''

    def __init__(self, df, input_plates, output_name='enz_scr',
                 replicates=2):
        self.__df = df
        self.__input_plates = input_plates
        self.__replicates = replicates
        GraphWriter.__init__(self, output_name)

    def _initialise(self):
        assay_ids = []

        for _, row in self.__df.iterrows():
            for _ in range(self.__replicates):
                part_ids = [part_id.strip()
                            for part_id in str(row['Part ID (s)']).split('+')]

                locations = row[['Lysate 1 location', 'Lysate 2 location']]

                self.__add_lysate(zip(part_ids, locations))

                self.__add_assay(part_ids, row['Substrate'], assay_ids)

    def __add_lysate(self, part_locs):
        '''Add lysate.'''
        for part_loc in part_locs:
            if part_loc[1]:
                if part_loc[1][0] in self.__input_plates:
                    input_plate = self.__input_plates[part_loc[1][0]]
                else:
                    input_plate = plate.Plate(part_loc[1][0])
                    self.__input_plates[part_loc[1][0]] = input_plate

                row, col = plate.get_indices(part_loc[1][1:])
                input_plate.set({'id': part_loc[0] + '_lys'}, row, col)

    def __add_assay(self, part_ids, substrate, assay_ids):
        '''Add assay.'''
        assay_id = tuple(part_ids + [substrate])
        assay_ids.append(assay_id)
        count = str(Counter(assay_ids)[assay_id])

        assay = self._add_vertex('_'.join(part_ids + [substrate, count]),
                                 {'is_reagent': False})

        substrate = self._add_vertex(substrate, {'is_reagent': False})

        self._add_edge(substrate, assay, {'Volume': 1.0})

        for part_id in part_ids:
            part = self._add_vertex(part_id + '_lys',
                                    {'is_reagent': False})

            self._add_edge(part, assay, {'Volume': 19.0 / len(part_ids)})


def _get_substrate_plates(wklst_dfs, substrates, name='substrate_treatment'):
    '''Get substrate plate.'''
    plate_dfs = []

    for count, wklst_df in enumerate(wklst_dfs):
        substrates_df = wklst_df[wklst_df['ComponentName'].isin(substrates)]
        substrates_df = substrates_df[['ComponentName',
                                       'DestinationPlateWell']]
        substrates_df.columns = ['id', 'well']
        plate_dfs.append(plate.from_table(substrates_df, name + '_' +
                                          str(count + 1)))

    return plate_dfs


def main(args):
    '''main method.'''
    recipe_df = pd.read_csv(args[0], dtype=str)

    dte = strftime("%y%m%d", gmtime())

    for name, group_df in recipe_df.groupby('Project'):
        out_dir_name = os.path.join(os.path.join('out', dte + args[1]), name)
        input_plates = pipeline.get_input_plates('data/plates')

        writers = [EnzymeScreenWriter(group_df, input_plates,
                                      dte + 'ENZ' + name[:1].upper())]

        pipeline.run(writers, input_plates, {'reagents': 'reagents'},
                     out_dir_name)

        wklst_dfs = worklist.format_worklist(out_dir_name)
        plts = _get_substrate_plates(wklst_dfs, group_df['Substrate'].unique())

        for plt in plts:
            plate_dir = os.path.join(out_dir_name, 'plates')

            if not os.path.exists(plate_dir):
                os.makedirs(plate_dir)

            plt.to_csv(plate_dir)


if __name__ == '__main__':
    main(sys.argv[1:])
