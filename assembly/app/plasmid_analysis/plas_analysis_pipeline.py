'''
AssemblyGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=wrong-import-order
from functools import partial
import os
import re
import sys
from time import gmtime, strftime

from assembly import pipeline, plate, worklist
from assembly.app.plasmid_analysis import colony_pcr, colony_qc
import pandas as pd


def _get_colony_dfs(dir_name):
    '''Get colony plates.'''
    col_dfs = []

    for dirpath, _, filenames in os.walk(os.path.abspath(dir_name)):
        for filename in filenames:
            filename = os.path.join(dirpath, filename)

            if filename.endswith('.csv'):
                col_dfs.append(pd.read_csv(filename))

    return col_dfs


def _get_colony_plates(colony_df, input_plates):
    '''Get colony plates.'''
    barcode_plates = _get_barcode_plates(input_plates)

    colony_df.drop_duplicates(subset=['DWPBarcode', 'DWPWell'],
                              keep='last', inplace=True)
    colony_df.rename(columns={'DWPWell': 'well',
                              'ColonyID': 'id'}, inplace=True)

    colony_df['actual_ice_id'] = \
        colony_df['id'].apply(lambda x: x.split('_')[0])

    colony_df['plate_idx'] = \
        colony_df['DWPBarcode'].apply(lambda x: int(x.split('_')[-1]))

    _get_barcodes_plates = partial(_get_barcodes, barcode_plates)
    colony_df = colony_df.apply(_get_barcodes_plates, axis=1)

    plate_idx = colony_df['plate_idx'].iloc[0]
    plate_name = colony_df['DWPBarcode'].iloc[0]
    plate_df = colony_df[['well', 'id']]
    plt = plate.from_table(plate_df, plate_name)
    colony_ids = plate_df.values.tolist()

    return plt, plate_idx, plate_name, colony_ids, colony_df


def _get_barcode_plates(input_plates):
    '''Get barcode plates.'''
    barcode_plates = {}

    for plt_id, plt in input_plates.items():
        search = re.search(r'(?<=ONT)(\d+)', plt_id)

        if search:
            barcode_plates[int(search.group(1))] = plt

    return barcode_plates


def _get_barcodes(barcode_plates, row):
    '''Get barcodes from colony picker table.'''
    properties = barcode_plates[row['plate_idx']].get_by_well(row['well'])
    row['forward'] = properties['forward']
    row['reverse'] = properties['reverse']
    return row


def _get_frag_anal_labels(plt, name, out_dir_name):
    '''Get fragment analysis plates.'''
    components = plt.get_all()
    labels = [[plate.get_idx(*plate.get_indices(pos), col_ord=True) + 1,
               pos, vals['id']]
              for pos, vals in components.items()]

    df = pd.DataFrame(labels)
    df.to_csv(os.path.join(out_dir_name, name + '_frag_anal.csv'),
              index=False, header=False)


def main(args):
    '''main method.'''
    colony_dfs = []
    dte = strftime("%y%m%d", gmtime())
    out_dir_name = os.path.join(args[3], dte + args[1])

    for colony_df in _get_colony_dfs(args[4]):
        # Ignore colony in H12 (if it exists).
        # This position is reserved for the fragment analyser ladder.
        colony_df.drop(colony_df[colony_df['DWPWell'] == 'H12'].index,
                       inplace=True)

        input_plates = pipeline.get_input_plates(args[0])

        # Parse colony pick output:
        colony_plate, plate_idx, plate_name, colony_ids, colony_df = \
            _get_colony_plates(colony_df, input_plates)

        input_plates[plate_name] = colony_plate

        plate_dir_name = os.path.join(out_dir_name, 'plate_' + str(plate_idx))

        #  Write PCR worklists:
        writers = [colony_pcr.ColonyPcrWriter(colony_ids, plate_idx,
                                              dte + 'COL' + args[1]),
                   colony_qc.ColonyQcWriter([colony_id[1]
                                             for colony_id in colony_ids],
                                            output_name=dte + 'FPT' + args[1])]
        pipeline.run(writers, input_plates, {'reagents': args[2]},
                     plate_dir_name)
        worklist.format_worklist(plate_dir_name)

        # Generate fragment analyse labels:
        colony_dfs.append(colony_df)
        _get_frag_anal_labels(colony_plate, plate_name, plate_dir_name)

    # Write NGS files:
    barcodes_df = pd.concat(colony_dfs, axis=0, ignore_index=True)

    barcodes_df.to_csv(os.path.join(
        plate_dir_name, 'barcodes.csv'), index=False)

    pd.DataFrame(barcodes_df['actual_ice_id'].unique()).to_csv(os.path.join(
        out_dir_name, 'ice_ids.txt'), index=False, header=False)


if __name__ == '__main__':
    main(sys.argv[1:])
