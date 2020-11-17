'''
AssemblyGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
import os
import sys
from time import gmtime, strftime

from assembly import pipeline, worklist
from assembly.app.lcr2 import lcr2, part_pcr, part_qc, part_dig, utils


def main(args):
    '''main method.'''
    assert len(args[3]) < 4

    ice_helper = utils.ICEHelper(args[0], args[1], args[2])

    if args[4] == 'True':
        plasmid_parts = ice_helper.get_plasmid_parts(args[5:])
    else:
        plasmid_parts = ice_helper.get_plasmid_parts_designs(args[5:])

    parts_ice = {ice_id: part_ice
                 for _, parts_map in plasmid_parts.items()
                 for ice_id, part_ice in parts_map.items()
                 if part_ice.get_parameter('Type') != 'DOMINO'}
    part_ids = parts_ice.keys()

    dte = strftime("%y%m%d", gmtime())

    part_vol = 1.0
    pcr_numbers, pcr_df = part_pcr.get_pcr_numbers(plasmid_parts, part_vol)

    writers = [part_pcr.SpecificPartPcrWriter(parts_ice, pcr_numbers,
                                              ice_helper,
                                              dte + 'PCR' + args[3]),
               part_dig.PartDigestWriter(part_ids, pcr_numbers,
                                         dte + 'DIG' + args[3]),
               part_qc.PartQcWriter(part_ids, pcr_numbers,
                                    dte + 'FPT' + args[3]),
               lcr2.DominoPoolWriter(plasmid_parts, dte + 'DOM' + args[3]),
               part_dig.PartPoolWriter(pcr_numbers,
                                       dte + 'POO' + args[3]),
               lcr2.LcrWriter(plasmid_parts, part_vol, dte + 'LCR' + args[3])]

    input_plates = pipeline.get_input_plates('data/plates')
    out_dir_name = os.path.join('out', dte + args[3])
    pipeline.run(writers, input_plates,
                 {'reagents': 'reagents'}, out_dir_name)

    worklist.format_worklist(out_dir_name)

    pcr_df.to_csv(os.path.join(out_dir_name, 'pcr_details.csv'))

    ice_helper.close()


if __name__ == '__main__':
    main(sys.argv[1:])
