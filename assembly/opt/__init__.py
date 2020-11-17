'''
AssemblyGenie (c) University of Manchester 2018
All rights reserved.
@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=wrong-import-order
from itertools import zip_longest
import math
import random

from assembly import plate
import pandas as pd


def score(df, channels=8):
    '''Score.'''
    min_score = (math.ceil(len(df.index) / 8) - 1) * 2

    scores = [_score_group(pd.DataFrame([list(row) for row in rows
                                         if row is not None],
                                        columns=df.columns))
              for rows in _grouper(channels, df.values)]

    return sum(scores) - min_score


def get_shuffled_wklst(num_wells):
    '''Get shuffled worklist.'''
    idx = list(range(0, num_wells))
    df = _get_wklst(idx, idx)
    return df.sample(frac=1)


def get_semirandom_wklst(num_wells, plate_size=96):
    '''Get semi-random worklist.'''
    idxs = random.sample(range(0, plate_size), num_wells)
    return _get_wklst(idxs, idxs)


def get_random_wklst(num_wells, plate_size=96):
    '''Get random worklist.'''
    src_idxs = random.sample(range(0, plate_size), num_wells)
    dest_idxs = random.sample(range(0, plate_size), num_wells)
    return _get_wklst(src_idxs, dest_idxs)


def _score_group(group_df):
    '''Score group.'''
    return _score_cols(group_df, ['src_row', 'src_col']) + \
        _score_cols(group_df, ['dest_row', 'dest_col'])


def _score_cols(group_df, columns):
    '''Score group by either src or dest (specified by columns).'''
    col_score = 0
    head_col = 0
    head_pos = [[row, head_col] for row in range(0, 8)]
    curr_pos = [list(pos) for pos in group_df[columns].values]

    while True:
        # Eliminate wells in correct position (i.e. pick up or dispense:
        head_curr = [list(pair)
                     if (pair[0] != pair[1]) or (pair[1] == [-1, -1])
                     else [pair[0], [-1, -1]]
                     for pair in zip(head_pos, curr_pos)]

        # If all wells have been picked up / dispensed:
        next_idx = float('NaN')

        for idx, pair in enumerate(head_curr):
            if pair[1] != [-1, -1]:
                next_idx = idx
                break

        if math.isnan(next_idx):
            break

        # Move col:
        head_col = head_curr[next_idx][1][1]
        col_score += abs(head_col - head_curr[next_idx][0][1])

        # Move row:
        head_row_offset = head_curr[next_idx][1][0] - head_curr[next_idx][0][0]
        col_score += abs(head_row_offset)
        head_pos = [[row + head_row_offset, head_col] for row, _ in head_pos]
        curr_pos[next_idx] = [-1, -1]

    return col_score


def _grouper(n, iterable, fillvalue=None):
    '''Collect data into fixed-length chunks or blocks.'''
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def _get_wklst(src_idxs, dest_idxs):
    '''Get worklist.'''
    data = []

    for src, dest in zip(src_idxs, dest_idxs):
        data.append(_get_well_details(src) + _get_well_details(dest))

    cols = ['src_well', 'src_idx', 'src_row', 'src_col',
            'dest_well', 'dest_idx', 'dest_row', 'dest_col']

    return pd.DataFrame(data, columns=cols)


def _get_well_details(well_idx):
    '''Get well details.'''
    row, col = plate.get_row_col(well_idx)

    return [plate.get_well_name(row, col), well_idx, row, col]
