'''
PathwayGenie (c) GeneGenie Bioinformatics Ltd. 2018

PathwayGenie is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-many-arguments
import RNA


def run(cmd, sequences, temp=37.0, dangles=2, energy_gap=None,
        bp_x=None, bp_y=None):
    '''Runs ViennaRNA.'''
    sequences = [str(seq) for seq in sequences]

    if cmd == 'mfe':
        return _mfe(sequences, temp, dangles)

    if cmd == 'subopt':
        return _subopt(sequences, energy_gap, temp, dangles)

    if cmd == 'energy':
        return _energy(sequences, bp_x, bp_y, temp, dangles)

    return None


def _mfe(sequences, temp, dangles):
    '''mfe.'''
    model = RNA.md()
    model.temperature = temp
    model.dangles = dangles
    result = RNA.fold_compound(sequences[0], model).mfe()
    bp_x, bp_y = _get_numbered_pairs(result[0])

    if bp_x and bp_y:
        return [result[1]], [bp_x], [bp_y]

    return [float('NaN')], [[]], [[]]


def _subopt(sequences, energy_gap, temp, dangles):
    '''subopt.'''
    model = RNA.md()
    model.temperature = temp
    model.dangles = dangles

    results = \
        RNA.fold_compound('&'.join(sequences), model).subopt(int(energy_gap))

    energies = []
    bp_xs = []
    bp_ys = []

    for result in results:
        bp_x, bp_y = _get_numbered_pairs(result.structure)

        if bp_x and bp_y:
            energies.append(result.energy)
            bp_xs.append(bp_x)
            bp_ys.append(bp_y)

    return energies, bp_xs, bp_ys


def _energy(sequences, bp_x, bp_y, temp, dangles):
    '''energy.'''
    model = RNA.md()
    model.temperature = temp
    model.dangles = dangles
    sequence = '&'.join(sequences)
    structure = _get_brackets([len(seq) for seq in sequences], bp_x, bp_y)
    return RNA.fold_compound(sequence, model).eval_structure(structure)


def _get_numbered_pairs(bracket_str):
    '''_get_numbered_pairs'''
    bracket_count = bracket_str.count(')')

    if not bracket_count:
        return [None, None]

    bp_x = []
    bp_y = [None for _ in range(bracket_count)]
    last_nt_x = []
    strand_num = 0

    for pos, letter in enumerate(bracket_str):
        if letter == '(':
            bp_x.append(pos - strand_num)
            last_nt_x.append(pos - strand_num)
        elif letter == ')':
            nt_x = last_nt_x.pop()
            nt_x_pos = bp_x.index(nt_x)
            bp_y[nt_x_pos] = pos - strand_num
        elif letter == '&':
            strand_num += 1

    return [[pos + 1 for pos in bp_x], [pos + 1 for pos in bp_y]]


def _get_brackets(seq_lens, bp_x, bp_y):
    '''_get_brackets'''
    bp_x = [pos - 1 for pos in bp_x]
    bp_y = [pos - 1 for pos in bp_y]
    brackets = []
    counter = 0

    for seq_len in seq_lens:
        for pos in range(counter, seq_len + counter):
            if pos in bp_x:
                brackets.append('(')
            elif pos in bp_y:
                brackets.append(')')
            else:
                brackets.append('.')
        counter += seq_len

    return ''.join(brackets)
