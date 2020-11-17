'''
AssemblyGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=too-few-public-methods
import pandas as pd


class Optimiser():
    '''Recipes optimiser.'''

    def __init__(self, ingredients):
        self.__df = pd.DataFrame()
        self.__reagents = []
        self.__intermediates = []
        self.__products = []
        self.__get_components(ingredients[0], ingredients[1], ingredients[2])
        self.__drop()

    def optimise(self):
        '''Optimise dataframe.'''
        while True:
            mask = None
            max_match_col = pd.Series(0.0, index=self.__df.index)

            for idx, col1 in enumerate(self.__df.columns):
                for col2 in self.__df[self.__df.columns[idx + 1:]]:
                    match_col = self.__df[col1] * \
                        (self.__df[col1] == self.__df[col2])

                    if sum(match_col.astype(bool)) > \
                            sum(max_match_col.astype(bool)):
                        mask = pd.DataFrame(0.0, columns=self.__df.columns,
                                            index=self.__df.index)
                        mask[col1] = match_col.values
                        max_match_col = match_col

                    if match_col.any() and (match_col == max_match_col).all():
                        mask[col2] = match_col.values

            if mask is None or max(mask.astype(bool).sum()) < 2:
                break

            self.__add_intermediate(mask, max_match_col)

    def get_matrix(self):
        '''Gets matrix.'''
        return self.__df

    def get_reagents(self):
        '''Gets reagents.'''
        return self.__reagents

    def save_matrix(self, outfile):
        '''Saves matrix as csv.'''
        self.__df.to_csv(outfile)

    def __get_components(self, comps, vol, is_reagent, dest=None):
        '''Gets components.'''
        if isinstance(comps, str):
            comp_id = comps
            self.__add_row_col(comp_id)
        else:
            comp_id = self.__get_intermediate_name(vol)

            self.__add_row_col(comp_id)

            for comp in comps:
                self.__get_components(comp[0], comp[1], comp[2], comp_id)

        if is_reagent and comp_id not in self.__reagents:
            self.__reagents.append(comp_id)

        if dest:
            self.__df[dest][comp_id] = vol

    def __drop(self):
        '''Drop empty columns and rows.'''
        self.__df = self.__df[self.__df.columns[(self.__df != 0).any()]]
        self.__df = self.__df[(self.__df.T != 0).any()]

    def __add_row_col(self, comp_id):
        '''Add new row and column.'''
        if self.__df.empty:
            self.__df[comp_id] = pd.Series([0.0], index=[comp_id])
        elif comp_id not in self.__df:
            # Add row:
            new_row = pd.Series([0.0] * len(self.__df.columns),
                                index=self.__df.columns)
            new_row.name = comp_id
            self.__df = self.__df.append(new_row)

            # Add col:
            self.__df[comp_id] = pd.Series([0.0] * len(self.__df.index),
                                           index=self.__df.index)

    def __add_intermediate(self, mask, max_match_col):
        '''Adds an intermediate.'''
        self.__df = self.__df - mask
        new_id = self.__get_intermediate_name()
        self.__df[new_id] = max_match_col
        new_row = mask.any().astype(float)
        new_row.name = new_id
        self.__df = self.__df.append(new_row)
        self.__df = self.__df.fillna(0)
        self.__drop()

    def __get_intermediate_name(self, intermediate=True):
        '''Get unique intermediate name.'''
        if intermediate:
            int_id = 'i' + str(len(self.__intermediates) + 1)
            self.__intermediates.append(int_id)
        else:
            int_id = 'p' + str(len(self.__products) + 1)
            self.__products.append(int_id)

        return int_id


def main():
    '''main method.'''
    ingredients = (
        (((((('K', 1, False), ('L', 2, False)), 3, False), ('D', 4, False),
           ('E', 5, False), ('F', 6, False)), 0, False),
         ((('D', 4, False), ('E', 5, False), ('F', 6, False),
           ('Y', 10, False)), 0, False),
         ((('D', 4, False), ('E', 5, False), ('F', 13, False),
           ('Z', 14, False)), 0, False)
         ), 0, False)

    optim = Optimiser(ingredients)
    optim.save_matrix('init.csv')
    optim.optimise()
    optim.save_matrix('optim.csv')


if __name__ == '__main__':
    main()
