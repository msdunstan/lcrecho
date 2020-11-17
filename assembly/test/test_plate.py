'''
AssemblyGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=too-many-public-methods
import random
import unittest

from assembly.plate import Plate


class TestPlate(unittest.TestCase):
    '''Test class for Plate class.'''

    @classmethod
    def setUpClass(cls):
        cls.__plate_col = Plate('col', col_ord=True)
        cls.__plate_row = Plate('row')

    def test_get_set(self):
        '''Tests get and set methods.'''
        self.__test_get_set(self.__plate_col)
        self.__test_get_set(self.__plate_row)

    def test_add_col(self):
        '''Tests add method.'''
        obj1 = {'id': 'OBJ1'}
        obj2 = {'id': 'OBJ2'}
        row = 3
        col = 5

        self.__plate_col.set(obj1, row, col)
        self.__plate_col.add(obj2)
        self.assertEqual(obj2, self.__plate_col.get(row, col + 1))

    def test_add_row(self):
        '''Tests add method.'''
        obj1 = {'id': 'OBJ1'}
        obj2 = {'id': 'OBJ2'}
        row = 3
        col = 5

        self.__plate_row.set(obj1, row, col)
        self.__plate_row.add(obj2)
        self.assertEqual(obj2, self.__plate_row.get(row + 1, col))

    def test_find(self):
        '''Tests find methods.'''
        self.__test_find(self.__plate_col)
        self.__test_find(self.__plate_row)

    def test_idx_row_col(self):
        '''Tests get_idx and get_row_col methods.'''
        self.__test_idx_row_col(self.__plate_col)
        self.__test_idx_row_col(self.__plate_row)

    def __test_get_set(self, plate):
        '''Tests get, set and get by well methods.'''
        obj = {'id': 'OBJ1'}
        row = 3
        col = 5

        plate.set(obj, row, col)
        self.assertEqual(obj, plate.get(row, col))
        self.assertEqual(obj, plate.get_by_well('D6'))

    def __test_find(self, plate):
        '''Tests find method.'''
        obj = {'id': 'OBJ1'}
        plate.set(obj, 3, 5)
        plate.set(obj, 6, 10)
        self.assertEqual(sorted(['D6', 'G11']),
                         sorted(self.__plate_col.find(obj)))

    def __test_idx_row_col(self, plate):
        '''Tests get_idx and get_row_col methods.'''
        rows, cols = plate.shape()

        for _ in range(0, 16):
            idx = random.randint(0, rows * cols)
            row, col = plate.get_row_col(idx)
            self.assertEqual(idx, plate.get_idx(row, col))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
