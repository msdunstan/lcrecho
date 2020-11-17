'''
AssemblyGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
# pylint: disable=useless-object-inheritance
# pylint: disable=wrong-import-order
from collections import defaultdict
import re

from synbiochem.utils.ice_utils import ICEClientFactory


class ICEHelper(object):
    '''Helper class for accessing ICE.'''

    def __init__(self, ice_url, ice_username, ice_password):
        self.__ice_factory = ICEClientFactory()
        self.__ice_client = self.__ice_factory.get_ice_client(ice_url,
                                                              ice_username,
                                                              ice_password)
        self.__ice_entries = {}

    def close(self):
        '''Close.'''
        self.__ice_factory.close()

    def get_plasmid_parts_designs(self, design_numbers):
        '''Get parts from design numbers.'''
        parts = {}

        for design_number in design_numbers:
            plasmid_ids = [entry['entryInfo']['partId']
                           for entry in
                           self.__ice_client.search_design(design_number)]
            parts.update(self.get_plasmid_parts(plasmid_ids))

        return parts

    def get_plasmid_parts(self, plasmid_ids, type_filter=None):
        '''Get parts from part plasmid ids.'''
        parts = defaultdict(dict)

        for plasmid_id in plasmid_ids:
            parts[plasmid_id] = {}

            for part_ice in self.__get_parts(plasmid_id):
                part_id = part_ice.get_ice_id()

                if part_id not in parts and \
                        (not type_filter or
                         (part_ice.get_parameter('Type') and
                          re.match(type_filter,
                                   part_ice.get_parameter('Type')))):
                    parts[plasmid_id][part_id] = part_ice

        assert len(plasmid_ids) == len(parts)

        return parts

    def get_ice_entry(self, ice_id):
        '''Get ICE entry.'''
        if ice_id not in self.__ice_entries:
            ice_entry = self.__ice_client.get_ice_entry(ice_id)
            self.__ice_entries[ice_id] = ice_entry

        return self.__ice_entries[ice_id]

    def __get_parts(self, plasmid_id):
        ''''Get parts from plasmid id.'''
        ice_entry = self.get_ice_entry(plasmid_id)

        part_ids = [part['id']
                    for part in ice_entry.get_metadata()['linkedParts']]

        return [self.get_ice_entry(part_id) for part_id in part_ids]
