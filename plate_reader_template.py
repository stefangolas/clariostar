# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 09:40:15 2023

@author: stefa
"""
import os
import csv

class GeneralPlateReader:
    
    def __init__(self, simulating = False):
        self.simulating = simulating
    
    def __enter__(self):
        pass
    
    def run_protocols(self, protocol_names, plate_id_1):
        if self.simulating:
            return [simulated_plate_data] * len(protocol_names)


class PlateData:

    header_divider = 'End_of_header'

    def __init__(self, path):
        try:
            path = os.path.abspath(path())
            self._path_getter = lambda: path
        except TypeError:
            path = os.path.abspath(path)
            self._path_getter = lambda: path
        except IOError:
            self._path_getter = path
        self._path = None
        self._text = None
        self._csvrows = None
        self._blockdata = None
        self._header_namespace = None

    @property
    def path(self):
        if self._path is None:
            try:
                self._path = self._path_getter()
            except IOError:
                pass
        return self._path

    @property
    def text(self):
        self._assert_file_exists()
        if self._text is None:
            with open(self._path) as f:
                self._text = str(f.read())
        return self._text

    @property
    def header(self):
        if self._header_namespace:
            return self._header_namespace
        header = {}
        header_str = str(self.csv_rows).split(PlateData.header_divider)[0]
        property_patterns = [('test_name', r"(?<=Testname:\s)(\S*)(?='\])"),
                             ('date', r"(?<=Date:\s)(\d*)/(\d*)/(\d*)\s"),
                             ('time', r"(?<=Time:\s)(\d*):(\d*):(\d*)\s(\SM)"),
                             ('num_channels', r"(?<=No.\sof\sChannels\s/\sMultichromatics:\s)(\d*)"),
                             ('num_cycles', r"(?<=No.\sof\sCycles:\s)(\d*)"),
                             ('configuration', r"(?<=Configuration:\s)(\S*)(?='\])"),
                             ('focal_height', r"(?<=Focal\sheight\s\[mm\]:\s)(\d*.\d*)"),
                             ('plate_ids', r"ID1:.*ID2:.*ID3:.*?(?='\])")]
        for prop, pattern in property_patterns:
            match = re.compile(pattern).search(header_str)
            header[prop] = match.group(0).strip() if match else None
        ids = header['plate_ids']
        if ids is not None:
            id_segments = ids.split('ID')[1:] # element 0 will always be ''
            header['plate_ids'] = tuple(seg[2:].strip() for seg in id_segments) # Remove id number and colon
        self._header_namespace = types.SimpleNamespace(**header)
        return self._header_namespace

    @property
    def csv_rows(self):
        self._assert_file_exists()
        if self._csvrows is None:
            with open(self._path) as f:
                self._csvrows = list(csv.reader(f))
        return self._csvrows

    @property
    def data_array(self):
        if self._blockdata is None:
            block_data_text = self.text.split(PlateData.header_divider)[1].strip()
            chromblocks = block_data_text.split('Chromatic:')
            if chromblocks.pop(0) != '':
                logging.error('THERE WAS A PARSE ERROR WITH ' + str(self._path))
                exit()
            lettermap = {a:i for i, a in enumerate('ABCDEFGHIJKLMNO')}
            chrom_data_list = []
            for chromblock in chromblocks:
                cycle_data_list = []
                cycleblocks = chromblock.split('Cycle:')[1:]
                for cycleblock in cycleblocks:
                    maxcol = 0
                    maxrow = 0
                    well_vals = np.zeros((6000, 6000))
                    for line in cycleblock.split('\n'):
                        try:
                            well_id, val_str = line.split(':')
                            col, row = PlateData.well_id_coords(well_id)
                            if col > maxcol:
                                maxcol = col
                            if row > maxrow:
                                maxrow = row
                            well_vals[col, row] = float(val_str)
                        except ValueError:
                            continue
                        except KeyError:
                            continue
                    cycle_data_list.append(well_vals[:maxcol + 1, :maxrow + 1])
                chrom_data_list.append(cycle_data_list)
            self._blockdata = np.array(chrom_data_list)
        return self._blockdata

    def value_at(self, col, row, chromatic=0, cycle=0):
        coords = (chromatic, cycle, col, row)
        if not all(0 <= idx < bound for idx, bound in zip(coords, self.data_array.shape)):
            raise ValueError('Reading position [' + ', '.join((c_name + '=' + str(coord) for c_name, coord
                    in zip(('chromatic', 'cycle', 'col', 'row'), coords))) +
                    '] out of bounds of plate data array with shape ' + str(self.data_array.shape))
        return self.data_array[chromatic, cycle, col, row]

    def _assert_file_exists(self):
        if self._path is None:
            try:
                self._path = self._path_getter()
            except IOError:
                raise PlateDataPathError('Path to PlateData file not resolved yet')
        if not os.path.isfile(self._path):
            raise PlateDataPathError('Path to PlateData file ' + str(self._path) + ' does not exist')

    def wait_for_file(self, timeout=float('inf')):
        start_time = time.time()
        while time.time() - start_time < timeout:
            try:
                self._assert_file_exists()
                return
            except PlateDataPathError:
                time.sleep(1)
        raise PlateDataPathError('PlateData timed out while waiting for ' + str(self._path) + ' to appear')

    def reload(self):
        self._path = self._text = self._csvrows = self._blockdata = self._header_namespace = None

    @staticmethod
    def parse_well_id(well_id):
        # convert a well identifier ('AB0025') to a tuple (row letters, col int) ('AB', 25)
        for i in range(len(well_id)):
            try:
                return (well_id[:i], int(well_id[i:]))
            except ValueError:
                continue
        else:
            raise ValueError('Well id ' + well_id + ' not correctly formatted')

    @staticmethod
    def fixed_width_well_id(well_id, width=3):
        # Pad the integer (column) portion of a well identifier ('B2') with zeros so that it is width characters long ('B02')
        row_letters, colnum = PlateData.parse_well_id(well_id)
        col_num_str = str(colnum)
        num_zeros = width - len(row_letter) - len(col_num_str)
        if num_zeros < 0:
            raise ValueError('Well id ' + well_id + ' cannot fit into width ' + str(width))
        return row_letter + '0'*num_zeros + col_num_str

    @staticmethod
    def well_id_coords(well_id):
        # Return zero-indexed col, row of well id ('F03' becomes (2, 5))
        letters, colnum = PlateData.parse_well_id(well_id)
        if not letters:
            raise ValueError('Row letter label cannot be empty')
        if len(letters) != 1:
            raise NotImplementedError('Well ids with more than one row letter not supported yet')
        row = string.ascii_uppercase.index(letters)
        col = colnum - 1 # wells are 1-indexed
        return col, row


sim_data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dummy_platedata.csv')
simulated_plate_data = PlateData(sim_data_path)
