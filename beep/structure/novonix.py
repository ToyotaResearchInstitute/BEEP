import pandas as pd
import os
from beep.structure.base import BEEPDatapath
from beep.conversion_schemas import NOVONIX_CONFIG
from beep import logger, VALIDATION_SCHEMA_DIR


class NovonixDatapath(BEEPDatapath):
    """A BEEPDatapath for ingesting and structuring Novonix data files.
    """

    @classmethod
    def from_file(cls, path):
        """Create a NovonixDatapath from a raw Novonix cycler file.

        Args:
            path (str, Pathlike): file path for novonix file.

        Returns:
            (NonovixDatapath)
        """
        #format raw data
        with open(path, "rb") as f:
            i = 1
            search_lines = 200
            header_starts_line = None
            while header_starts_line is None:
                line = f.readline()
                if b'Cycle Number' in line:
                    header_starts_line = i
                i += 1
                if i > search_lines:
                    raise LookupError("Unable to find the header line in first {} lines of file".format(search_lines))
        raw = pd.read_csv(path, header=None)
        raw.dropna(axis=0, how='all', inplace=True)
        data = raw.iloc[header_starts_line-1:]
        data = data[0].str.split(',', expand = True)
        headers = data.iloc[0]
        data = pd.DataFrame(data.values[1:], columns=headers, index=None)

        # format columns
        map = NOVONIX_CONFIG['data_columns']
        name_map = {i:map[i]['beep_name'] for i in map}
        type_map = {j:map[j]['data_type'] for j in map}
        data = data.astype(type_map)
        data.rename(name_map, axis="columns", inplace=True)

        # format capacity and energy
        # rest = data['step_type'] == '0'
        cc_charge = data['step_type_num'] == '1'
        cc_discharge = data['step_type_num'] == '2'
        cccv_charge = data['step_type_num'] == '7'
        cv_hold_discharge = data['step_type_num'] == '8'
        cccv_discharge = data['step_type_num'] == '9'
        cccv_hold_discharge = data['step_type_num'] == '10'

        data['charge_capacity'] = data[cc_charge | cccv_charge]['capacity'].astype('float')
        data['discharge_capacity'] = data[cc_discharge | cv_hold_discharge | cccv_discharge | cccv_hold_discharge][
            'capacity'].astype('float')
        data['charge_energy'] = data[cc_charge | cccv_charge]['energy'].astype('float')
        data['discharge_energy'] = data[cc_discharge | cv_hold_discharge | cccv_discharge | cccv_hold_discharge][
            'energy'].astype('float')

        # add step type #todo add schema
        step_map = {0: 'discharge',
                    1: 'charge',
                    2: 'discharge',
                    7: 'charge',
                    8: 'discharge',
                    9: 'discharge',
                    10: 'discharge'}
        data['step_type'] = data['step_type_num'].replace(step_map)
        data.fillna(0)

        #paths
        metadata = {}
        paths = {
            "raw": path,
            "metadata": path if metadata else None
        }
        # todo convert time
        # validation
        schema = os.path.join(VALIDATION_SCHEMA_DIR, "schema-novonix.yaml")

        return cls(data, metadata, paths=paths, schema=schema)

    # todo change base for more than 1 cycle index