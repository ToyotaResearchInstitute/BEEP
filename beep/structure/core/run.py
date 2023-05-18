import os
from typing import Optional, Union

import tqdm
import dask.bag as bag
import pandas as pd
import numpy as np

from monty.json import MSONable
from monty.serialization import loadfn, dumpfn
from dask.diagnostics import ProgressBar

from beep import logger
from beep.structure.diagnostic import DiagnosticConfig
from beep.structure.core.cycles_container import CyclesContainer
from beep.structure.core.util import (
    label_chg_state,
    get_max_paused_over_threshold,
    get_cv_stats,
    DFSelectorIndexError,
    CVStatsError,
    TQDM_STYLE_ARGS
)
from beep.structure.core.interpolate import interpolate_cycle, CONTAINER_CONFIG_DEFAULT
from beep.structure.core.validate import SimpleValidator

"""
Things assumed to make sure it all works correctly with checking:
    1. All cycles have a unique cycle index, and no cycles have more than one cycle index.
    2. All steps have a unique step index, and no steps have more than one step index.


# THINGS A RUN MUST IMPLEMENT
 - from_file (method)
 - conversion schema
    - must have 3 root keys: file_pattern, data_columns, data_types. 
        Optionally, can have metadata_fields
 - validation schema (class attribute)

"""


class Run(MSONable):
    """
    A Run object represents an entire cycler run as well as it's structured (interpolated)
    data. It is the top level object in the structured data hierarchy.

    A run object has its own config. This config mainly determines what columns
    will be kept and what data types they will possess.

    Args:
        raw_cycle_container (CyclesContainer): CyclesContainer object containing raw data
        structured_cycle_container (Optional[CyclesContainer], optional): CyclesContainer object containing structured data. Defaults to None.
        metadata (dict, optional): Dictionary of metadata. Defaults to None.
        schema (dict, optional): Dictionary to perform validation.
        paths (dict, optional): Dictionary of paths from which this object was derived. 
            Useful for keeping track of what Run file corresponds with what cycler run
            output file.
    """
    # Set the default datatypes for ingestion and raw data
    # to the same used for interpolated data (for simplicity)
    DEFAULT_DTYPES = CONTAINER_CONFIG_DEFAULT["dtypes"]

    # Basic LiFePO4 validation, as a backup/default
    DEFAULT_VALIDATION_SCHEMA = {
        'charge_capacity': {
            'schema': {
                'max': 2.0,
                'min': 0.0,
                'type': 'float'
            },
            'type': 'list'
        },
        'cycle_index': {
            'schema': {
                'min': 0,
                'max_at_least': 1,
                'type': 'integer'
            },
            'type': 'list'
        },
        'discharge_capacity': {
            'schema': {
                'max': 2.0,
                'min': 0.0,
                'type': 'float'
            },
            'type': 'list'
        },
        'temperature': {
            'schema': {
                'max': 80.0,
                'min': 20.0,
                'type': 'float'
            },
            'type': 'list'
        },
        'test_time': {
            'schema': {
                'type': 'float'
            },
            'type': 'list'
        },
        'voltage': {
            'schema': {
                'max': 3.8,
                'min': 0.0,
                'type': 'float'
            },
            'type': 'list'
        }
    }

    # Data types for all summaries (diag and regular)
    SUMMARY_DTYPES = {
        'cycle_index': 'int32',
        'discharge_capacity': 'float64',
        'charge_capacity': 'float64',
        'discharge_energy': 'float64',
        'charge_energy': 'float64',
        'dc_internal_resistance': 'float32',
        'temperature_maximum': 'float32',
        'temperature_average': 'float32',
        'temperature_minimum': 'float32',
        'date_time_iso': 'object',
        'energy_efficiency': 'float32',
        'charge_throughput': 'float32',
        'energy_throughput': 'float32',
        'charge_duration': 'float32',
        'time_temperature_integrated': 'float64',
        'paused': 'int32',
        'CV_time': 'float32',
        'CV_current': 'float32',
        'CV_capacity': 'float32',
        "coulombic_efficiency": "float64"
    }

    def __init__(
            self,
            raw_cycle_container: Optional[CyclesContainer] = None,
            structured_cycle_container: Optional[CyclesContainer] = None,
            diagnostic: Optional[DiagnosticConfig] = None,
            metadata: Optional[dict] = None,
            schema: Optional[dict] = None,
            paths: Optional[dict] = None,
            summary_cycles: Optional[pd.DataFrame] = None,
            summary_diagnostic: Optional[pd.DataFrame] = None,
    ):
        self.raw = raw_cycle_container
        self.structured = structured_cycle_container
        if diagnostic:
            # This is needed because the setter for diagnostic
            self.diagnostic = diagnostic
        else:
            self._diagnostic = None
        self.paths = paths if paths else {}
        self.schema = schema if schema else self.DEFAULT_VALIDATION_SCHEMA
        self.metadata = metadata if metadata else {}
        self.summary_cycles = summary_cycles
        self.summary_diagnostic = summary_diagnostic
    
    def __repr__(self) -> str:
        has_raw = True if self.raw else False
        has_structured = True if self.structured else False
        has_diagnostic = True if self.diagnostic else False
        from_path = self.paths.get("raw", "unknown")
        return f"{self.__class__.__name__} (" \
            f"raw={has_raw}, structured={has_structured}, diagnostic={has_diagnostic})"\
            f" from {from_path}"

    def validate(self):
        """
        Validate the run object against the validation schema.
        If a validation schema is not passed to __init__, a default is used.
        """
        logger.warning("Validation requires loading entire df into memory!")
        validator = SimpleValidator(self.schema)
        is_valid, reason = validator.validate(self.raw.cycles.data)
        return is_valid, reason

    def structure(self):
        pbar = ProgressBar(dt=1, width=10)
        pbar.register()
        cycles_interpolated = bag.from_sequence(
            bag.map(
                interpolate_cycle, 
                self.raw.cycles.items,

                # remaining kwargs are broadcast to all calls
                cconfig=self.raw.config
            ).compute()
        )
        cycles_interpolated.repartition(npartitions=cycles_interpolated.count().compute())
        cycles_interpolated.remove(lambda xdf: xdf is None)
        self.structured = CyclesContainer(cycles_interpolated)
        pbar.unregister()

    # Diagnostic config methods
    @property
    def diagnostic(self):
        return self._diagnostic

    @diagnostic.setter
    def diagnostic(self, diagnostic_config: DiagnosticConfig):

        if not isinstance(diagnostic_config, DiagnosticConfig):
            logger.warning(
                f"Diagnostic config passed does not inherit "
                "DiagnosticConfig, can cause downstream errors."
            )
        self._diagnostic = diagnostic_config

        for cycle in tqdm.tqdm(
            self.raw.cycles, 
            total=self.raw.cycles.items_length,
            desc="Updating cycle labels based on diagnostic config",
            **TQDM_STYLE_ARGS
        ):
            for step in cycle.steps:
                step.data["cycle_label"] = step.data["cycle_index"].apply(
                    lambda cix: diagnostic_config.type_by_ix.get(cix, "regular")
            )    

    @diagnostic.deleter
    def diagnostic(self):
        del self._diagnostic

    # Serialization methods
    # maybe do not need to do this for this class according to monty docstrings
    @classmethod
    def from_dict(cls, d: dict):
        """
        Create a Run object from a dictionary.
        """
        d.pop("@module")
        d.pop("@class")
        return cls(**d)

    def as_dict(self):
        """
        Convert a Run object to a dictionary.
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "raw_cycle_container": self.raw.as_dict(),
            "structured_cycle_container": self.structured.as_dict() if self.structured else None,
            "diagnostic_config": self.diagnostic.as_dict() if self.diagnostic else None,
            "metadata": self.metadata,
            "schema": self.schema,
            "paths": self.paths
        }

    # Convenience methods for loading and saving
    @classmethod
    def load(cls, path: Union[str, os.PathLike]):
        """
        Load a Run object from a file or list of files.
        """
        return loadfn(path)
    
    def save(self, path: Union[str, os.PathLike]):
        """
        Save a Run object to a file.
        """
        dumpfn(self, path)

    @classmethod
    def from_dataframe(
        cls, 
        df: pd.DataFrame, 
        **kwargs
    ):
        """
        Convenience method to create an unstructured Run object from a raw dataframe.
        """
        # Assign a per-cycle step index counter
        df.loc[:, "step_counter"] = 0
        for cycle_index in tqdm.tqdm(
            df.cycle_index.unique(),
            desc="Assigning step counter",
            **TQDM_STYLE_ARGS
            ):
            indices = df.loc[df.cycle_index == cycle_index].index
            step_index_list = df.step_index.loc[indices]
            shifted = step_index_list.ne(step_index_list.shift()).cumsum()
            df.loc[indices, "step_counter"] = shifted - 1

        # Assign an absolute step index counter
        compounded_counter = df.step_counter.astype(str) + "-" + df.cycle_index.astype(str)
        absolute_shifted = compounded_counter.ne(compounded_counter.shift()).cumsum()
        df["step_counter_absolute"] = absolute_shifted - 1

        # Assign step label if not known
        if "step_label" not in df.columns:
            df["step_label"] = None
            for sca, df_sca in tqdm.tqdm(df.groupby("step_counter_absolute"),
                                         desc="Determining charge/discharge steps",
                                         **TQDM_STYLE_ARGS):
                indices = df_sca.index
                df.loc[indices, "step_label"] = label_chg_state(df_sca)

        # Assign cycle label from diagnostic config
        df["cycle_index"] = df["cycle_index"].astype(cls.DEFAULT_DTYPES["cycle_index"])
        df["cycle_label"] = "regular"

        diagnostic = kwargs.get("diagnostic", None)
        if diagnostic:
            df["cycle_label"] = df["cycle_index"].apply(
                lambda cix: diagnostic.type_by_ix.get(cix, "regular")
            )    
        if "datum" not in df.columns:
            df["datum"] = df.index

        # Note this will not convert columns
        # not listed in the default dtypes

        dtypes = {c: dtype for c, dtype in cls.DEFAULT_DTYPES.items() if c in df.columns}
        df = df.astype(dtypes)
        raw = CyclesContainer.from_dataframe(df, tqdm_desc_suffix="(raw)")
        return cls(raw, **kwargs)

    def summarize_cycles(
            self,
            nominal_capacity=1.1,
            full_fast_charge=0.8,
            cycle_complete_discharge_ratio=0.97,
            cycle_complete_vmin=3.3,
            cycle_complete_vmax=3.3,
            error_threshold=1e6
    ):
        """
        Gets summary statistics for data according to cycle number. Summary data
        must be float or int type for compatibility with other methods

        Note: Tries to avoid loading large dfs into memory at any given time
        by using dask. Small dfs are permissible.

        Args:
            nominal_capacity (float): nominal capacity for summary stats
            full_fast_charge (float): full fast charge for summary stats
            cycle_complete_discharge_ratio (float): expected ratio
                discharge/charge at the end of any complete cycle
            cycle_complete_vmin (float): expected voltage minimum achieved
                in any complete cycle
            cycle_complete_vmax (float): expected voltage maximum achieved
                in any complete cycle
            error_threshold (float): threshold to consider the summary value
                an error (applied only to specific columns that should reset
                each cycle)

        Returns:
            (pandas.DataFrame): summary statistics by cycle.
        """
        agg = {
            "cycle_index": "first",
            "discharge_capacity": "max",
            "charge_capacity": "max",
            "discharge_energy": "max",
            "charge_energy": "max",
            "internal_resistance": "last",
            "date_time_iso": "first",
            "test_time": "first"
        }

        # Aggregate and format a potentially large dataframe
        raw_lazy_df = self.raw.cycles["regular"].data_lazy
        summary = raw_lazy_df. \
            groupby("cycle_index").agg(agg). \
            rename(columns={"internal_resistance": "dc_internal_resistance"})
        summary = summary.compute()

        summary["energy_efficiency"] = \
                summary["discharge_energy"] / summary["charge_energy"]
        summary.loc[
            ~np.isfinite(summary["energy_efficiency"]), "energy_efficiency"
        ] = np.NaN
        # This code is designed to remove erroneous energy values
        for col in ["discharge_energy", "charge_energy"]:
            summary.loc[summary[col].abs() > error_threshold, col] = np.NaN
        summary["charge_throughput"] = summary.charge_capacity.cumsum()
        summary["energy_throughput"] = summary.charge_energy.cumsum()

        # Computing charge durations
        # This method for computing charge start and end times implicitly
        # assumes that a cycle starts with a charge step and is then followed
        # by discharge step.
        charge_start_time = raw_lazy_df. \
            groupby("cycle_index")["date_time_iso"]. \
            agg("first").compute().to_frame()
        charge_finish_time = raw_lazy_df \
            [raw_lazy_df.charge_capacity >= nominal_capacity * full_fast_charge]. \
            groupby("cycle_index")["date_time_iso"]. \
            agg("first").compute().to_frame()

        # Left merge, since some cells might not reach desired levels of
        # charge_capacity and will have NaN for charge duration
        merged = charge_start_time.merge(
            charge_finish_time, on="cycle_index", how="left"
        )

        # Charge duration stored in seconds -
        # note that date_time_iso is only ~1sec resolution
        time_diff = np.subtract(
            pd.to_datetime(merged.date_time_iso_y, utc=True, errors="coerce"),
            pd.to_datetime(merged.date_time_iso_x, utc=True, errors="coerce"),
        )
        summary["charge_duration"] = np.round(
            time_diff / np.timedelta64(1, "s"), 2)

        # Compute time-temeprature integral, if available
        if "temperature" in raw_lazy_df.columns:
            # Compute time since start of cycle in minutes. This comes handy
            # for featurizing time-temperature integral
            raw_lazy_df["time_since_cycle_start"] = \
                raw_lazy_df["date_time_iso"].apply(pd.to_datetime) - \
                raw_lazy_df.groupby("cycle_index")["date_time_iso"].transform("first")

            raw_lazy_df["time_since_cycle_start"] = \
                (raw_lazy_df["time_since_cycle_start"] / np.timedelta64(1, "s")) / 60

            # Group by cycle index and integrate time-temperature
            summary["time_temperature_integrated"] = raw_lazy_df.groupby(
                "cycle_index").apply(
                lambda g: np.integrate.trapz(
                    g.temperature,
                    x=g.time_since_cycle_start
                ).compute()
            )
            raw_lazy_df.drop(columns=["time_since_cycle_start"]).compute()

        # Determine if any of the cycles has been paused
        summary["paused"] = raw_lazy_df.groupby("cycle_index").apply(
            get_max_paused_over_threshold,
            meta=pd.Series([], dtype=float)).compute()

        # Find CV step data
        cv_data = []
        for cyc in self.raw.cycles["regular"]:
            cix = cyc.cycle_index
            try:
                cv = get_cv_stats(cyc.steps["charge"].data)
            except (CVStatsError, DFSelectorIndexError):
                logger.debug(f"Cannot extract CV charge segment for cycle {cix}!")
                continue
            cv_data.append(cv)
        cv_summary = pd.DataFrame(cv_data).set_index("cycle_index")
        summary = summary.merge(
            cv_summary,
            how="outer",
            right_index=True,
            left_index=True
        ).set_index("cycle_index")

        summary = summary.astype(
            {c: v for c, v in self.SUMMARY_DTYPES.items() if c in summary.columns}
        )

        # Avoid returning empty summary dataframe for single cycle raw_data
        if summary.shape[0] == 1:
            return summary

        # Ensure final cycle has actually been completed; if not, exclude it.
        last_cycle = raw_lazy_df.cycle_index.max().compute().item()
        last_voltages = self.raw.cycles[last_cycle].data.voltage
        min_voltage_ok = last_voltages.min() < cycle_complete_vmin
        max_voltage_ok = last_voltages.max() > cycle_complete_vmax
        dchg_ratio_ok = (summary.iloc[[-1]])["discharge_capacity"].iloc[0] \
            > cycle_complete_discharge_ratio \
            * ((summary.iloc[[-1]])["charge_capacity"].iloc[0])

        if all([min_voltage_ok, max_voltage_ok, dchg_ratio_ok]):
            return summary
        else:
            return summary.iloc[:-1]