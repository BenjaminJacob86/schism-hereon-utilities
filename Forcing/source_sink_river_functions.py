#### Collection of functions for SCHISM to write msource, vsource and source_sink.in

import numpy as np
from typing import Optional, List

def create_source_sink(s, coords, name='source_sink.in',mindepth=3.0, elements=None,numcells=5):
    """
    Write source_sink index list for SCHISM.

    Parameters
    ----------
    s : schism_setup
        SCHISM object.
    coords : list of (lon, lat) tuples
        List of river coordinates, e.g. [(lon1, lat1), (lon2, lat2), ...].
    name : str, optional
        Output filename (default: 'source_sink.in').
    elements : list, optional
        Predefined list of SCHISM element IDs (must match length of coords).
    numcells : int
        number of grid elements to return as neighbours.
    return elements: 1 based list of elements select for source application    
    """
    s.init_element_tree(latlon=True)

    with open(name, 'w') as f:
        n_sources = len(coords) if elements is None else len(elements)
        f.write(f"{n_sources}   ! number of elements with sources\n")

        if elements is None:
            # Find nearest element for each coordinate
            elements=[]
            for i, (lon, lat) in enumerate(coords, 1):
                print(i,lon,lat)
                el_id = s.find_nearest_element(lon, lat,mindepth=mindepth,numcells=numcells)
                el_depth = s.element_depth[el_id]
                f.write(f"{el_id}    ! river{i}, element depth: {el_depth:.1f} m\n")
                elements.append(el_id)
        else:
            # Predefined list of element IDs
            for i, el_id in enumerate(elements, 1):
                el_depth = s.element_depth[el_id]
                f.write(f"{el_id}    ! river{i}, element depth: {el_depth:.1f} m\n")

        f.write('\n')
        f.write('0   ! number of elements with sinks\n')
        
    return elements            

def create_vsource(time_in_seconds: np.ndarray, Q: np.ndarray, name: str = "vsource.th") -> None:
    """
    Create time x source entries for a source file.

    Parameters
    ----------
    time_in_seconds : np.ndarray
        1D array of times (seconds since model start).
    Q : np.ndarray
        1D or 2D array of discharges. If 2D, shape should be (n_times, n_sources).
    name : str, optional
        Output filename, by default "vsource.th".
    """
    # Ensure Q is 2D
    Q = np.atleast_2d(Q)
    if Q.shape[0] != len(time_in_seconds):
        Q = Q.T  # adjust if given as (n_sources, n_times)

    # Stack time + discharge matrix
    M = np.column_stack((time_in_seconds, Q))

    # Format: first column int, rest floats with 2 decimals
    fmt = ["%d"] + ["%.2f"] * Q.shape[1]

    # Write file
    np.savetxt(name, M, fmt=fmt,  header = "!time (s), sources (m^3/s) at elem. #1,2 (see source_sink.in)")

def create_msource(
    time_in_seconds: np.ndarray,
    temperature: np.ndarray,
    salinity: np.ndarray,
    others: Optional[List[np.ndarray]] = None,
    name: str = "msource.th"
) -> None:
    """
    Create time x source entries for an msource file (T, S, and optional tracers).

    Parameters
    ----------
    time_in_seconds : np.ndarray
        1D array of times (seconds since model start).
    temperature : np.ndarray
        1D or 2D array of temperatures (n_times, n_sources).
    salinity : np.ndarray
        1D or 2D array of salinities (n_times, n_sources).
    others : list of np.ndarray, optional
        List of additional arrays (e.g., tracers). Each must be (n_times, n_sources).
    name : str, optional
        Output filename, by default "msource.th".
    -9999 for tracers to inject ambient values in tracer matirces
    """
    n_times = len(time_in_seconds)

    def ensure_2d(arr: np.ndarray, name: str) -> np.ndarray:
        arr = np.atleast_2d(arr)
        if arr.shape[0] != n_times:
            arr = arr.T
        if arr.shape[0] != n_times:
            raise ValueError(f"{name} must have {n_times} rows (one per time step).")
        return arr

    # Ensure T and S are consistent
    T = ensure_2d(temperature, "temperature")
    S = ensure_2d(salinity, "salinity")

    n_sources = T.shape[1]
    if S.shape[1] != n_sources:
        raise ValueError("Temperature and salinity must have the same number of sources (columns).")

    # Collect blocks
    blocks = [T, S]

    # Add optional others
    if others is not None:
        for k, arr in enumerate(others, start=1):
            arr = ensure_2d(arr, f"others[{k}]")
            if arr.shape[1] != n_sources:
                raise ValueError("All variables must have the same number of sources (columns).")
            blocks.append(arr)

    # Final stacked matrix
    M = np.column_stack([time_in_seconds] + blocks)

    # Format: int for time, floats for others
    fmt = ["%d"] + ["%.2f"] * (M.shape[1] - 1)

    # Build header
    header_parts = ["time (s)"]
    header_parts += [f"T at elem. #{i+1}" for i in range(n_sources)]
    header_parts += [f"S at elem. #{i+1}" for i in range(n_sources)]
    if others is not None:
        for k, arr in enumerate(others, start=1):
            header_parts += [f"tracer{k} at elem. #{i+1}" for i in range(n_sources)]
    header = "!" + ", ".join(header_parts) + " (see source_sink.in)"

    # Save file
    np.savetxt(name, M, fmt=fmt, header=header, comments="")





