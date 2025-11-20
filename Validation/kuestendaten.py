import re
import pandas as pd
from datetime import datetime
from io import StringIO

class KuestenDataFile:
    """
    Parser for ASCII time-series files downloaded from the Kuestendaten.de portal
    (WISKI coastal water monitoring system).

    This class reads a standard Kuestendaten text export and provides:

    - Metadata as a dictionary (`.metadata`)
    - Metadata fields also accessible as Python attributes (e.g. `a.East`, `a.Parameter`)
    - Time-series data as a pandas DataFrame (`.data`)
    - Automatic parsing of timestamps, values, status, and remark columns

    File format structure expected
    ------------------------------
    The input ASCII file must follow this pattern:

        Key : Value
        Key : Value
        ... (multiple header metadata lines)

        ================================================
        Datum/Uhrzeit   Parameter [unit]    Status   Bemerkung
        2013-05-02 10:00:00   0.5    240
        2013-05-02 10:05:00   0.5    240
        ...

    Metadata lines prior to the separator (`=====`) are extracted into `metadata`.
    Table lines after the separator are parsed into a pandas DataFrame.

    Attribute conversion
    --------------------
    Metadata keys are normalized so that spaces and special characters become `_`,
    allowing attribute access. For example:

        "WGS84 lat"  →  obj.WGS84_lat
        "Erhebende Organisation" → obj.Erhebende_Organisation

    Numeric metadata values are automatically converted to int/float when possible.

    Attributes
    ----------
    filepath : str
        Path to the source ASCII file.
    metadata : dict
        Dictionary of parsed metadata entries.
    data : pandas.DataFrame
        Parsed time-series data with columns: timestamp, value, status, remark.

    Methods
    -------
    summary()
        Print basic file summary: station, parameter, row count, etc.
    to_csv(path)
        Save parsed time-series to CSV.

    Usage example
    -------------
    >>> a = KuestenDataFile("ABBENFLETH_SPERRWERK!Salzgehalt.txt")

    # Access metadata
    >>> a.Stationsname
    'Abbenfleth Sperrwerk'

    >>> a.Parameter, a.Parameter_Einheit
    ('Salzgehalt', 'pmill')

    >>> a.East, a.North
    (532690, 5947800)

    # Access time-series data
    >>> a.data.head()
           timestamp  value  status remark
    0 2013-05-02 ...   0.5     240   NaN

    # Export
    >>> a.to_csv("abbenfleth_salinity.csv")

    Notes
    -----
    This reader targets WISKI-formatted text exports provided by
    https://www.kuestendaten.de and used in the German coastal water network.
    """

    def __init__(self, filepath):
        self.filepath = filepath
        self.metadata = {}
        self.data = None

        lines = self._read_file(filepath)
        self._parse(lines)
        self._attach_metadata_as_attributes()
        
    def _read_file(self, filepath):
        encodings = ["utf-8", "cp1252", "latin-1"]

        for enc in encodings:
            try:
                with open(filepath, "r", encoding=enc) as f:
                    return f.readlines()
            except UnicodeDecodeError:
                continue

        raise UnicodeDecodeError(
            "Failed to read file with encodings: " + ", ".join(encodings)
        )
        

    def _parse(self, lines):
        header_lines = []
        table_lines = []
        separator_found = False

        # Normalize line endings
        for line in lines:
            line = line.rstrip("\r\n")
            if line.strip().startswith("==="):
                separator_found = True
                continue

            if not separator_found:
                header_lines.append(line)
            else:
                table_lines.append(line)

        self._parse_metadata_header(header_lines)
        self._parse_table(table_lines)

    def _parse_metadata_header(self, lines):
        for line in lines:
            if ":" not in line:
                continue
            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()
            # normalize key to python attribute style
            key_attr = re.sub(r"[^0-9a-zA-Z_]", "_", key).strip("_")
            if key:
                self.metadata[key_attr] = value

    def _attach_metadata_as_attributes(self):
        """
        Attach metadata items as object attributes, e.g.:
        self.metadata["East"] → self.East
        """
        for key, value in self.metadata.items():
            # attempt numeric conversion where reasonable
            try:
                if "." in value:
                    val = float(value)
                else:
                    val = int(value)
            except:
                val = value
            setattr(self, key, val)

    def _parse_table(self, lines):
        lines = [line for line in lines if line.strip()]
        text = "\n".join(lines)

        df = pd.read_csv(
            StringIO(text),
            sep=r"\s{2,}|\t",
            engine="python"
        )

        df.columns = [c.strip() for c in df.columns]

        # Timestamp column (first)
        datetime_col = df.columns[0]
        df[datetime_col] = pd.to_datetime(df[datetime_col], errors='coerce')

        # Value column
        value_col = df.columns[1]

        df.rename(columns={
            datetime_col: "timestamp",
            value_col: "value",
            "Status": "status",
            "Bemerkung": "remark"
        }, inplace=True)

        self.data = df

    def summary(self):
        print("File:", self.filepath)
        print("Rows:", len(self.data))
        print("Columns:", list(self.data.columns))
        print("--- Station:", getattr(self, "Stationsname", None))
        print("Parameter:", getattr(self, "Parameter", None),
              getattr(self, "Parameter_Einheit", None))

    def to_csv(self, outfile):
        self.data.to_csv(outfile, index=False)



def _read_file(self, filepath):
    encodings = ["utf-8", "cp1252", "latin-1"]

    for enc in encodings:
        try:
            with open(filepath, "r", encoding=enc) as f:
                return f.readlines()
        except UnicodeDecodeError:
            continue

    raise UnicodeDecodeError(
        "Failed to read file with encodings: " + ", ".join(encodings)
    )
