#!/usr/bin/env python3
""" -- extract the counts matrix from IMvigor210CoreBiologies
mirror cds.RData into a plain CSV (entrez_id rows x SAM sample columns).

Uses a TOLERANT rdata converter (NotImplementedError -> None) because the cds
object tree contains promise/closure objects from the package namespace that
pure-python rdata cannot convert. R matrices with dimnames become xarray
DataArrays (values reshaped, dimnames as coordinates); we search the converted
object for shape (G, 348) with G > 20000.

Usage: cds_extract.py <file.RData> <out.csv>
"""
import sys
import numpy as np
import rdata
from rdata.conversion import SimpleConverter


class TolerantConverter(SimpleConverter):
    def _convert_next(self, obj):
        try:
            return super()._convert_next(obj)
        except NotImplementedError:
            return None
        except RecursionError:
            return None


def walk(obj, found, path="root", depth=0, seen=None):
    if depth > 20:
        return
    if seen is None:
        seen = set()
    if id(obj) in seen:
        return
    seen.add(id(obj))
    try:
        if hasattr(obj, "dims") and hasattr(obj, "values"):  # xarray.DataArray
            if len(obj.shape) == 2 and obj.shape[1] == 348 and obj.shape[0] > 20000:
                found.append((path, obj))
            return
        if isinstance(obj, np.ndarray) and obj.ndim == 2 and obj.shape[1] == 348 and obj.shape[0] > 20000:
            found.append((path, obj))
            return
    except Exception:
        pass
    if isinstance(obj, dict):
        for k, v in obj.items():
            walk(v, found, f"{path}.{k}", depth + 1, seen)
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            walk(v, found, f"{path}[{i}]", depth + 1, seen)
    else:
        # SimpleNamespace (S4 slots) and Mapping-like (REnvironment)
        if hasattr(obj, "__dict__"):
            for k, v in vars(obj).items():
                walk(v, found, f"{path}.{k}", depth + 1, seen)
        try:
            for k, v in obj.items():
                walk(v, found, f"{path}{{{k}}}", depth + 1, seen)
        except Exception:
            pass


def main():
    src, out = sys.argv[1], sys.argv[2]
    p = rdata.parser.parse_file(src)
    conv = TolerantConverter()
    c = conv.convert(p.object) if hasattr(conv, "convert") else conv._convert_next(p.object)
    # normalize: SimpleConverter.convert returns dict of top-level symbols
    print("parsed & converted:", src, file=sys.stderr)
    try:
        keys = list(c.keys())
        print("top-level keys:", keys, file=sys.stderr)
    except Exception:
        pass

    found = []
    walk(c, found)
    print(f"candidates shape (G,348): {len(found)}", file=sys.stderr)
    for path, obj in found[:6]:
        print("  ", path, obj.shape, type(obj).__name__, file=sys.stderr)
    if not found:
        print("NO (G,348) matrix found", file=sys.stderr)
        sys.exit(1)

    path, obj = max(found, key=lambda x: x[1].shape[0])
    G, S = obj.shape
    print("chosen:", path, (G, S), file=sys.stderr)

    mat = np.asarray(obj.values if hasattr(obj, "values") else obj)
    rn = None
    cn = None
    if hasattr(obj, "coords") and hasattr(obj, "dims"):
        try:
            rn = [str(x) for x in obj.coords[obj.dims[0]].values]
            cn = [str(x) for x in obj.coords[obj.dims[1]].values]
        except Exception:
            rn = cn = None
    if rn is None or cn is None or len(rn) != G or len(cn) != S:
        print("dimnames missing/incomplete; cannot label rows", file=sys.stderr)
        sys.exit(2)
    print("rownames:", len(rn), rn[:5], "... colnames:", len(cn), cn[:3], file=sys.stderr)

    if mat.dtype.kind == "f":
        mat = np.rint(mat)
    mat = mat.astype(np.int64)
    with open(out, "w") as f:
        f.write("entrez," + ",".join(cn) + "\n")
        for e, row in zip(rn, mat):
            f.write(str(e) + "," + ",".join(str(int(x)) for x in row) + "\n")
    print("wrote", out, f"({G} x {S})", file=sys.stderr)


if __name__ == "__main__":
    main()
